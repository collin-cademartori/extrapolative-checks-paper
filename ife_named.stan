functions {

  matrix first_diffs(matrix Y) {
    int T = rows(Y);
    int N = cols(Y);
    matrix[T, N] Ydiffs;

    Ydiffs[1, :] = Y[1, :];
    Ydiffs[2:T, :] = Y[2:T, :] - Y[1:(T-1), :];

    return(Ydiffs);
  }

  vector ar_process(vector innovations,
                    real autocor,
                    real scale,
                    real mean) {

    int T = num_elements(innovations);
    vector[T] ar_process;
    // Scale the innovations so the AR(1) has unit marginal variance (stationary;
    // paper Section 2).
    real scale_inno = sqrt(1 - square(autocor));

    ar_process[1] = innovations[1];
    for(t in 2:T) {
      ar_process[t] = (autocor * ar_process[t-1])
                      + (scale_inno * innovations[t]);
    }

    ar_process = mean + scale * ar_process;

    return ar_process;
  }

  // Pearson correlation between two equal-length vectors.
  real pearson_cor(vector a, vector b) {
    vector[num_elements(a)] a_centered = a - mean(a);
    vector[num_elements(b)] b_centered = b - mean(b);
    return dot_product(a_centered, b_centered)
           / sqrt(dot_self(a_centered) * dot_self(b_centered));
  }

}

data {

  int<lower=1> M_units;
  int<lower=1> T_times;
  int<lower=1, upper=M_units> K_latent;

  matrix[T_times, M_units] Y;

  real<lower=0> a_rho;
  real<lower=0> b_rho;

  int<lower=0, upper=1> fit_overall_scales;
  vector<lower=0>[M_units] sigma_data;

  int<lower=0, upper=1> nonstationary;
  int<lower=0, upper=1> unit_intercepts;
  int<lower=0, upper=1> factor_means;
  int<lower=0, upper=1> sample_posterior;

  int<lower=0, upper=T_times> num_treated;
  real<lower=0> gamma_scale;
  real gamma_loc;

  real<lower=0> tau_val;
  real<lower=0> m_tau;
  real<lower=0> s_tau;

}

transformed data {

  vector[T_times] times;
  for(t in 1:T_times) {
    times[t] = t;
  }

  matrix[T_times, M_units] Y_outcome;
  cov_matrix[T_times] errors_cov;
  cov_matrix[T_times] errors_precision;

  if(nonstationary) {
    Y_outcome = first_diffs(Y);

    // errors_cov is the differenced-error covariance divided by 2v, where v is the
    // iid level-error variance (square(tau*sigma[n])). Differencing iid level
    // errors gives d_1 = eps_1 (the level, variance v) and d_t = eps_t - eps_{t-1}
    // for t>1 (variance 2v, lag-1 covariance -v). Scaling this matrix by 2v (as
    // Y_prior and error_precisions do) reproduces that exactly: the [1,1]=0.5
    // corner encodes the level's half-variance; the interior is the MA(1)
    // correlation (1 on the diagonal, -0.5 off).
    errors_cov[1,1] = 0.5;
    errors_cov[1,2] = -0.5;
    errors_cov[2,1] = -0.5;
    errors_cov[1, 3:T_times] = rep_row_vector(0, T_times - 2);
    errors_cov[3:T_times, 1] = rep_vector(0, T_times - 2);
    for(s in 2:T_times) {
      for(t in 2:T_times) {
        errors_cov[s, t] = (s == t) ? 1 : ((abs(s-t) == 1) ? -0.5 : 0);
      }
    }
    errors_precision = inverse_spd(errors_cov);
  } else {
    Y_outcome = Y;

    errors_precision = identity_matrix(T_times);
    errors_cov = errors_precision;
  }

  // Treatment-effect prior precision. delta_raw enters Y_means on the model's
  // working (differenced, when nonstationary) scale, so it borrows the error
  // covariance structure -- the leading num_treated block of errors_cov, whose
  // level corner seeds the integrated delta = cumsum(delta_raw). Scaled to the
  // treated unit's error scale by sigma[1] in the prior below.
  matrix[num_treated, num_treated] effects_prec = inverse_spd(errors_cov[1:num_treated, 1:num_treated]);


}

parameters {

  vector<lower=0>[fit_overall_scales == 1 ? M_units : 0] sigma_raw;
  vector<lower=0>[tau_val > 0 ? 0 : 1] tau_param;

  matrix[T_times, K_latent] Phi_innovations;
  vector[factor_means == 1 ? K_latent : 0] Phi_means_param;
  vector<lower=0, upper=1>[K_latent] rho;

  // Lower-triangular loadings (cholesky_factor_cov) fix the factor rotation for
  // identifiability: unit n loads only on the first min(K, n) factors.
  cholesky_factor_cov[M_units, K_latent] Lambda;

  vector[unit_intercepts ? M_units : 0] gamma_raw;

  vector[num_treated] delta_raw;

}

transformed parameters {

  vector[M_units] sigma;
  if(fit_overall_scales == 1) {
    sigma = sigma_data .* sigma_raw;
  } else {
    sigma = sigma_data;
  }

  vector[K_latent] Phi_means;
  if(factor_means) {
    Phi_means = gamma_loc + gamma_scale * Phi_means_param;
  } else {
    Phi_means = rep_vector(0, K_latent);
  }

  real<lower=0> tau;
  if(tau_val > 0) {
    tau = tau_val;
  } else {
    tau = tau_param[1];
  }

  vector<lower=0>[M_units] error_precisions = square(inv(tau * sigma));
  if(nonstationary) {
    // Differenced errors have variance 2v (see errors_cov); halve the precision.
    error_precisions = 0.5 * error_precisions;
  }

  matrix[T_times, K_latent] Phi;
  for(k in 1:K_latent) {
    Phi[:,k] = ar_process(Phi_innovations[:,k], rho[k], 1, Phi_means[k]);
  }

  vector[M_units] gamma;
  if(unit_intercepts) {
    gamma = gamma_loc + gamma_scale * gamma_raw;
  } else {
    gamma = rep_vector(0, M_units);
  }

  matrix[T_times, M_units] Lambda_Phi = Phi * (Lambda');

  // Each unit's outcome mean is sigma[n] * (intercept + latent component), so the
  // intercept gamma[n] is scaled per-unit by that unit's overall scale sigma[n]
  // (commensurate with how the loadings and errors are scaled). In the
  // nonstationary (differenced) branch the constant intercept survives
  // differencing only in the first (level) element, and reappears at every time
  // once Y is reconstructed by cumulative_sum.
  matrix[T_times, M_units] Y_means_0;
  matrix[T_times, M_units] Y_means;
  if(nonstationary) {
    for(n in 1:M_units) {
      Y_means_0[:,n] = (sigma[n]) * Lambda_Phi[:,n];
      Y_means_0[1,n] += gamma[n];
    }
    Y_means = Y_means_0;
    if(num_treated > 0) {
      Y_means[(T_times - num_treated + 1):T_times, 1] += delta_raw;
    }
  } else {
    for(n in 1:M_units) {
      Y_means_0[:,n] = gamma[n] + (sigma[n] * Lambda_Phi[:,n]);
    }

    Y_means = Y_means_0;
    if(num_treated > 0) {
      Y_means[(T_times - num_treated + 1):T_times, 1] += delta_raw;
    }
  }

}

model {

  if(fit_overall_scales == 1) {
    sigma_raw ~ normal(0, 5);
  }

  if(tau_val == 0) {
    tau_param ~ normal(m_tau, s_tau);
  }

  to_vector(Phi_innovations) ~ std_normal();
  rho ~ beta(a_rho, b_rho);

  // Loadings prior (paper eq. loadings_prior). The per-unit scale sigma[n] is
  // applied downstream in Y_means, so the loadings as they enter the outcome
  // have scale sigma[n] / sqrt(min(K, n)) -- the paper's sigma_i / sqrt(min(K, i)).
  // The positive diagonal comes from the cholesky_factor_cov type, making the
  // diagonal prior a half-normal (the paper's normal_+).
  for(n in 1:M_units) {
    Lambda[n,1:min(K_latent, n)] ~ normal(0, 1 / sqrt(min(K_latent, n)));
  }

  gamma_raw ~ std_normal();
  Phi_means_param ~ std_normal();

  if(num_treated > 0) {
    delta_raw ~ multi_normal_prec(rep_vector(0, num_treated), square(inv(sigma[1])) * effects_prec);
  }

  if(sample_posterior) {
    for(n in 1:M_units) {
      Y_outcome[:,n] ~ multi_normal_prec(Y_means[:,n], error_precisions[n] * errors_precision);
    }
  }

}

generated quantities {

  matrix[T_times, M_units] Y_prior;

  if(nonstationary) {
    for(n in 1:M_units) {
      Y_prior[:,n] = cumulative_sum(multi_normal_rng(Y_means[:,n], 2 * square(tau * sigma[n]) * errors_cov));
    }
  } else {
    for(n in 1:M_units) {
      Y_prior[:,n] =  multi_normal_rng(Y_means[:,n], square(tau * sigma[n]) * errors_cov);
    }
  }

  matrix[T_times, M_units] Y_latent;
  if(nonstationary) {
    for(n in 1:M_units) {
      Y_latent[:,n] = cumulative_sum(Y_means[,n]);
    }
  } else {
    for(n in 1:M_units) {
      Y_latent[:,n] = Y_means[,n];
    }
  }

  matrix[T_times, M_units] Y_pred;
  for(n in 1:M_units) {
    Y_pred[:,n] = Y_means[,n] + multi_normal_rng(rep_vector(0, T_times), inv(error_precisions[n]) * errors_cov);
  }
  if(nonstationary) {
    for(n in 1:M_units) {
      Y_pred[:,n] = cumulative_sum(Y_pred[,n]);
    }
  }

  // Statistic S1 (paper Section 5): the absolute correlation of the first untreated
  // unit's predictive replicate with time, a summary of monotone trend.
  real time_cor_pred = abs(pearson_cor(Y_pred[:, 2], times));

  // Statistic S2 (paper Section 5): the association across untreated units
  // between (i) each unit's correlation with the treated unit and (ii) its
  // long-run location (mean level). Evaluated on the predictive replicate
  // Y_pred over the pre-treatment window only, so the treated unit's treatment
  // period does not distort the correlations.
  real loc_cor_pred;
  {
    int T_pre = T_times - num_treated;
    vector[M_units - 1] cor_with_treated;
    vector[M_units - 1] unit_location;
    for(n in 2:M_units) {
      cor_with_treated[n - 1] = pearson_cor(Y_pred[1:T_pre, n], Y_pred[1:T_pre, 1]);
      unit_location[n - 1] = mean(Y_pred[1:T_pre, n]);
    }
    loc_cor_pred = pearson_cor(cor_with_treated, unit_location);
  }

  // Model-implied absolute cross-sectional correlation between the treated unit
  // and each other unit: the factor covariance (Lambda_tr . Lambda_un) over the
  // product of total standard deviations (factor variance plus idiosyncratic
  // variance var = 1 / error_precisions).
  vector[M_units - 1] abs_cors;
  {
    row_vector[K_latent] Lambda_tr = Lambda[1,:];
    real var_tr = inv(error_precisions[1]);
    for(n in 2:M_units) {
      row_vector[K_latent] Lambda_un = Lambda[n,:];
      real var_un = inv(error_precisions[n]);
      abs_cors[n-1] = abs(dot_product(Lambda_tr, Lambda_un)) / sqrt((dot_self(Lambda_tr) + var_tr) * (dot_self(Lambda_un) + var_un));
    }
  }

  real mean_abs_diffs = mean(abs(to_vector(Y_latent - Y)));

  vector[num_treated] delta;
  if(nonstationary) {
    delta = cumulative_sum(delta_raw);
  } else {
    delta = delta_raw;
  }

}
