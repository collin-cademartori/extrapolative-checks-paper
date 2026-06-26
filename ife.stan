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
    real scale_inno = sqrt(1 - square(autocor));
    
    ar_process[1] = innovations[1];
    for(t in 2:T) {
      ar_process[t] = (autocor * ar_process[t-1]) 
                      + (scale_inno * innovations[t]);
    }

    ar_process = mean + scale * ar_process;

    return ar_process;
  }

  real ns_scale_factor(int max_T, real a, real b) {

    real factor = max_T;

    vector[max_T - 1] ls = linspaced_vector(max_T - 1, 0, max_T - 2);
    print(ls[max_T-1]);
  
    for(k in 1:(max_T - 1)) {
      factor += 2 * (max_T - k) * prod((a + ls[1:k]) ./ (a + b + ls[1:k]));
    }

    factor = sqrt(factor);

    return factor;
  }

}

data {

  int<lower=1> N_units;
  int<lower=1> T_times;
  int<lower=1, upper=N_units> K_latent;

  matrix[T_times, N_units] Y;

  real<lower=0> autocor_a;
  real<lower=0> autocor_b;

  int<lower=0, upper=1> fit_overall_scales;
  vector<lower=0>[N_units] overall_scales_0;

  int<lower=0, upper=1> nonstationary;
  int<lower=0, upper=1> unit_intercepts;
  int<lower=0, upper=1> sample_posterior;

  int<lower=1, upper=T_times> T_ref;

  int<lower=0, upper=T_times> num_treated;
  real<lower=0> intercept_scale;

  real<lower=0> err_scale_val;
  real<lower=0> err_scale_mean;
  real<lower=0> err_scale_sd;

  //vector[K_latent] autocors_fixed;

}

transformed data {

  matrix[T_times, N_units] Y_outcome;
  cov_matrix[T_times] errors_cov;
  cov_matrix[T_times] errors_precision;

  if(nonstationary) {
    Y_outcome = first_diffs(Y);

    errors_cov[1,1] = 1;
    errors_cov[1,2] = -1.0/sqrt(2);
    errors_cov[2,1] = -1.0/sqrt(2);
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

  matrix[num_treated, num_treated] effects_prec = inverse_spd(errors_cov[1:num_treated, 1:num_treated]);
  

}

parameters {

  vector<lower=0>[fit_overall_scales == 1 ? N_units : 0] overall_scales_param;
  vector<lower=0>[err_scale_val > 0 ? 0 : 1] err_scale_param;

  matrix[T_times, K_latent] factor_innovations;
  vector<lower=0, upper=0.98>[K_latent] autocors;

  cholesky_factor_cov[N_units, K_latent] loadings;

  vector[unit_intercepts ? N_units : 0] intercepts_if_modeled;

  vector[num_treated] treatment_effects_0;

}

transformed parameters {

  vector[N_units] overall_scales;
  if(fit_overall_scales == 1) {
    overall_scales = overall_scales_0 .* overall_scales_param;
  } else {
    overall_scales = overall_scales_0;
  }

  real<lower=0> err_scale;
  if(err_scale_val > 0) {
    err_scale = err_scale_val;
  } else {
    err_scale = err_scale_param[1];
  }

  vector<lower=0>[N_units] unit_error_precs = square(inv(err_scale * overall_scales));
  if(nonstationary) {
    unit_error_precs = 0.5 * unit_error_precs;
  }

  matrix[T_times, K_latent] factors;
  for(k in 1:K_latent) {
    factors[:,k] = ar_process(factor_innovations[:,k], autocors[k], 1, 0);
  }

  vector[N_units] intercepts;
  if(unit_intercepts) {
    intercepts = intercept_scale * intercepts_if_modeled;
  } else {
    intercepts = rep_vector(0, N_units);
  }

  matrix[T_times, N_units] latent_component = factors * (loadings');

  matrix[T_times, N_units] Y_means_0;
  matrix[T_times, N_units] Y_means;
  if(nonstationary) {
    for(n in 1:N_units) {
      Y_means_0[:,n] = (overall_scales[n]) * latent_component[:,n];
      Y_means_0[1,n] += overall_scales[n] * intercepts[n];
    }
    Y_means = Y_means_0;
    if(num_treated > 0) {
      Y_means[(T_times - num_treated + 1):T_times, 1] += treatment_effects_0;
    }
  } else {
    for(n in 1:N_units) {
      Y_means_0[:,n] = (overall_scales[n] * intercepts[n]) + (overall_scales[n] * latent_component[:,n]);
    }

    Y_means = Y_means_0;
    if(num_treated > 0) {
      Y_means[(T_times - num_treated + 1):T_times, 1] += treatment_effects_0;
    }
  }

}

model {

  if(fit_overall_scales == 1) {
    overall_scales_param ~ cauchy(2, 5);
  }

  if(err_scale_val == 0) {
    err_scale_param ~ normal(err_scale_mean, err_scale_sd);
  }

  to_vector(factor_innovations) ~ std_normal();
  autocors ~ beta(autocor_a, autocor_b);

  for(n in 1:N_units) {
    loadings[n,1:min(K_latent, n)] ~ normal(0, 1 / sqrt(min(K_latent, n)));
  }

  intercepts_if_modeled ~ std_normal();

  if(num_treated > 0) {
    treatment_effects_0 ~ multi_normal_prec(rep_vector(0, num_treated), square(inv(sd(Y[:,1]))) * effects_prec);
    //treatment_effects_0 ~ normal(0, 100);
  }

  if(sample_posterior) {
    for(n in 1:N_units) {
      Y_outcome[:,n] ~ multi_normal_prec(Y_means[:,n], unit_error_precs[n] * errors_precision);
    }
  }

}

generated quantities {

  matrix[T_times, N_units] Y_prior;

  if(nonstationary) {
    for(n in 1:N_units) {
      Y_prior[:,n] = cumulative_sum(multi_normal_rng(Y_means[:,n], 2 * square(err_scale * overall_scales[n]) * errors_cov));
    }
  } else {
    for(n in 1:N_units) {
      Y_prior[:,n] =  multi_normal_rng(Y_means[:,n], square(err_scale * overall_scales[n]) * errors_cov);
    }
  }

  matrix[T_times, N_units] Y_latent;
  if(nonstationary) {
    for(n in 1:N_units) {
      Y_latent[:,n] = cumulative_sum(Y_means[,n]);
    }
  } else {
    for(n in 1:N_units) {
      Y_latent[:,n] = Y_means[,n];
    }
  }

  matrix[T_times, N_units] Y_pred;
  for(n in 1:N_units) {
    Y_pred[:,n] = Y_means[,n] + multi_normal_rng(rep_vector(0, T_times), inv(unit_error_precs[n]) * errors_cov);
  }
  if(nonstationary) {
    for(n in 1:N_units) {
      Y_pred[:,n] = cumulative_sum(Y_pred[,n]);
    }
  }

  matrix[T_times, N_units] Y0_pred;
  for(n in 1:N_units) {
    Y0_pred[:,n] = Y_means_0[,n] + multi_normal_rng(rep_vector(0, T_times), inv(unit_error_precs[n]) * errors_cov);
  }
  if(nonstationary) {
    for(n in 1:N_units) {
      Y0_pred[:,n] = cumulative_sum(Y0_pred[,n]);
    }
  }

  vector[N_units - 1] abs_cors;
  {
    // vector[T_times] Yt = Y0_pred[:,1];
    row_vector[K_latent] Ltr = loadings[1,:];
    real vtr = inv(unit_error_precs[1]);
    for(n in 2:N_units) {
      row_vector[K_latent] Lun = loadings[n,:];
      real vun = inv(unit_error_precs[n]);
      abs_cors[n-1] = abs(dot_product(Ltr, Lun)) / sqrt((dot_self(Ltr) + vtr) * (dot_self(Lun) + vun));
      // vector[T_times] Yu = Y0_pred[:,n];
      // abs_cors[n - 1] = abs(dot_product(Yt, Yu) - (T_times * mean(Yt) * mean(Yu))) / (T_times * sd(Yt) * sd(Yu));
    }
  }

  real mean_abs_diffs = mean(abs(to_vector(Y_latent - Y)));

  vector[num_treated] treatment_effects;
  if(nonstationary) {
    treatment_effects = cumulative_sum(treatment_effects_0);
  } else {
    treatment_effects = treatment_effects_0;
  }

}
