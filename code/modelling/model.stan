
functions {

  vector cumulative_prod(vector x) {
    vector[size(x)] y = exp(cumulative_sum(log(x)));
    return y;
  }
  
  vector e_x(vector M_x, vector agegaps) {
    int n_ages = size(agegaps);
    vector[n_ages] a_x;
    vector[n_ages] q_x;
    vector[n_ages] p_x;
    vector[n_ages] l_x;
    vector[n_ages] d_x;
    vector[n_ages] L_x;
    vector[n_ages] T_x;
    vector[n_ages] ex;
    
    a_x = append_row(append_row([0.07 + 1.7*M_x[1], 1.5]', rep_vector(2.5, n_ages-3)), 1/M_x[n_ages]); 
    q_x = agegaps .* M_x ./ (1 + (agegaps - a_x) .* M_x);
    // spot fix for problematic age group
    q_x[n_ages-1] = 1 - exp(-agegaps[n_ages-1]*M_x[n_ages-1]);
    p_x = 1 - q_x;
    l_x = append_row(1, cumulative_prod(p_x[1:(n_ages-1)])); 
    d_x = l_x - append_row(l_x[2:n_ages], 0);
    L_x = agegaps .* append_row(l_x[2:n_ages], 0) + (a_x .* d_x);
    T_x = reverse(cumulative_sum(reverse(L_x)));
    ex = T_x ./ l_x;
    
    return ex;
  
  }
  
  array[] real rep_array_1d() {
    
    vector[3] x;
    x[1] = 1;
    x[2] = 2;
    x[3] = 3;
    
    return to_array_1d(rep_matrix(x, 4));
    
  }
  
  vector horseshoe(vector z, vector lambda, real tau, real c2) {
    int K = rows(z);
    vector[K] lambda2 = square(lambda);
    vector[K] lambda_tilde = sqrt(c2 * lambda2 ./ (c2 + tau^2 * lambda2));
    return z .* lambda_tilde * tau;
  }
  
}

data {
  // constants ----------
  int<lower=0> n_countries;
  int<lower=0> n_periods;
  int<lower=0> n_ages;
  
  // design matrix stuff -----
  int<lower=0> N_long;
  
  int<lower=0> n_groups;
  array[N_long] int<lower=1, upper=n_groups> group_long;
  array[n_groups, n_ages*n_periods] int<lower=1, upper=N_long> group_ind;
  
   // training data ----------
  int<lower=0> N_fam; // number of country-gender-year combinations
  array[N_fam, n_ages] int<lower=0> d_fam;
  matrix<lower=0>[N_fam, n_ages] P_fam;
  int<lower=0> n_valid;
  array[n_valid] int<lower=1, upper=N_fam> valid_rows;
  array[N_fam] int<lower=1, upper=n_countries> country_fam;
  array[N_fam] int<lower=1, upper=n_periods> year_fam;
  array[N_fam] int<lower=1, upper=2> gender_fam;
  array[N_fam] int<lower=1, upper=n_groups> group_fam;
  
  array[N_long] int<lower=1, upper=2> gender_long;
  array[N_long %/% 2] int<lower=1, upper=N_long> male_long_ind;
  array[N_long %/% 2] int<lower=1, upper=N_long> female_long_ind;
  
  array[n_groups] int <lower=1, upper=2> gender_group;
  
  int<lower=0, upper=N_fam> N_hmd;
  array[N_hmd, n_ages] int<lower=0> d_hmd;
  matrix<lower=0>[N_hmd, n_ages] P_hmd;
//  array[N_hmd] int<lower=1, upper=n_countries> country_hmd;
//  array[N_fam] int<lower=1, upper=n_periods> year_hmd;
  array[N_hmd] int<lower=1, upper=N_fam> ii; // correspondence between HMD and FamiLinx data
  
  // age/period spline objects for mortality rate
  int<lower=0> k_pc;
  matrix[k_pc, n_ages] V_pc_female;
  matrix[k_pc, n_ages] V_pc_male;
  
  int<lower=0> k_X_pc;
  matrix[N_fam, k_X_pc] X_pc;
  
  int<lower=0> k_Z_pc;
  matrix[N_fam, k_Z_pc] Z_pc;
  
  // age/period spline objects for adjustment factor
  int<lower=0> k_X;
  matrix[N_long, k_X] X;
  
  int<lower=0> k_Z1;
  int<lower=0> k_Z2;
  int<lower=0> k_Z3;
  
  matrix[N_long, k_Z1] Z1;
  matrix[N_long, k_Z2] Z2;
  matrix[N_long, k_Z3] Z3;
  
  real<lower=0> scale_random_spline;
  
 
  // for life expectancy calculations
  int<lower=0> n_groups_unk;
  vector<lower=0>[n_ages] n_agegap;
  array[n_groups_unk] int<lower=1, upper= N_fam> j; // rows for which to produce life expectancy
  
  // prior for overdispersion
  real<lower=0> phi_numer_hmd;
  real<lower=0> phi_numer_fam;
  
  // for the horseshoe prior
  real<lower=0> hs_df;  // local degrees of freedom
  real<lower=0> hs_df_global;  // global degrees of freedom
  real<lower=0> hs_df_slab;  // slab degrees of freedom
  real<lower=0> hs_scale_global;  // global prior scale
  real<lower=0> hs_scale_slab;  // slab prior scale
}

transformed data {
  real<lower=0> scale_mu = 10;
  real<lower=0> scale_alpha = 10;
  real<lower=0> scale_eta = 10;
  
  row_vector<lower=0>[k_pc] scale_beta_pc = [20, 1, 1, 1];

  vector[N_fam] gender_fam_vec = to_vector(gender_fam);
}

parameters {
  // mortality rate parameters -------
  array[k_pc] real<lower=0> sigma_beta_pc;
  
  vector[N_fam*n_ages] eps_raw; // deviations around SVD mortality
  vector<lower=0>[N_fam*n_ages] hs_local;
  real<lower=0> hs_global;
  real<lower=0> hs_slab;
  
  // adjustment parameters --------
  vector[n_groups] intercept_c; 
  real<lower=0> sigma_intercept_c;
  real mean_intercept;
  
  vector[n_countries] age0_effect_raw;
  real mean_age0_effect;
  real<lower=0> sigma_age0_effect;
  
  vector[n_countries] young_effect_raw;
  real mean_young_effect;
  real<lower=0> sigma_young_effect;
  
  // spline parameters for the adjustment factor
  array[2] vector[k_X] beta_X; 
  array[2] vector[k_Z1] alpha_Z1;
  array[2] vector[k_Z2] alpha_Z2;
  array[2] vector[k_Z3] alpha_Z3;
  
  array[2] real<lower=0> sigma_alpha1;
  array[2] real<lower=0> sigma_alpha2;
  array[2] real<lower=0> sigma_alpha3;
  
  array[n_groups] vector[k_Z1] zeta_Z1;
  array[n_groups] vector[k_Z2] zeta_Z2;
  array[n_groups] vector[k_Z3] zeta_Z3;
  
  array[2] real<lower=0> sigma_zeta1;
  array[2] real<lower=0> sigma_zeta2;
  array[2] real<lower=0> sigma_zeta3;
  
  // random slopes for adjustment factor
  array[n_groups] vector[k_X] gamma_X;
  vector<lower=0>[k_X] sigma_gamma_X;
  
  // spline parameters for mortality rate
  array[2, k_pc] real mean_intercept_pc;
  array[n_countries, 2, k_pc] real intercept_pc;
  array[k_pc] real<lower=0> sigma_intercept_pc;
  array[n_countries, 2, k_pc] vector[k_X_pc] beta_X_pc;
  array[n_countries, 2, k_pc] vector[k_Z_pc] alpha_Z_pc;
  array[k_pc] real<lower=0> sigma_alpha_pc;
  
  // overdispersion ---------
  vector<lower=0>[n_ages] phi_inv_hmd;
  vector<lower=0>[n_ages] phi_inv_fam;
  
}

transformed parameters {

  vector[N_fam*n_ages] eps;  
  matrix[N_fam, n_ages] log_mu;
  matrix[N_fam, k_pc] beta_pc_long;
  
  matrix[N_hmd, n_ages] theta_hmd;
  
  vector[N_fam] intercept;
  matrix[N_fam, n_ages] smooth;
  matrix[N_fam, n_ages] log_psi;
  matrix[N_fam, n_ages] theta_fam;
  
  vector[n_countries] age0_effect = sigma_age0_effect*age0_effect_raw + mean_age0_effect;
  vector[n_countries] young_effect = sigma_young_effect*young_effect_raw + mean_young_effect;
  
  // mortality rate parameters ---------
  for (i in 1:N_fam) {
    for (k in 1:k_pc) {
      beta_pc_long[i, k] = scale_beta_pc[k]*(mean_intercept_pc[gender_fam[i], k] + intercept_pc[country_fam[i], gender_fam[i], k])+ 
        X_pc[i] * beta_X_pc[country_fam[i], gender_fam[i], k] * sigma_beta_pc[k] + 
        Z_pc[i] * (alpha_Z_pc[country_fam[i], gender_fam[i], k] * sigma_alpha_pc[k]);
    }
  }
  
  eps = horseshoe(eps_raw, hs_local, hs_global, hs_scale_slab^2 * hs_slab);
  
  log_mu = diag_pre_multiply(gender_fam_vec-2.0, beta_pc_long * V_pc_female) + // if female
    diag_pre_multiply(gender_fam_vec-1.0, beta_pc_long * V_pc_male) + // if male
    to_matrix(eps, N_fam, n_ages); 
  
  // adjustment factor parameters --------
  intercept = mean_intercept + intercept_c[group_fam];

  // form full set 
  {
    vector[N_long] smooth_long;
    int g;

    for (i in 1:n_groups) {
      
      g = gender_group[i];
      
      smooth_long[group_ind[i]] = X[group_ind[i]] * (beta_X[g] + gamma_X[i] .* sigma_gamma_X) +
        Z1[group_ind[i]] * (alpha_Z1[g] * sigma_alpha1[g] + scale_random_spline*zeta_Z1[i]*sigma_zeta1[g]) +
        Z2[group_ind[i]] * (alpha_Z2[g] * sigma_alpha2[g] + scale_random_spline*zeta_Z2[i]*sigma_zeta2[g]) + 
        Z3[group_ind[i]] * (alpha_Z3[g] * sigma_alpha3[g] + scale_random_spline*zeta_Z3[i]*sigma_zeta3[g]);
    }
    
    smooth = to_matrix(smooth_long, N_fam, n_ages, 0); // last arg 0 for filling in row-major order
  }
  
    theta_hmd = log(P_hmd) + log_mu[ii]; // HMD mean
    
    log_psi = rep_matrix(intercept, n_ages) + 
      append_col(age0_effect[country_fam], rep_matrix(0, N_fam, n_ages-1)) + 
      append_col(rep_matrix(young_effect[country_fam], 4), rep_matrix(0, N_fam, n_ages-4)) +
      smooth;
      
    theta_fam = log(P_fam) + log_mu + log_psi ; // adjusted mortality rate
}

model {
  
  mean_intercept ~ normal(0, 3);
  
  intercept_c ~ normal(0, sigma_intercept_c);
  target += normal_lpdf(sigma_intercept_c | 0, 1) - log(0.5);
  
  mean_age0_effect ~ normal(0, 3);
  target += normal_lpdf(sigma_age0_effect | 0, 1) - log(0.5); 
  age0_effect_raw ~ std_normal();
  
  mean_young_effect ~ normal(0, 3);
  target += normal_lpdf(sigma_young_effect | 0, 1) - log(0.5);
  young_effect_raw ~ std_normal();
  
  // regularized horseshoe error parameters
  eps_raw ~ std_normal();
  target += student_t_lpdf(hs_local | hs_df, 0, 1) - rows(hs_local) * log(0.5);
  target += student_t_lpdf(hs_global | hs_df_global, 0, hs_scale_global)
    - 1 * log(0.5);
  target += inv_gamma_lpdf(hs_slab | 0.5 * hs_df_slab, 0.5 * hs_df_slab);
  
  //mortality rate smoothing parameters
  target += normal_lpdf(sigma_alpha_pc | 0, 1) - k_pc*log(0.5);
  target += normal_lpdf(sigma_beta_pc | 0 , 1) - k_pc*log(0.5);
  
  // adjustment factor smoothing parameters
  for (g in 1:2) {
    alpha_Z1[g] ~ std_normal();
    alpha_Z2[g] ~ std_normal();
    alpha_Z3[g] ~ std_normal();
    beta_X[g] ~ normal(0, 3);
  }
  
  target += normal_lpdf(sigma_alpha1 | 0, 1) - 2*log(0.5); 
  target += normal_lpdf(sigma_alpha2 | 0, 1) - 2*log(0.5); 
  target += normal_lpdf(sigma_alpha3 | 0, 1) - 2*log(0.5); 
  
  phi_inv_hmd ~ normal(0, 1);
  phi_inv_fam ~ normal(0, 1);
  
  for (i in 1:n_groups) {
    gamma_X[i] ~ std_normal();
    
    zeta_Z1[i] ~ std_normal();
    zeta_Z2[i] ~ std_normal();
    zeta_Z3[i] ~ std_normal();
  }
  
  target += normal_lpdf(sigma_gamma_X | 0, 1) - k_X*log(0.5); 
  
  sigma_zeta1 ~ std_normal();
  sigma_zeta2 ~ std_normal();
  sigma_zeta3 ~ std_normal();
  
  mean_intercept_pc[1] ~ std_normal();
  mean_intercept_pc[2] ~ std_normal();
  sigma_intercept_pc ~ std_normal(); 
  
  for (i in 1:n_countries) {
    for (k in 1:k_pc) {
      intercept_pc[i, 1, k] ~ normal(0, sigma_intercept_pc[k]);
      beta_X_pc[i, 1, k] ~ std_normal();
      alpha_Z_pc[i, 1, k] ~ std_normal(); 
      
      intercept_pc[i, 2, k] ~ normal(0, sigma_intercept_pc[k]);
      beta_X_pc[i, 2, k] ~ std_normal();
      alpha_Z_pc[i, 2, k] ~ std_normal(); 
    }
  }
  
  // likelihood
  to_array_1d(d_hmd) ~ neg_binomial_2_log(to_array_1d(theta_hmd'), to_array_1d(rep_matrix(phi_numer_hmd ./ phi_inv_hmd, N_hmd)));
  to_array_1d(d_fam[valid_rows,]) ~ neg_binomial_2_log(to_array_1d((theta_fam[valid_rows,])'), to_array_1d(rep_matrix(phi_numer_fam ./ phi_inv_fam, n_valid)));
  
}

generated quantities {
  
  vector[n_groups_unk] e0; 
  vector[n_groups_unk] e10;
  vector[n_groups_unk] e20;
  
  {
    vector[n_ages] e_x_temp;
    for (group in 1:n_groups_unk) {
      e_x_temp = e_x(exp(log_mu[j[group]])', n_agegap);
      
      e0[group] = e_x_temp[1];
      e10[group] = e_x_temp[4];
      e20[group] = e_x_temp[6];
    }
  }

}
