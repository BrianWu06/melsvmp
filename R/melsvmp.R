psi_alpha <- function(uis, mu_alpha, Sigma_alpha){
  lin <- -uis %*% mu_alpha
  quad <- 0.5 * rowSums((uis %*% Sigma_alpha) * uis)
  res <- exp(lin + quad)
  
  as.vector(res)
}

psi_tau <- function(wijs, mu_tau, Sigma_tau){
  lin <- -wijs %*% mu_tau
  quad <- 0.5 * rowSums((wijs %*% Sigma_tau) * wijs)
  res <- exp(lin + quad)
  
  as.vector(res)
}

psi_omega <- function(mu_omegas, Sigma_omegas){
  lin <- -mu_omegas
  quad <- 0.5 * Sigma_omegas
  res <- exp(lin + quad)
  
  as.vector(res)
}

beta_update <- function(X_mat, y_vec, id, uq_id, mu_nu, psi_taus, psi_omegas, beta_prior){
  psi_omega_ij <- psi_omegas[match(id, uq_id)]
  X_weights <- psi_taus * psi_omega_ij
  nuij <- mu_nu[match(id, uq_id)]
  sum_xx <- t(X_mat) %*% sweep(X_mat, 1, X_weights, "*")
  
  weight_resid <- (y_vec - nuij) * X_weights
  sum_yx <- t(X_mat) %*% weight_resid
  
  prior_mat <- diag(1 / beta_prior, ncol(X_mat))
  
  Sigma_beta_new <- solve(sum_xx + prior_mat)
  mu_beta_new <- Sigma_beta_new %*% sum_yx
  
  list(mu_beta_new = mu_beta_new, Sigma_beta_new = Sigma_beta_new)
}

nui_grad_hes <- function(idx, X_mat, y_vec, uq_id, id_list, mu_beta, psi_taus, psi_omegas, psi_alphas, mu_nu){
  idx_char <- as.character(uq_id[idx]) 
  cur_rows <- id_list[[idx_char]]
  
  yi <- y_vec[cur_rows]
  xi <- X_mat[cur_rows, ]
  psi_taui <- psi_taus[cur_rows]
  psi_omegai <- psi_omegas[idx]
  psi_alphai <- psi_alphas[idx]
  mu_nui <- mu_nu[idx]
  
  psi_i <- psi_taui * psi_omegai
  resid_i <- yi - xi %*% mu_beta - mu_nui
  
  lik_grad <- sum(psi_i * resid_i)
  prior_grad <- -psi_alphai * mu_nui
  tol_grad <- lik_grad + prior_grad
  
  lik_hess <- -sum(psi_i)
  prior_hess <- -psi_alphai
  tol_hess <- lik_hess + prior_hess
  
  c(grad = tol_grad, hessian = tol_hess)
}

nu_update <- function(X_mat, y_vec, uq_id, id_list, mu_beta, psi_taus, psi_omegas, psi_alphas, mu_nu){
  num_ids <- length(uq_id)
  nu_grad_hes <- sapply(1:num_ids, nui_grad_hes, X_mat = X_mat, y_vec = y_vec, uq_id = uq_id, id_list = id_list, 
                        mu_beta = mu_beta, psi_taus = psi_taus, psi_omegas = psi_omegas, 
                        psi_alphas = psi_alphas, mu_nu = mu_nu)
  nu_grad <- nu_grad_hes[1,]
  nu_hes <- nu_grad_hes[2,]
  
  Sigma_nu_new <- -1 / nu_hes
  mu_nu_new <- mu_nu + Sigma_nu_new * nu_grad
  
  list(mu_nu_new = mu_nu_new, Sigma_nu_new = Sigma_nu_new)
}

alpha_update <- function(uis, psi_alphas, mu_nu, Sigma_nu, mu_alpha, alpha_prior){
  nu_sq <- mu_nu^2 + Sigma_nu
  
  Sigma_alpha_w <- 0.5 * psi_alphas * nu_sq
  Sigma_alpha_sum <- t(uis) %*% (Sigma_alpha_w * uis)
  prior_mat <- diag(1 / alpha_prior, ncol(uis))
  Sigma_alpha_new <- solve(Sigma_alpha_sum + prior_mat)
  
  mu_alpha_w <- 0.5 * (psi_alphas * nu_sq - 1)
  mu_alpha_sum <- t(uis) %*% mu_alpha_w
  prior_grad <- -(1 / alpha_prior) * mu_alpha
  tol_grad <- mu_alpha_sum + prior_grad
  mu_alpha_new <- mu_alpha + Sigma_alpha_new %*% tol_grad
  
  list(mu_alpha_new = mu_alpha_new, Sigma_alpha_new = Sigma_alpha_new)
}

tau_update <- function(X_mat, y_vec, w_mat, id, uq_id, mu_nu, Sigma_nu, mu_beta, Sigma_beta, psi_taus, 
                           psi_omegas, mu_tau, tau_prior){
  mu_nuis <- mu_nu[match(id, uq_id)]
  Sigma_nuis <- Sigma_nu[match(id, uq_id)]
  psi_omega_ij <- psi_omegas[match(id, uq_id)]
  
  fixed <- y_vec - (X_mat %*% mu_beta) - mu_nuis
  var_beta <- rowSums((X_mat %*% Sigma_beta) * X_mat)
  err_sq <- fixed^2 + var_beta + Sigma_nuis
  hij <- as.vector(err_sq * psi_omega_ij)
  
  hess_w <- 0.5 * psi_taus * hij
  sum_ww <- t(w_mat) %*% sweep(w_mat, 1, hess_w, "*")
  prior_mat <- diag(1/tau_prior, ncol(w_mat))
  Sigma_tau_new <- solve(sum_ww + prior_mat)
  
  grad_w <- 0.5 * (psi_taus * hij - 1)
  lik_grad <- t(w_mat) %*% grad_w
  prior_grad <- -(1/tau_prior) * mu_tau
  tol_grad <- lik_grad + prior_grad
  mu_tau_new <- mu_tau + Sigma_tau_new %*% tol_grad
  
  list(mu_tau_new = mu_tau_new, Sigma_tau_new = Sigma_tau_new)
}

omega_grad_hess <- function(idx, X_mat, y_vec, w_mat, uq_id, id_list, mu_beta, Sigma_beta, mu_nu, Sigma_nu, 
                                psi_taus, mu_inv_sigma, psi_omegas, mu_omegas){
  idx_char <- as.character(uq_id[idx])
  cur_rows <- id_list[[idx_char]]
  
  yi <- y_vec[cur_rows]
  xi <- X_mat[cur_rows, ]
  mu_nui <- mu_nu[idx]
  Sigma_nui <- Sigma_nu[idx]
  mu_omegai <- mu_omegas[idx]
  psi_omegai <- psi_omegas[idx]
  psi_taui <- psi_taus[cur_rows]
  
  fixed_i <- yi - xi %*% mu_beta - mu_nui
  var_beta_i <- rowSums((xi %*% Sigma_beta) * xi)
  err_sq_i <- fixed_i^2 + var_beta_i + Sigma_nui
  hij <- err_sq_i * psi_taui
  
  lik_grad <- 0.5 * psi_omegai * sum(hij) - 0.5 * length(yi)
  prior_grad <- -mu_inv_sigma * mu_omegai
  tol_grad <- lik_grad + prior_grad
  
  lik_hess <- -0.5 * psi_omegai * sum(hij)
  prior_hess <- -mu_inv_sigma
  tol_hess <- lik_hess + prior_hess
  
  c(grad = tol_grad, hessian = tol_hess)
}

omega_update <- function(X_mat, y_vec, w_mat, uq_id, id_list, mu_beta, Sigma_beta, mu_nu, Sigma_nu, mu_inv_sigma, 
                             psi_taus, psi_omegas, mu_omegas){
  num_ids <- length(uq_id)
  omega_grad_hess <- sapply(1:num_ids, omega_grad_hess, X_mat = X_mat, y_vec = y_vec, w_mat = w_mat, 
                            uq_id = uq_id, id_list = id_list, mu_beta = mu_beta, Sigma_beta = Sigma_beta, 
                            mu_nu = mu_nu, Sigma_nu = Sigma_nu, psi_taus = psi_taus, mu_inv_sigma = mu_inv_sigma, 
                            psi_omegas = psi_omegas, mu_omegas = mu_omegas)
  
  omega_grad <- omega_grad_hess[1, ]
  omega_hess <- omega_grad_hess[2, ]
  
  Sigma_omega_new <- -1 / omega_hess
  mu_omega_new <- mu_omegas + Sigma_omega_new * omega_grad
  
  list(mu_omega_new = as.vector(mu_omega_new), Sigma_omega_new = as.vector(Sigma_omega_new))
}

sigma_inv_update <- function(num_ids, mu_omega, Sigma_omega, mu_inv_a){
  sigma_shape <- (num_ids+1)/2
  sum_omega_sq <- sum(Sigma_omega + mu_omega^2)
  rate_sigma <- mu_inv_a + 0.5 * sum_omega_sq
  
  mu_inv_sigma_new <- sigma_shape / rate_sigma
  mu_log_sigma <- log(rate_sigma) - digamma(sigma_shape)
  
  list(mu_inv_sigma_new = mu_inv_sigma_new, mu_log_sigma = mu_log_sigma)
}

a_inv_update <- function(mu_inv_sigma, A_omega){
  rate_a <- (1 / A_omega^2) + mu_inv_sigma
  mu_inv_a_new <- 1 / rate_a
  mu_log_a <- log(rate_a) - digamma(1)
  
  list(mu_inv_a_new = mu_inv_a_new, mu_log_a = mu_log_a)
}

elbo_cal <- function(X_mat, y_vec, u_mat, w_mat, mu_beta, Sigma_beta, mu_nu, Sigma_nu, 
                         mu_alpha, Sigma_alpha, mu_tau, Sigma_tau, mu_omega, Sigma_omega, 
                         mu_inv_sigma, mu_log_sigma, mu_inv_a, mu_log_a, 
                         psi_ij_tau, psi_i_alpha, psi_i_omega,
                         beta_prior, alpha_prior, tau_prior, A_prior,
                         id, uq_id){
  num_ids <- length(uq_id)
  
  mu_omega_expanded <- mu_omega[match(id, uq_id)]
  psi_omega_expanded <- psi_i_omega[match(id, uq_id)]
  mu_nuis <- mu_nu[match(id, uq_id)]
  Sigma_nuis <- Sigma_nu[match(id, uq_id)]
  
  fixed <- y_vec - (X_mat %*% mu_beta) - mu_nuis
  var_beta <- rowSums((X_mat %*% Sigma_beta) * X_mat)
  err_sq <- fixed^2 + var_beta + Sigma_nuis
  hij <- as.vector(err_sq * psi_ij_tau)
  
  log_var_mean_expanded <- as.vector(w_mat %*% mu_tau) + mu_omega_expanded
  log_p_y <- -0.5 * sum(log_var_mean_expanded) - 0.5 * sum(psi_omega_expanded * hij)
  
  log_p_beta <- -0.5 * (1 / beta_prior) * (sum(mu_beta^2) + sum(diag(Sigma_beta)))
  log_p_tau <- -0.5 * (1 / tau_prior) * (sum(mu_tau^2) + sum(diag(Sigma_tau)))
  log_p_alpha <- -0.5 * (1 / alpha_prior) * (sum(mu_alpha^2) + sum(diag(Sigma_alpha)))
  
  log_p_nu <- -0.5 * sum(u_mat %*% mu_alpha + psi_i_alpha * (mu_nu^2 + Sigma_nu))
  log_p_omega <- -0.5 * mu_inv_sigma * sum(mu_omega^2 + Sigma_omega) 
  
  log_p_sigma <- -0.5 * mu_log_a - 1.5 * mu_log_sigma - mu_inv_a * mu_inv_sigma
  log_p_a <- -1.5 * mu_log_a - mu_inv_a * (1 / A_prior^2)
  
  log_q_beta <- 0.5 * determinant(Sigma_beta, logarithm = TRUE)$modulus
  log_q_tau <- 0.5 * determinant(Sigma_tau, logarithm = TRUE)$modulus
  log_q_alpha <- 0.5 * determinant(Sigma_alpha, logarithm = TRUE)$modulus
  log_q_nu <- 0.5 * sum(log(Sigma_nu))
  log_q_omega <- 0.5 * sum(log(Sigma_omega)) 
  
  rate_sigma <- mu_inv_a + 0.5 * sum(Sigma_omega + mu_omega^2) 
  log_q_sigma <- -lgamma((num_ids+1)/2) + ((num_ids+1)/2)*log(rate_sigma)
  
  rate_a <- (1 / A_prior^2) + mu_inv_sigma
  log_q_a <- -lgamma(1) + log(rate_a)
  
  p_terms <- log_p_y + log_p_beta + log_p_tau + log_p_alpha + log_p_nu + log_p_omega + log_p_sigma + log_p_a
  q_terms <- log_q_beta + log_q_alpha + log_q_tau + log_q_nu + log_q_omega + log_q_sigma + log_q_a
  
  elbo <- p_terms - q_terms
  elbo
}

san_est <- function(X_matrix, y_vector, u_matrix, w_matrix, ids, unique_ids, id_indices,
                    mu_beta_q, Sigma_beta_q, mu_nu_i_q, Sigma_nu_i_q,
                    mu_alpha_q, Sigma_alpha_q, mu_tau_q, Sigma_tau_q,
                    mu_omega_q, Sigma_omega_q, mu_inv_sigma_omega_q,
                    beta_prior, alpha_prior, tau_prior) {
  
  # Get parameter dimensions
  p_beta <- ncol(X_matrix)
  p_alpha <- ncol(u_matrix)
  p_tau <- ncol(w_matrix)
  p_total <- p_beta + p_alpha + p_tau
  num_ids <- length(unique_ids)
  
  # Initialize the Bread (A) and Meat (B) matrices
  A_hat <- matrix(0, p_total, p_total)
  B_hat <- matrix(0, p_total, p_total)
  
  # Pre-calculate all converged psi expectations
  psi_i_alpha_all <- psi_alpha(u_matrix, mu_alpha_q, Sigma_alpha_q)
  psi_ij_tau_all <- psi_tau(w_matrix, mu_tau_q, Sigma_tau_q)
  psi_i_omega_all <- psi_omega(mu_omega_q, Sigma_omega_q)
  
  # Loop over each subject (i) to build A_i and B_i
  for (i in 1:num_ids) {
    
    # --- 1. Get subject-specific data and params ---
    idx_char <- as.character(unique_ids[i])
    cur_rows <- id_indices[[idx_char]]
    
    # Skip subjects with no data (can happen in bootstrap)
    if (length(cur_rows) == 0) {
      next
    }
    
    xi <- X_matrix[cur_rows, , drop = FALSE]
    yi <- y_vector[cur_rows]
    ui <- u_matrix[i, , drop = FALSE]
    wi <- w_matrix[cur_rows, , drop = FALSE]
    
    mu_nui <- mu_nu_i_q[i]
    Sigma_nui <- Sigma_nu_i_q[i]
    mu_omegai <- mu_omega_q[i]
    Sigma_omegai <- Sigma_omega_q[i]
    
    # --- 2. Calculate subject-specific expectations ---
    psi_alphai <- psi_i_alpha_all[i]
    psi_taui <- psi_ij_tau_all[cur_rows]
    psi_omegai <- psi_i_omega_all[i]
    
    Lambda_i <- psi_taui * psi_omegai
    
    # --- THIS IS THE FIX ---
    # Force resid_i to be a vector, not an n_i x 1 matrix
    resid_i <- as.vector(yi - xi %*% mu_beta_q - mu_nui)
    # --- END FIX ---
    
    H_i_terms <- resid_i^2 + rowSums((xi %*% Sigma_beta_q) * xi) + Sigma_nui
    nu_sq_i <- Sigma_nui + mu_nui^2
    E_inv_sigma_omega <- mu_inv_sigma_omega_q
    
    # --- 3. Calculate G_i (Gradient) for the "Meat" ---
    # (vector[n_i] * vector[n_i]) * matrix[n_i, p_beta] -> matrix[n_i, p_beta]
    G_beta_i <- colSums((Lambda_i * resid_i) * xi, na.rm = TRUE)
    G_alpha_i <- 0.5 * (psi_alphai * nu_sq_i - 1) * as.vector(ui)
    G_tau_i <- 0.5 * colSums((Lambda_i * H_i_terms - 1) * wi, na.rm = TRUE)
    
    G_i <- c(G_beta_i, G_alpha_i, G_tau_i)
    B_hat <- B_hat + G_i %*% t(G_i)
    
    # --- 4. Calculate H_i (Adjusted Hessian) for the "Bread" ---
    
    # H_local: D^2 v_i w.r.t. psi_i (local params) [2x2]
    H_vv <- -sum(Lambda_i, na.rm = TRUE) - psi_alphai
    H_oo <- -0.5 * sum(Lambda_i * H_i_terms, na.rm = TRUE) - E_inv_sigma_omega
    H_vo <- -sum(Lambda_i * resid_i, na.rm = TRUE)
    H_local <- matrix(c(H_vv, H_vo, H_vo, H_oo), 2, 2)
    
    # Check for singularity
    if (abs(det(H_local)) < 1e-10) {
      next # Skip this subject if H_local is singular
    }
    H_local_inv <- solve(H_local)
    
    # H_global: D^2 v_i w.r.t. theta (global params) [p_total x p_total]
    H_bb <- -t(xi) %*% (Lambda_i * xi)
    H_aa <- -0.5 * psi_alphai * nu_sq_i * (t(ui) %*% ui)
    H_tt <- -0.5 * t(wi) %*% (Lambda_i * H_i_terms * wi)
    H_bt <- -t(xi) %*% ((Lambda_i * resid_i) * wi) # (vector * matrix)
    H_tb <- t(H_bt)
    
    # Zero blocks
    H_ba <- matrix(0, p_beta, p_alpha)
    H_ab <- matrix(0, p_alpha, p_beta)
    H_at <- matrix(0, p_alpha, p_tau)
    H_ta <- matrix(0, p_tau, p_alpha)
    
    H_global <- rbind(cbind(H_bb, H_ba, H_bt),
                      cbind(H_ab, H_aa, H_at),
                      cbind(H_tb, H_ta, H_tt))
    
    # H_cross: D^2 v_i w.r.t. theta and psi_i [p_total x 2]
    H_bv <- -colSums(Lambda_i * xi, na.rm = TRUE)
    H_bo <- -colSums((Lambda_i * resid_i) * xi, na.rm = TRUE)
    
    H_av <- as.vector(psi_alphai * mu_nui * ui)
    H_ao <- rep(0, p_alpha)
    
    H_tv <- -colSums((Lambda_i * resid_i) * wi, na.rm = TRUE)
    H_to <- -0.5 * colSums((Lambda_i * H_i_terms) * wi, na.rm = TRUE)
    
    H_cross <- rbind(cbind(H_bv, H_bo),
                     cbind(H_av, H_ao),
                     cbind(H_tv, H_to))
    
    # H_i = H_global - H_cross %*% solve(H_local) %*% t(H_cross)
    H_i <- H_global - H_cross %*% H_local_inv %*% t(H_cross)
    
    A_hat <- A_hat + H_i
  }
  
  # --- 5. Add Priors to the Bread ---
  H_prior_beta <- diag(-1 / beta_prior, p_beta)
  H_prior_alpha <- diag(-1 / alpha_prior, p_alpha)
  H_prior_tau <- diag(-1 / tau_prior, p_tau)
  
  # Use Matrix::bdiag to create a block-diagonal matrix
  H_prior <- as.matrix(Matrix::bdiag(H_prior_beta, H_prior_alpha, H_prior_tau))
  
  A_hat_final <- A_hat + H_prior
  
  # --- 6. Assemble the Sandwich ---
  # Use solve() with error handling
  A_hat_inv <- tryCatch(
    solve(A_hat_final),
    error = function(e) {
      warning("Hessian matrix (Bread) is singular. Sandwich estimator failed. \n", e)
      return(matrix(NA, p_total, p_total))
    }
  )
  
  if (any(is.na(A_hat_inv))) {
    V_hat <- matrix(NA, p_total, p_total)
  } else {
    V_hat <- A_hat_inv %*% B_hat %*% A_hat_inv
  }
  
  return(V_hat)
}

mels_vmp_fitter <- function(X_matrix, y_vector, w_matrix, u_matrix, ids, max_iter = 1000, tol = 1e-6, 
                                  beta_prior = 1e+4, alpha_prior = 1e+4, tau_prior = 1e+4, A_prior = 1e+4, 
                                  verbose = TRUE){
  unique_ids <- unique(ids)
  num_ids <- length(unique_ids)
  id_indices <- split(1:length(ids), ids)
  
  lmer_df <- data.frame(y = as.vector(y_vector), id = ids)
  covariates_df <- as.data.frame(X_matrix[, -1, drop = FALSE])
  lmer_df <- cbind(lmer_df, covariates_df)
  
  covariate_names <- colnames(covariates_df)
  lmer_formula <- as.formula(paste("y ~", paste(covariate_names, collapse = " + "), "+ (1 | id)"))
  
  lme_fit <- lme4::lmer(lmer_formula, data = lmer_df)
  beta_init <- lme4::fixef(lme_fit)
  
  ranef_df <- lme4::ranef(lme_fit)$id
  nu_init <- as.vector(ranef_df[match(as.character(unique_ids), rownames(ranef_df)), "(Intercept)"])
  
  resid_init <- residuals(lme_fit)
  
  p_beta <- ncol(X_matrix); mu_beta_q <- beta_init; Sigma_beta_q <- as.matrix(vcov(lme_fit))
  mu_nu_i_q <- nu_init; Sigma_nu_i_q <- rep(0.1, num_ids)
  p_alpha <- ncol(u_matrix); mu_alpha_q <- rep(0, p_alpha); Sigma_alpha_q <- diag(1, p_alpha)
  p_tau <- ncol(w_matrix); lm_tau <- lm(log(resid_init^2 + 1e-4) ~ w_matrix - 1)
  mu_tau_q <- coef(lm_tau); Sigma_tau_q <- as.matrix(vcov(lm_tau))
  mu_omega_q <- tapply(residuals(lm_tau), ids, mean)[as.character(unique_ids)]
  Sigma_omega_q <- rep(0.1, num_ids)
  mu_inv_sigma_omega_q <- 1; mu_log_sigma_omega_q <- 0
  mu_inv_a_omega_q <- 1; mu_log_a_omega_q <- 0
  
  psi_i_alpha <- psi_alpha(u_matrix, mu_alpha_q, Sigma_alpha_q)
  psi_ij_tau <- psi_tau(w_matrix, mu_tau_q, Sigma_tau_q)
  psi_i_omega <- psi_omega(mu_omega_q, Sigma_omega_q)
  
  elbo_history <- c()
  for (i in 1:max_iter) {
    beta_new <- beta_update(X_matrix, y_vector, ids, unique_ids, mu_nu_i_q, 
                                psi_ij_tau, psi_i_omega, beta_prior)
    mu_beta_q <- beta_new$mu_beta_new; Sigma_beta_q <- beta_new$Sigma_beta_new
    
    nu_new <- nu_update(X_matrix, y_vector, unique_ids, id_indices, mu_beta_q, 
                            psi_ij_tau, psi_i_omega, psi_i_alpha, mu_nu_i_q)
    mu_nu_i_q <- nu_new$mu_nu_new; Sigma_nu_i_q <- nu_new$Sigma_nu_new
    
    alpha_new <- alpha_update(u_matrix, psi_i_alpha, mu_nu_i_q, Sigma_nu_i_q, mu_alpha_q, alpha_prior)
    mu_alpha_q <- alpha_new$mu_alpha_new; Sigma_alpha_q <- alpha_new$Sigma_alpha_new
    psi_i_alpha <- psi_alpha(u_matrix, mu_alpha_q, Sigma_alpha_q)
    
    tau_new <- tau_update(X_matrix, y_vector, w_matrix, ids, unique_ids, mu_nu_i_q, Sigma_nu_i_q,
                              mu_beta_q, Sigma_beta_q, psi_ij_tau, psi_i_omega, mu_tau_q, tau_prior)
    mu_tau_q <- tau_new$mu_tau_new; Sigma_tau_q <- tau_new$Sigma_tau_new
    psi_ij_tau <- psi_tau(w_matrix, mu_tau_q, Sigma_tau_q)
    
    omega_new <- omega_update(X_matrix, y_vector, w_matrix, unique_ids, id_indices, mu_beta_q, 
                                  Sigma_beta_q, mu_nu_i_q, Sigma_nu_i_q,
                                  mu_inv_sigma_omega_q, psi_ij_tau, psi_i_omega, mu_omega_q)
    mu_omega_q <- omega_new$mu_omega_new; Sigma_omega_q <- omega_new$Sigma_omega_new
    psi_i_omega <- psi_omega(mu_omega_q, Sigma_omega_q)
    
    sigma_omega_new <- sigma_inv_update(num_ids, mu_omega_q, Sigma_omega_q, mu_inv_a_omega_q)
    mu_inv_sigma_omega_q <- sigma_omega_new$mu_inv_sigma_new; mu_log_sigma_omega_q <- sigma_omega_new$mu_log_sigma
    
    a_omega_new <- a_inv_update(mu_inv_sigma_omega_q, A_prior)
    mu_inv_a_omega_q <- a_omega_new$mu_inv_a_new; mu_log_a_omega_q <- a_omega_new$mu_log_a
    
    elbo <- elbo_cal(X_matrix, y_vector, u_matrix, w_matrix, mu_beta_q, Sigma_beta_q, mu_nu_i_q, Sigma_nu_i_q,
                         mu_alpha_q, Sigma_alpha_q, mu_tau_q, Sigma_tau_q, mu_omega_q, Sigma_omega_q,
                         mu_inv_sigma_omega_q, mu_log_sigma_omega_q, mu_inv_a_omega_q, mu_log_a_omega_q,
                         psi_ij_tau, psi_i_alpha, psi_i_omega,
                         beta_prior, alpha_prior, tau_prior, A_prior,
                         ids, unique_ids)
    elbo_history <- c(elbo_history, elbo)
    
    if ((i > 1) && (abs(elbo_history[i] - elbo_history[i-1]) < tol * abs(elbo_history[i]))) {
      if (verbose) { cat("CAVI converges at iteration ", i, "\n") }
      break
    }
    if (i == max_iter && verbose) {
      cat("Algorithm reached max iterations (", i, ") without converging.\n")
    }
  }
  
  # Call the new function with all converged parameters
  sand_cov_mat <- san_est(X_matrix, y_vector, u_matrix, w_matrix, ids, unique_ids, id_indices,
                          mu_beta_q, Sigma_beta_q, mu_nu_i_q, Sigma_nu_i_q,
                          mu_alpha_q, Sigma_alpha_q, mu_tau_q, Sigma_tau_q,
                          mu_omega_q, Sigma_omega_q, mu_inv_sigma_omega_q,
                          beta_prior, alpha_prior, tau_prior) # <--- ADDED
  
  # Partition the full sandwich matrix into blocks
  p_total <- p_beta + p_alpha + p_tau
  beta_sand_cov <- sand_cov_mat[1:p_beta, 1:p_beta, drop = FALSE] # <--- ADDED
  alpha_sand_cov <- sand_cov_mat[(p_beta + 1):(p_beta + p_alpha), (p_beta + 1):(p_beta + p_alpha), drop = FALSE] # <--- ADDED
  tau_sand_cov <- sand_cov_mat[(p_beta + p_alpha + 1):p_total, (p_beta + p_alpha + 1):p_total, drop = FALSE] # <--- ADDED
  
  sigma_omega_rate <- mu_inv_a_omega_q + 0.5 * sum(Sigma_omega_q + mu_omega_q^2)
  sigma_omega_A <- (num_ids + 1) / 2
  sigma_omega_mean <- sigma_omega_rate / ((num_ids + 1) / 2 - 1)
  if (sigma_omega_A > 2) {
    var_sigma_omega_sq <- sigma_omega_rate^2 / ((sigma_omega_A - 1)^2 * (sigma_omega_A - 2))
    # Apply Delta Method: Var(sqrt(X)) approx = Var(X) / (4 * E[X])
    var_sigma_omega_approx <- var_sigma_omega_sq / (4 * sigma_omega_mean)
    se_sigma_omega_approx <- sqrt(var_sigma_omega_approx)
  } else {
    se_sigma_omega_approx <- NA # Variance is undefined
  }
  
  rownames(beta_new$mu_beta_new) <- colnames(X_matrix)
  rownames(alpha_new$mu_alpha_new) <- colnames(u_matrix)
  names(tau_new$mu_tau_new) <- colnames(w_matrix)
  
  results <- list(
    beta = list(params = mu_beta_q, cov_mat = beta_sand_cov, vmp_cov_mat = Sigma_beta_q), # <--- CHANGED
    alpha = list(params = mu_alpha_q, cov_mat = alpha_sand_cov, vmp_cov_mat = Sigma_alpha_q), # <--- CHANGED
    tau = list(params = mu_tau_q, cov_mat = tau_sand_cov, vmp_cov_mat = Sigma_tau_q), # <--- CHANGED
    omega = list(std_dev = sqrt(sigma_omega_mean), approx_se = se_sigma_omega_approx), 
    r = NULL,
    elbo_history = elbo_history, 
    iterations = i
  )
  
  return(results)
}


#' Using a Variational Message Passing (VMP) algorithm to fit a Mixed-Effects Location Scale Model.
#' 
#' @details
#' This function fits a mixed-effects model with structured variances.
#' The model is defined as:
#' 
#' \deqn{y_{ij} = x_{ij}^\top \beta + \nu_i + \varepsilon_{ij}}
#' \deqn{\nu_i \sim \mathcal{N}(0, \exp{(U_i^\top \alpha)})}
#' \deqn{\varepsilon_{ij} \sim \mathcal{N}(0, \exp{(W_{ij}^\top \tau + \omega_i)})}
#' \deqn{\omega_i \sim \mathcal{N}(0, \sigma_\omega^2)}
#' 
#' Where \eqn{\beta} and \eqn{\alpha} are fixed effects for the mean and
#' between-subject variance, and \eqn{\tau} are fixed effects for the 
#' within-subject variance.
#'
#' @param y The name of the response variable in 'data'.
#' @param beta_formula A formula for the mean model (fixed effects).
#' @param alpha_formula A formula for the between-subject variance.
#' @param tau_formula A formula for the within-subject variance.
#' @param id The name of the subject ID variable in 'data'.
#' @param data A data.frame containing all variables.
#' @param ... Additional arguments passed to the fitter function (e.g., 
#'   `max_iter`, `tol`).
#'
#' @param y The name of the response variable in 'data'.
#' @param beta_formula A formula for the mean model (fixed effects).
#' @param alpha_formula A formula for the between-subject variance.
#' @param tau_formula A formula for the within-subject variance.
#' @param id The name of the subject ID variable in 'data'.
#' @param data A data.frame containing all variables.
#' @param ... Additional arguments passed to the fitter function.
#'
#' @return An object of class 'mels_vmp'.
#' 
#' @seealso \code{\link{summary.mels_vmp}}, \code{\link{bootstrap_mels_vmp}}
#' 
#' @export 
#' @example /examples/mels_vmp_example.R
#' @importFrom lme4 lmer fixef ranef
#' @importFrom stats as.formula coef lm model.matrix na.omit pnorm residuals vcov
#'
mels_vmp <- function(y, beta_formula, alpha_formula, tau_formula, id, data, ...){
  st_time <- Sys.time()
  all_vars <- c(y, id, all.vars(beta_formula), all.vars(alpha_formula), all.vars(tau_formula))
  all_vars <- unique(all_vars)
  all_vars <- all_vars[all_vars %in% colnames(data)]
  clean_data <- na.omit(data[, all_vars])
  
  y_vector <- as.matrix(clean_data[[y]])
  ids_vector <- clean_data[[id]]
  X_matrix <- model.matrix(beta_formula, data = clean_data)
  w_matrix <- model.matrix(tau_formula, data = clean_data)
  
  subject_data <- clean_data[!duplicated(clean_data[[id]]), ]
  u_matrix <- model.matrix(alpha_formula, subject_data)
  
  results <- mels_vmp_fitter(X_matrix = X_matrix, y_vector = y_vector, w_matrix = w_matrix, 
                                   u_matrix = u_matrix, ids = ids_vector, ...) 
  results$call <- match.call()
  results$data <- data
  
  ed_time <- Sys.time()
  runtime <- format(round(ed_time - st_time, 2))
  results$runtime <- runtime
  
  class(results) <- "mels_vmp"
  return(results)
}

#' Summarize a mels_vmp Model Fit
#'
#' Prints a formatted summary of a fitted Mixed-Effects Location Scale model,
#' including parameter estimates, robust sandwich standard errors, z-values,
#' and p-values.
#'
#' @details
#' The summary is organized into sections:
#' \itemize{
#'   \item \strong{Mean Model Parameters (beta):} Estimates for the fixed effects in the mean model.
#'   \item \strong{Between-Subject Variance Parameters (alpha):} Estimates for the fixed effects in the between-subject variance model.
#'   \item \strong{Within-Subject Variance Parameters (tau):} Estimates for the fixed effects in the within-subject variance model.
#'   \item \strong{Random Effect Standard Deviation (omega):} Estimate for the random effect standard deviation.
#'   \item \strong{Convergence Details:} Number of iterations and total runtime.
#' }
#' All standard errors and confidence intervals are based on the robust 
#' sandwich estimator.
#'
#' @param object An object of class `mels_vmp`.
#' @param ... Additional arguments (currently unused).
#'
#' @return
#' Invisibly returns the original `mels_vmp` object.
#'
#' @seealso \code{\link{mels_vmp}}
#' @exportS3Method summary mels_vmp
#'
summary.mels_vmp <- function(object, ...){
  cat("## VMP for MELS (with Robust Sandwich Errors) ##\n")
  cat("--------------------------------------------------------\n")
  
  cat("--- Mean Model Parameters (beta) ---\n")
  beta_se <- sqrt(diag(object$beta$cov_mat))
  beta_z <- object$beta$params / beta_se
  beta_p <- 2 * pnorm(-abs(beta_z))
  beta_df <- data.frame(
    Estimate = object$beta$params,
    'Robust SE' = beta_se,
    'CI.Lower' = object$beta$params - 1.96 * beta_se,  # <--- ADDED
    'CI.Upper' = object$beta$params + 1.96 * beta_se,  # <--- ADDED
    'z value' = beta_z,
    'p-value' = beta_p,
    check.names = FALSE
  )
  print(round(beta_df, 4))
  cat("\n")
  
  cat("--- Between-Subject Variance Parameters (alpha) ---\n")
  alpha_se <- sqrt(diag(object$alpha$cov_mat))
  alpha_z <- object$alpha$params / alpha_se
  alpha_p <- 2 * pnorm(-abs(alpha_z))
  alpha_df <- data.frame(
    Estimate = object$alpha$params,
    'Robust SE' = alpha_se,
    'CI.Lower' = object$alpha$params - 1.96 * alpha_se, # <--- ADDED
    'CI.Upper' = object$alpha$params + 1.96 * alpha_se, # <--- ADDED
    'z value' = alpha_z,
    'p-value' = alpha_p,
    check.names = FALSE
  )
  print(round(alpha_df, 4))
  cat("\n")
  
  cat("--- Within-Subject Variance Parameters (tau) ---\n")
  tau_se <- sqrt(diag(object$tau$cov_mat))
  tau_z <- object$tau$params / tau_se
  tau_p <- 2 * pnorm(-abs(tau_z))
  tau_df <- data.frame(
    Estimate = object$tau$params,
    'Robust SE' = tau_se,
    'CI.Lower' = object$tau$params - 1.96 * tau_se,  # <--- ADDED
    'CI.Upper' = object$tau$params + 1.96 * tau_se,  # <--- ADDED
    'z value' = tau_z,
    'p-value' = tau_p,
    check.names = FALSE
  )
  print(round(tau_df, 4))
  
  cat("\n")
  cat("--- Random Effect Standard Deviation (omega) ---\n")
  omega_est <- object$omega$std_dev
  omega_se <- object$omega$approx_se
  
  # Calculate CI, but ensure it's not negative (since std dev > 0)
  omega_ci_lower <- pmax(0, omega_est - 1.96 * omega_se) 
  omega_ci_upper <- omega_est + 1.96 * omega_se
  
  omega_df <- data.frame(
    Estimate = omega_est,
    'Approx. SE' = omega_se, # Labelled as Approx. SE (Delta Method)
    'CI.Lower' = omega_ci_lower,
    'CI.Upper' = omega_ci_upper,
    check.names = FALSE
  )
  rownames(omega_df) <- "omega_std_dev"
  print(round(omega_df, 4))
  
  cat("-------------------------------------------------------\n")
  cat("Convergence Details:\n")
  cat(paste0("  Algorithm converged in ", object$iterations, " iterations.\n"))
  cat(paste0("  Total Runtime: ", object$runtime, " \n"))
  
  invisible(object)
}

#' Non-parametric Bootstrap for mels_vmp Models
#'
#' Performs a non-parametric bootstrap (resampling subjects with replacement)
#' to estimate standard errors and confidence intervals for a mels_vmp model.
#'
#' @details
#' This function implements a non-parametric bootstrap by resampling subjects
#' (as defined by the `id` variable) with replacement. For each of the `B`
#' replicates, it refits the `mels_vmp` model on the resampled data.
#' It can be run in parallel by setting `parallel = TRUE`.
#'
#' @param model_object An object of class 'mels_vmp'.
#' @param B Number of bootstrap replicates.
#' @param parallel Logical. If `TRUE` (the default), run in parallel.
#' @param cores Number of cores to use for parallel execution.
#' @param ... Additional arguments (currently unused).
#'
#' @return
#' An object of class `mels_vmp_bootstrap`. This is a list containing:
#' \itemize{
#'   \item `beta`: A list with original estimates, bootstrap SE, and percentile CIs.
#'   \item `alpha`: A list with original estimates, bootstrap SE, and percentile CIs.
#'   \item `tau`: A list with original estimates, bootstrap SE, and percentile CIs.
#'   \item `omega`: A list with original estimates, bootstrap SE, and percentile CIs.
#'   \item `n_reps`: The number of successful replicates.
#'   \item `runtime`: The total runtime.
#' }
#'
#' @seealso \code{\link{mels_vmp}}, \code{\link{summary.mels_vmp_bootstrap}}
#'
#'
#' @export
#' @importFrom dplyr left_join
#' @importFrom foreach foreach %dopar%
#' @importFrom doParallel registerDoParallel
#' @importFrom parallel makeCluster stopCluster
#' @importFrom stats sd quantile
#'
bootstrap_mels_vmp <- function(model_object, B = 1000, 
                               parallel = TRUE, cores = 5, ...) {
  
  st_time <- Sys.time()
  
  original_call <- model_object$call
  original_data <- model_object$data
  id_col_name <- as.character(original_call$id)
  
  unique_ids <- unique(original_data[[id_col_name]])
  
  beta_results <- list()
  alpha_results <- list()
  tau_results <- list()
  omega_results <- list()
  
  if (parallel) {
    cl <- parallel::makeCluster(cores)
    doParallel::registerDoParallel(cl)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    
    cat(paste0("Starting PARALLEL bootstrap with ", B, " replicates on ", cores, " cores...\n"))
    
    functions_to_export <- c(
      "mels_vmp", "mels_vmp_fitter", "mels_vmp_fitter", "psi_alpha", 
      "psi_tau", "psi_omega", "beta_update", "nui_grad_hes", 
      "nu_update", "alpha_update", "tau_update", "omega_grad_hess", 
      "omega_update", "sigma_inv_update", "a_inv_update", "elbo_cal", "san_est"
    )
    
    results_list <- foreach(
      b = 1:B, .packages = c("dplyr", "lme4"), .export = functions_to_export, .errorhandling = "pass"
    ) %dopar% {
      resampled_ids <- sample(unique_ids, size = length(unique_ids), replace = TRUE)
      resampled_df <- data.frame(new_id = 1:length(unique_ids))
      resampled_df[[id_col_name]] <- resampled_ids
      bootstrap_data <- dplyr::left_join(resampled_df, original_data, by = id_col_name, relationship = "many-to-many")
      
      boot_call <- original_call
      boot_call$data <- bootstrap_data
      boot_call$id <- "new_id"
      boot_call$verbose <- FALSE
      boot_fit <- eval(boot_call)
      
      if (!is.null(boot_fit) && !inherits(boot_fit, "error")) {
        list(beta = as.vector(boot_fit$beta$params), alpha = as.vector(boot_fit$alpha$params),
             tau = as.vector(boot_fit$tau$params), omega = boot_fit$omega$std_dev)
      } else { NULL }
    }
    
    results_list <- results_list[!sapply(results_list, is.null)]
    beta_results <- lapply(results_list, `[[`, "beta")
    alpha_results <- lapply(results_list, `[[`, "alpha")
    tau_results <- lapply(results_list, `[[`, "tau")
    omega_results <- lapply(results_list, `[[`, "omega")
    
  } else {
    cat(paste0("Starting SEQUENTIAL bootstrap with ", B, " replicates...\n"))
    
    for (b in 1:B) {
      if (b %% 100 == 0) { cat(paste0("  Running replicate: ", b, " of ", B, "\n")) }
      
      resampled_ids <- sample(unique_ids, size = length(unique_ids), replace = TRUE)
      resampled_df <- data.frame(new_id = 1:length(unique_ids))
      resampled_df[[id_col_name]] <- resampled_ids
      bootstrap_data <- dplyr::left_join(resampled_df, original_data, by = id_col_name, relationship = "many-to-many")
      
      boot_call <- original_call
      boot_call$data <- bootstrap_data
      boot_call$id <- "new_id"
      boot_call$verbose <- FALSE
      boot_fit <- eval(boot_call)
      
      if (!is.null(boot_fit) && !inherits(boot_fit, "error")) {
        beta_results[[b]] <- as.vector(boot_fit$beta$params)
        alpha_results[[b]] <- as.vector(boot_fit$alpha$params)
        tau_results[[b]] <- as.vector(boot_fit$tau$params)
        omega_results[[b]] <- boot_fit$omega$std_dev
      }
    }
  }
  
  summarize_boots <- function(results_list, original_estimates) {
    results_list <- results_list[!sapply(results_list, is.null)]
    if (length(results_list) == 0) {
      warning("All bootstrap replicates failed.", call. = FALSE)
      return(list(Estimate=original_estimates, Boot.SE=NA, CI.Lower=NA, CI.Upper=NA))
    }
    estimates_mat <- do.call(rbind, results_list)
    se <- apply(estimates_mat, 2, sd, na.rm = TRUE)
    ci <- apply(estimates_mat, 2, quantile, probs = c(0.025, 0.975), na.rm = TRUE)
    
    list(Estimate = original_estimates, Boot.SE = se, CI.Lower = ci[1, ], CI.Upper = ci[2, ])
  }
  
  beta_summary <- summarize_boots(beta_results, model_object$beta$params)
  alpha_summary <- summarize_boots(alpha_results, model_object$alpha$params)
  tau_summary <- summarize_boots(tau_results, model_object$tau$params)
  omega_summary <- summarize_boots(omega_results, model_object$omega$std_dev)
  
  ed_time <- Sys.time()
  runtime <- format(round(ed_time - st_time, 2))
  
  output <- list(
    beta = beta_summary, alpha = alpha_summary, tau = tau_summary, omega = omega_summary,
    n_reps = length(beta_results[!sapply(beta_results, is.null)]), runtime = runtime
  )
  
  class(output) <- "mels_vmp_bootstrap"
  return(output)
}

#' Summarize a mels_vmp Bootstrap
#'
#' Prints a formatted summary of the non-parametric bootstrap results,
#' including the original estimates, bootstrap standard errors, and
#' 95% percentile confidence intervals.
#'
#' @details
#' The summary is printed to the console and organized into sections:
#' \itemize{
#'   \item Header: Reports the number of successful replicates and total runtime.
#'   \item Mean Model Parameters (beta): Estimates, Bootstrap SE, and CIs.
#'   \item Between-Subject Variance Parameters (alpha): Estimates, Bootstrap SE, and CIs.
#'   \item Within-Subject Variance Parameters (tau): Estimates, Bootstrap SE, and CIs.
#'   \item Random Effect Standard Deviation (omega): Estimates, Bootstrap SE, and CIs.
#' }
#'
#' @param object An object of class `mels_vmp_bootstrap`.
#' @param ... Additional arguments (currently unused).
#'
#' @return
#' Invisibly returns the original `mels_vmp_bootstrap` object.
#'
#' @seealso \code{\link{bootstrap_mels_vmp}}
#' @exportS3Method summary mels_vmp_bootstrap
#'
summary.mels_vmp_bootstrap <- function(object, ...) {
  # Header
  cat("## Bootstrap Summary for MELS Model ##\n")
  cat("--------------------------------------\n")
  cat(paste0("Successful replicates: ", object$n_reps, "\n"))
  cat(paste0("Total runtime: ", object$runtime, "\n\n"))
  
  # Mean Model Parameters (beta)
  cat("--- Mean Model Parameters (beta) ---\n")
  beta_df <- data.frame(object$beta)
  rownames(beta_df) <- row.names(object$beta$Estimate)
  print(round(beta_df, 4))
  cat("\n")
  
  # Alpha parameters
  cat("--- Between-Subject Variance Parameters (alpha) ---\n")
  alpha_df <- data.frame(object$alpha)
  rownames(alpha_df) <- row.names(object$alpha$Estimate)
  print(round(alpha_df, 4))
  cat("\n")
  
  # Tau parameters
  cat("--- Within-Subject Variance Parameters (tau) ---\n")
  tau_df <- data.frame(object$tau)
  rownames(tau_df) <- row.names(object$tau$Estimate)
  print(round(tau_df, 4))
  cat("\n")
  
  # Omega parameter
  cat("--- Random Effect Standard Deviation (omega) ---\n")
  omega_df <- data.frame(object$omega)
  rownames(omega_df) <- "omega_std_dev"
  print(round(omega_df, 4))
  cat("--------------------------------------\n")
  
  invisible(object)
}

