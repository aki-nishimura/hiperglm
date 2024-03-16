calc_loglik <- function(model, reg_coef, ...) {
  UseMethod("calc_loglik")
}

calc_loglink_deriv <- function(model, ...) {
  UseMethod("calc_loglink_deriv")
}

calc_grad <- function(model, reg_coef, ...) {
  loglink_grad <- calc_loglink_deriv(model, reg_coef, order = 1)
  grad <- t(model$design) %*% loglink_grad
  grad <- as.vector(grad)
  return(grad)
}

calc_loglik.linear_model <- function(model, reg_coef, ...) {
  args <- list(...)
  noise_var <- if ("noise_var" %in% names(args)) args$noise_var else 1

  design <- model$design
  outcome <- model$outcome

  predicted_val <- design %*% reg_coef
  loglik <- -0.5 * sum((outcome - predicted_val)^2) / noise_var
  return(loglik)
}

calc_loglik.logit_model <- function(model, reg_coef, ...) {
  design <- model$design
  outcome <- model$outcome

  if (is.list(outcome)) {
    n_success <- outcome$n_success
    n_trial <- outcome$n_trial
  } else {
    n_success <- outcome
    n_trial <- rep(1, length(n_success)) # Assume binary outcome
  }

  logit_prob <- design %*% reg_coef
  loglik <- sum(n_success * logit_prob - n_trial * log(1 + exp(logit_prob)))
  return(loglik)
}

calc_loglink_deriv.linear_model <- function(model, reg_coef, ...) {
  args <- list(...)
  noise_var <- if ("noise_var" %in% names(args)) args$noise_var else 1
  predicted_val <- model$design %*% reg_coef
  deriv <- (model$outcome - predicted_val) / noise_var
  deriv <- as.vector(deriv)
  return(deriv)
}

calc_loglink_deriv.logit_model <- function(model, reg_coef, ...) {
  args <- list(...)
  design <- model$design
  outcome <- model$outcome
  order <- if ("order" %in% names(args)) args$order else 1

  if (is.list(outcome)) {
    n_success <- outcome$n_success
    n_trial <- outcome$n_trial
  } else {
    n_success <- outcome
    n_trial <- rep(1, length(n_success)) # Assume binary outcome
  }

  logit_prob <- as.vector(design %*% reg_coef)
  predicted_prob <- 1 / (1 + exp(-logit_prob))
  if (order == 1) {
    deriv <- n_success - n_trial * predicted_prob
  } else if (order == 2) {
    deriv <- n_trial * predicted_prob * (1 - predicted_prob)
  } else {
    stop("3rd+ order derivative calculations are not supported")
  }
  deriv <- as.vector(deriv)
  return(deriv)
}

calc_loglik.poisson_model <- function(model, reg_coef, ...) {
  design <- model$design
  outcome <- model$outcome

  log_mu <- design %*% reg_coef
  loglik <- sum(outcome * log_mu - exp(log_mu) - lgamma(outcome + 1))
  return(loglik)
}

calc_loglink_deriv.poisson_model <- function(model, reg_coef, ...) {
  args <- list(...)
  design <- model$design
  outcome <- model$outcome
  order <- if ("order" %in% names(args)) args$order else 1

  log_mu <- as.vector(design %*% reg_coef)
  predicted_outcome <- exp(log_mu)
  if (order == 1) {
    deriv <- outcome - predicted_outcome
  } else if (order == 2) {
    deriv <- predicted_outcome
  } else {
    stop("3rd+ order derivative calculations are not supported")
  }
  deriv <- as.vector(deriv)
  return(deriv)
}

calc_logit_hessian <- function(model, reg_coef) {
  weight <- calc_loglink_deriv(model, reg_coef, order = 2)
  hess <- -t(model$design) %*% (outer(weight, rep(1, ncol(model$design))) * model$design)
  return(hess)
}

calc_logit_hessian_inverse <- function(model, reg_coef) {
  weight <- calc_loglink_deriv(model, reg_coef, order = 2)
  sqrt_weighted_design <- outer(sqrt(weight), rep(1, ncol(model$design))) * model$design
  R <- qr_wrapper(sqrt_weighted_design)$R
  inverse <- -invert_gram_mat_from_qr(R)
  return(inverse)
}
