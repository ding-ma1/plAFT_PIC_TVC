# forked from plAFT_RC_TVC_Gaussian_1_autosmooth.R on 20220406
# based on "plAFT_RC_TVC_Gaussian.R", with different strategy of calculating the knots and sigmas of Gaussian basis functions
# search "comment" and "under construction" for special comments
plAFT_GIC_TVC_Gaussian_1_autosmooth <- function(tmat, Xmat, Zmat, beta_initial, gamma_initial, theta_initial, smooth_initial = 1e-6, autosmooth_no = 5, numIntKnt = 7, maxiter = 10000, threshold = 1e-6, knots_min_quantile = 1, knots_max_quantile = 99, sig_coverage = 0.6, Hes_addition = 0, knots_fix_threshold = 1e-4, knots_fix_iter_no = 10000, frequent_knots_update = TRUE, eta = 3e-1, X_option = 2, knots_option = "equal_space", beta_Hessian_option = "modified", gamma_Hessian_option = "modified", draw_plot = TRUE) {
  # tmat includes two columns: 1. termination time, 
  #   2. censoring indicator (censor_n=1 for event time, censor_n=0 for right-censored time).
  #
  # Xmat includes columns: 1. individual id, 2-. time-invariant covariate(s) X.
  #
  # Zmat includes columns: 1. individual id, 2. time interval index, 
  #   time intervals (two columns, 3. interval_L and 4. interval_R), 
  #   5-. time-varying covariate(s) Z.
  # browser()
  # handle missing arguments
  if (missing(Xmat)) {stop("Please include the time-invariant covariates (Xmat) in the input arguments.\n")}
  if (missing(tmat)) {stop("Please specify follow up time and censoring indicator (tmat).\n")}
  if (missing(Zmat) & !missing(gamma_initial)) {
    stop("No need to specify initial value for gamma if time-varying covariates (Zmat) are not included.\n")
  }
  if (!missing(Zmat) & missing(gamma_initial)) {
    stop("Please specify intial value for gamma (gamma_initial).\n")
  }
  
  # library(splines2)
  
  # time[is.infinite(time) | time==0] = NA
  
  # information from tmat matrix
  tmat = as.matrix(tmat)
  rownames(tmat) = NULL
  time = tmat[, 1:2, drop = FALSE] # n by 2 matrix
  rownames(time) = 1:(nrow(time))
  censor_n = tmat[, 3, drop = FALSE] # n by 1 matrix
  time[censor_n==2L, 2] = NA # for right censoring, time[,2]==0
  time[censor_n==3L, 1] = 0 # for left censoring, time[,1]==0
  if (any(time[censor_n!=2, 1] > time[censor_n!=2, 2])) { # check time[,1]<=time[,2] except right-censored cases
    stop(which(time[censor_n!=2, 1] > time[censor_n!=2, 2]), " th row(s) in tmat: input error.\n")
  }
  
  # information from Xmat matrix
  Xmat = as.matrix(Xmat)
  rownames(Xmat) = NULL
  # id_X = Xmat[, 1] # n vector
  X = Xmat[, -1, drop = FALSE] # n by p matrix, 20200330 the negative doesn't work with class as data.frame if drop = FALSE
  n = dim(X)[1]
  p = dim(X)[2] # number of time-invariant covariates
  if (X_option==1) {X = X - matrix(1, n, 1) %*% colMeans(X, na.rm = TRUE)} # standardisation by subtracting the mean of each column
  if (X_option==2) {X = X} # without standardisation
  
  n_E = sum(censor_n==1) # n_E denotes number of event observations
  n_RC = sum(censor_n==2) # n_RC denotes number of right-censored observations
  n_LC = sum(censor_n==3) # n_LC denotes number of left-censored observations
  n_IC = sum(censor_n==4) # n_IC denotes number of interval-censored observations
  
  if (n_E==0) {time_E = NA}
  else {time_E = time[censor_n==1L, 2, drop = FALSE]}
  if (n_RC==0) {time_RC = NA}
  else {time_RC = time[censor_n==2L, 1, drop = FALSE]}
  if (n_LC==0) {time_LC = NA}
  else {time_LC = time[censor_n==3L, 2, drop = FALSE]}
  if (n_IC==0) {time_IC_L = time_IC_R = NA}
  else {
    time_IC_L = time[censor_n==4L, 1, drop = FALSE]
    time_IC_R = time[censor_n==4L, 2, drop = FALSE]
  }
  time_bind = rbind(time_E, time_RC, time_LC, time_IC_L, time_IC_R)
  time_bind = time_bind[!is.na(time_bind), , drop = FALSE]
  time_bind = time_bind[order(as.numeric(rownames(time_bind))), , drop = FALSE]
  
  if (n_E==0) {X_E = matrix(0, 1, p)}
  else {X_E = X[censor_n==1L, , drop = FALSE]} # n_E by p
  if (n_RC==0) {X_RC = matrix(0, 1, p)}
  else {X_RC = X[censor_n==2L, , drop = FALSE]} # n_RC by p
  if (n_LC==0) {X_LC = matrix(0, 1, p)}
  else {X_LC = X[censor_n==3L, , drop = FALSE]} # n_LC by p
  if (n_IC==0) {X_IC = matrix(0, 1, p)}
  else {X_IC = X[censor_n==4L, , drop = FALSE]} # n_IC by p
  
  # information from Zmat matrix
  if (missing(Zmat) | missing(gamma_initial)) {dim_Zmat = c(0, 0)} # 20210126 added to accommodate (missing(Zmat) & missing(gamma_initial))
  else {dim_Zmat = dim(Zmat)} # 20210126 added to accommodate (missing(Zmat) & missing(gamma_initial))
  
  # 20201112 added (dim(Zmat)[1]>(dim(Xmat)[1]+1)) condition to avoid computation issue when none or only one individual took Z treatment
  # if (!missing(Zmat) & !missing(gamma_initial) & (dim(Zmat)[1]>(dim(Xmat)[1]+1))) { #20201112
  if (!missing(Zmat) & !missing(gamma_initial) & (dim_Zmat[1]>n+1)) { # 20210126 modified to accommodate (missing(Zmat) & missing(gamma_initial))
    Zmat = as.matrix(Zmat)
    rownames(Zmat) = NULL
    id_Z = Zmat[, 1] # N vector, N is greater than or equal to n
    # interval_index = Zmat[, 2] # N vector
    # ni = tapply(interval_index, id_Z, FUN=which.max) # n vector, change FUN from max to which.max on 20220318, robust
    ni = diff(c(0, which(c(diff(id_Z)!=0, TRUE)))) # n vector, 20220318, robust, not needing interval_index any more
    interval_length = Zmat[, 4, drop = FALSE] - Zmat[, 3, drop = FALSE] # N by 1 matrix
    Z = Zmat[, -c(1:4), drop = FALSE] # N by q matrix
    N = dim(Z)[1]
    q = dim(Z)[2] # number of time-varying covariates
    Z_ni = Z[c(diff(id_Z)!=0, TRUE), , drop = FALSE] # n by q matrix
    if (n_E==0) {Z_ni_E = matrix(0, 1, q)}
    else {Z_ni_E = Z_ni[censor_n==1L, , drop = FALSE]} # n_E by q matrix
    # for full Hessian of gamma
    # if (gamma_Hessian_option=="full") {
    ni_E = ni[censor_n==1L] # n_E vector
    ni_RC = ni[censor_n==2L] # n_RC vector
    ni_LC = ni[censor_n==3L] # n_LC vector
    ni_IC = ni[censor_n==4L] # n_IC vector
    
    censor_N = rep(censor_n, times = ni) # N vector
    
    N_E = sum(censor_N==1L)
    N_RC = sum(censor_N==2L)
    N_LC = sum(censor_N==3L)
    N_IC = sum(censor_N==4L)
    
    id_Z_E = id_Z[censor_N==1L] # N_E vector
    id_Z_RC = id_Z[censor_N==2L] # N_RC vector
    id_Z_LC = id_Z[censor_N==3L] # N_LC vector
    id_Z_IC = id_Z[censor_N==4L] # N_IC vector
    
    # interval_length_E = Zmat[censor_N==1L, 4, drop = FALSE] - Zmat[censor_N==1L, 3, drop = FALSE] # N_E by 1 matrix
    interval_length_E = interval_length[censor_N==1L, , drop = FALSE] # N_E by 1 matrix
    interval_length_RC = interval_length[censor_N==2L, , drop = FALSE] # N_RC by 1 matrix
    interval_length_LC = interval_length[censor_N==3L, , drop = FALSE] # N_LC by 1 matrix
    interval_length_IC = interval_length[censor_N==4L, , drop = FALSE] # N_IC by 1 matrix
    
    if (N_E==0) {Z_E = matrix(0, 1, q)}
    else {Z_E = Z[censor_N==1L, , drop = FALSE]} # N_E by 1 matrix
    if (N_RC==0) {Z_RC = matrix(0, 1, q)}
    else {Z_RC = Z[censor_N==2L, , drop = FALSE]} # N_RC by 1 matrix
    if (N_LC==0) {Z_LC = matrix(0, 1, q)}
    else {Z_LC = Z[censor_N==3L, , drop = FALSE]} # N_LC by 1 matrix
    if (N_IC==0) {Z_IC = matrix(0, 1, q)}
    else {Z_IC = Z[censor_N==4L, , drop = FALSE]} # N_IC by 1 matrix
    # }
    if (N_IC==0) {
      id_Z_IC_subset = numeric(0)
      ni_IC_subset = numeric(0)
      interval_length_IC_subset = NA
      Z_IC_subset = matrix(0, 1, q)
      N_IC_subset = 0
      # X_IC_subset = matrix(0, 1, p)
      # n_IC_subset = 0
    }
    else {
      Z_IC_subset_indicator = (Zmat[censor_N==4L, 4, drop = FALSE] <= rep(time_IC_L, times = ni_IC))
      id_Z_IC_subset = id_Z_IC[Z_IC_subset_indicator==TRUE] # < N_IC vector, get the id_Z's of interval-censored entries till the left censored time points
      # interval_index_IC = Zmat[censor_N==4L, 2] # N_IC vector
      # interval_index_IC_subset = interval_index_IC[Z_IC_subset_indicator == TRUE]
      # ni_IC_subset = as.numeric(tapply(interval_index_IC_subset, id_Z_IC_subset, FUN=max)) # n_IC vector
      ni_IC_subset = diff(c(0, which(c(diff(id_Z_IC_subset)!=0, TRUE))))
      interval_length_IC_subset = interval_length_IC[Z_IC_subset_indicator == TRUE, , drop = FALSE]
      Z_IC_subset = Z_IC[Z_IC_subset_indicator == TRUE, , drop = FALSE]
      N_IC_subset = sum(Z_IC_subset_indicator)
      # X_IC_subset = X[unique(id_Z_IC_subset), , drop = FALSE]
      # n_IC_subset = length(ni_IC_subset)
    }
  }
  # browser()
  
  # Gaussian basis part, necessary functions to calculate the derivative and integral of Gaussian basis
  Gaussian_basis <- function(x, mean, sd) {
    if (nrow(x)==1) { # x must be matrix class (i.e. length(dim(x))==2)
      if (length(sd)==1) {
        return(matrix(sapply(X = 1:(numIntKnt+2), FUN = function(i) {dnorm(x, mean[i], sd) * sd*sqrt(2*pi)}), nrow = 1))
      }
      else {return(matrix(sapply(X = 1:(numIntKnt+2), FUN = function(i) {dnorm(x, mean[i], sd[i]) * sd[i]*sqrt(2*pi)}), nrow = 1))}
    }
    else {
      if (length(sd)==1) {
        return(sapply(X = 1:(numIntKnt+2), FUN = function(i) {dnorm(x, mean[i], sd) * sd*sqrt(2*pi)}))
      }
      else {return(sapply(X = 1:(numIntKnt+2), FUN = function(i) {dnorm(x, mean[i], sd[i]) * sd[i]*sqrt(2*pi)}))}
    }
  }

  Gaussian_basis_integ1 <- function(q, mean, sd) {
    if (nrow(q)==1) { # q must be matrix class (i.e. length(dim(q))==2)
      if (length(sd)==1) {
        return(matrix(sapply(X = 1:(numIntKnt+2), FUN = function(i) {(pnorm(q, mean[i], sd) - matrix(1, NROW(q), 1)%*%pnorm(0, mean[i], sd)) * sd*sqrt(2*pi)}), nrow = 1))
      }
      else {return(matrix(sapply(X = 1:(numIntKnt+2), FUN = function(i) {(pnorm(q, mean[i], sd[i]) - matrix(1, NROW(q), 1)%*%pnorm(0, mean[i], sd[i])) * sd[i]*sqrt(2*pi)}), nrow = 1))}
    }
    else {
      if (length(sd)==1) {
        return(sapply(X = 1:(numIntKnt+2), FUN = function(i) {(pnorm(q, mean[i], sd) - matrix(1, NROW(q), 1)%*%pnorm(0, mean[i], sd)) * sd*sqrt(2*pi)}))
      }
      else {return(sapply(X = 1:(numIntKnt+2), FUN = function(i) {(pnorm(q, mean[i], sd[i]) - matrix(1, NROW(q), 1)%*%pnorm(0, mean[i], sd[i])) * sd[i]*sqrt(2*pi)}))}
    }
  }

  Gaussian_basis_deriv1 <- function(x, mean, sd) {
    if (nrow(x)==1) { # x must be matrix class (i.e. length(dim(x))==2)
      if (length(sd)==1) {
        return(matrix(sapply(X = 1:(numIntKnt+2), FUN = function(i) {-(x-mean[i])/sd^2 * dnorm(x, mean[i], sd) * sd*sqrt(2*pi)}), nrow = 1))
      }
      else {return(matrix(sapply(X = 1:(numIntKnt+2), FUN = function(i) {-(x-mean[i])/sd[i]^2 * dnorm(x, mean[i], sd[i]) * sd[i]*sqrt(2*pi)}), nrow = 1))}
    }
    else {
      if (length(sd)==1) {
        return(sapply(X = 1:(numIntKnt+2), FUN = function(i) {-(x-mean[i])/sd^2 * dnorm(x, mean[i], sd) * sd*sqrt(2*pi)}))
      }
      else {return(sapply(X = 1:(numIntKnt+2), FUN = function(i) {-(x-mean[i])/sd[i]^2 * dnorm(x, mean[i], sd[i]) * sd[i]*sqrt(2*pi)}))}
    }
  }

  Gaussian_basis_deriv2 <- function(x, mean, sd) {
    if (nrow(x)==1) { # x must be matrix class (i.e. length(dim(x))==2)
      if (length(sd)==1) {
        return(matrix(sapply(X = 1:(numIntKnt+2), FUN = function(i) {((x-mean[i])^2/sd^4 - 1/sd^2) * dnorm(x, mean[i], sd) * sd*sqrt(2*pi)}), nrow = 1))
      }
      else {return(matrix(sapply(X = 1:(numIntKnt+2), FUN = function(i) {((x-mean[i])^2/sd[i]^4 - 1/sd[i]^2) * dnorm(x, mean[i], sd[i]) * sd[i]*sqrt(2*pi)}), nrow = 1))}
    }
    else {
      if (length(sd)==1) {
        return(sapply(X = 1:(numIntKnt+2), FUN = function(i) {((x-mean[i])^2/sd^4 - 1/sd^2) * dnorm(x, mean[i], sd) * sd*sqrt(2*pi)}))
      }
      else {return(sapply(X = 1:(numIntKnt+2), FUN = function(i) {((x-mean[i])^2/sd[i]^4 - 1/sd[i]^2) * dnorm(x, mean[i], sd[i]) * sd[i]*sqrt(2*pi)}))}
    }
  }
  
  
  # # inverse_quadratic_basis
  # Gaussian_basis <- function(x, mean, sd) {
  #   if (nrow(x)==1) {
  #     if (length(sd)==1) {
  #       return(matrix(sapply(X = 1:(numIntKnt+2), FUN = function(i) {1/(1+((x-mean[i])/sd)^2)}), nrow = 1))
  #     }
  #     else {return(matrix(sapply(X = 1:(numIntKnt+2), FUN = function(i) {1/(1+((x-mean[i])/sd[i])^2)}), nrow = 1))}
  #   }
  #   else {
  #     if (length(sd)==1) {
  #       return(sapply(X = 1:(numIntKnt+2), FUN = function(i) {1/(1+((x-mean[i])/sd)^2)}))
  #     }
  #     else {return(sapply(X = 1:(numIntKnt+2), FUN = function(i) {1/(1+((x-mean[i])/sd[i])^2)}))}
  #   }
  # }
  # # inverse_quadratic_basis_integ1
  # Gaussian_basis_integ1 <- function(q, mean, sd) {
  #   if (nrow(q)==1) {
  #     if (length(sd)==1) {
  #       return(matrix(sapply(X = 1:(numIntKnt+2), FUN = function(i) {atan(abs(q-mean[i])/sd)}), nrow = 1))
  #     }
  #     else {return(matrix(sapply(X = 1:(numIntKnt+2), FUN = function(i) {atan(abs(q-mean[i])/sd[i])}), nrow = 1))}
  #   }
  #   else {
  #     if (length(sd)==1) {
  #       return(sapply(X = 1:(numIntKnt+2), FUN = function(i) {atan(abs(q-mean[i])/sd)}))
  #     }
  #     else {return(sapply(X = 1:(numIntKnt+2), FUN = function(i) {atan(abs(q-mean[i])/sd[i])}))}
  #   }
  # }
  # # inverse_quadratic_basis_deriv1
  # Gaussian_basis_deriv1 <- function(x, mean, sd) {
  #   if (nrow(x)==1) {
  #     if (length(sd)==1) {
  #       return(matrix(sapply(X = 1:(numIntKnt+2), FUN = function(i) {-2*(x-mean[i])/sd/(1+((x-mean[i])/sd)^2)^2}), nrow = 1))
  #     }
  #     else {return(matrix(sapply(X = 1:(numIntKnt+2), FUN = function(i) {-2*(x-mean[i])/sd[i]/(1+((x-mean[i])/sd[i])^2)^2}), nrow = 1))}
  #   }
  #   else {
  #     if (length(sd)==1) {
  #       return(sapply(X = 1:(numIntKnt+2), FUN = function(i) {-2*(x-mean[i])/sd/(1+((x-mean[i])/sd)^2)^2}))
  #     }
  #     else {return(sapply(X = 1:(numIntKnt+2), FUN = function(i) {-2*(x-mean[i])/sd[i]/(1+((x-mean[i])/sd[i])^2)^2}))}
  #   }
  # }
  # # inverse_quadratic_basis_deriv2
  # Gaussian_basis_deriv2 <- function(x, mean, sd) {
  #   if (nrow(x)==1) {
  #     if (length(sd)==1) {
  #       return(matrix(sapply(X = 1:(numIntKnt+2), FUN = function(i) {(-2/sd*(1+((x-mean[i])/sd)^2) + 8*((x-mean[i])/sd)^2)/(1+((x-mean[i])/sd)^2)^3}), nrow = 1))
  #     }
  #     else {return(matrix(sapply(X = 1:(numIntKnt+2), FUN = function(i) {(-2/sd[i]*(1+((x-mean[i])/sd[i])^2) + 8*((x-mean[i])/sd[i])^2)/(1+((x-mean[i])/sd[i])^2)^3}), nrow = 1))}
  #   }
  #   else {
  #     if (length(sd)==1) {
  #       return(sapply(X = 1:(numIntKnt+2), FUN = function(i) {(-2/sd*(1+((x-mean[i])/sd)^2) + 8*((x-mean[i])/sd)^2)/(1+((x-mean[i])/sd)^2)^3}))
  #     }
  #     else {return(sapply(X = 1:(numIntKnt+2), FUN = function(i) {(-2/sd[i]*(1+((x-mean[i])/sd[i])^2) + 8*((x-mean[i])/sd[i])^2)/(1+((x-mean[i])/sd[i])^2)^3}))}
  #   }
  # }
  
  
  # outer loop
  smooth_iter = matrix(NA, autosmooth_no, 1)
  smooth = smooth_initial
  df_old = 0
  for (j in 1:autosmooth_no) {
    cat(j, "th outer loop begins. The current smoothing parameter is", smooth, ".\n")
    smooth_iter[j] = smooth
    
    # inner loop
    
    # dgrSp = ordSp - 1
    # numSp = numIntKnt + ordSp
    
    # 0. initial value of beta, gamma and theta
    # library(survival)
    # iniFit = survreg(Surv(time = time[, 1], time2 = time[, 2], type = "interval2") ~ X, dist="weibull")
    # beta_old = matrix(iniFit$coef[2:(p+1)], ncol = 1)
    # beta_old = matrix(0, p, 1) # alternative initial value for beta
    if (missing(beta_initial)) {
      beta_old = matrix(0L, p, 1)
      cat("beta_initial is not provided. Set beta_initial to", beta_old, ".\n")
    }
    else {beta_old = matrix(beta_initial)} # p by 1 matrix
    # theta_old = matrix(1L, numSp, 1) # numSp by 1 matrix
    if (missing(theta_initial)) {theta_old = matrix(1L, numIntKnt+2, 1)} # Gaussian basis part, 1 less knots compare to spline basis
    else {theta_old = theta_initial}
    
    # outputs of this function (faster than the alternative assignment)
    grad_beta_iter = matrix(NA, maxiter, p)
    Hes_beta_iter = matrix(NA, maxiter, p^2)
    det_Hes_beta_iter = matrix(NA, maxiter, 1)
    eigenvl_Hes_beta_iter = matrix(NA, maxiter, p)
    likelihood_beta_iter_before_ls = matrix(NA, maxiter, 1)
    iter_Newton_beta = matrix(NA, maxiter, 1)
    likelihood_beta_iter_after_ls = matrix(NA, maxiter, 1)
    ts_range_beta_iter = matrix(NA, maxiter, 2) # ts_bind_range after beta update
    beta_iter = matrix(NA, maxiter, p)
    
    # grad_theta_iter = matrix(NA, maxiter, numSp)
    deno_iter = matrix(NA, maxiter, numIntKnt+2) # 20220502
    multiplier_iter = matrix(NA, maxiter, numIntKnt+2) # 20220502
    grad_theta_iter = matrix(NA, maxiter, numIntKnt+2) # Gaussian basis part, 1 less theta's compare to spline basis
    penlike_theta_iter_before_ls = matrix(NA, maxiter, 1)
    iter_MI_theta = matrix(NA, maxiter, 1)
    penlike_theta_iter_after_ls = matrix(NA, maxiter, 1)
    ts_range_theta_iter = matrix(NA, maxiter, 2) # ts_bind_range after theta update
    # theta_iter = matrix(NA, maxiter, numSp)
    theta_iter = matrix(NA, maxiter, numIntKnt+2)# Gaussian basis part, 1 less theta's compare to spline basis
    # Alternative assignment
    # likelihood_beta_iter_after_ls = NULL ## example: likelihood_beta_iter_after_ls = c(likelihood_beta_iter_after_ls, likelihood_beta_iter_1)
    # iter_Newton_beta = NULL
    # penlike_theta_iter_after_ls = NULL
    # iter_MI_theta = NULL
    # grad_beta_iter = NULL
    # ts_range_iter = NULL
    # beta_iter = NULL
    # theta_iter = NULL
    knots_iter = matrix(NA, maxiter, numIntKnt+2)
    
    if (!missing(Zmat) & !missing(gamma_initial) & (dim_Zmat[1]<=n+1)) {
      warning("There is NO time-varying covariates in this dataset!\n")
      # gamma_old = matrix(gamma_initial) # comment this line in simulations considering time-varying covariates
      # gamma_new = gamma_old + 1 # comment this line in simulations considering time-varying covariates
    }
    
    if (!missing(Zmat) & !missing(gamma_initial) & (dim_Zmat[1]>n+1)) {
      gamma_old = matrix(gamma_initial) # q by 1
      gamma_new = gamma_old + 1
      
      grad_gamma_iter = matrix(NA, maxiter, q)
      Hes_gamma_iter = matrix(NA, maxiter, q^2)
      det_Hes_gamma_iter = matrix(NA, maxiter, 1)
      eigenvl_Hes_gamma_iter = matrix(NA, maxiter, q)
      likelihood_gamma_iter_before_ls = matrix(NA, maxiter, 1)
      iter_Newton_gamma = matrix(NA, maxiter, 1)
      likelihood_gamma_iter_after_ls = matrix(NA, maxiter, 1)
      ts_range_gamma_iter = matrix(NA, maxiter, 2) # ts_bind_range after gamma update
      gamma_iter = matrix(NA, maxiter, q)
    }
    else {
      gamma_old = 0
      gamma_new = gamma_old + knots_fix_threshold
    }
    
    # 0. accelerated time
    if (!missing(Zmat) & !missing(gamma_initial) & (dim_Zmat[1]>n+1)) {
      if (n_E==0) {ts_E = NA}
      else {ts_E = exp(-X_E%*%beta_old) * as.matrix(tapply(exp(-Z_E%*%gamma_old)*interval_length_E, id_Z_E, FUN = sum))}
      if (n_RC==0) {ts_RC = NA}
      else {ts_RC = exp(-X_RC%*%beta_old) * as.matrix(tapply(exp(-Z_RC%*%gamma_old)*interval_length_RC, id_Z_RC, FUN = sum))}
      if (n_LC==0) {ts_LC = NA}
      else {ts_LC = exp(-X_LC%*%beta_old) * as.matrix(tapply(exp(-Z_LC%*%gamma_old)*interval_length_LC, id_Z_LC, FUN = sum))}
      if (n_IC==0) {ts_IC_L = ts_IC_R = NA}
      else {
        ts_IC_L = exp(-X_IC%*%beta_old) * as.matrix(tapply(exp(-Z_IC_subset%*%gamma_old)*interval_length_IC_subset, id_Z_IC_subset, FUN = sum)) # special requirement
        ts_IC_R = exp(-X_IC%*%beta_old) * as.matrix(tapply(exp(-Z_IC%*%gamma_old)*interval_length_IC, id_Z_IC, FUN = sum))
      }
    }
    else {
      if (n_E==0) {ts_E = NA}
      else {ts_E = time_E*exp(-X_E%*%beta_old)} # n_E by 1 matrix
      if (n_RC==0) {ts_RC = NA}
      else {ts_RC = time_RC*exp(-X_RC%*%beta_old)} # n_RC by 1 matrix
      if (n_LC==0) {ts_LC = NA}
      else {ts_LC = time_LC*exp(-X_LC%*%beta_old)} # n_LC by 1 matrix
      if (n_IC==0) {ts_IC_L = ts_IC_R = NA}
      else {
        ts_IC_L = time_IC_L*exp(-X_IC%*%beta_old) # n_IC by 1 matrix
        ts_IC_R = time_IC_R*exp(-X_IC%*%beta_old) # n_IC by 1 matrix
      }
    }
    ts_bind = rbind(ts_E, ts_RC, ts_LC, ts_IC_L, ts_IC_R) # ? not sure if this is suitable for "percentile" knots
    ts_bind = ts_bind[!is.na(ts_bind), , drop = FALSE]
    ts_bind = ts_bind[order(as.numeric(rownames(ts_bind))), , drop = FALSE]
    # ts_E = as.matrix(ts[censor_n==1])
    # ts_E = ts[censor_n==1, , drop = FALSE]
    
    # 0. knots
    # # knots option 1
    # if (knots_option==1) {
    #   bryKnt = c(min(ts), max(ts)+1e-40)
    #   IntKnt = quantile(sort(ts)[2:(n-1)], seq(5, 95, length.out = numIntKnt)/100, type=1)
    # }
    # # knots option 2
    # if (knots_option==2) {
    #   Knt = quantile(sort(ts), seq(0, 1, length.out = numSp-dgrSp+1), type=1)
    #   bryKnt = c(Knt[1], Knt[numSp-dgrSp+1]+1e-40)
    #   IntKnt = Knt[2:(numSp-dgrSp)]
    # }
    # # knots option 3
    # if (knots_option==3) {
    #   ts_min = min(ts)
    #   ts_max = max(ts)
    #   bryKnt = c(ts_min, ts_max+1e-40)
    #   IntKnt = seq(ts_min, ts_max, length.out = numSp-dgrSp+1)[2:(numSp-dgrSp)]
    # }
    # bryKnt_initial = bryKnt
    # IntKnt_initial = IntKnt
    # Gaussian basis part
    if (knots_option=="equal_space") {
      # bryKnt = range(ts)
      bryKnt = range(ts_bind)
      bin_width = (bryKnt[2] - bryKnt[1]) / (numIntKnt+1)
      gknots = seq(bryKnt[1], bryKnt[2], bin_width)
      sig = (2/3) * bin_width
      # gauss = pnorm(0, gknots, sig)*sqrt(2*pi*sig^2)
      # gauss = pnorm(bryKnt[1], gknots, sig)*sqrt(2*pi*sig^2) # which one to choose?
    }
    else if (knots_option=="percentile") {
      # different strategy of calculating knots and sigmas
      gknots = quantile(sort(ts_bind), seq(knots_min_quantile, knots_max_quantile, length.out = numIntKnt+2)/100, type=1)
      dist = sapply(gknots, function(x) {abs(x - ts_bind)})
      sig = apply(dist, 2, function(x) {quantile(x, sig_coverage) / 2})
    }
    else if (knots_option=="percentile_1") { # test on 20220525
      ts_bind_1 = rbind(ts_E, ts_RC, ts_LC, ts_IC_R)
      ts_bind_1 = ts_bind_1[!is.na(ts_bind_1), , drop = FALSE]
      ts_bind_1 = ts_bind_1[order(as.numeric(rownames(ts_bind_1))), , drop = FALSE]
      gknots = quantile(sort(ts_bind_1), seq(knots_min_quantile, knots_max_quantile, length.out = numIntKnt+2)/100, type=1)
      dist = sapply(gknots, function(x) {abs(x - ts_bind_1)})
      sig = apply(dist, 2, function(x) {quantile(x, sig_coverage) / 2})
    }
    else {stop("knots_option is either percentile or equal_space.\n")}
    
    
    # # Gaussian basis part
    # # bryKnt = range(ts)
    # # basis_coverage = 0.3
    # ts_sorted = sort(ts)
    # basis_coverage = 1 / (numSp - 1)
    # sig = rep(NA, numSp - 1)
    # gknots = quantile(ts_sorted, seq(5, 95, length.out = numSp - 1)/100, type=1)
    # sig[1] = (quantile(ts_sorted, basis_coverage)-gknots[1])/1.96
    # # bin_width = (bryKnt[2] - bryKnt[1]) / (numIntKnt+1)
    # # gknots = seq(bryKnt[1], bryKnt[2], bin_width)
    # # gknots = quantile(sort(ts), seq(0, 1, length.out = numSp-dgrSp+1), type=1)
    # # sig = (2/3) * bin_width
    # # gauss = pnorm(0, gknots, sig)*sqrt(2*pi*sig^2)
    # 
    # sum(ts>gknots[1]-1.96*sig[1] & ts<gknots[1]+1.96*sig[1])/n
    # 
    # browser()
    
    
    
    
    # 0. spline and cumulative spline
    # if (n_E==0) {phi_E = matrix(0, 1, numSp)}
    # else {
    #   phi_E = mSpline(ts_E, knots=IntKnt, degree=dgrSp, intercept=TRUE, Boundary.knots=bryKnt)
    # }
    # cphi = iSpline(ts, knots=IntKnt, degree=dgrSp, intercept=TRUE, Boundary.knots=bryKnt)
    
    # 0. hazard and cumulative hazard
    # if (n_E==0) {h0ts_E = h0tss_E = 1}
    # else {
    #   h0ts_E = phi_E%*%theta_old
    #   h0tss_E = h0ts_E
    #   h0tss_E[h0tss_E<1e-40] = 1e-40
    # }
    # ch0ts = cphi%*%theta_old
    
    # Gaussian basis part
    # h0ts = Gaussian_basis(x = ts, mean = gknots, sd = sig)%*%theta_old
    # ch0ts = Gaussian_basis_integ1(q = ts, mean = gknots, sd = sig)%*%theta_old
    
    if (n_E==0) {
      h0ts_E = h0tss_E = 1
      ch0ts_E = 0
    }
    else {
      h0ts_E = Gaussian_basis(x = ts_E, mean = gknots, sd = sig)%*%theta_old
      h0tss_E = h0ts_E
      h0tss_E[h0tss_E<1e-40] = 1e-40
      ch0ts_E = Gaussian_basis_integ1(q = ts_E, mean = gknots, sd = sig)%*%theta_old
    }
    if (n_RC==0) {ch0ts_RC = 0}
    else {ch0ts_RC = Gaussian_basis_integ1(q = ts_RC, mean = gknots, sd = sig)%*%theta_old}
    if (n_LC==0) {
      S0ts_LC = 0
      S0ts_LC_diff = 1
    }
    else {
      S0ts_LC = exp(-Gaussian_basis_integ1(q = ts_LC, mean = gknots, sd = sig)%*%theta_old)
      S0ts_LC_diff = 1 - S0ts_LC
      S0ts_LC_diff[S0ts_LC_diff<1e-40] = 1e-40 # 20220514
    }
    if (n_IC==0) {
      S0ts_IC_L = S0ts_IC_R = 0
      S0ts_IC_diff = 1
    }
    else {
      S0ts_IC_L = exp(-Gaussian_basis_integ1(q = ts_IC_L, mean = gknots, sd = sig)%*%theta_old)
      S0ts_IC_R = exp(-Gaussian_basis_integ1(q = ts_IC_R, mean = gknots, sd = sig)%*%theta_old)
      S0ts_IC_diff = S0ts_IC_L - S0ts_IC_R
      S0ts_IC_diff[S0ts_IC_diff<1e-40] = 1e-40 # 20220514
    }
    # browser()
    # 0. log-likelihood
    if (!missing(Zmat) & !missing(gamma_initial) & (dim_Zmat[1]>n+1)) {
      like_beta_old = sum(log(h0tss_E) - X_E%*%beta_old - Z_ni_E%*%gamma_old - ch0ts_E) - sum(ch0ts_RC) + sum(log(S0ts_LC_diff)) + sum(log(S0ts_IC_diff))
    }
    else {like_beta_old = sum(log(h0tss_E) - X_E%*%beta_old - ch0ts_E) - sum(ch0ts_RC) + sum(log(S0ts_LC_diff)) + sum(log(S0ts_IC_diff))}
    
    ####################### main iterative loop #########################
    # tryCatch({ #20181129
    for(k in 1L:maxiter) {
      ########################## Newton step for beta ##########################
      # cat(k, "th iteration.\n")
      # 0.5. differential spline for gradient and modified Hessian
      # if (n_E==0) {dphi_E = matrix(0, 1, numSp)}
      # else {dphi_E = mSpline(ts_E, knots=IntKnt, degree=dgrSp, intercept=TRUE, Boundary.knots=bryKnt, derivs=1)}
      # phi = mSpline(ts, knots=IntKnt, degree=dgrSp, intercept=TRUE, Boundary.knots=bryKnt)
      
      # 0.5. segmentations for gradient and modified Hessian
      # dh0ts_E = dphi_E%*%theta_old
      # h0ts_E = phi_E%*%theta_old
      # h0tss_E = h0ts_E
      # h0tss_E[h0tss_E<1e-40] = 1e-40
      # h0ts = phi%*%theta_old
      
      # Gaussian basis part
      # dh0ts = Gaussian_basis_deriv1(x = ts, mean = gknots, sd = sig)%*%theta_old
      # h0ts = Gaussian_basis(x = ts, mean = gknots, sd = sig)%*%theta_old
      
      # dh0ts = matrix(0, n, 1)
      # h0ts = matrix(0, n, 1)
      # for (i in 1:n) { # will improve using apply function
      #   dh0ts[i] = t(dnorm_deriv1(ts[i], gknots, sig)*sqrt(2*pi*sig^2))%*%theta_old
      #   h0ts[i] = t(dnorm(ts[i], gknots, sig)*sqrt(2*pi*sig^2))%*%theta_old
      # }
      if (n_E==0) {
        dh0ts_E = 0
        h0ts_E = h0tss_E = 1
        if (beta_Hessian_option=="full") {ddh0ts_E = 0}
      }
      else {
        dh0ts_E = Gaussian_basis_deriv1(x = ts_E, mean = gknots, sd = sig)%*%theta_old
        h0ts_E = Gaussian_basis(x = ts_E, mean = gknots, sd = sig)%*%theta_old
        h0tss_E = h0ts_E
        h0tss_E[h0tss_E<1e-40] = 1e-40
        if (beta_Hessian_option=="full") {ddh0ts_E = Gaussian_basis_deriv2(x = ts_E, mean = gknots, sd = sig)%*%theta_old}
      }
      if (n_RC==0) {
        h0ts_RC = 0
        if (beta_Hessian_option=="full") {dh0ts_RC = 0}
      }
      else {
        h0ts_RC = Gaussian_basis(x = ts_RC, mean = gknots, sd = sig)%*%theta_old
        if (beta_Hessian_option=="full") {dh0ts_RC = Gaussian_basis_deriv1(x = ts_RC, mean = gknots, sd = sig)%*%theta_old}
      }
      if (n_LC==0) {
        h0ts_LC = S0ts_LC = 0
        S0ts_LC_diff = 1
        if (beta_Hessian_option=="full") {dh0ts_LC = 0}
      }
      else {
        h0ts_LC = Gaussian_basis(x = ts_LC, mean = gknots, sd = sig)%*%theta_old
        S0ts_LC = exp(-Gaussian_basis_integ1(q = ts_LC, mean = gknots, sd = sig)%*%theta_old)
        S0ts_LC_diff = 1 - S0ts_LC
        S0ts_LC_diff[S0ts_LC_diff<1e-40] = 1e-40 # 20220514
        if (beta_Hessian_option=="full") {dh0ts_LC = Gaussian_basis_deriv1(x = ts_LC, mean = gknots, sd = sig)%*%theta_old}
      }
      if (n_IC==0) {
        h0ts_IC_L = S0ts_IC_L = h0ts_IC_R = S0ts_IC_R = 0
        S0ts_IC_diff = 1
        if (beta_Hessian_option=="full") {dh0ts_IC_L = dh0ts_IC_R = 0}
      }
      else {
        h0ts_IC_L = Gaussian_basis(x = ts_IC_L, mean = gknots, sd = sig)%*%theta_old
        S0ts_IC_L = exp(-Gaussian_basis_integ1(q = ts_IC_L, mean = gknots, sd = sig)%*%theta_old)
        h0ts_IC_R = Gaussian_basis(x = ts_IC_R, mean = gknots, sd = sig)%*%theta_old
        S0ts_IC_R = exp(-Gaussian_basis_integ1(q = ts_IC_R, mean = gknots, sd = sig)%*%theta_old)
        S0ts_IC_diff = S0ts_IC_L - S0ts_IC_R
        S0ts_IC_diff[S0ts_IC_diff<1e-40] = 1e-40 # 20220514
        if (beta_Hessian_option=="full") {
          dh0ts_IC_L = Gaussian_basis_deriv1(x = ts_IC_L, mean = gknots, sd = sig)%*%theta_old
          dh0ts_IC_R = Gaussian_basis_deriv1(x = ts_IC_R, mean = gknots, sd = sig)%*%theta_old
        }
      }
      
      # 0.5. gradient
      if (n_E==0) {quan_E_grad_beta = matrix(0, 1, 1)}
      else {quan_E_grad_beta = (dh0ts_E/h0tss_E - h0ts_E)*ts_E + 1}
      if (n_RC==0) {quan_RC_grad_beta = matrix(0, 1, 1)}
      else {quan_RC_grad_beta = -h0ts_RC*ts_RC}
      if (n_LC==0) {quan_LC_grad_beta = matrix(0, 1, 1)}
      else {quan_LC_grad_beta = S0ts_LC*h0ts_LC*ts_LC / S0ts_LC_diff}
      if (n_IC==0) {quan_IC_grad_beta = matrix(0, 1, 1)}
      else {quan_IC_grad_beta = (-S0ts_IC_L*h0ts_IC_L*ts_IC_L + S0ts_IC_R*h0ts_IC_R*ts_IC_R) / S0ts_IC_diff}
      # quan_grad_beta = h0ts*ts
      # grad_beta = t(-X_E)%*%quan_E_grad_beta - t(-X)%*%quan_grad_beta
      grad_beta = t(-X_E)%*%quan_E_grad_beta +
                  t(-X_RC)%*%quan_RC_grad_beta +
                  t(-X_LC)%*%quan_LC_grad_beta +
                  t(-X_IC)%*%quan_IC_grad_beta
      grad_beta_iter[k, ] = t(grad_beta)
      
      # # 0.5. differential spline for Hessian
      # ddphi_E = mSpline(ts_E, knots=IntKnt, degree=dgrSp, intercept=TRUE, Boundary.knots=bryKnt, derivs=2)
      # dphi_RC = mSpline(ts_RC, knots=IntKnt, degree=dgrSp, intercept=TRUE, Boundary.knots=bryKnt, derivs=1)
      # dphi_LC = mSpline(ts_LC, knots=IntKnt, degree=dgrSp, intercept=TRUE, Boundary.knots=bryKnt, derivs=1)
      # dphi_IC_L = mSpline(ts_IC_L, knots=IntKnt, degree=dgrSp, intercept=TRUE, Boundary.knots=bryKnt, derivs=1)
      # dphi_IC_R = mSpline(ts_IC_R, knots=IntKnt, degree=dgrSp, intercept=TRUE, Boundary.knots=bryKnt, derivs=1)
      # # 0.5. supplement segmentations for Hessian
      # ddh0ts_E = ddphi_E%*%theta_old
      # dh0ts_RC = dphi_RC%*%theta_old
      # dh0ts_LC = dphi_LC%*%theta_old
      # dh0ts_IC_L = dphi_IC_L%*%theta_old
      # dh0ts_IC_R = dphi_IC_R%*%theta_old
      # 0.5. Hessian
      # quan_E_Hes = (ddh0ts_E*(ts_E)^2*h0ts_E-(dh0ts_E*ts_E)^2)/(h0tss_E)^2 - dh0ts_E*(ts_E)^2 - h0ts_E*ts_E
      # quan_RC_Hes = dh0ts_RC*(ts_RC)^2 - h0ts_RC*ts_RC
      # quan_LC_Hes = exp(-ch0ts_LC)*ts_LC*(-h0ts_LC^2*ts_LC + dh0ts_LC*ts_LC + h0ts_LC)/(1-exp(-ch0ts_LC)) - (exp(-ch0ts_LC)*h0ts_LC*ts_LC/(1-exp(-ch0ts_LC)))^2
      # quan_IC_Hes = (exp(-ch0ts_IC_L)*ts_IC_L*(h0ts_IC_L^2*ts_IC_L - dh0ts_IC_L*ts_IC_L - h0ts_IC_L) + exp(-ch0ts_IC_R)*ts_IC_R*(-h0ts_IC_R^2*ts_IC_R+dh0ts_IC_R*ts_IC_R+h0ts_IC_R))/(exp(-ch0ts_IC_L)-exp(-ch0ts_IC_R)) - ((-exp(-ch0ts_IC_L)*h0ts_IC_L*ts_IC_L+exp(-ch0ts_IC_R)*h0ts_IC_R*ts_IC_R)/(exp(-ch0ts_IC_L)-exp(-ch0ts_IC_R)))^2
      # 
      # Hes = t(X_E)%*%diag(x = as.vector(quan_E_Hes), nrow = n_E, ncol = n_E)%*%X_E + t(X_RC)%*%diag(x = as.vector(quan_RC_Hes), nrow = n_RC, ncol = n_RC)%*%X_RC + t(X_LC)%*%diag(x= as.vector(quan_LC_Hes), nrow = n_LC, ncol = n_LC)%*%X_LC + t(X_IC)%*%diag(x = as.vector(quan_IC_Hes), nrow = n_IC, ncol = n_IC)%*%X_IC
      
      if (beta_Hessian_option=="full") {
        # 0.5. full Hessian
        if (n_E==0) {quan_E_Hes_beta = matrix(0, 1, 1)}
        else {quan_E_Hes_beta = (ddh0ts_E / h0tss_E - (dh0ts_E / h0tss_E)^2 - dh0ts_E) * ts_E^2 + (dh0ts_E / h0tss_E - h0ts_E) * ts_E}
        if(n_RC==0) {quan_RC_Hes_beta = matrix(0, 1, 1)}
        else {quan_RC_Hes_beta = - dh0ts_RC * ts_RC^2 - h0ts_RC * ts_RC}
        if(n_LC==0) {quan_LC_Hes_beta = matrix(0, 1, 1)}
        else {quan_LC_Hes_beta = S0ts_LC*(-h0ts_LC^2 + dh0ts_LC)*ts_LC^2 / S0ts_LC_diff - (S0ts_LC*h0ts_LC*ts_LC / S0ts_LC_diff)^2 + S0ts_LC*h0ts_LC*ts_LC / S0ts_LC_diff}
        if(n_IC==0) {quan_IC_Hes_beta = matrix(0, 1, 1)}
        else {
          quan_IC_Hes_beta = S0ts_IC_L*(h0ts_IC_L^2 - dh0ts_IC_L)*ts_IC_L^2 / S0ts_IC_diff - (S0ts_IC_L*h0ts_IC_L*ts_IC_L / S0ts_IC_diff)^2 - S0ts_IC_L*h0ts_IC_L*ts_IC_L / S0ts_IC_diff +
                             S0ts_IC_R*(-h0ts_IC_R^2 + dh0ts_IC_R)*ts_IC_R^2 / S0ts_IC_diff - (S0ts_IC_R*h0ts_IC_R*ts_IC_R / S0ts_IC_diff)^2 + S0ts_IC_R*h0ts_IC_R*ts_IC_R / S0ts_IC_diff +
                             2*(S0ts_IC_R*h0ts_IC_R*ts_IC_R * S0ts_IC_L*h0ts_IC_L*ts_IC_L) / S0ts_IC_diff^2
        }
        # Hes_beta = t(-X_E)%*%(quan_E_Hes_beta%*%matrix(1, 1, p)*(-X_E)) - t(-X)%*%(quan_Hes_beta%*%matrix(1, 1, p)*(-X))
        # Hes_beta = t(-X_E)%*%(as.numeric(quan_E_Hes_beta)*(-X_E)) - t(-X)%*%(as.numeric(quan_Hes_beta)*(-X)) # 20220117 faster alternative of above line, use R recycling rule, may cause trouble
        Hes_beta = t(-X_E)%*%(as.numeric(quan_E_Hes_beta)*(-X_E)) + t(-X_RC)%*%(as.numeric(quan_RC_Hes_beta)*(-X_RC)) + t(-X_LC)%*%(as.numeric(quan_LC_Hes_beta)*(-X_LC)) + t(-X_IC)%*%(as.numeric(quan_IC_Hes_beta)*(-X_IC))
      }
      else if (beta_Hessian_option=="modified") {
        # 0.5. modified Hessian
        if (n_E==0) {quan_E_Hes_beta = matrix(0, 1, 1)}
        else {quan_E_Hes_beta = (dh0ts_E / h0tss_E * ts_E)^2 + h0ts_E * ts_E} # only pick the negative definite terms from the full Hessian and make them positive
        if(n_RC==0) {quan_RC_Hes_beta = matrix(0, 1, 1)}
        else {quan_RC_Hes_beta = h0ts_RC * ts_RC}
        if(n_LC==0) {quan_LC_Hes_beta = matrix(0, 1, 1)}
        else {
          quan_LC_Hes_beta = S0ts_LC*(h0ts_LC*ts_LC)^2 / S0ts_LC_diff + (S0ts_LC*h0ts_LC*ts_LC / S0ts_LC_diff)^2 + h0ts_LC*ts_LC*(S0ts_LC/ S0ts_LC_diff)^2
          # # alternative of above line, slightly faster for X_1, slightly faster for X_3
          # quan_LC_Hes_beta = S0ts_LC*(h0ts_LC*ts_LC)^2 / S0ts_LC_diff + (S0ts_LC*h0ts_LC*ts_LC / S0ts_LC_diff)^2
        }
        
        if(n_IC==0) {quan_IC_Hes_beta = matrix(0, 1, 1)}
        else {
          quan_IC_Hes_beta = S0ts_IC_L*S0ts_IC_R*(h0ts_IC_L*ts_IC_L / S0ts_IC_diff)^2 + (S0ts_IC_L*h0ts_IC_L*ts_IC_L / S0ts_IC_diff)^2 + S0ts_IC_L*h0ts_IC_L*ts_IC_L / S0ts_IC_diff +
                             S0ts_IC_R*(h0ts_IC_R*ts_IC_R)^2 / S0ts_IC_diff + (S0ts_IC_R*h0ts_IC_R*ts_IC_R / S0ts_IC_diff)^2 + h0ts_IC_R*ts_IC_R*(S0ts_IC_R / S0ts_IC_diff)^2
          # # alternative of above line, slower for X_1, slightly faster for X_3
          # quan_IC_Hes_beta = (S0ts_IC_L*h0ts_IC_L*ts_IC_L / S0ts_IC_diff)^2 + S0ts_IC_L*h0ts_IC_L*ts_IC_L / S0ts_IC_diff +
          #                    S0ts_IC_R*(h0ts_IC_R*ts_IC_L)^2 / S0ts_IC_diff + (S0ts_IC_R*h0ts_IC_R*ts_IC_R / S0ts_IC_diff)^2
        }
        # Hes_beta = t(X_E)%*%diag(x = as.vector(quan_E_Hes_beta), nrow = dim(X_E)[1], ncol = dim(X_E)[1])%*%X_E + t(X)%*%diag(x = as.vector(h0ts*ts), nrow = dim(X)[1], ncol = dim(X)[1])%*%X
        # Hes_beta = t(-X_E)%*%(quan_E_Hes_beta%*%matrix(1, 1, p)*(-X_E)) + t(-X)%*%(quan_Hes_beta%*%matrix(1, 1, p)*(-X)) # faster alternative, but cause tiny difference in some elements
        Hes_beta = t(-X_E)%*%(as.numeric(quan_E_Hes_beta)*(-X_E)) +
                   t(-X_RC)%*%(as.numeric(quan_RC_Hes_beta)*(-X_RC)) +
                   t(-X_LC)%*%(as.numeric(quan_LC_Hes_beta)*(-X_LC)) +
                   t(-X_IC)%*%(as.numeric(quan_IC_Hes_beta)*(-X_IC)) # 20220117 faster alternative of above line, use R recycling rule, may cause trouble
      }
      else if (beta_Hessian_option=="modified_1") {
        # 0.5. modified Hessian
        if (n_E==0) {quan_E_Hes_beta = matrix(0, 1, 1)}
        else {quan_E_Hes_beta = (dh0ts_E / h0tss_E * ts_E)^2 + h0ts_E * ts_E} # only pick the negative definite terms from the full Hessian and make them positive
        if(n_RC==0) {quan_RC_Hes_beta = matrix(0, 1, 1)}
        else {quan_RC_Hes_beta = h0ts_RC * ts_RC}
        if(n_LC==0) {quan_LC_Hes_beta = matrix(0, 1, 1)}
        else {
          # quan_LC_Hes_beta = S0ts_LC*(h0ts_LC*ts_LC)^2 / S0ts_LC_diff + (S0ts_LC*h0ts_LC*ts_LC / S0ts_LC_diff)^2 + h0ts_LC*ts_LC*(S0ts_LC/ S0ts_LC_diff)^2
          # alternative of above line, slightly faster for X_1, slightly faster for X_3
          quan_LC_Hes_beta = S0ts_LC*(h0ts_LC*ts_LC)^2 / S0ts_LC_diff + (S0ts_LC*h0ts_LC*ts_LC / S0ts_LC_diff)^2
        }
        
        if(n_IC==0) {quan_IC_Hes_beta = matrix(0, 1, 1)}
        else {
          # quan_IC_Hes_beta = S0ts_IC_L*S0ts_IC_R*(h0ts_IC_L*ts_IC_L / S0ts_IC_diff)^2 + (S0ts_IC_L*h0ts_IC_L*ts_IC_L / S0ts_IC_diff)^2 + S0ts_IC_L*h0ts_IC_L*ts_IC_L / S0ts_IC_diff +
          #                    S0ts_IC_R*(h0ts_IC_R*ts_IC_R)^2 / S0ts_IC_diff + (S0ts_IC_R*h0ts_IC_R*ts_IC_R / S0ts_IC_diff)^2 + h0ts_IC_R*ts_IC_R*(S0ts_IC_R / S0ts_IC_diff)^2
          # alternative of above line, slower for X_1, slightly faster for X_3
          quan_IC_Hes_beta = (S0ts_IC_L*h0ts_IC_L*ts_IC_L / S0ts_IC_diff)^2 + S0ts_IC_L*h0ts_IC_L*ts_IC_L / S0ts_IC_diff +
                             S0ts_IC_R*(h0ts_IC_R*ts_IC_L)^2 / S0ts_IC_diff + (S0ts_IC_R*h0ts_IC_R*ts_IC_R / S0ts_IC_diff)^2
        }
        # Hes_beta = t(X_E)%*%diag(x = as.vector(quan_E_Hes_beta), nrow = dim(X_E)[1], ncol = dim(X_E)[1])%*%X_E + t(X)%*%diag(x = as.vector(h0ts*ts), nrow = dim(X)[1], ncol = dim(X)[1])%*%X
        # Hes_beta = t(-X_E)%*%(quan_E_Hes_beta%*%matrix(1, 1, p)*(-X_E)) + t(-X)%*%(quan_Hes_beta%*%matrix(1, 1, p)*(-X)) # faster alternative, but cause tiny difference in some elements
        Hes_beta = t(-X_E)%*%(as.numeric(quan_E_Hes_beta)*(-X_E)) +
                   t(-X_RC)%*%(as.numeric(quan_RC_Hes_beta)*(-X_RC)) +
                   t(-X_LC)%*%(as.numeric(quan_LC_Hes_beta)*(-X_LC)) +
                   t(-X_IC)%*%(as.numeric(quan_IC_Hes_beta)*(-X_IC)) # 20220117 faster alternative of above line, use R recycling rule, may cause trouble
      }
      else {stop("beta_Hessian_option is either full or modified.\n")}
      
      Hes_beta_iter[k, ] = c(t(Hes_beta)) # convert the Hes_beta into vector by row, and save
      det_Hes_beta_iter[k] = det(Hes_beta)
      eigenvl_Hes_beta_iter[k, ] = eigen(Hes_beta)$values
      Hes_beta = Hes_beta + diag(Hes_addition, p) # 20190812 to avoid singular in solve()
      # eigen(Hes_beta)
      # 0.5. update beta temporarily
      beta_inc = solve(Hes_beta, grad_beta)
      beta_new = beta_old + beta_inc
      
      # 0.5. accelerated time
      if (!missing(Zmat) & !missing(gamma_initial) & (dim_Zmat[1]>n+1)) {
        if (n_E==0) {ts_E = NA}
        else {ts_E = exp(-X_E%*%beta_new) * as.matrix(tapply(exp(-Z_E%*%gamma_old)*interval_length_E, id_Z_E, FUN = sum))}
        if (n_RC==0) {ts_RC = NA}
        else {ts_RC = exp(-X_RC%*%beta_new) * as.matrix(tapply(exp(-Z_RC%*%gamma_old)*interval_length_RC, id_Z_RC, FUN = sum))}
        if (n_LC==0) {ts_LC = NA}
        else {ts_LC = exp(-X_LC%*%beta_new) * as.matrix(tapply(exp(-Z_LC%*%gamma_old)*interval_length_LC, id_Z_LC, FUN = sum))}
        if (n_IC==0) {ts_IC_L = ts_IC_R = NA}
        else {
          ts_IC_L = exp(-X_IC%*%beta_new) * as.matrix(tapply(exp(-Z_IC_subset%*%gamma_old)*interval_length_IC_subset, id_Z_IC_subset, FUN = sum)) # special requirement
          ts_IC_R = exp(-X_IC%*%beta_new) * as.matrix(tapply(exp(-Z_IC%*%gamma_old)*interval_length_IC, id_Z_IC, FUN = sum))
        }
      }
      else {
        if (n_E==0) {ts_E = NA}
        else {ts_E = time_E*exp(-X_E%*%beta_new)} # n_E by 1 matrix
        if (n_RC==0) {ts_RC = NA}
        else {ts_RC = time_RC*exp(-X_RC%*%beta_new)} # n_RC by 1 matrix
        if (n_LC==0) {ts_LC = NA}
        else {ts_LC = time_LC*exp(-X_LC%*%beta_new)} # n_LC by 1 matrix
        if (n_IC==0) {ts_IC_L = ts_IC_R = NA}
        else {
          ts_IC_L = time_IC_L*exp(-X_IC%*%beta_new) # n_IC by 1 matrix
          ts_IC_R = time_IC_R*exp(-X_IC%*%beta_new) # n_IC by 1 matrix
        }
      }
      ts_bind = rbind(ts_E, ts_RC, ts_LC, ts_IC_L, ts_IC_R) # ? not sure if this is suitable for "percentile" knots
      ts_bind = ts_bind[!is.na(ts_bind), , drop = FALSE]
      ts_bind = ts_bind[order(as.numeric(rownames(ts_bind))), , drop = FALSE]
      
      # 0.5. knots
      # knots option 1
      # if (knots_option==1) {
      #   bryKnt = c(min(ts), max(ts)+1e-40)
      #   IntKnt = quantile(sort(ts)[2:(n-1)], seq(5, 95, length.out = numIntKnt)/100, type=1)
      # }
      # # knots option 2
      # if (knots_option==2) {
      #   Knt = quantile(sort(ts), seq(0, 1, length.out = numSp-dgrSp+1), type=1)
      #   bryKnt = c(Knt[1], Knt[numSp-dgrSp+1]+1e-40)
      #   IntKnt = Knt[2:(numSp-dgrSp)]
      # }
      # # knots option 3
      # if (knots_option==3) {
      #   ts_min = min(ts)
      #   ts_max = max(ts)
      #   bryKnt = c(ts_min, ts_max+1e-40)
      #   IntKnt = seq(ts_min, ts_max, length.out = numSp-dgrSp+1)[2:(numSp-dgrSp)]
      # }
      
      # if (all(abs(beta_new-beta_old)<1e-3)) {browser()}
      # Gaussian basis part
      # if (all(abs(beta_new-beta_old)<1e-4) & abs(gamma_new-gamma_old)<1e-4) {browser()}
      if ((all(abs(beta_new-beta_old)>knots_fix_threshold) | all(abs(gamma_new-gamma_old)>knots_fix_threshold)) & k<=knots_fix_iter_no & frequent_knots_update==TRUE) { # 20210630 # 20220208 add frequent_knots_update check
        if (knots_option=="equal_space") {
          bryKnt = range(ts_bind)
          bin_width = (bryKnt[2] - bryKnt[1]) / (numIntKnt+1)
          gknots = seq(bryKnt[1], bryKnt[2], bin_width)
          sig = (2/3) * bin_width
        }
        else if (knots_option=="percentile") {
          # different strategy of calculating knots and sigmas
          gknots = quantile(sort(ts_bind), seq(knots_min_quantile, knots_max_quantile, length.out = numIntKnt+2)/100, type=1)
          dist = sapply(gknots, function(x) {abs(x - ts_bind)})
          sig = apply(dist, 2, function(x) {quantile(x, sig_coverage) / 2})
        }
        else if (knots_option=="percentile_1") { # test on 20220525
          ts_bind_1 = rbind(ts_E, ts_RC, ts_LC, ts_IC_R)
          ts_bind_1 = ts_bind_1[!is.na(ts_bind_1), , drop = FALSE]
          ts_bind_1 = ts_bind_1[order(as.numeric(rownames(ts_bind_1))), , drop = FALSE]
          gknots = quantile(sort(ts_bind_1), seq(knots_min_quantile, knots_max_quantile, length.out = numIntKnt+2)/100, type=1)
          dist = sapply(gknots, function(x) {abs(x - ts_bind_1)})
          sig = apply(dist, 2, function(x) {quantile(x, sig_coverage) / 2})
        }
      }
      
      
      # gauss = pnorm(0, gknots, sig)*sqrt(2*pi*sig^2)
      # gauss = pnorm(bryKnt[1], gknots, sig)*sqrt(2*pi*sig^2) # which one to choose?
      
      # 0.5. spline and cumulative spline
      # if (n_E==0) {phi_E = matrix(0, 1, numSp)}
      # else {
      #   phi_E = mSpline(ts_E, knots=IntKnt, degree=dgrSp, intercept=TRUE, Boundary.knots=bryKnt)
      # }
      # cphi = iSpline(ts, knots=IntKnt, degree=dgrSp, intercept=TRUE, Boundary.knots=bryKnt)
      
      # 0.5. hazard and cumulative hazard
      # if (n_E==0) {h0ts_E = h0tss_E = 1}
      # else {
      #   h0ts_E = phi_E%*%theta_old
      #   h0tss_E = h0ts_E
      #   h0tss_E[h0tss_E<1e-40] = 1e-40
      # }
      # ch0ts = cphi%*%theta_old
      
      # Gaussian basis part
      if (n_E==0) {
        h0ts_E = h0tss_E = 1
        ch0ts_E = 0
      }
      else {
        h0ts_E = Gaussian_basis(x = ts_E, mean = gknots, sd = sig)%*%theta_old
        h0tss_E = h0ts_E
        h0tss_E[h0tss_E<1e-40] = 1e-40
        ch0ts_E = Gaussian_basis_integ1(q = ts_E, mean = gknots, sd = sig)%*%theta_old
      }
      if (n_RC==0) {ch0ts_RC = 0}
      else {ch0ts_RC = Gaussian_basis_integ1(q = ts_RC, mean = gknots, sd = sig)%*%theta_old}
      if (n_LC==0) {
        S0ts_LC = 0
        S0ts_LC_diff = 1
      }
      else {
        S0ts_LC = exp(-Gaussian_basis_integ1(q = ts_LC, mean = gknots, sd = sig)%*%theta_old)
        S0ts_LC_diff = 1 - S0ts_LC
        S0ts_LC_diff[S0ts_LC_diff<1e-40] = 1e-40 # 20220514
      }
      if (n_IC==0) {
        S0ts_IC_L = S0ts_IC_R = 0
        S0ts_IC_diff = 1
      }
      else {
        S0ts_IC_L = exp(-Gaussian_basis_integ1(q = ts_IC_L, mean = gknots, sd = sig)%*%theta_old)
        S0ts_IC_R = exp(-Gaussian_basis_integ1(q = ts_IC_R, mean = gknots, sd = sig)%*%theta_old)
        S0ts_IC_diff = S0ts_IC_L - S0ts_IC_R
        S0ts_IC_diff[S0ts_IC_diff<1e-40] = 1e-40 # 20220514
      }
      
      # 0. log-likelihood
      if (!missing(Zmat) & !missing(gamma_initial) & (dim_Zmat[1]>n+1)) {
        like_beta_new = sum(log(h0tss_E) - X_E%*%beta_new - Z_ni_E%*%gamma_old - ch0ts_E) - sum(ch0ts_RC) + sum(log(S0ts_LC_diff)) + sum(log(S0ts_IC_diff))
      }
      else {like_beta_new = sum(log(h0tss_E) - X_E%*%beta_new - ch0ts_E) - sum(ch0ts_RC) + sum(log(S0ts_LC_diff)) + sum(log(S0ts_IC_diff))}
      
      likelihood_beta_iter_before_ls[k] = like_beta_new
      # 1. backtracking line search for beta
      alpha_N_beta = 1
      iteration_beta = 0L
      # message("for loop No. is ", k)
      while (like_beta_new <= like_beta_old) {
        iteration_beta = iteration_beta + 1L
        # message("like in while loop equals ", like)
        if (alpha_N_beta >= 1e-2) {alpha_N_beta = alpha_N_beta*0.6} # original step size multiplier 0.6, alternative 0.1
        else if (alpha_N_beta < 1e-2 & alpha_N_beta >= 1e-5) {alpha_N_beta = alpha_N_beta*5e-2} # original step size multiplier 5e-2, alternative 1e-3
        else if (alpha_N_beta < 1e-5 & alpha_N_beta >= 1e-20) {alpha_N_beta = alpha_N_beta*1e-5} # original step size multiplier 1e-5, alternative 1e-8
        # else if (alpha_N_beta<1e-20 & alpha_N_beta>=1e-100) {alpha_N_beta=alpha_N_beta*1e-10}
        # if (any(abs(grad_beta)>0.1)) {beta_new = beta_old + 0.01}
        else {break}
        # else if (all(abs(grad_beta)<0.01)) {break}
        
        # 1. update beta with search direction
        beta_new = beta_old + alpha_N_beta * beta_inc # bug of version 0.1.0 fixed
        
        # 1. accelerated time
        if (!missing(Zmat) & !missing(gamma_initial) & (dim_Zmat[1]>n+1)) {
          if (n_E==0) {ts_E = NA}
          else {ts_E = exp(-X_E%*%beta_new) * as.matrix(tapply(exp(-Z_E%*%gamma_old)*interval_length_E, id_Z_E, FUN = sum))}
          if (n_RC==0) {ts_RC = NA}
          else {ts_RC = exp(-X_RC%*%beta_new) * as.matrix(tapply(exp(-Z_RC%*%gamma_old)*interval_length_RC, id_Z_RC, FUN = sum))}
          if (n_LC==0) {ts_LC = NA}
          else {ts_LC = exp(-X_LC%*%beta_new) * as.matrix(tapply(exp(-Z_LC%*%gamma_old)*interval_length_LC, id_Z_LC, FUN = sum))}
          if (n_IC==0) {ts_IC_L = ts_IC_R = NA}
          else {
            ts_IC_L = exp(-X_IC%*%beta_new) * as.matrix(tapply(exp(-Z_IC_subset%*%gamma_old)*interval_length_IC_subset, id_Z_IC_subset, FUN = sum)) # special requirement
            ts_IC_R = exp(-X_IC%*%beta_new) * as.matrix(tapply(exp(-Z_IC%*%gamma_old)*interval_length_IC, id_Z_IC, FUN = sum))
          }
        }
        else {
          if (n_E==0) {ts_E = NA}
          else {ts_E = time_E*exp(-X_E%*%beta_new)} # n_E by 1 matrix
          if (n_RC==0) {ts_RC = NA}
          else {ts_RC = time_RC*exp(-X_RC%*%beta_new)} # n_RC by 1 matrix
          if (n_LC==0) {ts_LC = NA}
          else {ts_LC = time_LC*exp(-X_LC%*%beta_new)} # n_LC by 1 matrix
          if (n_IC==0) {ts_IC_L = ts_IC_R = NA}
          else {
            ts_IC_L = time_IC_L*exp(-X_IC%*%beta_new) # n_IC by 1 matrix
            ts_IC_R = time_IC_R*exp(-X_IC%*%beta_new) # n_IC by 1 matrix
          }
        }
        ts_bind = rbind(ts_E, ts_RC, ts_LC, ts_IC_L, ts_IC_R) # ? not sure if this is suitable for "percentile" knots
        ts_bind = ts_bind[!is.na(ts_bind), , drop = FALSE]
        ts_bind = ts_bind[order(as.numeric(rownames(ts_bind))), , drop = FALSE]
        
        # 1. knots
        # # knots option 1
        # if (knots_option==1) {
        #   bryKnt = c(min(ts), max(ts)+1e-40)
        #   IntKnt = quantile(sort(ts)[2:(n-1)], seq(5, 95, length.out = numIntKnt)/100, type=1)
        # }
        # # knots option 2
        # if (knots_option==2) {
        #   Knt = quantile(sort(ts), seq(0, 1, length.out = numSp-dgrSp+1), type=1)
        #   bryKnt = c(Knt[1], Knt[numSp-dgrSp+1]+1e-40)
        #   IntKnt = Knt[2:(numSp-dgrSp)]
        # }
        # # knots option 3
        # if (knots_option==3) {
        #   ts_min = min(ts)
        #   ts_max = max(ts)
        #   bryKnt = c(ts_min, ts_max+1e-40)
        #   IntKnt = seq(ts_min, ts_max, length.out = numSp-dgrSp+1)[2:(numSp-dgrSp)]
        # }
        
        # Gaussian basis part
        if ((all(abs(beta_new-beta_old)>knots_fix_threshold) | all(abs(gamma_new-gamma_old)>knots_fix_threshold)) & k<=knots_fix_iter_no & frequent_knots_update==TRUE) { # 20210630 # 20220208 add frequent_knots_update check
          if (knots_option=="equal_space") {
            bryKnt = range(ts_bind)
            bin_width = (bryKnt[2] - bryKnt[1]) / (numIntKnt+1)
            gknots = seq(bryKnt[1], bryKnt[2], bin_width)
            sig = (2/3) * bin_width
          }
          else if (knots_option=="percentile") {
            # different strategy of calculating knots and sigmas
            gknots = quantile(sort(ts_bind), seq(knots_min_quantile, knots_max_quantile, length.out = numIntKnt+2)/100, type=1)
            dist = sapply(gknots, function(x) {abs(x - ts_bind)})
            sig = apply(dist, 2, function(x) {quantile(x, sig_coverage) / 2})
          }
          else if (knots_option=="percentile_1") { # test on 20220525
            ts_bind_1 = rbind(ts_E, ts_RC, ts_LC, ts_IC_R)
            ts_bind_1 = ts_bind_1[!is.na(ts_bind_1), , drop = FALSE]
            ts_bind_1 = ts_bind_1[order(as.numeric(rownames(ts_bind_1))), , drop = FALSE]
            gknots = quantile(sort(ts_bind_1), seq(knots_min_quantile, knots_max_quantile, length.out = numIntKnt+2)/100, type=1)
            dist = sapply(gknots, function(x) {abs(x - ts_bind_1)})
            sig = apply(dist, 2, function(x) {quantile(x, sig_coverage) / 2})
          }
        }
        
        # gauss = pnorm(0, gknots, sig)*sqrt(2*pi*sig^2)
        # gauss = pnorm(bryKnt[1], gknots, sig)*sqrt(2*pi*sig^2) # which one to choose?
        
        # 1. spline and cumulative spline
        # if (n_E==0) {phi_E = matrix(0, 1, numSp)}
        # else {
        #   phi_E = mSpline(ts_E, knots=IntKnt, degree=dgrSp, intercept=TRUE, Boundary.knots=bryKnt)
        # }
        # cphi = iSpline(ts, knots=IntKnt, degree=dgrSp, intercept=TRUE, Boundary.knots=bryKnt)
        
        # 1. hazard and cumulative hazard
        # if (n_E==0) {h0ts_E = h0tss_E = 1}
        # else {
        #   h0ts_E = phi_E%*%theta_old
        #   h0tss_E = h0ts_E
        #   h0tss_E[h0tss_E<1e-40] = 1e-40
        # }
        # ch0ts = cphi%*%theta_old
        
        # Gaussian basis part
        if (n_E==0) {
          h0ts_E = h0tss_E = 1
          ch0ts_E = 0
        }
        else {
          h0ts_E = Gaussian_basis(x = ts_E, mean = gknots, sd = sig)%*%theta_old
          h0tss_E = h0ts_E
          h0tss_E[h0tss_E<1e-40] = 1e-40
          ch0ts_E = Gaussian_basis_integ1(q = ts_E, mean = gknots, sd = sig)%*%theta_old
        }
        if (n_RC==0) {ch0ts_RC = 0}
        else {ch0ts_RC = Gaussian_basis_integ1(q = ts_RC, mean = gknots, sd = sig)%*%theta_old}
        if (n_LC==0) {
          S0ts_LC = 0
          S0ts_LC_diff = 1
        }
        else {
          S0ts_LC = exp(-Gaussian_basis_integ1(q = ts_LC, mean = gknots, sd = sig)%*%theta_old)
          S0ts_LC_diff = 1 - S0ts_LC
          S0ts_LC_diff[S0ts_LC_diff<1e-40] = 1e-40 # 20220514
        }
        if (n_IC==0) {
          S0ts_IC_L = S0ts_IC_R = 0
          S0ts_IC_diff = 1
        }
        else {
          S0ts_IC_L = exp(-Gaussian_basis_integ1(q = ts_IC_L, mean = gknots, sd = sig)%*%theta_old)
          S0ts_IC_R = exp(-Gaussian_basis_integ1(q = ts_IC_R, mean = gknots, sd = sig)%*%theta_old)
          S0ts_IC_diff = S0ts_IC_L - S0ts_IC_R
          S0ts_IC_diff[S0ts_IC_diff<1e-40] = 1e-40 # 20220514
        }
        
        # 0. log-likelihood
        if (!missing(Zmat) & !missing(gamma_initial) & (dim_Zmat[1]>n+1)) {
          like_beta_new = sum(log(h0tss_E) - X_E%*%beta_new - Z_ni_E%*%gamma_old - ch0ts_E) - sum(ch0ts_RC) + sum(log(S0ts_LC_diff)) + sum(log(S0ts_IC_diff))
        }
        else {like_beta_new = sum(log(h0tss_E) - X_E%*%beta_new - ch0ts_E) - sum(ch0ts_RC) + sum(log(S0ts_LC_diff)) + sum(log(S0ts_IC_diff))}
      }
      # # one way to escape from local maximum 20210608
      # if (any(abs(grad_beta)>=0.01) & iteration_beta == 7) {beta_new = beta_old - 1e-1*abs(beta_old)*sign(grad_beta)}
      
      ts_range_beta_iter[k, ] = range(ts_bind)
      likelihood_beta_iter_after_ls[k] = like_beta_new
      iter_Newton_beta[k] = iteration_beta
      beta_iter[k, ] = beta_new
      
      ########################## Newton step for gamma ##########################
      if (!missing(Zmat) & !missing(gamma_initial) & (dim_Zmat[1]>n+1)) {
        like_gamma_old = like_beta_new
        # 0.5. differential spline for gradient and modified Hessian
        # if (n_E==0) {dphi_E = matrix(0, 1, numSp)}
        # else {dphi_E = mSpline(ts_E, knots=IntKnt, degree=dgrSp, intercept=TRUE, Boundary.knots=bryKnt, derivs=1)}
        # phi = mSpline(ts, knots=IntKnt, degree=dgrSp, intercept=TRUE, Boundary.knots=bryKnt)
        
        # 0.5. segmentations for gradient and modified Hessian
        # dh0ts_E = dphi_E%*%theta_old
        # h0ts_E = phi_E%*%theta_old
        # h0tss_E = h0ts_E
        # h0tss_E[h0tss_E<1e-40] = 1e-40 # to avoid 0 denominator
        # h0ts = phi%*%theta_old
        
        # Gaussian basis part
        # dh0ts = Gaussian_basis_deriv1(x = ts, mean = gknots, sd = sig)%*%theta_old
        # h0ts = Gaussian_basis(x = ts, mean = gknots, sd = sig)%*%theta_old
        # dh0ts = matrix(0, n, 1)
        # h0ts = matrix(0, n, 1)
        # for (i in 1:n) { # will improve using apply function
        #   dh0ts[i] = t(dnorm_deriv1(ts[i], gknots, sig)*sqrt(2*pi*sig^2))%*%theta_old
        #   h0ts[i] = t(dnorm(ts[i], gknots, sig)*sqrt(2*pi*sig^2))%*%theta_old
        # }
        
        if (n_E==0) {
          dh0ts_E = 0
          h0ts_E = h0tss_E = 1
          if (gamma_Hessian_option=="full") {ddh0ts_E = 0}
        }
        else {
          dh0ts_E = Gaussian_basis_deriv1(x = ts_E, mean = gknots, sd = sig)%*%theta_old
          h0ts_E = Gaussian_basis(x = ts_E, mean = gknots, sd = sig)%*%theta_old
          h0tss_E = h0ts_E
          h0tss_E[h0tss_E<1e-40] = 1e-40
          if (gamma_Hessian_option=="full") {ddh0ts_E = Gaussian_basis_deriv2(x = ts_E, mean = gknots, sd = sig)%*%theta_old} # for full Hessian of gamma}
        }
        if (n_RC==0) {
          h0ts_RC = 0
          if (gamma_Hessian_option=="full") {dh0ts_RC = 0}
        }
        else {
          h0ts_RC = Gaussian_basis(x = ts_RC, mean = gknots, sd = sig)%*%theta_old
          if (gamma_Hessian_option=="full") {dh0ts_RC = Gaussian_basis_deriv1(x = ts_RC, mean = gknots, sd = sig)%*%theta_old}
        }
        if (n_LC==0) {
          h0ts_LC = S0ts_LC = 0
          S0ts_LC_diff = 1
          if (gamma_Hessian_option=="full") {dh0ts_LC = 0}
        }
        else {
          h0ts_LC = Gaussian_basis(x = ts_LC, mean = gknots, sd = sig)%*%theta_old
          S0ts_LC = exp(-Gaussian_basis_integ1(q = ts_LC, mean = gknots, sd = sig)%*%theta_old)
          S0ts_LC_diff = 1 - S0ts_LC
          S0ts_LC_diff[S0ts_LC_diff<1e-40] = 1e-40 # 20220514
          if (gamma_Hessian_option=="full") {dh0ts_LC = Gaussian_basis_deriv1(x = ts_LC, mean = gknots, sd = sig)%*%theta_old}
        }
        if (n_IC==0) {
          h0ts_IC_L = S0ts_IC_L = h0ts_IC_R = S0ts_IC_R = 0
          S0ts_IC_diff = 1
          if (gamma_Hessian_option=="full") {dh0ts_IC_L = dh0ts_IC_R = 0}
        }
        else {
          h0ts_IC_L = Gaussian_basis(x = ts_IC_L, mean = gknots, sd = sig)%*%theta_old
          S0ts_IC_L = exp(-Gaussian_basis_integ1(q = ts_IC_L, mean = gknots, sd = sig)%*%theta_old)
          h0ts_IC_R = Gaussian_basis(x = ts_IC_R, mean = gknots, sd = sig)%*%theta_old
          S0ts_IC_R = exp(-Gaussian_basis_integ1(q = ts_IC_R, mean = gknots, sd = sig)%*%theta_old)
          S0ts_IC_diff = S0ts_IC_L - S0ts_IC_R
          S0ts_IC_diff[S0ts_IC_diff<1e-40] = 1e-40 # 20220514
          if (gamma_Hessian_option=="full") {
            dh0ts_IC_L = Gaussian_basis_deriv1(x = ts_IC_L, mean = gknots, sd = sig)%*%theta_old
            dh0ts_IC_R = Gaussian_basis_deriv1(x = ts_IC_R, mean = gknots, sd = sig)%*%theta_old
          }
        }
        
        # 0.5. first derivative of t_star w.r.t. gamma
        if (q==1) {
          # dts_gamma = (exp(-X%*%beta_new)%*%matrix(1, 1, q)) * as.matrix(tapply((-Z)*((exp(-Z%*%gamma_old)*interval_length)%*%matrix(1, 1, q)), id_Z, FUN = sum)) # n by 1
          # dts_gamma = as.numeric(exp(-X%*%beta_new)) * as.matrix(tapply((-Z)*as.numeric(exp(-Z%*%gamma_old)*interval_length), id_Z, FUN = sum)) # 20220119 faster alternative of above line, use R recycling rule, may cause trouble
          # dts_gamma = as.numeric(exp(-X%*%beta_new)) * apply(as.matrix((-Z)*as.numeric(exp(-Z%*%gamma_old)*interval_length)), 2, function(x) tapply(x, id_Z, sum)) # 20220119 slower alternative of above line, but works for q > 0
          if (n_E==0) {dts_E_gamma = matrix(0, 1, q)}
          else {dts_E_gamma = as.numeric(exp(-X_E%*%beta_new)) * as.matrix(tapply((-Z_E)*as.numeric(exp(-Z_E%*%gamma_old)*interval_length_E), id_Z_E, FUN = sum))} # n_E by q
          if (n_RC==0) {dts_RC_gamma = matrix(0, 1, q)}
          else {dts_RC_gamma = as.numeric(exp(-X_RC%*%beta_new)) * as.matrix(tapply((-Z_RC)*as.numeric(exp(-Z_RC%*%gamma_old)*interval_length_RC), id_Z_RC, FUN = sum))}
          if (n_LC==0) {dts_LC_gamma = matrix(0, 1, q)}
          else {dts_LC_gamma = as.numeric(exp(-X_LC%*%beta_new)) * as.matrix(tapply((-Z_LC)*as.numeric(exp(-Z_LC%*%gamma_old)*interval_length_LC), id_Z_LC, FUN = sum))}
          if (n_IC==0) {dts_IC_L_gamma = dts_IC_R_gamma = matrix(0, 1, q)}
          else {
            dts_IC_L_gamma = as.numeric(exp(-X_IC%*%beta_new)) * as.matrix(tapply((-Z_IC_subset)*as.numeric(exp(-Z_IC_subset%*%gamma_old)*interval_length_IC_subset), id_Z_IC_subset, FUN = sum))
            dts_IC_R_gamma = as.numeric(exp(-X_IC%*%beta_new)) * as.matrix(tapply((-Z_IC)*as.numeric(exp(-Z_IC%*%gamma_old)*interval_length_IC), id_Z_IC, FUN = sum))
          }
        }
        else {
          # dts_gamma = (exp(-X%*%beta_new)%*%matrix(1, 1, q)) * as.matrix(tapply((-Z)*((exp(-Z%*%gamma_old)*interval_length)%*%matrix(1, 1, q)), id_Z, FUN = colSums)) # n by q
          # dts_gamma = as.numeric(exp(-X%*%beta_new)) * apply((-Z)*as.numeric(exp(-Z%*%gamma_old)*interval_length), 2, function(x) tapply(x, id_Z, sum))
          if (n_E==0) {dts_E_gamma = matrix(0, 1, q)}
          else {dts_E_gamma = as.numeric(exp(-X_E%*%beta_new)) * apply((-Z_E)*as.numeric(exp(-Z_E%*%gamma_old)*interval_length_E), 2, function(x) tapply(x, id_Z_E, sum))}
          if (n_RC==0) {dts_RC_gamma = matrix(0, 1, q)}
          else {dts_RC_gamma = as.numeric(exp(-X_RC%*%beta_new)) * apply((-Z_RC)*as.numeric(exp(-Z_RC%*%gamma_old)*interval_length_RC), 2, function(x) tapply(x, id_Z_RC, sum))}
          if (n_LC==0) {dts_LC_gamma = matrix(0, 1, q)}
          else {dts_LC_gamma = as.numeric(exp(-X_LC%*%beta_new)) * apply((-Z_LC)*as.numeric(exp(-Z_LC%*%gamma_old)*interval_length_LC), 2, function(x) tapply(x, id_Z_LC, sum))}
          if (n_IC==0) {dts_IC_L_gamma = dts_IC_R_gamma = matrix(0, 1, q)}
          else {
            dts_IC_L_gamma = as.numeric(exp(-X_IC%*%beta_new)) * apply((-Z_IC_subset)*as.numeric(exp(-Z_IC_subset%*%gamma_old)*interval_length_IC_subset), 2, function(x) tapply(x, id_Z_IC_subset, sum))
            dts_IC_R_gamma = as.numeric(exp(-X_IC%*%beta_new)) * apply((-Z_IC)*as.numeric(exp(-Z_IC%*%gamma_old)*interval_length_IC), 2, function(x) tapply(x, id_Z_IC, sum))
          }
        }
        
        # 0.5. gradient
        if (n_E==0) {quan_E_grad_gamma = matrix(0, 1, 1)}
        else {quan_E_grad_gamma = dh0ts_E/h0tss_E - h0ts_E} # n_E by q
        if(n_RC==0) {quan_RC_grad_gamma = matrix(0, 1, 1)}
        else {quan_RC_grad_gamma = -h0ts_RC}
        if (n_LC==0) {quan_LC_grad_gamma = matrix(0, 1, 1)}
        else {quan_LC_grad_gamma = S0ts_LC*h0ts_LC / S0ts_LC_diff}
        if (n_IC==0) {quan_IC_L_grad_gamma = quan_IC_R_grad_gamma = matrix(0, 1, 1)}
        else {
          quan_IC_L_grad_gamma = -S0ts_IC_L*h0ts_IC_L / S0ts_IC_diff
          quan_IC_R_grad_gamma = S0ts_IC_R*h0ts_IC_R / S0ts_IC_diff
        }
        
        # grad_gamma = matrix(colSums(quan_E_grad_gamma*dts_E_gamma - Z_ni_E) - colSums(h0ts%*%matrix(1, 1, q)*dts_gamma)) # q by 1
        # grad_gamma = matrix(t(dts_E_gamma)%*%quan_E_grad_gamma - colSums(Z_ni_E) - t(dts_gamma)%*%h0ts) # q by 1
        # ? negative grad_gamma
        grad_gamma = t(dts_E_gamma)%*%quan_E_grad_gamma - colSums(Z_ni_E) +
                     t(dts_RC_gamma)%*%quan_RC_grad_gamma + 
                     t(dts_LC_gamma)%*%quan_LC_grad_gamma + 
                     t(dts_IC_L_gamma)%*%quan_IC_L_grad_gamma + t(dts_IC_R_gamma)%*%quan_IC_R_grad_gamma # q by 1
        
        grad_gamma_iter[k, ] = t(grad_gamma)
        
        # # 0.5. differential spline for Hessian
        # ddphi_E = mSpline(ts_E, knots=IntKnt, degree=dgrSp, intercept=TRUE, Boundary.knots=bryKnt, derivs=2)
        # dphi_RC = mSpline(ts_RC, knots=IntKnt, degree=dgrSp, intercept=TRUE, Boundary.knots=bryKnt, derivs=1)
        # dphi_LC = mSpline(ts_LC, knots=IntKnt, degree=dgrSp, intercept=TRUE, Boundary.knots=bryKnt, derivs=1)
        # dphi_IC_L = mSpline(ts_IC_L, knots=IntKnt, degree=dgrSp, intercept=TRUE, Boundary.knots=bryKnt, derivs=1)
        # dphi_IC_R = mSpline(ts_IC_R, knots=IntKnt, degree=dgrSp, intercept=TRUE, Boundary.knots=bryKnt, derivs=1)
        # # 0.5. supplement segmentations for Hessian
        # ddh0ts_E = ddphi_E%*%theta_old
        # dh0ts_RC = dphi_RC%*%theta_old
        # dh0ts_LC = dphi_LC%*%theta_old
        # dh0ts_IC_L = dphi_IC_L%*%theta_old
        # dh0ts_IC_R = dphi_IC_R%*%theta_old
        # 0.5. Hessian
        # quan_E_Hes = (ddh0ts_E*(ts_E)^2*h0ts_E-(dh0ts_E*ts_E)^2)/(h0tss_E)^2 - dh0ts_E*(ts_E)^2 - h0ts_E*ts_E
        # quan_RC_Hes = dh0ts_RC*(ts_RC)^2 - h0ts_RC*ts_RC
        # quan_LC_Hes = exp(-ch0ts_LC)*ts_LC*(-h0ts_LC^2*ts_LC + dh0ts_LC*ts_LC + h0ts_LC)/(1-exp(-ch0ts_LC)) - (exp(-ch0ts_LC)*h0ts_LC*ts_LC/(1-exp(-ch0ts_LC)))^2
        # quan_IC_Hes = (exp(-ch0ts_IC_L)*ts_IC_L*(h0ts_IC_L^2*ts_IC_L - dh0ts_IC_L*ts_IC_L - h0ts_IC_L) + exp(-ch0ts_IC_R)*ts_IC_R*(-h0ts_IC_R^2*ts_IC_R+dh0ts_IC_R*ts_IC_R+h0ts_IC_R))/(exp(-ch0ts_IC_L)-exp(-ch0ts_IC_R)) - ((-exp(-ch0ts_IC_L)*h0ts_IC_L*ts_IC_L+exp(-ch0ts_IC_R)*h0ts_IC_R*ts_IC_R)/(exp(-ch0ts_IC_L)-exp(-ch0ts_IC_R)))^2
        # 
        # Hes = t(X_E)%*%diag(x = as.vector(quan_E_Hes), nrow = n_E, ncol = n_E)%*%X_E + t(X_RC)%*%diag(x = as.vector(quan_RC_Hes), nrow = n_RC, ncol = n_RC)%*%X_RC + t(X_LC)%*%diag(x= as.vector(quan_LC_Hes), nrow = n_LC, ncol = n_LC)%*%X_LC + t(X_IC)%*%diag(x = as.vector(quan_IC_Hes), nrow = n_IC, ncol = n_IC)%*%X_IC
        
        if (gamma_Hessian_option=="full") {
          # 0.5. full Hessian
          if (n_E==0) {quan_E_Hes_gamma_1 = quan_E_Hes_gamma_2 = matrix(0, 1, 1)}
          else {
            quan_E_Hes_gamma_1 = ddh0ts_E / h0tss_E - (dh0ts_E / h0tss_E)^2 - dh0ts_E # n_E by 1
            quan_E_Hes_gamma_2 = rep((dh0ts_E / h0tss_E - h0ts_E) * exp(-X_E%*%beta_new), ni_E)*(exp(-Z_E%*%gamma_old)*interval_length_E)
          }
          if(n_RC==0) {quan_RC_Hes_gamma_1 = quan_RC_Hes_gamma_2 = matrix(0, 1, 1)}
          else {
            quan_RC_Hes_gamma_1 = - dh0ts_RC
            quan_RC_Hes_gamma_2 = rep(- h0ts_RC * exp(-X_RC%*%beta_new), ni_RC)*(exp(-Z_RC%*%gamma_old)*interval_length_RC)
          }
          if(n_LC==0) {quan_LC_Hes_gamma_1 = quan_LC_Hes_gamma_2 = matrix(0, 1, 1)}
          else {
            quan_LC_Hes_gamma_1 = S0ts_LC*(-h0ts_LC^2 + dh0ts_LC) / S0ts_LC_diff - (S0ts_LC*h0ts_LC / S0ts_LC_diff)^2
            quan_LC_Hes_gamma_2 = rep(S0ts_LC*h0ts_LC / S0ts_LC_diff * exp(-X_LC%*%beta_new), ni_LC)*(exp(-Z_LC%*%gamma_old)*interval_length_LC)
          }
          if(n_IC==0) {quan_IC_L_Hes_gamma_1 = quan_IC_L_Hes_gamma_2 = quan_IC_L_Hes_gamma_3 = quan_IC_R_Hes_gamma_1 = quan_IC_R_Hes_gamma_2 = quan_IC_R_Hes_gamma_3 = matrix(0, 1, 1)}
          else {
            quan_IC_L_Hes_gamma_1 = S0ts_IC_L*(h0ts_IC_L^2 - dh0ts_IC_L) / S0ts_IC_diff - (S0ts_IC_L*h0ts_IC_L / S0ts_IC_diff)^2
            quan_IC_L_Hes_gamma_2 = rep(- S0ts_IC_L*h0ts_IC_L / S0ts_IC_diff * exp(-X_IC%*%beta_new), ni_IC_subset)*(exp(-Z_IC_subset%*%gamma_old)*interval_length_IC_subset)
            quan_IC_L_Hes_gamma_3 = (S0ts_IC_L*h0ts_IC_L * S0ts_IC_R*h0ts_IC_R) / S0ts_IC_diff^2
            quan_IC_R_Hes_gamma_1 = S0ts_IC_R*(-h0ts_IC_R^2 + dh0ts_IC_R) / S0ts_IC_diff - (S0ts_IC_R*h0ts_IC_R / S0ts_IC_diff)^2
            quan_IC_R_Hes_gamma_2 = rep(S0ts_IC_R*h0ts_IC_R / S0ts_IC_diff * exp(-X_IC%*%beta_new), ni_IC)*(exp(-Z_IC%*%gamma_old)*interval_length_IC)
            quan_IC_R_Hes_gamma_3 = (S0ts_IC_R*h0ts_IC_R * S0ts_IC_L*h0ts_IC_L) / S0ts_IC_diff^2
          }
          
          # Hes_gamma = t(dts_E_gamma)%*%diag(x = as.vector(quan_E_Hes_gamma), nrow = dim(X_E)[1], ncol = dim(X_E)[1])%*%dts_E_gamma + t(Z)%*%diag(x = as.vector(quan_Hes_gamma), nrow = dim(Z)[1], ncol = dim(Z)[1])%*%(Z)
          # Hes_gamma = t(dts_E_gamma)%*%(quan_E_Hes_gamma_1%*%matrix(1, 1, q)*dts_E_gamma) +
          #             t(Z_E)%*%(quan_E_Hes_gamma_2%*%matrix(1, 1, q)*Z_E) -
          #             t(dts_gamma)%*%(quan_Hes_gamma_1%*%matrix(1, 1, q)*dts_gamma) -
          #             t(Z)%*%(quan_Hes_gamma_2%*%matrix(1, 1, q)*Z) # faster alternative
          # Hes_gamma = t(dts_E_gamma)%*%(as.numeric(quan_E_Hes_gamma_1)*dts_E_gamma) +
          #   t(Z_E)%*%(as.numeric(quan_E_Hes_gamma_2)*Z_E) -
          #   t(dts_gamma)%*%(as.numeric(quan_Hes_gamma_1)*dts_gamma) -
          #   t(Z)%*%(as.numeric(quan_Hes_gamma_2)*Z) # 20220117 faster alternative of above line, use R recycling rule, may cause trouble
          Hes_gamma = t(dts_E_gamma)%*%(as.numeric(quan_E_Hes_gamma_1)*dts_E_gamma) +
                      t(-Z_E)%*%(as.numeric(quan_E_Hes_gamma_2)*(-Z_E)) +
                      t(dts_RC_gamma)%*%(as.numeric(quan_RC_Hes_gamma_1)*dts_RC_gamma) +
                      t(-Z_RC)%*%(as.numeric(quan_RC_Hes_gamma_2)*(-Z_RC)) +
                      t(dts_LC_gamma)%*%(as.numeric(quan_LC_Hes_gamma_1)*dts_LC_gamma) +
                      t(-Z_LC)%*%(as.numeric(quan_LC_Hes_gamma_2)*(-Z_LC)) +
                      t(dts_IC_L_gamma)%*%(as.numeric(quan_IC_L_Hes_gamma_1)*dts_IC_L_gamma) +
                      t(-Z_IC_subset)%*%(as.numeric(quan_IC_L_Hes_gamma_2)*(-Z_IC_subset)) +
                      t(dts_IC_L_gamma)%*%(as.numeric(quan_IC_L_Hes_gamma_3)*dts_IC_R_gamma) +
                      t(dts_IC_R_gamma)%*%(as.numeric(quan_IC_R_Hes_gamma_1)*dts_IC_R_gamma) +
                      t(-Z_IC)%*%(as.numeric(quan_IC_R_Hes_gamma_2)*(-Z_IC)) +
                      t(dts_IC_R_gamma)%*%(as.numeric(quan_IC_R_Hes_gamma_3)*dts_IC_L_gamma)
        }
        else if (gamma_Hessian_option=="modified") {
          # 0.5. modified Hessian
          if (n_E==0) {quan_E_Hes_gamma_1 = quan_E_Hes_gamma_2 = matrix(0, 1, 1)}
          else {
            quan_E_Hes_gamma_1 = (dh0ts_E / h0tss_E)^2 # n_E by 1
            quan_E_Hes_gamma_2 = rep(h0ts_E * exp(-X_E%*%beta_new), ni_E)*(exp(-Z_E%*%gamma_old)*interval_length_E)
          }
          if(n_RC==0) {quan_RC_Hes_gamma_2 = matrix(0, 1, 1)}
          else {
            quan_RC_Hes_gamma_2 = rep(h0ts_RC * exp(-X_RC%*%beta_new), ni_RC)*(exp(-Z_RC%*%gamma_old)*interval_length_RC)
          }
          if(n_LC==0) {quan_LC_Hes_gamma_1 = quan_LC_Hes_gamma_2 = matrix(0, 1, 1)}
          else {
            quan_LC_Hes_gamma_1 = S0ts_LC*h0ts_LC^2 / S0ts_LC_diff + (S0ts_LC*h0ts_LC / S0ts_LC_diff)^2
            quan_LC_Hes_gamma_2 = rep(h0ts_LC*(S0ts_LC/ S0ts_LC_diff)^2 * exp(-X_LC%*%beta_new), ni_LC)*(exp(-Z_LC%*%gamma_old)*interval_length_LC)
          }
          if(n_IC==0) {quan_IC_L_Hes_gamma_1 = quan_IC_L_Hes_gamma_2 = quan_IC_R_Hes_gamma_1 = quan_IC_R_Hes_gamma_2 = matrix(0, 1, 1)}
          else {
            quan_IC_L_Hes_gamma_1 = S0ts_IC_L*S0ts_IC_R*(h0ts_IC_L / S0ts_IC_diff)^2 + (S0ts_IC_L*h0ts_IC_L / S0ts_IC_diff)^2
            quan_IC_L_Hes_gamma_2 = rep(S0ts_IC_L*h0ts_IC_L / S0ts_IC_diff * exp(-X_IC%*%beta_new), ni_IC_subset)*(exp(-Z_IC_subset%*%gamma_old)*interval_length_IC_subset)
            quan_IC_R_Hes_gamma_1 = S0ts_IC_R*h0ts_IC_R^2 / S0ts_IC_diff + (S0ts_IC_R*h0ts_IC_R / S0ts_IC_diff)^2
            quan_IC_R_Hes_gamma_2 = rep(h0ts_IC_R*(S0ts_IC_R / S0ts_IC_diff)^2 * exp(-X_IC%*%beta_new), ni_IC)*(exp(-Z_IC%*%gamma_old)*interval_length_IC)
          }
          
          Hes_gamma = t(dts_E_gamma)%*%(as.numeric(quan_E_Hes_gamma_1)*dts_E_gamma) +
                      t(-Z_E)%*%(as.numeric(quan_E_Hes_gamma_2)*(-Z_E)) +
                      t(-Z_RC)%*%(as.numeric(quan_RC_Hes_gamma_2)*(-Z_RC)) +
                      t(dts_LC_gamma)%*%(as.numeric(quan_LC_Hes_gamma_1)*dts_LC_gamma) +
                      t(-Z_LC)%*%(as.numeric(quan_LC_Hes_gamma_2)*(-Z_LC)) +
                      t(dts_IC_L_gamma)%*%(as.numeric(quan_IC_L_Hes_gamma_1)*dts_IC_L_gamma) +
                      t(-Z_IC_subset)%*%(as.numeric(quan_IC_L_Hes_gamma_2)*(-Z_IC_subset)) +
                      t(dts_IC_R_gamma)%*%(as.numeric(quan_IC_R_Hes_gamma_1)*dts_IC_R_gamma) +
                      t(-Z_IC)%*%(as.numeric(quan_IC_R_Hes_gamma_2)*(-Z_IC)) # 20220117 faster alternative of above line, use R recycling rule, may cause trouble
        }
        else if (gamma_Hessian_option=="modified_1") { # 20220509
          # 0.5. modified Hessian
          if (n_E==0) {quan_E_Hes_gamma_1 = quan_E_Hes_gamma_2 = matrix(0, 1, 1)}
          else {
            quan_E_Hes_gamma_1 = (dh0ts_E / h0tss_E)^2 # n_E by 1
            quan_E_Hes_gamma_2 = rep(h0ts_E * exp(-X_E%*%beta_new), ni_E)*(exp(-Z_E%*%gamma_old)*interval_length_E)
          }
          if(n_RC==0) {quan_RC_Hes_gamma_2 = matrix(0, 1, 1)}
          else {
            quan_RC_Hes_gamma_2 = rep(h0ts_RC * exp(-X_RC%*%beta_new), ni_RC)*(exp(-Z_RC%*%gamma_old)*interval_length_RC)
          }
          if(n_LC==0) {quan_LC_Hes_gamma_1 = quan_LC_Hes_gamma_2 = matrix(0, 1, 1)}
          else {
            quan_LC_Hes_gamma_1 = S0ts_LC*h0ts_LC^2 / S0ts_LC_diff + (S0ts_LC*h0ts_LC / S0ts_LC_diff)^2
          }
          if(n_IC==0) {quan_IC_L_Hes_gamma_1 = quan_IC_L_Hes_gamma_2 = quan_IC_R_Hes_gamma_1 = quan_IC_R_Hes_gamma_2 = matrix(0, 1, 1)}
          else {
            quan_IC_L_Hes_gamma_1 = (S0ts_IC_L*h0ts_IC_L / S0ts_IC_diff)^2
            quan_IC_L_Hes_gamma_2 = rep(S0ts_IC_L*h0ts_IC_L / S0ts_IC_diff * exp(-X_IC%*%beta_new), ni_IC_subset)*(exp(-Z_IC_subset%*%gamma_old)*interval_length_IC_subset)
            quan_IC_R_Hes_gamma_1 = S0ts_IC_R*h0ts_IC_R^2 / S0ts_IC_diff + (S0ts_IC_R*h0ts_IC_R / S0ts_IC_diff)^2
          }
          
          Hes_gamma = t(dts_E_gamma)%*%(as.numeric(quan_E_Hes_gamma_1)*dts_E_gamma) +
                      t(-Z_E)%*%(as.numeric(quan_E_Hes_gamma_2)*(-Z_E)) +
                      t(-Z_RC)%*%(as.numeric(quan_RC_Hes_gamma_2)*(-Z_RC)) +
                      t(dts_LC_gamma)%*%(as.numeric(quan_LC_Hes_gamma_1)*dts_LC_gamma) +
                      t(dts_IC_L_gamma)%*%(as.numeric(quan_IC_L_Hes_gamma_1)*dts_IC_L_gamma) +
                      t(-Z_IC_subset)%*%(as.numeric(quan_IC_L_Hes_gamma_2)*(-Z_IC_subset)) +
                      t(dts_IC_R_gamma)%*%(as.numeric(quan_IC_R_Hes_gamma_1)*dts_IC_R_gamma)
        }
        else {stop("gamma_Hessian_option is either full or modified.\n")}
        
        Hes_gamma_iter[k, ] = c(t(Hes_gamma)) # convert the Hes_gamma into vector by row, and save
        det_Hes_gamma_iter[k] = det(Hes_gamma)
        eigenvl_Hes_gamma_iter[k, ] = eigen(Hes_gamma)$values
        Hes_gamma = Hes_gamma + diag(Hes_addition, q) # to avoid singular in solve()
        # eigen(Hes_gamma)
        
        # 0.5. update gamma temporarily
        gamma_inc = solve(Hes_gamma, grad_gamma)
        gamma_new = gamma_old + gamma_inc
        
        # 0.5. accelerated time
        if (n_E==0) {ts_E = NA}
        else {ts_E = exp(-X_E%*%beta_new) * as.matrix(tapply(exp(-Z_E%*%gamma_new)*interval_length_E, id_Z_E, FUN = sum))}
        if (n_RC==0) {ts_RC = NA}
        else {ts_RC = exp(-X_RC%*%beta_new) * as.matrix(tapply(exp(-Z_RC%*%gamma_new)*interval_length_RC, id_Z_RC, FUN = sum))}
        if (n_LC==0) {ts_LC = NA}
        else {ts_LC = exp(-X_LC%*%beta_new) * as.matrix(tapply(exp(-Z_LC%*%gamma_new)*interval_length_LC, id_Z_LC, FUN = sum))}
        if (n_IC==0) {ts_IC_L = ts_IC_R = NA}
        else {
          ts_IC_L = exp(-X_IC%*%beta_new) * as.matrix(tapply(exp(-Z_IC_subset%*%gamma_new)*interval_length_IC_subset, id_Z_IC_subset, FUN = sum)) # special requirement
          ts_IC_R = exp(-X_IC%*%beta_new) * as.matrix(tapply(exp(-Z_IC%*%gamma_new)*interval_length_IC, id_Z_IC, FUN = sum))
        }
        ts_bind = rbind(ts_E, ts_RC, ts_LC, ts_IC_L, ts_IC_R) # ? not sure if this is suitable for "percentile" knots
        ts_bind = ts_bind[!is.na(ts_bind), , drop = FALSE]
        ts_bind = ts_bind[order(as.numeric(rownames(ts_bind))), , drop = FALSE]
        
        # 0.5. knots
        # # knots option 1
        # if (knots_option==1) {
        #   bryKnt = c(min(ts), max(ts)+1e-40)
        #   IntKnt = quantile(sort(ts)[2:(n-1)], seq(5, 95, length.out = numIntKnt)/100, type=1)
        # }
        # # knots option 2
        # if (knots_option==2) {
        #   Knt = quantile(sort(ts), seq(0, 1, length.out = numSp-dgrSp+1), type=1)
        #   bryKnt = c(Knt[1], Knt[numSp-dgrSp+1]+1e-40)
        #   IntKnt = Knt[2:(numSp-dgrSp)]
        # }
        # # knots option 3
        # if (knots_option==3) {
        #   ts_min = min(ts)
        #   ts_max = max(ts)
        #   bryKnt = c(ts_min, ts_max+1e-40)
        #   IntKnt = seq(ts_min, ts_max, length.out = numSp-dgrSp+1)[2:(numSp-dgrSp)]
        # }
        
        # Gaussian basis part
        if ((all(abs(beta_new-beta_old)>knots_fix_threshold) | all(abs(gamma_new-gamma_old)>knots_fix_threshold)) & k<=knots_fix_iter_no & frequent_knots_update==TRUE) { # 20210630 # 20220208 add frequent_knots_update check
          if (knots_option=="equal_space") {
            bryKnt = range(ts_bind)
            bin_width = (bryKnt[2] - bryKnt[1]) / (numIntKnt+1)
            gknots = seq(bryKnt[1], bryKnt[2], bin_width)
            sig = (2/3) * bin_width
          }
          else if (knots_option=="percentile") {
            # different strategy of calculating knots and sigmas
            gknots = quantile(sort(ts_bind), seq(knots_min_quantile, knots_max_quantile, length.out = numIntKnt+2)/100, type=1)
            dist = sapply(gknots, function(x) {abs(x - ts_bind)})
            sig = apply(dist, 2, function(x) {quantile(x, sig_coverage) / 2})
          }
          else if (knots_option=="percentile_1") { # test on 20220525
            ts_bind_1 = rbind(ts_E, ts_RC, ts_LC, ts_IC_R)
            ts_bind_1 = ts_bind_1[!is.na(ts_bind_1), , drop = FALSE]
            ts_bind_1 = ts_bind_1[order(as.numeric(rownames(ts_bind_1))), , drop = FALSE]
            gknots = quantile(sort(ts_bind_1), seq(knots_min_quantile, knots_max_quantile, length.out = numIntKnt+2)/100, type=1)
            dist = sapply(gknots, function(x) {abs(x - ts_bind_1)})
            sig = apply(dist, 2, function(x) {quantile(x, sig_coverage) / 2})
          }
        }
        
        # gauss = pnorm(0, gknots, sig)*sqrt(2*pi*sig^2)
        # gauss = pnorm(bryKnt[1], gknots, sig)*sqrt(2*pi*sig^2) # which one to choose?
        
        # 0.5. spline and cumulative spline
        # if (n_E==0) {phi_E = matrix(0, 1, numSp)}
        # else {
        #   phi_E = mSpline(ts_E, knots=IntKnt, degree=dgrSp, intercept=TRUE, Boundary.knots=bryKnt)
        # }
        # cphi = iSpline(ts, knots=IntKnt, degree=dgrSp, intercept=TRUE, Boundary.knots=bryKnt)
        
        # 0.5. hazard and cumulative hazard
        # if (n_E==0) {h0ts_E = h0tss_E = 1}
        # else {
        #   h0ts_E = phi_E%*%theta_old
        #   h0tss_E = h0ts_E
        #   h0tss_E[h0tss_E<1e-40] = 1e-40
        # }
        # ch0ts = cphi%*%theta_old
        
        # Gaussian basis part
        if (n_E==0) {
          h0ts_E = h0tss_E = 1
          ch0ts_E = 0
        }
        else {
          h0ts_E = Gaussian_basis(x = ts_E, mean = gknots, sd = sig)%*%theta_old
          h0tss_E = h0ts_E
          h0tss_E[h0tss_E<1e-40] = 1e-40
          ch0ts_E = Gaussian_basis_integ1(q = ts_E, mean = gknots, sd = sig)%*%theta_old
        }
        if (n_RC==0) {ch0ts_RC = 0}
        else {ch0ts_RC = Gaussian_basis_integ1(q = ts_RC, mean = gknots, sd = sig)%*%theta_old}
        if (n_LC==0) {
          S0ts_LC = 0
          S0ts_LC_diff = 1
        }
        else {
          S0ts_LC = exp(-Gaussian_basis_integ1(q = ts_LC, mean = gknots, sd = sig)%*%theta_old)
          S0ts_LC_diff = 1 - S0ts_LC
          S0ts_LC_diff[S0ts_LC_diff<1e-40] = 1e-40 # 20220514
        }
        if (n_IC==0) {
          S0ts_IC_L = S0ts_IC_R = 0
          S0ts_IC_diff = 1
        }
        else {
          S0ts_IC_L = exp(-Gaussian_basis_integ1(q = ts_IC_L, mean = gknots, sd = sig)%*%theta_old)
          S0ts_IC_R = exp(-Gaussian_basis_integ1(q = ts_IC_R, mean = gknots, sd = sig)%*%theta_old)
          S0ts_IC_diff = S0ts_IC_L - S0ts_IC_R
          S0ts_IC_diff[S0ts_IC_diff<1e-40] = 1e-40 # 20220514
        }
        
        # 0.5. log-likelihood
        like_gamma_new = sum(log(h0tss_E) - X_E%*%beta_new - Z_ni_E%*%gamma_new - ch0ts_E) - sum(ch0ts_RC) + sum(log(S0ts_LC_diff)) + sum(log(S0ts_IC_diff))
        
        likelihood_gamma_iter_before_ls[k] = like_gamma_new
        # 1. backtracking line search for gamma
        alpha_N_gamma = 1
        iteration_gamma = 0L
        # message("for loop No. is ", k)
        while (like_gamma_new <= like_gamma_old) {
          iteration_gamma = iteration_gamma + 1L
          # if (iteration_gamma==7) {browser()}
          # message("like in while loop equals ", like)
          if (alpha_N_gamma >= 1e-2) {alpha_N_gamma = alpha_N_gamma*0.6} # original step size multiplier 0.6
          else if (alpha_N_gamma < 1e-2 & alpha_N_gamma >= 1e-5) {alpha_N_gamma = alpha_N_gamma*5e-2} # original step size multiplier 5e-2
          else if (alpha_N_gamma < 1e-5 & alpha_N_gamma >= 1e-20) {alpha_N_gamma = alpha_N_gamma*1e-5} # original step size multiplier 1e-5
          # else if (alpha_N_gamma<1e-20 & alpha_N_gamma>=1e-100) {alpha_N_beta=alpha_N_gamma*1e-10}
          else {break}
          
          # iteration_gamma==1, alpha_N_gamma == 1e-1
          # iteration_gamma==2, alpha_N_gamma == 1e-2
          # iteration_gamma==3, alpha_N_gamma == 1e-3
          # iteration_gamma==4, alpha_N_gamma == 1e-6
          # iteration_gamma==5, alpha_N_gamma == 1e-14
          # iteration_gamma==6, alpha_N_gamma == 1e-22
          
          # 1. update gamma with search direction
          gamma_new = gamma_old + alpha_N_gamma * gamma_inc # bug of version 0.1.0 fixed
          
          # 1. accelerated time
          if (n_E==0) {ts_E = NA}
          else {ts_E = exp(-X_E%*%beta_new) * as.matrix(tapply(exp(-Z_E%*%gamma_new)*interval_length_E, id_Z_E, FUN = sum))}
          if (n_RC==0) {ts_RC = NA}
          else {ts_RC = exp(-X_RC%*%beta_new) * as.matrix(tapply(exp(-Z_RC%*%gamma_new)*interval_length_RC, id_Z_RC, FUN = sum))}
          if (n_LC==0) {ts_LC = NA}
          else {ts_LC = exp(-X_LC%*%beta_new) * as.matrix(tapply(exp(-Z_LC%*%gamma_new)*interval_length_LC, id_Z_LC, FUN = sum))}
          if (n_IC==0) {ts_IC_L = ts_IC_R = NA}
          else {
            ts_IC_L = exp(-X_IC%*%beta_new) * as.matrix(tapply(exp(-Z_IC_subset%*%gamma_new)*interval_length_IC_subset, id_Z_IC_subset, FUN = sum)) # special requirement
            ts_IC_R = exp(-X_IC%*%beta_new) * as.matrix(tapply(exp(-Z_IC%*%gamma_new)*interval_length_IC, id_Z_IC, FUN = sum))
          }
          ts_bind = rbind(ts_E, ts_RC, ts_LC, ts_IC_L, ts_IC_R) # ? not sure if this is suitable for "percentile" knots
          ts_bind = ts_bind[!is.na(ts_bind), , drop = FALSE]
          ts_bind = ts_bind[order(as.numeric(rownames(ts_bind))), , drop = FALSE]
          
          # 1. knots
          # # knots option 1
          # if (knots_option==1) {
          #   bryKnt = c(min(ts), max(ts)+1e-40)
          #   IntKnt = quantile(sort(ts)[2:(n-1)], seq(5, 95, length.out = numIntKnt)/100, type=1)
          # }
          # # knots option 2
          # if (knots_option==2) {
          #   Knt = quantile(sort(ts), seq(0, 1, length.out = numSp-dgrSp+1), type=1)
          #   bryKnt = c(Knt[1], Knt[numSp-dgrSp+1]+1e-40)
          #   IntKnt = Knt[2:(numSp-dgrSp)]
          # }
          # # knots option 3
          # if (knots_option==3) {
          #   ts_min = min(ts)
          #   ts_max = max(ts)
          #   bryKnt = c(ts_min, ts_max+1e-40)
          #   IntKnt = seq(ts_min, ts_max, length.out = numSp-dgrSp+1)[2:(numSp-dgrSp)]
          # }
          
          # Gaussian basis part
          if ((all(abs(beta_new-beta_old)>knots_fix_threshold) | all(abs(gamma_new-gamma_old)>knots_fix_threshold)) & k<=knots_fix_iter_no & frequent_knots_update==TRUE) { # 20210630 # 20220208 add frequent_knots_update check
            if (knots_option=="equal_space") {
              bryKnt = range(ts_bind)
              bin_width = (bryKnt[2] - bryKnt[1]) / (numIntKnt+1)
              gknots = seq(bryKnt[1], bryKnt[2], bin_width)
              sig = (2/3) * bin_width
            }
            else if (knots_option=="percentile") {
              # different strategy of calculating knots and sigmas
              gknots = quantile(sort(ts_bind), seq(knots_min_quantile, knots_max_quantile, length.out = numIntKnt+2)/100, type=1)
              dist = sapply(gknots, function(x) {abs(x - ts_bind)})
              sig = apply(dist, 2, function(x) {quantile(x, sig_coverage) / 2})
            }
            else if (knots_option=="percentile_1") { # test on 20220525
              ts_bind_1 = rbind(ts_E, ts_RC, ts_LC, ts_IC_R)
              ts_bind_1 = ts_bind_1[!is.na(ts_bind_1), , drop = FALSE]
              ts_bind_1 = ts_bind_1[order(as.numeric(rownames(ts_bind_1))), , drop = FALSE]
              gknots = quantile(sort(ts_bind_1), seq(knots_min_quantile, knots_max_quantile, length.out = numIntKnt+2)/100, type=1)
              dist = sapply(gknots, function(x) {abs(x - ts_bind_1)})
              sig = apply(dist, 2, function(x) {quantile(x, sig_coverage) / 2})
            }
          }
          
          # gauss = pnorm(0, gknots, sig)*sqrt(2*pi*sig^2)
          # gauss = pnorm(bryKnt[1], gknots, sig)*sqrt(2*pi*sig^2) # which one to choose?
          
          # 1. spline and cumulative spline
          # if (n_E==0) {phi_E = matrix(0, 1, numSp)}
          # else {
          #   phi_E = mSpline(ts_E, knots=IntKnt, degree=dgrSp, intercept=TRUE, Boundary.knots=bryKnt)
          # }
          # cphi = iSpline(ts, knots=IntKnt, degree=dgrSp, intercept=TRUE, Boundary.knots=bryKnt)
          
          # 1. hazard and cumulative hazard
          # if (n_E==0) {h0ts_E = h0tss_E = 1}
          # else {
          #   h0ts_E = phi_E%*%theta_old
          #   h0tss_E = h0ts_E
          #   h0tss_E[h0tss_E<1e-40] = 1e-40
          # }
          # ch0ts = cphi%*%theta_old
          
          # Gaussian basis part
          if (n_E==0) {
            h0ts_E = h0tss_E = 1
            ch0ts_E = 0
          }
          else {
            h0ts_E = Gaussian_basis(x = ts_E, mean = gknots, sd = sig)%*%theta_old
            h0tss_E = h0ts_E
            h0tss_E[h0tss_E<1e-40] = 1e-40
            ch0ts_E = Gaussian_basis_integ1(q = ts_E, mean = gknots, sd = sig)%*%theta_old
          }
          if (n_RC==0) {ch0ts_RC = 0}
          else {ch0ts_RC = Gaussian_basis_integ1(q = ts_RC, mean = gknots, sd = sig)%*%theta_old}
          if (n_LC==0) {
            S0ts_LC = 0
            S0ts_LC_diff = 1
          }
          else {
            S0ts_LC = exp(-Gaussian_basis_integ1(q = ts_LC, mean = gknots, sd = sig)%*%theta_old)
            S0ts_LC_diff = 1 - S0ts_LC
            S0ts_LC_diff[S0ts_LC_diff<1e-40] = 1e-40 # 20220514
          }
          if (n_IC==0) {
            S0ts_IC_L = S0ts_IC_R = 0
            S0ts_IC_diff = 1
          }
          else {
            S0ts_IC_L = exp(-Gaussian_basis_integ1(q = ts_IC_L, mean = gknots, sd = sig)%*%theta_old)
            S0ts_IC_R = exp(-Gaussian_basis_integ1(q = ts_IC_R, mean = gknots, sd = sig)%*%theta_old)
            S0ts_IC_diff = S0ts_IC_L - S0ts_IC_R
            S0ts_IC_diff[S0ts_IC_diff<1e-40] = 1e-40 # 20220514
          }
          
          # 1. log-likelihood
          like_gamma_new = sum(log(h0tss_E) - X_E%*%beta_new - Z_ni_E%*%gamma_new - ch0ts_E) - sum(ch0ts_RC) + sum(log(S0ts_LC_diff)) + sum(log(S0ts_IC_diff))
        }
        
        # # one way to escape from local maximum 20210608
        # if (any(abs(grad_gamma)>=0.01) & iteration_gamma == 7) {gamma_new = gamma_old - 1e-1*abs(gamma_old)*sign(grad_gamma)}
        
        ts_range_gamma_iter[k, ] = range(ts_bind)
        likelihood_gamma_iter_after_ls[k] = like_gamma_new
        iter_Newton_gamma[k] = iteration_gamma
        gamma_iter[k, ] = gamma_new
        
        like_theta_old = like_gamma_new
      }
      else {like_theta_old = like_beta_new}
      
      
      ################ roughness penalty matrix ##################
      # 1. penalty equals theta_transpose %*% R %*% theta, R equals to phi_sd %*% phi_sd_transpose
      # R=matrix(0, nrow=numSp, ncol=numSp) 
      # xknots = c(rep(bryKnt[1], ordSp), IntKnt, rep(bryKnt[2], ordSp))
      # # browser()
      # for (ii in 1:numSp){
      #   for (jj in ii:numSp){
      #     if (jj - ii<ordSp){
      #       kntset = xknots[xknots>=xknots[jj] & xknots<=xknots[ii+ordSp]]
      #       kntsum = 0
      #       for (kk in 1:(length(kntset)-1)){
      #         kntsum = kntsum + mSpline(kntset[kk], knots=IntKnt, degree=dgrSp, intercept=T, Boundary.knots=bryKnt, 
      #                                   derivs=dgrSp)[ii]*mSpline(kntset[kk], knots=IntKnt, degree=dgrSp, intercept=T, 
      #                                                             Boundary.knots=bryKnt, derivs=dgrSp)[jj]*(kntset[kk+1]-kntset[kk])
      #       }
      #       R[ii, jj] = kntsum
      #     }
      #   }
      # }
      # R[lower.tri(R, diag = FALSE)] = t(R)[lower.tri(R, diag = FALSE)]
      
      
      # try alternative
      
      # R=matrix(0, nrow=numSp, ncol=numSp) 
      # xknots = c(rep(bryKnt[1], ordSp), IntKnt, rep(bryKnt[2], ordSp))
      # for (ii in 1:numSp){
      #   for (jj in ii:numSp){
      #     if (jj - ii<ordSp){
      #       kntset = xknots[xknots>=xknots[jj] & xknots<=xknots[ii+ordSp]]
      #       ddphi_kntset = mSpline(kntset, knots=IntKnt, degree=dgrSp, intercept=T, Boundary.knots=bryKnt, derivs=dgrSp)
      #       kntsum = 0
      #       for (kk in 1:(length(kntset)-1)){
      #         kntsum = kntsum + ddphi_kntset[kk, ii]*ddphi_kntset[kk, jj]*(kntset[kk+1]-kntset[kk])
      #       }
      #       R[ii, jj] = kntsum
      #     }
      #   }
      # }
      # R[lower.tri(R, diag = FALSE)] = t(R)[lower.tri(R, diag = FALSE)]
      # browser()
      if ((all(abs(beta_new-beta_old)>knots_fix_threshold) | all(abs(gamma_new-gamma_old)>knots_fix_threshold)) & k<=knots_fix_iter_no) {
        # Ding's optimised version
        num_dots = 5*n
        ts_bind_range = range(ts_bind)
        bin_edges = seq(ts_bind_range[1], ts_bind_range[2], length.out=num_dots+1)
        bin_middle_points = as.matrix((bin_edges[2:(num_dots+1)] + bin_edges[1:num_dots])/2)
        ddphi_bin_middle_points = Gaussian_basis_deriv2(x = bin_middle_points, mean = gknots, sd = sig)
        # R = matrix(0, numIntKnt+2, numIntKnt+2)
        # for (u in 1:(numIntKnt+2)) {
        #   for (r in u:(numIntKnt+2)) {
        #     R[u, r] = (ts_bind_range[2]-ts_bind_range[1])/num_dots * sum(ddphi_bin_middle_points[, u]*ddphi_bin_middle_points[, r])
        #   }
        # }
        # R[lower.tri(R, diag = FALSE)] = t(R)[lower.tri(R, diag = FALSE)]
        R = t(ddphi_bin_middle_points)%*%ddphi_bin_middle_points * (ts_bind_range[2]-ts_bind_range[1])/num_dots # 20211102
      }
      
      # ddphi_bind = mSpline(ts, knots=IntKnt, degree=dgrSp, intercept=TRUE, Boundary.knots=bryKnt, derivs = 2)
      # int_ddphi = matrix(NA, length(ts)-1, 1)
      # for(interval in 1:(length(ts)-1)) {
      #   int_ddphi[interval] = ddphi_bind[interval, ]%*%t(ddphi_bind[interval, ])*(ts[interval+1]-ts[interval])
      # }
      # browser()
      # penalty = t(theta_old)%*%(h0ts_bind)%*%theta_old
      
      ########################## MI step for theta ##########################
      # 1. penalised likelihood
      # if (!missing(Zmat) & !missing(gamma_initial) & (dim_Zmat[1]>n+1)) {like_theta_old = like_gamma_new}
      # else {like_theta_old = like_beta_new}
      plike_old = like_theta_old - smooth*t(theta_old)%*%R%*%theta_old
      
      # 1.5. terms of theta update algorithm
      # if (n_E==0) {quan_E_grad_theta_1 = matrix(0, 1, numSp)}
      # else {quan_E_grad_theta_1 = phi_E/(h0tss_E%*%matrix(1, 1, numSp))}
      # quan_deno = cphi
      
      
      if (n_E==0) {quan_E_grad_theta_1 = quan_E_grad_theta_2 = matrix(0, 1, numIntKnt+2)}
      else {
        phi_E = Gaussian_basis(x = ts_E, mean = gknots, sd = sig)
        h0ts_E = phi_E%*%theta_old
        h0tss_E = h0ts_E
        h0tss_E[h0tss_E<1e-40] = 1e-40
        # quan_E_grad_theta_1 = phi_E/(h0tss_E%*%matrix(1, 1, numIntKnt+2))
        quan_E_grad_theta_1 = phi_E / as.numeric(h0tss_E) # 20220117 faster alternative of above line, use R recycling rule, may cause trouble
        quan_E_grad_theta_2 = Gaussian_basis_integ1(q = ts_E, mean = gknots, sd = sig)
      }
      if(n_RC==0) {quan_RC_grad_theta = matrix(0, 1, numIntKnt+2)}
      else {quan_RC_grad_theta = Gaussian_basis_integ1(q = ts_RC, mean = gknots, sd = sig)}
      if(n_LC==0) {quan_LC_grad_theta = matrix(0, 1, numIntKnt+2)}
      else {
        cphi_LC = Gaussian_basis_integ1(q = ts_LC, mean = gknots, sd = sig)
        S0ts_LC = exp(-cphi_LC%*%theta_old)
        S0ts_LC_diff = 1 - S0ts_LC
        S0ts_LC_diff[S0ts_LC_diff<1e-40] = 1e-40 # 20220514
        quan_LC_grad_theta = cphi_LC * as.numeric(S0ts_LC / S0ts_LC_diff)
      }
      if(n_IC==0) {quan_IC_R_grad_theta = quan_IC_L_grad_theta = matrix(0, 1, numIntKnt+2)}
      else {
        cphi_IC_L = Gaussian_basis_integ1(q = ts_IC_L, mean = gknots, sd = sig)
        S0ts_IC_L = exp(-cphi_IC_L%*%theta_old)
        cphi_IC_R = Gaussian_basis_integ1(q = ts_IC_R, mean = gknots, sd = sig)
        S0ts_IC_R = exp(-cphi_IC_R%*%theta_old)
        S0ts_IC_diff = S0ts_IC_L - S0ts_IC_R
        S0ts_IC_diff[S0ts_IC_diff<1e-40] = 1e-40 # 20220514
        quan_IC_R_grad_theta = cphi_IC_R * as.numeric(S0ts_IC_R / S0ts_IC_diff)
        quan_IC_L_grad_theta = cphi_IC_L * as.numeric(S0ts_IC_L / S0ts_IC_diff)
      }
      
      # nume = colSums(quan_E_grad_theta_1) + colSums(quan_LC_grad_theta) + colSums(quan_IC_R_grad_theta) - 2*smooth*R%*%theta_old*(R%*%theta_old<=0) + eta
      # deno = colSums(quan_E_grad_theta_2) + colSums(quan_RC_grad_theta) + colSums(quan_IC_L_grad_theta) + 2*smooth*R%*%theta_old*(R%*%theta_old>0) + eta
      
      nume = colSums(quan_E_grad_theta_1) + colSums(quan_LC_grad_theta) + colSums(quan_IC_R_grad_theta) - 2*smooth*R%*%theta_old*(R%*%theta_old<=0) # 20220512
      deno = colSums(quan_E_grad_theta_2) + colSums(quan_RC_grad_theta) + colSums(quan_IC_L_grad_theta) + 2*smooth*R%*%theta_old*(R%*%theta_old>0) # 20220512
      
      # 1.5. update theta temporarily
      # theta_new = theta_old*(nume/deno)
      # theta_inc = theta_new-theta_old
      grad_theta = nume - deno
      # browser()
      
      # grad_theta = colSums(quan_E_grad_theta_1 - quan_E_grad_theta_2) + colSums(-quan_RC_grad_theta) + colSums(quan_LC_grad_theta) + colSums(-quan_IC_L_grad_theta+quan_IC_R_grad_theta) - 2*smooth*R%*%theta_old
      # if (k==knots_fix_iter_no) {theta_old[grad_theta > 0] = 1} # reset theta 20220420
      # deno = colSums(quan_E_grad_theta_2) + colSums(quan_RC_grad_theta) + colSums(quan_IC_L_grad_theta) + 2*smooth*R%*%theta_old*(R%*%theta_old>0) + eta
      
      deno[deno<1e-40] = 1e-40 # 20220512
      A_mat = theta_old / deno
      theta_inc = A_mat*grad_theta
      theta_new = theta_old + theta_inc
      # theta_new[theta_new<1e-8]=0
      
      # grad_theta = nume - deno
      # theta_old[grad_theta > 0] = 1 # reset theta 20220420
      # if (any(grad_theta > 0)) {
      #   theta_inc = grad_theta / deno
      # }
      # else {theta_inc = (theta_old*grad_theta) / deno}
      # theta_new = theta_old + theta_inc
      deno_iter[k, ] = t(deno)
      multiplier_iter[k, ] = t(nume/deno)
      grad_theta_iter[k, ] = t(grad_theta)
      
      # 1.5. hazard and cumulative hazard
      # if (n_E==0) {h0ts_E = h0tss_E = 1}
      # else {
      #   h0ts_E = phi_E%*%theta_new
      #   h0tss_E = h0ts_E
      #   h0tss_E[h0tss_E<1e-40] = 1e-40
      # }
      # ch0ts = cphi%*%theta_new
      
      # Gaussian basis part
      if (n_E==0) {
        h0ts_E = h0tss_E = 1
        ch0ts_E = 0
      }
      else {
        h0ts_E = Gaussian_basis(x = ts_E, mean = gknots, sd = sig)%*%theta_new
        h0tss_E = h0ts_E
        h0tss_E[h0tss_E<1e-40] = 1e-40
        ch0ts_E = Gaussian_basis_integ1(q = ts_E, mean = gknots, sd = sig)%*%theta_new
      }
      if (n_RC==0) {ch0ts_RC = 0}
      else {ch0ts_RC = Gaussian_basis_integ1(q = ts_RC, mean = gknots, sd = sig)%*%theta_new}
      if (n_LC==0) {
        S0ts_LC = 0
        S0ts_LC_diff = 1
      }
      else {
        S0ts_LC = exp(-Gaussian_basis_integ1(q = ts_LC, mean = gknots, sd = sig)%*%theta_new)
        S0ts_LC_diff = 1 - S0ts_LC
        S0ts_LC_diff[S0ts_LC_diff<1e-40] = 1e-40 # 20220514
      }
      if (n_IC==0) {
        S0ts_IC_L = S0ts_IC_R = 0
        S0ts_IC_diff = 1
      }
      else {
        S0ts_IC_L = exp(-Gaussian_basis_integ1(q = ts_IC_L, mean = gknots, sd = sig)%*%theta_new)
        S0ts_IC_R = exp(-Gaussian_basis_integ1(q = ts_IC_R, mean = gknots, sd = sig)%*%theta_new)
        S0ts_IC_diff = S0ts_IC_L - S0ts_IC_R
        S0ts_IC_diff[S0ts_IC_diff<1e-40] = 1e-40 # 20220514
      }
      
      # 1.5. log-likelihood and penalised log-likelihood
      if (!missing(Zmat) & !missing(gamma_initial) & (dim_Zmat[1]>n+1)) {
        like_theta_new = sum(log(h0tss_E) - X_E%*%beta_new - Z_ni_E%*%gamma_new - ch0ts_E) - sum(ch0ts_RC) + sum(log(S0ts_LC_diff)) + sum(log(S0ts_IC_diff))
      }
      else {like_theta_new = sum(log(h0tss_E) - X_E%*%beta_new - ch0ts_E) - sum(ch0ts_RC) + sum(log(S0ts_LC_diff)) + sum(log(S0ts_IC_diff))}
      plike_new = like_theta_new - smooth*t(theta_new)%*%R%*%theta_new
      
      penlike_theta_iter_before_ls[k] = plike_new
      # 2. backtracking line search for theta
      alpha_MI = 1
      iteration_theta = 0L
      while (plike_new <= plike_old) {
        iteration_theta = iteration_theta + 1L
        if (alpha_MI >= 1e-2) {alpha_MI = alpha_MI*0.6}
        else if (alpha_MI < 1e-2 & alpha_MI >= 1e-5) {alpha_MI = alpha_MI*5e-2}
        else if (alpha_MI < 1e-5 & alpha_MI >= 1e-20) {alpha_MI = alpha_MI*1e-5}
        # else if (alpha_MI<1e-20 & alpha_MI>=1e-100) {alpha_MI=alpha_MI*1e-10}
        else {break}
        
        # 2. update theta with search direction
        theta_new = theta_old + alpha_MI * theta_inc
        
        # Gaussian basis part
        if (n_E==0) {
          h0ts_E = h0tss_E = 1
          ch0ts_E = 0
        }
        else {
          h0ts_E = Gaussian_basis(x = ts_E, mean = gknots, sd = sig)%*%theta_new
          h0tss_E = h0ts_E
          h0tss_E[h0tss_E<1e-40] = 1e-40
          ch0ts_E = Gaussian_basis_integ1(q = ts_E, mean = gknots, sd = sig)%*%theta_new
        }
        if (n_RC==0) {ch0ts_RC = 0}
        else {ch0ts_RC = Gaussian_basis_integ1(q = ts_RC, mean = gknots, sd = sig)%*%theta_new}
        if (n_LC==0) {
          S0ts_LC = 0
          S0ts_LC_diff = 1
        }
        else {
          S0ts_LC = exp(-Gaussian_basis_integ1(q = ts_LC, mean = gknots, sd = sig)%*%theta_new)
          S0ts_LC_diff = 1 - S0ts_LC
          S0ts_LC_diff[S0ts_LC_diff<1e-40] = 1e-40 # 20220514
        }
        if (n_IC==0) {
          S0ts_IC_L = S0ts_IC_R = 0
          S0ts_IC_diff = 1
        }
        else {
          S0ts_IC_L = exp(-Gaussian_basis_integ1(q = ts_IC_L, mean = gknots, sd = sig)%*%theta_new)
          S0ts_IC_R = exp(-Gaussian_basis_integ1(q = ts_IC_R, mean = gknots, sd = sig)%*%theta_new)
          S0ts_IC_diff = S0ts_IC_L - S0ts_IC_R
          S0ts_IC_diff[S0ts_IC_diff<1e-40] = 1e-40 # 20220514
        }
        
        # 2. log-likelihood and penalised log-likelihood
        if (!missing(Zmat) & !missing(gamma_initial) & (dim_Zmat[1]>n+1)) {
          like_theta_new = sum(log(h0tss_E) - X_E%*%beta_new - Z_ni_E%*%gamma_new - ch0ts_E) - sum(ch0ts_RC) + sum(log(S0ts_LC_diff)) + sum(log(S0ts_IC_diff))
        }
        else {like_theta_new = sum(log(h0tss_E) - X_E%*%beta_new - ch0ts_E) - sum(ch0ts_RC) + sum(log(S0ts_LC_diff)) + sum(log(S0ts_IC_diff))}
        plike_new = like_theta_new - smooth*t(theta_new)%*%R%*%theta_new
      }
      
      penlike_theta_iter_after_ls[k] = plike_new
      iter_MI_theta[k] = iteration_theta
      theta_iter[k, ] = theta_new
      ts_range_theta_iter[k, ] = range(ts_bind)
      
      knots_iter[k, ] = gknots
      
      iteration_no = k
      
      if ((all(abs(beta_new-beta_old)>knots_fix_threshold) | all(abs(gamma_new-gamma_old)>knots_fix_threshold)) & k<=knots_fix_iter_no & frequent_knots_update==FALSE) { # 20210630 # 20220208 add frequent_knots_update check
        if (knots_option=="equal_space") {
          bryKnt = range(ts_bind)
          bin_width = (bryKnt[2] - bryKnt[1]) / (numIntKnt+1)
          gknots = seq(bryKnt[1], bryKnt[2], bin_width)
          sig = (2/3) * bin_width
        }
        else if (knots_option=="percentile") {
          # different strategy of calculating knots and sigmas
          gknots = quantile(sort(ts_bind), seq(knots_min_quantile, knots_max_quantile, length.out = numIntKnt+2)/100, type=1)
          dist = sapply(gknots, function(x) {abs(x - ts_bind)})
          sig = apply(dist, 2, function(x) {quantile(x, sig_coverage) / 2})
        }
        else if (knots_option=="percentile_1") { # test on 20220525
          ts_bind_1 = rbind(ts_E, ts_RC, ts_LC, ts_IC_R)
          ts_bind_1 = ts_bind_1[!is.na(ts_bind_1), , drop = FALSE]
          ts_bind_1 = ts_bind_1[order(as.numeric(rownames(ts_bind_1))), , drop = FALSE]
          gknots = quantile(sort(ts_bind_1), seq(knots_min_quantile, knots_max_quantile, length.out = numIntKnt+2)/100, type=1)
          dist = sapply(gknots, function(x) {abs(x - ts_bind_1)})
          sig = apply(dist, 2, function(x) {quantile(x, sig_coverage) / 2})
        }
      }
      ########################## stopping criteria ##########################
      if (!missing(Zmat) & !missing(gamma_initial) & (dim_Zmat[1]>n+1)) {
        beta_hat = beta_new
        gamma_hat = gamma_new
        theta_hat = theta_new
        if (all(abs(beta_new-beta_old)<threshold) & all(abs(gamma_new-gamma_old)<threshold) & all(abs(theta_new-theta_old)<threshold)) {
        # if (all(abs(beta_new-beta_old)<threshold) & all(abs(gamma_new-gamma_old)<threshold)) { # 20220523 test, remove the limitation from theta as there are active constraints
          cat("The stopping criteria are met and the results are obtained after", iteration_no, "iterations.\n")
          break
        }
        else if (iteration_no==maxiter) {
          warning("The stopping criteria are NOT met after ", iteration_no, " itreations!\n")
        }
        else {
          beta_old = beta_new
          gamma_old = gamma_new
          if (k==knots_fix_iter_no) {theta_old = matrix(1, numIntKnt+2, 1)} # 20220208 inspired by Jun's modification
          else {theta_old = theta_new}
          # theta_old = theta_new
          like_beta_old = like_theta_new
          # print(paste("The penalised likelihood after iteration No.", k, " is ", plike_new))
        }
      }
      else {
        beta_hat = beta_new
        theta_hat = theta_new
        if (all(abs(beta_new-beta_old)<threshold) & all(abs(theta_new-theta_old)<threshold)) {
        # if (all(abs(beta_new-beta_old)<threshold)) { # 20220523 test, remove the limitation from theta as there are active constraints
          cat("The stopping criteria are met and the results are obtained after", iteration_no, "iterations.\n")
          break
        }
        else if (iteration_no==maxiter) {
          warning("The stopping criteria are NOT met after ", iteration_no, " itreations!\n")
        }
        else {
          beta_old = beta_new
          if (k==knots_fix_iter_no) {theta_old = matrix(1, numIntKnt+2, 1)} # 20220208 inspired by Jun's modification
          else {theta_old = theta_new}
          # theta_old = theta_new
          like_beta_old = like_theta_new
          # print(paste("The penalised likelihood after iteration No.", k, " is ", plike_new))
        }
      }
      # test if nonparametric estimates of AFT model are noisy with only 5 iterations
      # if (k==5){break}
    }
    
    # }, error=function(e){ #20181129
    #   theta_new=matrix(NA, numSp, 1)
    #   beta_new=matrix(NA, p, 1)
    # })
    
    ########################## estimation results ##########################
    time_range = rbind(range(time_bind), range(ts_bind))
    rownames(time_range) <- c("original time", "accelerated time")
    colnames(time_range) <- c("min", "max")
    ts_range_beta_iter = ts_range_beta_iter[1:k, ]
    ts_range_theta_iter = ts_range_theta_iter[1:k, ]
    knots_iter = knots_iter[1:k, ]
    if (any(abs(diff(knots_iter[, 1]))<=1e-12 & abs(diff(knots_iter[, ncol(knots_iter)]))<=1e-12)) {
      if(identical(max(which(abs(diff(knots_iter[, 1]))>1e-12 | abs(diff(knots_iter[, ncol(knots_iter)]))>1e-12))+1,-Inf)) {
        cat("The knots location stops changing since the 1st iteration.\n")
      }
      else {
        cat("The knots location stops changing since the", max(which(abs(diff(knots_iter[, 1]))>1e-12 | abs(diff(knots_iter[, ncol(knots_iter)]))>1e-12))+1, "th iteration.\n")
        # cat("The knots location stops changing since the", max(which(abs(diff(knots_iter[, 1]))>.Machine$double.eps | abs(diff(knots_iter[, ncol(knots_iter)]))>.Machine$double.eps))+1, "th iteration.\n")
      }
    }
    multiple_ts_time = diff(range(ts_bind))/diff(range(time, na.rm=TRUE))
    if (multiple_ts_time>1) {
      warning("This is a decelerated failure time case! The range of the decelerated time is ", multiple_ts_time, " times that of the original time!\n")
    }
    
    grad_beta_iter = grad_beta_iter[1:k, ]
    Hes_beta_iter = Hes_beta_iter[1:k, ]
    det_Hes_beta_iter = det_Hes_beta_iter[1:k]
    eigenvl_Hes_beta_iter = eigenvl_Hes_beta_iter[1:k, ]
    beta_iter = beta_iter[1:k, ]
    
    deno_iter = deno_iter[1:k, ]
    multiplier_iter = multiplier_iter[1:k, ]
    grad_theta_iter = grad_theta_iter[1:k, ]
    theta_iter = theta_iter[1:k, ]
    
    if (!missing(Zmat) & !missing(gamma_initial) & (dim_Zmat[1]>n+1)) {
      cvg = cbind(likelihood_beta_iter_before_ls, iter_Newton_beta, likelihood_beta_iter_after_ls, likelihood_gamma_iter_before_ls, iter_Newton_gamma, likelihood_gamma_iter_after_ls, penlike_theta_iter_before_ls, iter_MI_theta, penlike_theta_iter_after_ls)[1:k, ]
      colnames(cvg) = c("like_beta_iter_before_ls", "iters_ls_beta", "like_beta_iter_after_ls", "like_gamma_iter_before_ls", "iters_ls_gamma", "like_gamma_iter_after_ls", "penlike_theta_iter_before_ls", "iters_ls_theta", "penlike_theta_iter_after_ls")
      grad_gamma_iter = grad_gamma_iter[1:k, ]
      Hes_gamma_iter = Hes_gamma_iter[1:k, ]
      det_Hes_gamma_iter = det_Hes_gamma_iter[1:k]
      eigenvl_Hes_gamma_iter = eigenvl_Hes_gamma_iter[1:k, ]
      ts_range_gamma_iter = ts_range_gamma_iter[1:k, ]
      gamma_iter = gamma_iter[1:k, ]
    }
    else {
      cvg = cbind(likelihood_beta_iter_before_ls, iter_Newton_beta, likelihood_beta_iter_after_ls, penlike_theta_iter_before_ls, iter_MI_theta, penlike_theta_iter_after_ls)[1:k, ]
      colnames(cvg) = c("like_beta_iter_before_ls", "iters_ls_beta", "like_beta_iter_after_ls", "penlike_theta_iter_before_ls", "iters_ls_theta", "penlike_theta_iter_after_ls")
    }
    
    # browser()
    ########################## Asymptotic Normality ##########################
    # elements for Hessian calculation
    if (n_E==0) {dts_E_beta = matrix(0, 1, p)}
    else {dts_E_beta = - X_E * as.numeric(ts_E)}
    if (n_RC==0) {dts_RC_beta = matrix(0, 1, p)}
    else {dts_RC_beta = - X_RC * as.numeric(ts_RC)}
    if (n_LC==0) {dts_LC_beta = matrix(0, 1, p)}
    else {dts_LC_beta = - X_LC * as.numeric(ts_LC)}
    if (n_IC==0) {dts_IC_L_beta = dts_IC_R_beta = matrix(0, 1, p)}
    else {
      dts_IC_L_beta = - X_IC * as.numeric(ts_IC_L)
      dts_IC_R_beta = - X_IC * as.numeric(ts_IC_R)
    }
    
    if (n_E==0) {
      dphi_E = phi_E = 0
      ddh0ts_E = dh0ts_E = 0
      h0ts_E = h0tss_E = 1
    }
    else {
      ddh0ts_E = Gaussian_basis_deriv2(x = ts_E, mean = gknots, sd = sig)%*%theta_hat
      dphi_E = Gaussian_basis_deriv1(x = ts_E, mean = gknots, sd = sig)
      dh0ts_E = dphi_E%*%theta_hat
      phi_E = Gaussian_basis(x = ts_E, mean = gknots, sd = sig)
      h0ts_E = phi_E%*%theta_hat
      h0tss_E = h0ts_E
      h0tss_E[h0tss_E<1e-40] = 1e-40
    }
    if (n_RC==0) {
      phi_RC = 0
      dh0ts_RC = h0ts_RC = 0
    }
    else {
      dh0ts_RC = Gaussian_basis_deriv1(x = ts_RC, mean = gknots, sd = sig)%*%theta_hat
      phi_RC = Gaussian_basis(x = ts_RC, mean = gknots, sd = sig)
      h0ts_RC = phi_RC%*%theta_hat
    }
    if (n_LC==0) {
      phi_LC = cphi_LC = 0
      dh0ts_LC = h0ts_LC = S0ts_LC = 0
      S0ts_LC_diff = 1
    }
    else {
      dh0ts_LC = Gaussian_basis_deriv1(x = ts_LC, mean = gknots, sd = sig)%*%theta_hat
      phi_LC = Gaussian_basis(x = ts_LC, mean = gknots, sd = sig)
      h0ts_LC = phi_LC%*%theta_hat
      cphi_LC = Gaussian_basis_integ1(q = ts_LC, mean = gknots, sd = sig)
      S0ts_LC = exp(-cphi_LC%*%theta_hat)
      S0ts_LC_diff = 1 - S0ts_LC
      S0ts_LC_diff[S0ts_LC_diff<1e-40] = 1e-40 # 20220514
    }
    if (n_IC==0) {
      phi_IC_L = cphi_IC_L = phi_IC_R = cphi_IC_R = 0
      dh0ts_IC_L = h0ts_IC_L = S0ts_IC_L = dh0ts_IC_R = h0ts_IC_R = S0ts_IC_R = 0
      S0ts_IC_diff = 1 
    }
    else {
      dh0ts_IC_L = Gaussian_basis_deriv1(x = ts_IC_L, mean = gknots, sd = sig)%*%theta_hat
      phi_IC_L = Gaussian_basis(x = ts_IC_L, mean = gknots, sd = sig)
      h0ts_IC_L = phi_IC_L%*%theta_hat
      cphi_IC_L = Gaussian_basis_integ1(q = ts_IC_L, mean = gknots, sd = sig)
      S0ts_IC_L = exp(-cphi_IC_L%*%theta_hat)
      dh0ts_IC_R = Gaussian_basis_deriv1(x = ts_IC_R, mean = gknots, sd = sig)%*%theta_hat
      phi_IC_R = Gaussian_basis(x = ts_IC_R, mean = gknots, sd = sig)
      h0ts_IC_R = phi_IC_R%*%theta_hat
      cphi_IC_R = Gaussian_basis_integ1(q = ts_IC_R, mean = gknots, sd = sig)
      S0ts_IC_R = exp(-cphi_IC_R%*%theta_hat)
      S0ts_IC_diff = S0ts_IC_L - S0ts_IC_R
      S0ts_IC_diff[S0ts_IC_diff<1e-40] = 1e-40 # 20220514
    }
    
    # Hessian calculation, 3 Full Hessians, 3 Cross Hessians
    # full Hessian beta
    if (n_E==0) {quan_E_Hessian_beta = matrix(0, 1, 1)}
    else {quan_E_Hessian_beta = (ddh0ts_E / h0tss_E - (dh0ts_E / h0tss_E)^2 - dh0ts_E) * ts_E^2 + (dh0ts_E / h0tss_E - h0ts_E) * ts_E}
    if(n_RC==0) {quan_RC_Hessian_beta = matrix(0, 1, 1)}
    else {quan_RC_Hessian_beta = - dh0ts_RC * ts_RC^2 - h0ts_RC * ts_RC}
    if(n_LC==0) {quan_LC_Hessian_beta = matrix(0, 1, 1)}
    else {quan_LC_Hessian_beta = S0ts_LC*(-h0ts_LC^2 + dh0ts_LC)*ts_LC^2 / S0ts_LC_diff - (S0ts_LC*h0ts_LC*ts_LC / S0ts_LC_diff)^2 + S0ts_LC*h0ts_LC*ts_LC / S0ts_LC_diff}
    if(n_IC==0) {quan_IC_Hessian_beta = matrix(0, 1, 1)}
    else {
      quan_IC_Hessian_beta = S0ts_IC_L*(h0ts_IC_L^2 - dh0ts_IC_L)*ts_IC_L^2 / S0ts_IC_diff - (S0ts_IC_L*h0ts_IC_L*ts_IC_L / S0ts_IC_diff)^2 - S0ts_IC_L*h0ts_IC_L*ts_IC_L / S0ts_IC_diff +
                             S0ts_IC_R*(-h0ts_IC_R^2 + dh0ts_IC_R)*ts_IC_R^2 / S0ts_IC_diff - (S0ts_IC_R*h0ts_IC_R*ts_IC_R / S0ts_IC_diff)^2 + S0ts_IC_R*h0ts_IC_R*ts_IC_R / S0ts_IC_diff +
                             2*(S0ts_IC_R*h0ts_IC_R*ts_IC_R * S0ts_IC_L*h0ts_IC_L*ts_IC_L) / S0ts_IC_diff^2
    }
    Hessian_beta = t(-X_E)%*%(as.numeric(quan_E_Hessian_beta)*(-X_E)) +
                   t(-X_RC)%*%(as.numeric(quan_RC_Hessian_beta)*(-X_RC)) +
                   t(-X_LC)%*%(as.numeric(quan_LC_Hessian_beta)*(-X_LC)) +
                   t(-X_IC)%*%(as.numeric(quan_IC_Hessian_beta)*(-X_IC))
    
    # cross Hessian beta theta
    if (n_E==0) {quan_E_Hessian_beta_theta = matrix(0, 1, numIntKnt+2)}
    else {quan_E_Hessian_beta_theta = dphi_E / as.numeric(h0tss_E) - phi_E * as.numeric(dh0ts_E / h0tss_E^2 + 1)}
    if (n_RC==0) {quan_RC_Hessian_beta_theta = matrix(0, 1, numIntKnt+2)}
    else {quan_RC_Hessian_beta_theta = - phi_RC}
    if (n_LC==0) {quan_LC_Hessian_beta_theta = matrix(0, 1, numIntKnt+2)}
    else {quan_LC_Hessian_beta_theta = (cphi_LC * as.numeric(- h0ts_LC) + phi_LC) * as.numeric(S0ts_LC / S0ts_LC_diff) - cphi_LC * as.numeric(h0ts_LC * (S0ts_LC / S0ts_LC_diff)^2)}
    if (n_IC==0) {quan_IC_L_Hessian_beta_theta = quan_IC_R_Hessian_beta_theta = matrix(0, 1, numIntKnt+2)}
    else {
      quan_IC_L_Hessian_beta_theta = (cphi_IC_L * as.numeric(h0ts_IC_L) - phi_IC_L) * as.numeric(S0ts_IC_L / S0ts_IC_diff) + (cphi_IC_L * as.numeric(-S0ts_IC_L) + cphi_IC_R * as.numeric(S0ts_IC_R)) * as.numeric(S0ts_IC_L*h0ts_IC_L/S0ts_IC_diff^2)
      quan_IC_R_Hessian_beta_theta = (cphi_IC_R * as.numeric(-h0ts_IC_R) + phi_IC_R) * as.numeric(S0ts_IC_R / S0ts_IC_diff) + (cphi_IC_L * as.numeric(S0ts_IC_L) - cphi_IC_R * as.numeric(S0ts_IC_R)) * as.numeric(S0ts_IC_R*h0ts_IC_R/S0ts_IC_diff^2)
    }
    
    Hessian_beta_theta = t(dts_E_beta)%*%quan_E_Hessian_beta_theta +
                         t(dts_RC_beta)%*%quan_RC_Hessian_beta_theta +
                         t(dts_LC_beta)%*%quan_LC_Hessian_beta_theta +
                         t(dts_IC_L_beta)%*%quan_IC_L_Hessian_beta_theta +
                         t(dts_IC_R_beta)%*%quan_IC_R_Hessian_beta_theta
    
    # full Hessian theta
    if (n_E==0) {quan_E_Hessian_theta_no_penalty = matrix(0, numIntKnt+2, numIntKnt+2)}
    else {quan_E_Hessian_theta_no_penalty = - t(phi_E / as.numeric(h0tss_E^2))%*%phi_E}
    if (n_LC==0) {quan_LC_Hessian_theta_no_penalty = matrix(0, numIntKnt+2, numIntKnt+2)}
    else {quan_LC_Hessian_theta_no_penalty = - t(cphi_LC * as.numeric(S0ts_LC / S0ts_LC_diff^2))%*%cphi_LC}
    if (n_IC==0) {quan_IC_Hessian_theta_no_penalty = matrix(0, numIntKnt+2, numIntKnt+2)}
    else {quan_IC_Hessian_theta_no_penalty = - t((cphi_IC_L - cphi_IC_R) * as.numeric(S0ts_IC_L*S0ts_IC_R / S0ts_IC_diff^2))%*%(cphi_IC_L - cphi_IC_R)}
    
    Hessian_theta_no_penalty = quan_E_Hessian_theta_no_penalty +
                               quan_LC_Hessian_theta_no_penalty +
                               quan_IC_Hessian_theta_no_penalty
    # Hessian_theta = Hessian_theta_no_penalty - 2*smooth*R # (numIntKnt+2) by (numIntKnt+2)
    
    # full Hessian gamma
    if (!missing(Zmat) & !missing(gamma_initial) & (dim_Zmat[1]>n+1)) {
      if (q==1) {
        if (n_E==0) {dts_E_gamma = matrix(0, 1, q)}
        else {dts_E_gamma = as.numeric(exp(-X_E%*%beta_hat)) * as.matrix(tapply((-Z_E)*as.numeric(exp(-Z_E%*%gamma_hat)*interval_length_E), id_Z_E, FUN = sum))} # n_E by q
        if (n_RC==0) {dts_RC_gamma = matrix(0, 1, q)}
        else {dts_RC_gamma = as.numeric(exp(-X_RC%*%beta_hat)) * as.matrix(tapply((-Z_RC)*as.numeric(exp(-Z_RC%*%gamma_hat)*interval_length_RC), id_Z_RC, FUN = sum))}
        if (n_LC==0) {dts_LC_gamma = matrix(0, 1, q)}
        else {dts_LC_gamma = as.numeric(exp(-X_LC%*%beta_hat)) * as.matrix(tapply((-Z_LC)*as.numeric(exp(-Z_LC%*%gamma_hat)*interval_length_LC), id_Z_LC, FUN = sum))}
        if (n_IC==0) {dts_IC_L_gamma = dts_IC_R_gamma = matrix(0, 1, q)}
        else {
          dts_IC_L_gamma = as.numeric(exp(-X_IC%*%beta_hat)) * as.matrix(tapply((-Z_IC_subset)*as.numeric(exp(-Z_IC_subset%*%gamma_hat)*interval_length_IC_subset), id_Z_IC_subset, FUN = sum))
          dts_IC_R_gamma = as.numeric(exp(-X_IC%*%beta_hat)) * as.matrix(tapply((-Z_IC)*as.numeric(exp(-Z_IC%*%gamma_hat)*interval_length_IC), id_Z_IC, FUN = sum))
        }
      }
      else {
        if (n_E==0) {dts_E_gamma = matrix(0, 1, q)}
        else {dts_E_gamma = as.numeric(exp(-X_E%*%beta_hat)) * apply((-Z_E)*as.numeric(exp(-Z_E%*%gamma_hat)*interval_length_E), 2, function(x) tapply(x, id_Z_E, sum))}
        if (n_RC==0) {dts_RC_gamma = matrix(0, 1, q)}
        else {dts_RC_gamma = as.numeric(exp(-X_RC%*%beta_hat)) * apply((-Z_RC)*as.numeric(exp(-Z_RC%*%gamma_hat)*interval_length_RC), 2, function(x) tapply(x, id_Z_RC, sum))}
        if (n_LC==0) {dts_LC_gamma = matrix(0, 1, q)}
        else {dts_LC_gamma = as.numeric(exp(-X_LC%*%beta_hat)) * apply((-Z_LC)*as.numeric(exp(-Z_LC%*%gamma_hat)*interval_length_LC), 2, function(x) tapply(x, id_Z_LC, sum))}
        if (n_IC==0) {dts_IC_L_gamma = dts_IC_R_gamma = matrix(0, 1, q)}
        else {
          dts_IC_L_gamma = as.numeric(exp(-X_IC%*%beta_hat)) * apply((-Z_IC_subset)*as.numeric(exp(-Z_IC_subset%*%gamma_hat)*interval_length_IC_subset), 2, function(x) tapply(x, id_Z_IC_subset, sum))
          dts_IC_R_gamma = as.numeric(exp(-X_IC%*%beta_hat)) * apply((-Z_IC)*as.numeric(exp(-Z_IC%*%gamma_hat)*interval_length_IC), 2, function(x) tapply(x, id_Z_IC, sum))
        }
      }
      
      if (n_E==0) {quan_E_Hessian_gamma_1 = quan_E_Hessian_gamma_2 = matrix(0, 1, 1)}
      else {
        quan_E_Hessian_gamma_1 = ddh0ts_E / h0tss_E - (dh0ts_E / h0tss_E)^2 - dh0ts_E # n_E by 1
        quan_E_Hessian_gamma_2 = rep((dh0ts_E / h0tss_E - h0ts_E) * exp(-X_E%*%beta_hat), ni_E)*(exp(-Z_E%*%gamma_hat)*interval_length_E)
      }
      if(n_RC==0) {quan_RC_Hessian_gamma_1 = quan_RC_Hessian_gamma_2 = matrix(0, 1, 1)}
      else {
        quan_RC_Hessian_gamma_1 = - dh0ts_RC
        quan_RC_Hessian_gamma_2 = rep(- h0ts_RC * exp(-X_RC%*%beta_hat), ni_RC)*(exp(-Z_RC%*%gamma_hat)*interval_length_RC)
      }
      if(n_LC==0) {quan_LC_Hessian_gamma_1 = quan_LC_Hessian_gamma_2 = matrix(0, 1, 1)}
      else {
        quan_LC_Hessian_gamma_1 = S0ts_LC*(-h0ts_LC^2 + dh0ts_LC) / S0ts_LC_diff - (S0ts_LC*h0ts_LC / S0ts_LC_diff)^2
        quan_LC_Hessian_gamma_2 = rep(S0ts_LC*h0ts_LC / S0ts_LC_diff * exp(-X_LC%*%beta_hat), ni_LC)*(exp(-Z_LC%*%gamma_hat)*interval_length_LC)
      }
      if(n_IC==0) {quan_IC_L_Hessian_gamma_1 = quan_IC_L_Hessian_gamma_2 = quan_IC_L_Hessian_gamma_3 = quan_IC_R_Hessian_gamma_1 = quan_IC_R_Hessian_gamma_2 = quan_IC_R_Hessian_gamma_3 = matrix(0, 1, 1)}
      else {
        quan_IC_L_Hessian_gamma_1 = S0ts_IC_L*(h0ts_IC_L^2 - dh0ts_IC_L) / S0ts_IC_diff - (S0ts_IC_L*h0ts_IC_L / S0ts_IC_diff)^2
        quan_IC_L_Hessian_gamma_2 = rep(- S0ts_IC_L*h0ts_IC_L / S0ts_IC_diff * exp(-X_IC%*%beta_hat), ni_IC_subset)*(exp(-Z_IC_subset%*%gamma_hat)*interval_length_IC_subset)
        quan_IC_L_Hessian_gamma_3 = (S0ts_IC_L*h0ts_IC_L * S0ts_IC_R*h0ts_IC_R) / S0ts_IC_diff^2
        quan_IC_R_Hessian_gamma_1 = S0ts_IC_R*(-h0ts_IC_R^2 + dh0ts_IC_R) / S0ts_IC_diff - (S0ts_IC_R*h0ts_IC_R / S0ts_IC_diff)^2
        quan_IC_R_Hessian_gamma_2 = rep(S0ts_IC_R*h0ts_IC_R / S0ts_IC_diff * exp(-X_IC%*%beta_hat), ni_IC)*(exp(-Z_IC%*%gamma_hat)*interval_length_IC)
        quan_IC_R_Hessian_gamma_3 = (S0ts_IC_R*h0ts_IC_R * S0ts_IC_L*h0ts_IC_L) / S0ts_IC_diff^2
      }
      
      Hessian_gamma = t(dts_E_gamma)%*%(as.numeric(quan_E_Hessian_gamma_1)*dts_E_gamma) +
                      t(-Z_E)%*%(as.numeric(quan_E_Hessian_gamma_2)*(-Z_E)) +
                      t(dts_RC_gamma)%*%(as.numeric(quan_RC_Hessian_gamma_1)*dts_RC_gamma) +
                      t(-Z_RC)%*%(as.numeric(quan_RC_Hessian_gamma_2)*(-Z_RC)) +
                      t(dts_LC_gamma)%*%(as.numeric(quan_LC_Hessian_gamma_1)*dts_LC_gamma) +
                      t(-Z_LC)%*%(as.numeric(quan_LC_Hessian_gamma_2)*(-Z_LC)) +
                      t(dts_IC_L_gamma)%*%(as.numeric(quan_IC_L_Hessian_gamma_1)*dts_IC_L_gamma) +
                      t(-Z_IC_subset)%*%(as.numeric(quan_IC_L_Hessian_gamma_2)*(-Z_IC_subset)) +
                      t(dts_IC_L_gamma)%*%(as.numeric(quan_IC_L_Hessian_gamma_3)*dts_IC_R_gamma) +
                      t(dts_IC_R_gamma)%*%(as.numeric(quan_IC_R_Hessian_gamma_1)*dts_IC_R_gamma) +
                      t(-Z_IC)%*%(as.numeric(quan_IC_R_Hessian_gamma_2)*(-Z_IC)) +
                      t(dts_IC_R_gamma)%*%(as.numeric(quan_IC_R_Hessian_gamma_3)*dts_IC_L_gamma)
      
      # cross Hessian beta gamma
      if (n_E==0) {X_E_N = matrix(0, 1, p)}
      else {X_E_N = X_E[rep(seq_len(n_E), ni_E), , drop = FALSE]}
      if (n_RC==0) {X_RC_N = matrix(0, 1, p)}
      else {X_RC_N = X_RC[rep(seq_len(n_RC), ni_RC), , drop = FALSE]}
      if (n_LC==0) {X_LC_N = matrix(0, 1, p)}
      else {X_LC_N = X_LC[rep(seq_len(n_LC), ni_LC), , drop = FALSE]}
      if (n_IC==0) {X_IC_N = X_IC_subset_N = matrix(0, 1, p)}
      else {
        X_IC_N = X_IC[rep(seq_len(n_IC), ni_IC), , drop = FALSE]
        X_IC_subset_N = X_IC[rep(seq_len(n_IC), ni_IC_subset), , drop = FALSE]
      }
      
      Hessian_beta_gamma = t(dts_E_beta)%*%(as.numeric(quan_E_Hessian_gamma_1)*dts_E_gamma) +
                           t(-X_E_N)%*%(as.numeric(quan_E_Hessian_gamma_2)*(-Z_E)) +
                           t(dts_RC_beta)%*%(as.numeric(quan_RC_Hessian_gamma_1)*dts_RC_gamma) +
                           t(-X_RC_N)%*%(as.numeric(quan_RC_Hessian_gamma_2)*(-Z_RC)) +
                           t(dts_LC_beta)%*%(as.numeric(quan_LC_Hessian_gamma_1)*dts_LC_gamma) +
                           t(-X_LC_N)%*%(as.numeric(quan_LC_Hessian_gamma_2)*(-Z_LC)) +
                           t(dts_IC_L_beta)%*%(as.numeric(quan_IC_L_Hessian_gamma_1)*dts_IC_L_gamma) +
                           t(-X_IC_subset_N)%*%(as.numeric(quan_IC_L_Hessian_gamma_2)*(-Z_IC_subset)) +
                           t(dts_IC_L_beta)%*%(as.numeric(quan_IC_L_Hessian_gamma_3)*dts_IC_R_gamma) +
                           t(dts_IC_R_beta)%*%(as.numeric(quan_IC_R_Hessian_gamma_1)*dts_IC_R_gamma) +
                           t(-X_IC_N)%*%(as.numeric(quan_IC_R_Hessian_gamma_2)*(-Z_IC)) +
                           t(dts_IC_R_beta)%*%(as.numeric(quan_IC_R_Hessian_gamma_3)*dts_IC_L_gamma)
      
      # cross Hessian gamma theta
      Hessian_gamma_theta = t(dts_E_gamma)%*%quan_E_Hessian_beta_theta +
                            t(dts_RC_gamma)%*%quan_RC_Hessian_beta_theta +
                            t(dts_LC_gamma)%*%quan_LC_Hessian_beta_theta +
                            t(dts_IC_L_gamma)%*%quan_IC_L_Hessian_beta_theta +
                            t(dts_IC_R_gamma)%*%quan_IC_R_Hessian_beta_theta
      
      # full Hessian
      F_mat = -rbind(cbind(Hessian_beta, Hessian_beta_gamma, Hessian_beta_theta), 
                     cbind(t(Hessian_beta_gamma), Hessian_gamma, Hessian_gamma_theta), 
                     cbind(t(Hessian_beta_theta), t(Hessian_gamma_theta), Hessian_theta_no_penalty)) # (p+q+numIntKnt+2) by (p+q+numIntKnt+2)
      Q_mat = matrix(0, nrow = p+q+numIntKnt+2, ncol = p+q+numIntKnt+2)
      Q_mat[(p+q+1):(p+q+numIntKnt+2), (p+q+1):(p+q+numIntKnt+2)] = 2*smooth*R
      F_mat_penalised = F_mat + Q_mat
      # F_mat_penalised = -rbind(cbind(Hessian_beta, Hessian_beta_gamma, Hessian_beta_theta), 
      #                          cbind(t(Hessian_beta_gamma), Hessian_gamma, Hessian_gamma_theta), 
      #                          cbind(t(Hessian_beta_theta), t(Hessian_gamma_theta), Hessian_theta)) # (p+q+numIntKnt+2) by (p+q+numIntKnt+2)
      
      # if (any(grad_theta<(-1e-1))) { # if else to deal with 0 active constraints
      if (any(grad_theta<0 & theta_hat<1e-2)) { # 20220213
        # U_index_active_constraints = which(grad_theta<(-1e-1))+p+q
        U_index_active_constraints = which(grad_theta<0 & theta_hat<1e-2)+p+q # 20220213
        # U_index_active_constraints = c(4, 5)
        U_mat = diag(1L, nrow = p+q+numIntKnt+2)[, -U_index_active_constraints]
      }
      else {U_mat = diag(1L, nrow = p+q+numIntKnt+2)}
    }
    else {
      # full Hessian
      F_mat = -rbind(cbind(Hessian_beta, Hessian_beta_theta), 
                     cbind(t(Hessian_beta_theta), Hessian_theta_no_penalty)) # (p+numIntKnt+2) by (p+numIntKnt+2)
      Q_mat = matrix(0, nrow = p+numIntKnt+2, ncol = p+numIntKnt+2)
      Q_mat[(p+1):(p+numIntKnt+2), (p+1):(p+numIntKnt+2)] = 2*smooth*R
      F_mat_penalised = F_mat + Q_mat
      # F_mat_penalised = -rbind(cbind(Hessian_beta, Hessian_beta_theta), 
      #                          cbind(t(Hessian_beta_theta), Hessian_theta)) # (p+numIntKnt+2) by (p+numIntKnt+2)
      # if (any(grad_theta<(-1e-1))) {
      if (any(grad_theta<0 & theta_hat<1e-2)) { # 20220213
        # U_index_active_constraints = which(grad_theta<(-1e-1))+p
        U_index_active_constraints = which(grad_theta<0 & theta_hat<1e-2)+p # 20220213
        U_mat = diag(1L, nrow = p+numIntKnt+2)[, -U_index_active_constraints]
      }
      else {U_mat = diag(1L, nrow = p+numIntKnt+2)}
    }
    
    # variance-covariance matrix calculation
    B_mat_inv = U_mat%*%solve(t(U_mat)%*%F_mat_penalised%*%U_mat)%*%t(U_mat)
    cov_mat = t(B_mat_inv)%*%F_mat%*%B_mat_inv
    if (any(diag(cov_mat)<0)) {cov_mat = B_mat_inv} # 20220208
    # table of estimates
    var_diag = diag(cov_mat)
    asym_var_beta = as.matrix(var_diag[1:p])
    asym_std_beta = sqrt(asym_var_beta)
    Z_score_beta = beta_hat / sqrt(asym_var_beta)
    p_value_beta = 2*pnorm(abs(Z_score_beta), lower.tail = FALSE)
    # significance_beta = rep(TRUE, p)
    # significance_beta[p_value_beta>=0.05] = FALSE
    if (!missing(Zmat) & !missing(gamma_initial) & (dim_Zmat[1]>n+1)) {
      asym_var_gamma = as.matrix(var_diag[(p+1):(p+q)])
      asym_std_gamma = sqrt(asym_var_gamma)
      Z_score_gamma = gamma_hat / sqrt(asym_var_gamma)
      p_value_gamma = 2*pnorm(abs(Z_score_gamma), lower.tail = FALSE)
      # significance_gamma = rep(TRUE, 1)
      # significance_gamma[p_value_gamma>=0.05] = FALSE
      asym_var_theta = as.matrix(var_diag[(p+q+1):(p+q+numIntKnt+2)])
      asym_std_theta = sqrt(asym_var_theta)
      Z_score_theta = theta_hat / sqrt(asym_var_theta)
      p_value_theta = 2*pnorm(abs(Z_score_theta), lower.tail = FALSE)
      # significance_theta = rep(TRUE, numIntKnt+2)
      # significance_theta[p_value_theta>=0.05] = FALSE
      table_coef = rbind(cbind(beta_hat, asym_std_beta, Z_score_beta, p_value_beta), 
                        cbind(gamma_hat, asym_std_gamma, Z_score_gamma, p_value_gamma), 
                        cbind(theta_hat, asym_std_theta, Z_score_theta, p_value_theta))
      rownames(table_coef)[(p+q+1):(p+q+numIntKnt+2)] <- c(sprintf("basis_function_%d", 1:(numIntKnt+2)))
    }
    else {
      asym_var_theta = as.matrix(var_diag[(p+1):(p+numIntKnt+2)])
      asym_std_theta = sqrt(asym_var_theta)
      Z_score_theta = theta_hat / sqrt(asym_var_theta)
      p_value_theta = 2*pnorm(abs(Z_score_theta), lower.tail = FALSE)
      # significance_theta = rep(TRUE, numIntKnt+2)
      # significance_theta[p_value_theta>=0.05] = FALSE
      table_coef = rbind(cbind(beta_hat, asym_std_beta, Z_score_beta, p_value_beta), 
                        cbind(theta_hat, asym_std_theta, Z_score_theta, p_value_theta))
      rownames(table_coef)[(p+1):(p+numIntKnt+2)] <- c(sprintf("basis_function_%d", 1:(numIntKnt+2)))
    }
    colnames(table_coef) <- c("coef", "se(coef)", "z", "p(>|z|)")
    # colnames(table_coef) <- c("est", "asym_std", "Z", "p", "significant? (<0.05)")
    
    # browser()
    df = numIntKnt+2 - sum(diag(solve(t(U_mat)%*%F_mat_penalised%*%U_mat)%*%(t(U_mat)%*%Q_mat%*%U_mat)))
    if (df < 0) {df = 1}
    smooth_new = as.numeric(df / (2 * t(theta_hat)%*%R%*%theta_hat))
    
    # # Asymptotics
    # dbeta=grad
    # dtheta=nume-deno
    # 
    # vv=as.numeric(exp(colMeans(Xmat)%*%beta_hat))
    # if (any(abs(theta_hat)<1e-5)) {
    #   id0theta=intersect(which(abs(theta_hat)<1e-5), which(dbeta<1e-5))
    #   # id0theta=which(abs(theta_hat)<1e-5)
    # }
    # else{id0theta = c()}
    # idtheta=setdiff(1:numSp, id0theta)
    # 
    # 
    # ddphi= mSpline(ts, knots=IntKnt, degree=dgrSp, intercept=TRUE, Boundary.knots=bryKnt, derivs=2)
    # ddh0ts=as.matrix(ddphi%*%theta_hat)
    # # 2nd derivative beta, p by p matrix
    # pl_beta=t(X)%*% ((( censor_n*((ddh0ts*ts^2+dh0ts*ts)/h0tss-(dh0ts*ts/h0tss)^2) - (ddh0ts*ts^2+dh0ts*ts) )%*%matrix(1, 1, p))*X)
    # # 2nd derivative theta, m by m matrix
    # l_theta=-t(h0tsM*phi)%*%((1/h0tss)%*%matrix(1, 1, numSp)*phi)
    # l_theta_s=l_theta[idtheta, idtheta]
    # pl_theta=l_theta-smooth*R
    # pl_theta_s=pl_theta[idtheta, idtheta]
    # # 2nd derivative beta and theta, p by m matrix
    # pl_beta_theta= -t(X)%*% ( (censor_n*ts/h0tss)%*%matrix(1, 1, numSp)*dphi - (censor_n*ts*dh0ts/h0tss^2)%*%matrix(1, 1, numSp)*phi  -  (ts%*%matrix(1, 1, numSp))*phi)
    # pl_beta_theta_s = pl_beta_theta[, idtheta]
    # 
    # # F-1G
    # # d2PLs=rbind(cbind(pl_beta, pl_beta_theta_s), cbind(t(pl_beta_theta_s), pl_theta_s))
    # # dfHs=solve(d2PLs, rbind(cbind(pl_beta, pl_beta_theta_s), cbind(t(pl_beta_theta_s), l_theta_s)))
    # # Asym_V_s = solve(d2PLs, dfHs);
    # # F_matrix_inv=solve(-rbind(cbind(pl_beta, pl_beta_theta), cbind(t(pl_beta_theta), pl_theta)) )
    # # G_matrix=-rbind(cbind(pl_beta, pl_beta_theta), cbind(t(pl_beta_theta), l_theta))
    # # Asym_V=F_matrix_inv%*%G_matrix%*%t(F_matrix_inv)
    # F_matrix_inv=solve(rbind(cbind(pl_beta, pl_beta_theta_s), cbind(t(pl_beta_theta_s), pl_theta_s)) )
    # G_matrix=-rbind(cbind(pl_beta, pl_beta_theta_s), cbind(t(pl_beta_theta_s), l_theta_s))
    # Asym_V_s=F_matrix_inv%*%G_matrix%*%t(F_matrix_inv)
    # # browser()
    # B11=Asym_V_s[1:p, 1:p]
    # B12=matrix(0, p, numSp)
    # B12[, idtheta]=Asym_V_s[1:p, (p+1):(p+length(idtheta))]
    # B22=matrix(0, numSp, numSp)
    # B22[idtheta, idtheta]=Asym_V_s[(p+1):(p+length(idtheta)), (p+1):(p+length(idtheta))]
    # 
    # Asym_V=rbind(cbind(B11, B12), cbind(t(B12), B22))
    # Asym_V_diag=diag(Asym_V)
    # 
    # Asym_V_beta=Asym_V_diag[1:p]
    # Asym_V_theta=Asym_V_diag[-(1:p)]
    # # Asym_v_h0=vv^2*diag(phi%*%diag(Asym_V_theta)%*%t(phi))
    # 
    # ASD_beta=sqrt(Asym_V_beta)
    # ASD_theta=sqrt(Asym_V_theta)
    # 
    # 
    # # S0=exp(-ch0ts*vv)
    # 
    # time_axis=seq(min(time), max(time), length.out=500)
    # phi_axis=mSpline(time_axis, knots = IntKnt, degree = dgrSp, intercept=TRUE)
    # Asym_v_h0=vv^2*diag(phi_axis%*%diag(Asym_V_theta)%*%t(phi_axis))
    # # Asym_v_h0=diag(phi_axis%*%diag(Asym_V_theta)%*%t(phi_axis))
    # 
    # ASD_h0=sqrt(Asym_v_h0)
    # h0_hat_plot=vv*phi_axis%*%theta_hat
    # h0_hat_plot=phi_axis%*%theta_hat
    # h0_hat_plot_a_upper=h0_hat_plot+1.96*ASD_h0
    # h0_hat_plot_a_lower=h0_hat_plot-1.96*ASD_h0
    # 
    # plot(time_axis, h0_hat_plot, type="l", main=bquote(paste("MPL approach with ", lambda==.(smooth))), xlab="Time", ylab="Estimated baseline hazard")
    # lines(time_axis, h0_hat_plot_a_upper, type="l", col="blue")
    # lines(time_axis, h0_hat_plot_a_lower, type="l", col="red")
    
    
    smooth_hat = smooth_new
    if (abs(df - df_old)<1) {
      break
    }
    else {
      smooth = smooth_new
      df_old = df
      beta_initial = beta_hat
      if (!missing(Zmat) & !missing(gamma_initial) & (dim_Zmat[1]>n+1)) {gamma_initial = gamma_hat}
    }
  }
  smooth_iter = smooth_iter[1:j, ]
  
  
  # browser()
  ################################# plots #################################
  if (draw_plot==TRUE) {
    # hazard plot
    # plots of final accelerated time
    # plot_hazard = cbind(ts, Gaussian_basis(x = ts, mean = gknots, sd = sig)%*%theta_hat)
    # plot_hazard = plot_hazard[order(plot_hazard[, 1]), ] # sort according to ascending order of ts
    # plot(plot_hazard[, 1], plot_hazard[, 2], xlim = c(0, max(plot_hazard[, 1])), ylim = c(0, max(plot_hazard[, 2])), main = "hazard plot", xlab = "accelerated/decelerated time", ylab = "hazard")
    plot(ts_bind, Gaussian_basis(x = ts_bind, mean = gknots, sd = sig)%*%theta_hat, xlim = c(0, max(ts_bind)), main = "hazard plot", xlab = "accelerated/decelerated time", ylab = "hazard")
    # smooth curve connecting data points
    lines(seq(0, max(ts_bind), length.out=max(1000, 10*n)), Gaussian_basis(x = as.matrix(seq(0, max(ts_bind), length.out=max(1000, 10*n))), mean = gknots, sd = sig)%*%theta_hat)
    lines(seq(0, max(ts_bind), length.out=max(1000, 10*n)), 3*(seq(0, max(ts_bind), length.out=max(1000, 10*n)))^2, lty=2)
    # plot(time_bind, Gaussian_basis(x = ts_bind, mean = gknots, sd = sig)%*%theta_hat, xlim = c(0, max(time_bind)), main = "hazard plot", xlab = "time", ylab = "hazard")
    # lines(seq(0, max(time_bind), length.out=max(1000, 10*n)), Gaussian_basis(x = seq(0, max(ts_bind), length.out=max(1000, 10*n)), mean = gknots, sd = sig)%*%theta_hat)
    
    # survival plot
    # plots of final accelerated time
    # plot_survival = cbind(ts, exp(-Gaussian_basis_integ1(q = ts, mean = gknots, sd = sig)%*%theta_hat))
    # plot_survival = plot_survival[order(plot_survival[, 1]), ] # sort according to ascending order of ts
    # plot(plot_survival[, 1], plot_survival[, 2], xlim = c(0, max(plot_hazard[, 1])), ylim = c(0, 1), main = "Survival plot", xlab = "accelerated/decelerated time", ylab = "Survival")
    plot(ts_bind, exp(-Gaussian_basis_integ1(q = ts_bind, mean = gknots, sd = sig)%*%theta_hat), xlim = c(0, max(ts_bind)), ylim = c(0, 1), main = "Survival plot", xlab = "accelerated/decelerated time", ylab = "Survival")
    # smooth curve connecting data points
    lines(seq(0, max(ts_bind), length.out=max(1000, 10*n)), exp(-Gaussian_basis_integ1(q = as.matrix(seq(0, max(ts_bind), length.out=max(1000, 10*n))), mean = gknots, sd = sig)%*%theta_hat))
    
    # below plots for original time may be incorrect
    if (knots_option=="equal_space") {
      bryKnt_0 = range(time_bind)
      bin_width_0 = (bryKnt_0[2] - bryKnt_0[1]) / (numIntKnt+1)
      gknots_0 = seq(bryKnt_0[1], bryKnt_0[2], bin_width_0)
      sig_0 = (2/3) * bin_width_0
    }
    else if (knots_option=="percentile") {
      gknots_0 = quantile(sort(time_bind), seq(knots_min_quantile, knots_max_quantile, length.out = numIntKnt+2)/100, type=1)
      dist_0 = sapply(gknots_0, function(x) {abs(x - time_bind)})
      sig_0 = apply(dist_0, 2, function(x) {quantile(x, sig_coverage) / 2})
    }
    else if (knots_option=="percentile_1") { # test on 20220525
      time_bind_1 = rbind(time_E, time_RC, time_LC, time_IC_R)
      time_bind_1 = time_bind_1[!is.na(time_bind_1), , drop = FALSE]
      time_bind_1 = time_bind_1[order(as.numeric(rownames(time_bind_1))), , drop = FALSE]
      gknots_0 = quantile(sort(time_bind_1), seq(knots_min_quantile, knots_max_quantile, length.out = numIntKnt+2)/100, type=1)
      dist_0 = sapply(gknots_0, function(x) {abs(x - time_bind_1)})
      sig_0 = apply(dist_0, 2, function(x) {quantile(x, sig_coverage) / 2})
    }
    # baseline hazard plot
    # plots of original time
    # plot_hazard_0 = cbind(time, Gaussian_basis(x = time, mean = gknots_0, sd = sig_0)%*%theta_hat)
    # plot_hazard_0 = plot_hazard_0[order(plot_hazard_0[, 1]), ] # sort according to ascending order of time
    # plot(plot_hazard_0[, 1], plot_hazard_0[, 2], xlim = c(0, max(plot_hazard_0[, 1])), ylim = c(0, max(plot_hazard_0[, 2])), main = "baseline hazard plot", xlab = "time", ylab = "baseline hazard")
    plot(time_bind, Gaussian_basis(x = time_bind, mean = gknots_0, sd = sig_0)%*%theta_hat, xlim = c(0, max(time_bind)), main = "baseline hazard plot", xlab = "time", ylab = "baseline hazard",ylim=c(0,1))
    # smooth curve connecting data points
    lines(seq(0, max(time_bind), length.out=max(1000, 10*n)), Gaussian_basis(x = as.matrix(seq(0, max(time_bind), length.out=max(1000, 10*n))), mean = gknots_0, sd = sig_0)%*%theta_hat)
    # lines(seq(0, max(time_bind), length.out=max(1000, 10*n)), 3*(seq(0, max(time_bind), length.out=max(1000, 10*n)))^2)
    lines(seq(0, max(time_bind), length.out=max(1000, 10*n)), 3*(seq(0, max(time_bind), length.out=max(1000, 10*n)))^2, lty=2)
    # baseline survival plot
    # plots of original time
    # plot_survival_0 = cbind(time, exp(-Gaussian_basis_integ1(q = time, mean = gknots_0, sd = sig_0)%*%theta_hat))
    # plot_survival_0 = plot_survival_0[order(plot_survival_0[, 1]), ] # sort according to ascending order of time
    # plot(plot_survival_0[, 1], plot_survival_0[, 2], xlim = c(0, max(plot_hazard_0[, 1])), ylim = c(0, 1), main = "baseline Survival plot", xlab = "time", ylab = "baseline Survival")
    plot(time_bind, exp(-Gaussian_basis_integ1(q = time_bind, mean = gknots_0, sd = sig_0)%*%theta_hat), xlim = c(0, max(time_bind)), ylim = c(0, 1), main = "baseline Survival plot", xlab = "time", ylab = "baseline Survival")
    # smooth curve connecting data points
    lines(seq(0, max(time_bind), length.out=max(1000, 10*n)), exp(-Gaussian_basis_integ1(q = as.matrix(seq(0, max(time_bind), length.out=max(1000, 10*n))), mean = gknots_0, sd = sig_0)%*%theta_hat))
  }
  
  # browser()
  ############################### Output ###############################
  if (!missing(Zmat) & !missing(gamma_initial) & (dim_Zmat[1]>n+1)) {
    return(list(beta_hat=beta_hat, gamma_hat=gamma_hat, theta_hat=theta_hat, 
                asym_var_beta=asym_var_beta, asym_var_gamma=asym_var_gamma, asym_var_theta=asym_var_theta, 
                asym_std_beta=asym_std_beta, asym_std_gamma=asym_std_gamma, asym_std_theta=asym_std_theta, 
                Z_score_beta=Z_score_beta, Z_score_gamma=Z_score_gamma, Z_score_theta=Z_score_theta, 
                p_value_beta=p_value_beta, p_value_gamma=p_value_gamma, p_value_theta=p_value_theta, 
                table_coef=table_coef, 
                smooth_hat=smooth_hat, smooth_iter=smooth_iter, 
                cvg=cvg, time_range=time_range, multiple_ts_time=multiple_ts_time, knots_iter=knots_iter, 
                ts_range_beta_iter=ts_range_beta_iter, ts_range_gamma_iter=ts_range_gamma_iter, ts_range_theta_iter=ts_range_theta_iter, 
                iteration_no=iteration_no, 
                grad_beta_iter=grad_beta_iter, grad_gamma_iter=grad_gamma_iter, grad_theta_iter=grad_theta_iter, multiplier_iter = multiplier_iter, deno_iter = deno_iter,
                Hes_beta_iter=Hes_beta_iter, det_Hes_beta_iter=det_Hes_beta_iter, eigenvl_Hes_beta_iter=eigenvl_Hes_beta_iter, 
                Hes_gamma_iter=Hes_gamma_iter, det_Hes_gamma_iter=det_Hes_gamma_iter, eigenvl_Hes_gamma_iter=eigenvl_Hes_gamma_iter, 
                Hes_beta_final=Hes_beta, Hes_gamma_final=Hes_gamma, 
                beta_iter=beta_iter, gamma_iter=gamma_iter, theta_iter=theta_iter)) #20181129
  }
  else {
    return(list(beta_hat=beta_hat, theta_hat=theta_hat, 
                asym_var_beta=asym_var_beta, asym_var_theta=asym_var_theta, 
                asym_std_beta=asym_std_beta, asym_std_theta=asym_std_theta, 
                Z_score_beta=Z_score_beta, Z_score_theta=Z_score_theta, 
                p_value_beta=p_value_beta, p_value_theta=p_value_theta, 
                table_coef=table_coef, 
                smooth_hat=smooth_hat, smooth_iter=smooth_iter, 
                cvg=cvg, time_range=time_range, multiple_ts_time=multiple_ts_time, knots_iter=knots_iter, 
                ts_range_beta_iter=ts_range_beta_iter, ts_range_theta_iter=ts_range_theta_iter, 
                iteration_no=iteration_no, 
                grad_beta_iter=grad_beta_iter, grad_theta_iter=grad_theta_iter, multiplier_iter = multiplier_iter, deno_iter = deno_iter,
                Hes_beta_iter=Hes_beta_iter, det_Hes_beta_iter=det_Hes_beta_iter, eigenvl_Hes_beta_iter=eigenvl_Hes_beta_iter, 
                Hes_beta_final=Hes_beta, 
                beta_iter=beta_iter, theta_iter=theta_iter)) #20181129
  }
  # return(list(beta=beta_hat, theta=theta_hat, beta_ASD=ASD_beta, h0_ASD=ASD_h0, h0=h0_hat_plot, theta_ASD=ASD_theta, h0_lower=h0_hat_plot_a_lower, h0_upper=h0_hat_plot_a_upper, beta_gradient=dbeta, theta_gradient=dtheta, penlike=penlike_theta_iter_after_ls))  #, "tstar"=t))
  # message("No. of iterations to convergence is ", k)
}
