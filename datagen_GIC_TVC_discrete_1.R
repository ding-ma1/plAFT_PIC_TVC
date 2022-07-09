# forked from datagen_RC_TVC_discrete_1.R on 20220406
# a fix of datagen_RC_TVC_discrete.R 20210915
datagen_GIC_TVC_discrete_1 <- function(n = 100, beta = c(1, -1), gamma = -0.1, t_treat_option = 1, t_treat_min = 0, t_treat_max = 0.5, pi_E, term_L, term_R, datagen_option = "GIC_TVC") {
  # 20201110 updated to incorporate a separate regulator for t_treat, i.e. t_treat_min != c_min and t_treat_max != c_max
  # if(missing(c_max)) {stop("Please specify a maximum value (c_max) in order to generate right-censored time")}
  
  #################### generate time-fixed covariate matrix Xmat ####################
  id_X = 1:n
  X_1 = rbinom(n, 1, 0.5)
  # X_2 = rlnorm(n)
  # X_2 = runif(n, min = 0, max = 3)
  X_2 = runif(n, min = 0, max = 1) # 20220530
  
  # X = cbind(X_1, X_2, X_3)
  X = cbind(X_1, X_2)
  Xmat = cbind(id_X, X)
  
  ############################ generate time matrix tmat ############################
  S = runif(n, min = 0, max = 1) # survival
  
  if (datagen_option=="GIC_TVC") {
    if (t_treat_option == 1) {t_treat = runif(n, min = t_treat_min, max = t_treat_max)} # treatment time (uniform distribution)
    else if (t_treat_option == 2) {t_treat = rep(t_treat_max, n)} # treatment time (all same and fixed)
    # # slower version of calculating t using for loop
    # # Z_i(t) == 0 if t < t_treat[i], Z_i(t) == 1 if t >= t_treat[i]
    # t = rep(NA, times = n)
    # for (i in 1:n) {
    #   if (-log(S[i]) < (exp(-X[i, ]%*%beta) * t_treat[i])^3) {t[i] = exp(X[i, ]%*%beta) * (-log(S[i]))^(1/3)}
    #   else {t[i] = exp(X[i, ]%*%beta + gamma) * (-log(S[i]) - (exp(-X[i, ]%*%beta) * t_treat[i])^3 + (exp(-X[i, ]%*%beta - gamma) * t_treat[i])^3)^(1/3)}
    # }
    # faster alternative of calculating t using matrix multiplication
    less_or_more = (-log(S) < (exp(-X%*%beta) * t_treat)^3)
  }
  
  t = exp(X%*%beta) * (-log(S))^(1/3)
  
  if (datagen_option=="GIC_TVC") {
    # # Jun's version
    # t_alternative = exp(X%*%beta + matrix(1, n, 1)%*%gamma) * (-log(S) - (exp(-X%*%beta) * t_treat)^3 + (exp(-X%*%beta - matrix(1, n, 1)%*%gamma) * t_treat)^3)^(1/3)
    # alternative version
    t_alternative = exp(X%*%beta + matrix(1, n, 1)%*%gamma) * (-log(S))^(1/3) - exp(matrix(1, n, 1)%*%gamma) * t_treat + t_treat
    t[less_or_more==FALSE] = t_alternative[less_or_more==FALSE]
  }
  t = as.numeric(t)
  
  # below changed to GIC version
  p_E = runif(n, min = 0, max = 1)
  
  p_L = runif(n, min = 0, max = 1)
  p_R = runif(n, min = p_L, max = 1)
  c_L = p_L * term_L
  c_R = p_R * term_R
  
  y_L = matrix(NA, n, 1)
  y_R = matrix(NA, n, 1)
  indicator_n = matrix(NA, n, 1)
  # event time
  event_cases = (p_E < pi_E) # TRUE for event, FALSE for non-event
  y_L[event_cases] = t[event_cases]
  y_R[event_cases] = t[event_cases]
  indicator_n[event_cases] = 1L
  e_prop = sum(event_cases) / n
  # right-censored time
  rc_cases = (p_E >= pi_E & t > c_R) # TRUE for rc, FALSE for non-rc
  y_L[rc_cases] = c_R[rc_cases]
  y_R[rc_cases] = NA
  indicator_n[rc_cases] = 2L
  rc_prop = sum(rc_cases) / n
  # left-censored time
  lc_cases = (p_E >= pi_E & t < c_L) # TRUE for lc, FALSE for non-lc
  y_L[lc_cases] = 0
  y_R[lc_cases] = c_L[lc_cases]
  indicator_n[lc_cases] = 3L
  lc_prop = sum(lc_cases) / n
  # interval-censored time
  ic_cases = (p_E >= pi_E & t >= c_L & t <= c_R) # TRUE for ic, FALSE for non-ic
  y_L[ic_cases] = c_L[ic_cases]
  y_R[ic_cases] = c_R[ic_cases]
  indicator_n[ic_cases] = 4L
  ic_prop = sum(ic_cases) / n
  
  censor_prop = rc_prop + lc_prop + ic_prop
  tmat = cbind(y_L, y_R, indicator_n)
  
  range_t = range(t)
  range_y = range(y_L, y_R, na.rm = TRUE)
  
  # y_upperbounds = as.numeric(pmax(y_L, y_R, na.rm = TRUE))
  if (datagen_option=="GIC_TVC") {
    
    # { # this doesn't work properly, must contain the information of t_i_L in the Zmat
    #   id_Z = c(rbind(id_X, id_X))
    #   interval_index = rep(1:2, times = n)
    #   Z = rep(c(0, 1), times = n)
    #   
    #   Zmat_deletion_mark_interval_1 = numeric(n)
    #   Zmat_deletion_mark_interval_2 = numeric(n)
    #   # if 0 < y_upperbounds <= t_treat, need to mark the rows of second intervals to delete later
    #   Zmat_deletion_mark_interval_2[y_upperbounds>0 & y_upperbounds<=t_treat] = 1
    #   # if 0 < y_upperbounds <= t_treat, need to replace the relevant t_treat's by y_upperbounds
    #   t_treat[y_upperbounds>0 & y_upperbounds<=t_treat] = y_upperbounds[y_upperbounds>0 & y_upperbounds<=t_treat]
    #   
    #   Zmat_deletion_interval_indicator = c(rbind(Zmat_deletion_mark_interval_1, Zmat_deletion_mark_interval_2))
    #   interval_L = c(rbind(numeric(n), t_treat))
    #   interval_R = c(rbind(t_treat, y_upperbounds))
    #   
    #   Zmat = cbind(id_Z, interval_index, interval_L, interval_R, Z)
    #   Zmat = Zmat[Zmat_deletion_interval_indicator==0, ]
    #   
    #   indicator_N = rep(indicator_n, times = tapply(Zmat[, 2], Zmat[, 1], FUN = max))
    # }
    
    # older and slower, but works as expected
    Zmat = NULL
    indicator_N = NULL
    for (i in 1:n) {
      # event and right-censored cases, two situations
      if (indicator_n[i] == 1 | indicator_n[i] == 2) {
        if (t_treat[i] < y_L[i]) {
          Zmat = rbind(Zmat, c(i, 1, 0, t_treat[i], 0), c(i, 2, t_treat[i], y_L[i], 1))
          indicator_N = c(indicator_N, indicator_n[i], indicator_n[i])
        }
        else {
          Zmat = rbind(Zmat, c(i, 1, 0, y_L[i], 0))
          indicator_N = c(indicator_N, indicator_n[i])
        }
      }
      # left-censored cases, two situations
      else if (indicator_n[i] == 3) {
        if (t_treat[i] < y_R[i]) {
          Zmat = rbind(Zmat, c(i, 1, 0, t_treat[i], 0), c(i, 2, t_treat[i], y_R[i], 1))
          indicator_N = c(indicator_N, indicator_n[i], indicator_n[i])
        }
        else {
          Zmat = rbind(Zmat, c(i, 1, 0, y_R[i], 0))
          indicator_N = c(indicator_N, indicator_n[i])
        }
      }
      # interval-censored cases, four situations
      else {
        if (t_treat[i] < y_L[i]) {
          Zmat = rbind(Zmat, c(i, 1, 0, t_treat[i], 0), c(i, 2, t_treat[i], y_L[i], 1), c(i, 3, y_L[i], y_R[i], 1))
          indicator_N = c(indicator_N, indicator_n[i], indicator_n[i], indicator_n[i])
        }
        else if (t_treat[i] == y_L[i]) {
          Zmat = rbind(Zmat, c(i, 1, 0, t_treat[i], 0), c(i, 2, t_treat[i], y_R[i], 1))
          indicator_N = c(indicator_N, indicator_n[i], indicator_n[i])
        }
        else if (t_treat[i] < y_R[i]) {
          # Zmat = rbind(Zmat, c(i, 1, 0, y_L[i], 0), c(i, 2, y_L[i], t_treat[i], 1), c(i, 3, t_treat[i], y_R[i], 1))
          Zmat = rbind(Zmat, c(i, 1, 0, y_L[i], 0), c(i, 2, y_L[i], t_treat[i], 0), c(i, 3, t_treat[i], y_R[i], 1)) # corrected on 20220519, stuck for more than 2 weeks
          indicator_N = c(indicator_N, indicator_n[i], indicator_n[i], indicator_n[i])
        }
        else {
          Zmat = rbind(Zmat, c(i, 1, 0, y_L[i], 0), c(i, 2, y_L[i], y_R[i], 0))
          indicator_N = c(indicator_N, indicator_n[i], indicator_n[i])
        }
      }
    }
    colnames(Zmat) <- c("id_Z","interval_index","interval_L","interval_R","Z(t)")
    # browser()
    treat_prop = sum(Zmat[, 2]==2)/n
    event_treated_prop = sum(indicator_N==1 & Zmat[, 5]==1) / n
    
    if (event_treated_prop==0) {warning("The generated dataset doesn't contain any event cases with time varying covariates! Please rerun the code and generate a different dataset!\n")}
    return(list(tmat=tmat, Xmat=Xmat, Zmat=Zmat, indicator_N=indicator_N, range_t=range_t, range_y=range_y, e_prop=e_prop, rc_prop=rc_prop, lc_prop=lc_prop, ic_prop=ic_prop, censor_prop=censor_prop, treat_prop=treat_prop, event_treated_prop=event_treated_prop))
  }
  else {
    return(list(tmat=tmat, Xmat=Xmat, range_t=range_t, range_y=range_y, e_prop=e_prop, rc_prop=rc_prop, lc_prop=lc_prop, ic_prop=ic_prop, censor_prop=censor_prop))
  }
}
