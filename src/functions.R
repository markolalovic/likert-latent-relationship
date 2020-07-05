##########################################################################################
### PACKAGES

# only for plots and tables
# sn (skew normal package) is not necessary (no faster than using dnorm, pnorm and also sn has problems with -Inf and +Inf)
ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) {
    #cat("installing package ", new.pkg)
    install.packages(new.pkg, dependencies = TRUE)
  }
  sapply(pkg, require, character.only = TRUE)
}
packages <- c("xtable", "gridExtra", "ggplot2", "latex2exp", "foreach", "sn")
ipak(packages)

### constants
labels_K_5 <- c(TeX("r_{1}"), TeX("r_{2}"), TeX("r_{3}"), TeX("r_{4}"), TeX("r_{5}"))
labels_K_10 <- c(TeX("r_{1}"), TeX("r_{2}"), TeX("r_{3}"), TeX("r_{4}"), TeX("r_{5}"), TeX("r_{6}"), TeX("r_{7}"), TeX("r_{8}"), TeX("r_{9}"), TeX("r_{10}"))

##########################################################################################
### CORE FUNCTIONS

define_problem <- function(alpha=2,
                          delta=0.2,
                          lambda=0.5,
                          eta=0,
                          disc_type="OD",
                          K=5,
                          lbound=-3, ubound=3,
                          trace=FALSE,
                          text_size=20,
                          dashed_size=1.2,
                          r_increase=1.2) {

  ##############################################################################################################
  # DEF_PROBLEM
  # disc_type \in {OD, EW}
  # alpha \in [-10, 10]
  # [lbound, ubound] = domain for plot of latent variable (e.g. -3, 3)
  # delta = difference in means betweeen X1 and X2
  # lambda = ratio of standard deviations of X1 and X2
  # eta = parameter of gaussian noise on thresholds xk

  # SETUP:
  # X ~ SN(alpha);              E[X] = mu;      sd(X) = sigma (e.g. 0, 1 for alpha = 0)
  # X_0 ~ X - mu;               E[X_0] = 0;     sd(X_0) = sigma           # hypothetical neutral group
  # X_1 ~ omega_1 * X_0 + xi_1; E[X_1] = xi_1;  sd(X_1) = omega_1 * sigma # group 1 (e. g. guys)
  # X_2 ~ omega_2 * X_0 + xi_2; E[X_2] = xi_2;  sd(X_2) = omega_2 * sigma # group 2 (e. g. girls)

  # such that:
  # xi_2 = xi_1 + delta * sigma # we multiply delta by sigma so results are comparable
  # omega_2 = lambda * omega_1

  ##############################################################################################################
  # X
  mu <- meanSN(alpha)
  sigma <- sqrt(varSN(alpha))


  # xi_1 and omega_1 are arbitrary, we choose them so that it's similar to real world problem we presented
  xi_1 <- -1*sigma
  omega_1 <- 1

  ##############################################################################################################
  # X_0 AKA hypothetical neutral group
  # calculate thresholds for X_0
  sl <- simulate_Likert(disc_type, K, 0, 1, alpha)

  # add noise to thresholds
  if (eta == 0) { # no noise
    xk <- c(-Inf, sl$xk, Inf)
    xkk <- round(xk, 2)
  } else {
    xk <- sl$xk + rnorm(K-1, 0, eta)
    xk <- c(-Inf, sort(xk), Inf) # we sort thresholds after adding noise, because it can happen that order get's mixed up (if noise is very high)
    xkk <- round(xk, 2)
  }

  fX_0 <- function(x) { dsn(x + mu, 0, 1, alpha) } # shift for mu
  pk_0 <- get_pk_of_xk(xk[2:K], fX_0)
  mu_X0_hat <- get_mean_of_pk(pk_0)
  sd_X0_hat <- sqrt(get_var_of_pk(pk_0))

  if (trace) {
    x <- seq(lbound, ubound, length.out = 1000)
    p10 <- ggplot(data.frame(x=x, y=fX_0(x)), aes(x = x, y = y)) +
      geom_line(col="blue") +
      xlab("x") + ylab("Density") +
      theme(text = element_text(size=text_size),
            axis.text.x = element_text(size=text_size),
            axis.text.y = element_text(size=text_size))
    p10 <- p10 + ggtitle(TeX("$X_{0} \\sim N(0, 1)$"))
    for (k in 2:K) {
      p10 <- p10 + geom_vline(xintercept = xkk[k], col="red", linetype="dashed", size=dashed_size)
    }
    p10 <- p10 + geom_area(mapping = aes(y = ifelse(x > xkk[2] + 10^-2 & x <= xkk[3] - 10^-3, y, 0)), fill = "grey")
    p10 <- p10 + geom_vline(xintercept = 0, col="blue", linetype="dashed", size=dashed_size)

    p20 <- ggplot(data.frame(k = 1:K, prob=pk_0), aes(x = k)) +
      geom_point(aes(y=prob), col="blue") +
      geom_linerange(aes(ymax=prob), ymin=0, col="blue") +
      ylab('Probability') + ylim(0, (max(pk_0))) +
      scale_x_continuous(breaks = 1:K, labels = labels_K_5) +
      theme(text = element_text(size=text_size),
            axis.title.x=element_blank(),
            axis.text.y = element_text(size=text_size),
            axis.text.x = element_text(size=r_increase * text_size))
    p20 <- p20 + ggtitle(TeX("$\\hat{X}_{0}$"))
    p20 <- p20 + geom_vline(xintercept = round(mu_X0_hat, 2), col="blue", linetype="dashed", size=dashed_size)
  }

  ##############################################################################################################
  # X_1 AKA group 1
  fX_1 <- function(x) { fX_0((x - xi_1)/omega_1)*(1/omega_1) } # shift and scale X_0
  pk_1 <- get_pk_of_xk(xk[2:K], fX_1)
  mu_X1_hat <- get_mean_of_pk(pk_1)
  sd_X1_hat <- sqrt(get_var_of_pk(pk_1))

  if (trace) {
    p11 <- ggplot(data.frame(x=x, y=fX_1(x)), aes(x = x, y = y)) +
      geom_line(col="blue") +
      xlab("x") + ylab("Density") +
      theme(text = element_text(size=text_size),
            axis.text.x = element_text(size=text_size),
            axis.text.y = element_text(size=text_size))
    p11 <- p11 + ggtitle(TeX("$X_{1} \\sim X_{0} - 1$"))
    for (k in 2:K) {
      p11 <- p11 + geom_vline(xintercept = xkk[k], col="red", linetype="dashed", size=dashed_size)
    }
    p11 <- p11 + geom_area(mapping = aes(y = ifelse(x > xkk[2] + 10^-2 & x <= xkk[3] - 10^-3, y, 0)), fill = "grey")
    p11 <- p11 + geom_vline(xintercept = xi_1, col="blue", linetype="dashed", size=dashed_size)

    p21 <- ggplot(data.frame(k = 1:K, prob=pk_1), aes(x = k)) +
      geom_point(aes(y=prob), col="blue") +
      geom_linerange(aes(ymax=prob), ymin=0, col="blue") +
      ylab('Probability') + ylim(0, (max(pk_1))) +
      scale_x_continuous(breaks = 1:K, labels = labels_K_5) +
      theme(text = element_text(size=text_size),
            axis.title.x=element_blank(),
            axis.text.y = element_text(size=text_size),
            axis.text.x = element_text(size=r_increase * text_size))
    p21 <- p21 + ggtitle(TeX("$\\hat{X}_{1}$"))
    p21 <- p21 + geom_vline(xintercept = round(mu_X1_hat, 2), col="blue", linetype="dashed", size=dashed_size)
  }

  ##############################################################################################################
  # X_2 AKA group 2
  xi_2 <- xi_1 + delta*sigma
  omega_2 <- omega_1 * lambda
  fX_2 <- function(x) { fX_0((x - xi_2)/omega_2)*(1/omega_2) } # scale and shift X_0
  pk_2 <- get_pk_of_xk(xk[2:K], fX_2)
  mu_X2_hat <- get_mean_of_pk(pk_2)
  sd_X2_hat <- sqrt(get_var_of_pk(pk_2))

  if (trace) {
    p12 <- ggplot(data.frame(x=x, y=fX_2(x)), aes(x = x, y = y)) +
      geom_line(col="blue") +
      xlab("x") + ylab("Density") +
      theme(text = element_text(size=text_size),
            axis.text.x = element_text(size=text_size),
            axis.text.y = element_text(size=text_size))
    p12 <- p12 + ggtitle(TeX("$X_{2} \\sim 0.5 \\cdot X_{0} - 1$"))
    for (k in 2:K) {
      p12 <- p12 + geom_vline(xintercept = xkk[k], col="red", linetype="dashed", size=dashed_size)
    }
    p12 <- p12 + geom_area(mapping = aes(y = ifelse(x > xkk[2] + 10^-2 & x <= xkk[3] - 10^-3, y, 0)), fill = "grey")
    p12 <- p12 + geom_vline(xintercept = xi_2, col="blue", linetype="dashed", size=dashed_size)

    p22 <- ggplot(data.frame(k = 1:K, prob=pk_2), aes(x = k)) +
    geom_point(aes(y=prob), col="blue") +
    geom_linerange(aes(ymax=prob), ymin=0, col="blue") +
    ylab('Probability') + ylim(0, (max(pk_2))) +
      scale_x_continuous(breaks = 1:K, labels = labels_K_5) +
      theme(text = element_text(size=text_size),
            axis.title.x=element_blank(),
            axis.text.y = element_text(size=text_size),
            axis.text.x = element_text(size=r_increase * text_size))
    p22 <- p22 + ggtitle(TeX("$\\hat{X}_{2}$"))
    p22 <- p22 + geom_vline(xintercept = round(mu_X2_hat, 2), col="blue", linetype="dashed", size=dashed_size)
  }

  if (trace) {
    grid.arrange(p10, p20, p11, p21, p12, p22)
  } else {
    p10 <- NA; p20 <- NA
    p11 <- NA; p21 <- NA
    p12 <- NA; p22 <- NA
  }

  return(list("xk"=xk, # return thresholds
              "alpha"=alpha, "mu"=mu,
              "pk_0"=pk_0, "pk_1"=pk_1, "pk_2"=pk_2, # return probabilities
              "xi_1"=xi_1, "omega_1"=omega_1, # return parameters
              "xi_2"=xi_2, "omega_2"=omega_2,
              "p10"=p10, "p20"=p20, # return plots
              "p11"=p11, "p21"=p21,
              "p12"=p12, "p22"=p22,
              "mu_X1_hat"=mu_X1_hat, "sd_X1_hat"=sd_X1_hat, # return means and sd's of hats
              "mu_X2_hat"=mu_X2_hat, "sd_X2_hat"=sd_X2_hat))
}

# function to get probability of rejection for small simulation test
get_p_vals_mat_sequential <- function(alpha=3, # shape parameter of distributions (could be generalized and so each group has it's own shape)
                                      eta=0.05, # noise parameter for thresholds
                                      n=1000, # size of samples (could be generalized so each group has it's own size)
                                      nsim=100, # number of simulations
                                      disc_type="OD", # can be EW or could be generalized to sth else
                                      rec_type="OD" # we can try to reconstruct using also EW, but get worse estimation of omega
) {
  p_vals_mat <- matrix(NA, nrow=4, ncol=2)
  colnames(p_vals_mat) <- c("T-Test", "RD-Test")
  rownames(p_vals_mat) <- c("0.0, 1", "0.0, 0.5", "0.2, 1", "0.2, 0.5")

  export_functions <- c("define_problem", "simulate_Likert", "newtons_method", "estimate_parameters",
                        "apply_Lloyd_Max", "get_MSE_fX_pk", "get_new_xk", "get_new_rk", "get_pk_of_xk",
                        "get_mean_of_pk", "get_var_of_pk", "rLikert", "meanSN", "varSN", "deltaSN",
                        "scale_shift_SN_values", "scale_shift_SN_breaks", "get_p_val_t_test", "RD_test")
  #export_packages <- c("sn")

  for (delta_lambda in rownames(p_vals_mat)) {
    delta_lambda_split <- strsplit(delta_lambda, ", ")[[1]]
    delta <- as.numeric(delta_lambda_split[1])
    lambda <- as.numeric(delta_lambda_split[2])

    # TODO: could combine in 1 for loop
    # TODO: could have different discretization in each group

    # t-test on coded responses
    rejected <- 0
    for(i in 1:nsim) {
      dp <- define_problem(alpha, delta, lambda, eta, disc_type)
      x1 <- rLikert(n, dp$pk_1)
      x2 <- rLikert(n, dp$pk_2)
      p_val <- t.test(x1, x2)$p.value
      rejected <- rejected + as.numeric(p_val < 0.05)
    }
    p_vals_mat[delta_lambda, 1] <- rejected/nsim

    # RD test
    rejected <- 0
    for(i in 1:nsim) {
      dp <- define_problem(alpha, delta, lambda, eta, disc_type)
      testing <- list("yes"=FALSE, "alpha"=alpha, "pk1"=dp$pk_1, "pk2"=dp$pk_2)
      x1 <- rLikert(n, dp$pk_1)
      x2 <- rLikert(n, dp$pk_2)
      p_val <- RD_test(x1, x2, testing=testing)
      rejected <- rejected + as.numeric(p_val < 0.05)
    }
    p_vals_mat[delta_lambda, 2] <- rejected/nsim

  }
  return(p_vals_mat)
}

# function to get probability of rejection for small simulation test in parallel
get_p_vals_mat_parallel <- function(alpha=3, # shape parameter of distributions (could be generalized and so each group has it's own shape)
                                    eta=0.05, # noise parameter for thresholds
                                    n=1000, # size of samples (could be generalized so each group has it's own size)
                                    nsim=100, # number of simulations
                                    disc_type="OD", # can be EW or could be generalized to sth else
                                    rec_type="OD" # we can try to reconstruct using also EW, but get worse estimation of omega
) {
  p_vals_mat <- matrix(NA, nrow=4, ncol=2)
  colnames(p_vals_mat) <- c("T-Test", "RD-Test")
  rownames(p_vals_mat) <- c("0.0, 1", "0.0, 0.5", "0.2, 1", "0.2, 0.5")

  export_functions <- c("define_problem", "simulate_Likert", "newtons_method", "estimate_parameters",
                        "apply_Lloyd_Max", "get_MSE_fX_pk", "get_new_xk", "get_new_rk", "get_pk_of_xk",
                        "get_mean_of_pk", "get_var_of_pk", "rLikert", "meanSN", "varSN", "deltaSN",
                        "scale_shift_SN_values", "scale_shift_SN_breaks", "get_p_val_t_test", "RD_test",
                        "is_rejected", "permutation.test", "bootstrap.diff.of.means", "get_p_val_t_test")

  export_packages <- c("sn")

  for (delta_lambda in rownames(p_vals_mat)) {
    delta_lambda_split <- strsplit(delta_lambda, ", ")[[1]]
    delta <- as.numeric(delta_lambda_split[1])
    lambda <- as.numeric(delta_lambda_split[2])

    for (test_type in colnames(p_vals_mat)) {
      cl <- parallel::makeCluster(4)
      doParallel::registerDoParallel(cl)

      rejected <- foreach(i=1:nsim,
                          .combine = 'c',
                          .packages = export_packages,
                          .export = export_functions) %dopar% {
        dp <- define_problem(alpha, delta, lambda, eta, disc_type)
        testing <- list("yes"=FALSE, "alpha"=alpha, "pk1"=dp$pk_1, "pk2"=dp$pk_2)

        x1 <- rLikert(n, dp$pk_1)
        x2 <- rLikert(n, dp$pk_2)
        is_rejected(x1, x2, test_type, testing)
      }
      p_vals_mat[delta_lambda, test_type] <- sum(rejected)/nsim

      parallel::stopCluster(cl)
    }

  }
  return(p_vals_mat)
}

# function to get probability of rejection for big simulation test
get_p_vals_mat_sequential_big <- function(alpha=3,
                                          delta=0.2,
                                          lambda=0.5,
                                          eta=0.1,
                                          n=2000,
                                          nsim=100) {
  p_vals_mat <- matrix(NA, nrow=3, ncol=5)
  colnames(p_vals_mat) <- c("T-Test", "RD-Test", "MWW", "Bootstr.", "Permut.")
  rownames(p_vals_mat) <- c("Latent", "Likert-EW", "Likert-OD")

  for (var_type in rownames(p_vals_mat)) {
    for (test_type in colnames(p_vals_mat)) {
      if ((var_type == "Latent") & (test_type == "RD-Test")) {
        p_vals_mat[var_type, test_type] <- -1 # RD-Test only works with Likert data
      } else {
        rejected <- 0
        for(i in 1:nsim) {
          if (var_type == "Latent") {
            dp <- define_problem(alpha, delta, lambda) # we only want latent parameters so no need to specify eta or disc_type
            testing <- list("yes"=FALSE, "alpha"=alpha, "pk1"=dp$pk_1, "pk2"=dp$pk_2)
            x1 <- rsn(n, dp$xi_1 - dp$omega_1 * dp$mu, dp$omega_1, dp$alpha)
            x2 <- rsn(n, dp$xi_2 - dp$omega_2 * dp$mu, dp$omega_2, dp$alpha)
          } else if (var_type == "Likert-EW") {
            dp <- define_problem(alpha, delta, lambda, eta, "EW")
            testing <- list("yes"=FALSE, "alpha"=alpha, "pk1"=dp$pk_1, "pk2"=dp$pk_2)
            x1 <- rLikert(n, dp$pk_1)
            x2 <- rLikert(n, dp$pk_2)
          } else if (var_type == "Likert-OD") {
            dp <- define_problem(alpha, delta, lambda, eta, "OD")
            testing <- list("yes"=FALSE, "alpha"=alpha, "pk1"=dp$pk_1, "pk2"=dp$pk_2)
            x1 <- rLikert(n, dp$pk_1)
            x2 <- rLikert(n, dp$pk_2)
          }
          rejected <- rejected + is_rejected(x1, x2, test_type, testing)
        }
        p_vals_mat[var_type, test_type] <- rejected/nsim
      }
    }
  }
  return(p_vals_mat)
}

# function to get probability of rejection for big simulation test in parallel
get_p_vals_mat_parallel_big <- function(alpha=3,
                                          delta=0.2,
                                          lambda=0.5,
                                          eta=0.1,
                                          n=2000,
                                          nsim=100) {
  export_functions <- c("define_problem", "simulate_Likert", "newtons_method", "estimate_parameters",
                        "apply_Lloyd_Max", "get_MSE_fX_pk", "get_new_xk", "get_new_rk", "get_pk_of_xk",
                        "get_mean_of_pk", "get_var_of_pk", "rLikert", "meanSN", "varSN", "deltaSN",
                        "scale_shift_SN_values", "scale_shift_SN_breaks", "get_p_val_t_test", "RD_test",
                        "is_rejected", "permutation.test", "bootstrap.diff.of.means", "get_p_val_t_test")

  export_packages <- c("sn")

  p_vals_mat <- matrix(NA, nrow=3, ncol=5)
  colnames(p_vals_mat) <- c("T-Test", "RD-Test", "MWW", "Bootstr.", "Permut.")
  rownames(p_vals_mat) <- c("Latent", "Likert-EW", "Likert-OD")
  options(warn=-1) # to suppress:
  # Warning message:
  #   In e$fun(obj, substitute(ex), parent.frame(), e$data) :
  #   already exporting variable(s): define_problem, simulate_Likert, newtons_method,...
  # Explanation: we are exporting functions each time in foreach loop and R complains it's already exported, same results as sequential way
  for (var_type in rownames(p_vals_mat)) {
    for (test_type in colnames(p_vals_mat)) {
      if ((var_type == "Latent") & (test_type == "RD-Test")) {
        p_vals_mat[var_type, test_type] <- -1 # RD-Test only works with Likert data
      } else {
        cl <- parallel::makeCluster(4)
        doParallel::registerDoParallel(cl)
        rejected <- foreach(i=1:nsim,
                            .packages = export_packages,
                            .export = export_functions,
                            .combine = 'c') %dopar% {
                              #rejected <- 0
                              #for(i in 1:nsim) {
                              if (var_type == "Latent") {
                                dp <- define_problem(alpha, delta, lambda) # we only want latent parameters so no need to specify eta or disc_type
                                testing <- list("yes"=FALSE, "alpha"=alpha, "pk1"=dp$pk_1, "pk2"=dp$pk_2)
                                x1 <- rsn(n, dp$xi_1 - dp$omega_1 * dp$mu, dp$omega_1, dp$alpha)
                                x2 <- rsn(n, dp$xi_2 - dp$omega_2 * dp$mu, dp$omega_2, dp$alpha)
                              } else if (var_type == "Likert-EW") {
                                dp <- define_problem(alpha, delta, lambda, eta, "EW")
                                testing <- list("yes"=FALSE, "alpha"=alpha, "pk1"=dp$pk_1, "pk2"=dp$pk_2)
                                x1 <- rLikert(n, dp$pk_1)
                                x2 <- rLikert(n, dp$pk_2)
                              } else if (var_type == "Likert-OD") {
                                dp <- define_problem(alpha, delta, lambda, eta, "OD")
                                testing <- list("yes"=FALSE, "alpha"=alpha, "pk1"=dp$pk_1, "pk2"=dp$pk_2)
                                x1 <- rLikert(n, dp$pk_1)
                                x2 <- rLikert(n, dp$pk_2)
                              }
                              is_rejected(x1, x2, test_type, testing)
                            }
        p_vals_mat[var_type, test_type] <- sum(rejected)/nsim
        parallel::stopCluster(cl)
      }
    }
  }
  options(warn=0)
  return(p_vals_mat)
}

# simulate_Likert(disc_type, K, 0, 1, alpha)
simulate_Likert <- function(disc_type, K, xi, omega, alpha,
                            plot_result=FALSE) {
  if (!(disc_type %in% c("OD", "EW"))) {
    print("Wrong disc_type")
    return(1)
  }
  # X ~ SN(\alpha)
  # X_0 ~ X - mu
  # f_X_0(x) = f_X(x + mu)
  #fX <- function(x, alpha) {2 * dnorm(x + meanSN(alpha)) * pnorm(alpha * (x + meanSN(alpha)))}
  fX <- function(x) { dsn(x + meanSN(alpha), 0, 1, alpha) } # same with package sn, probably not even faster

  if (disc_type == "OD") {
    # we could use precalculated xk and pk for K \in 3:10 and alpha \in -10:10 so we don't need to run Lloyd-Max algorithm each time
    sol <- apply_Lloyd_Max(fX, K)
    xk_est <- sol$xk_est
    pk_est <- sol$pk_est
  } else if (disc_type == "EW") {
    a <- -3*sqrt(varSN(alpha)) # no meanSN(alpha) because we are working with X_0
    b <-  3*sqrt(varSN(alpha))
    xk_est <- seq(a, b, length.out = K + 1) # K levels means K + 1 breaks
    xk_est <- xk_est[2:K] # we must extend lower and upper bound to domain bounds to get PMF
    pk_est <- get_pk_of_xk(xk_est, fX)
  } else if (disc_type == "RD") {
    xk_est <- get_random_breaks_discrete(fX_name, prior, K, noise=noise)
    pk_est <- get_pk_of_xk(xk_est, fX)
  }
  xk <- scale_shift_SN_breaks(xk_est, omega, xi, alpha) # scale and shift breaks for normal skew case
  pk <- get_pk_of_xk(xk, fX)

  if (plot_result) {
    xk <- c(-Inf, xk, Inf)
    x <- seq(-2, 3.5, length.out = 1000)
    p1 <- ggplot(data.frame(x=x, y=dsn(x, 0, 1, alpha)), aes(x=x, y=y)) +
      geom_line(col="blue") +
      xlab("x") + ylab("Density")
    p1 <- p1 + ggtitle(TeX(paste(disc_type, "cut points with $\\xi =$", xi, "and $\\omega =$", omega)))
    for (k in 2:K) {
      p1 <- p1 + geom_vline(xintercept = xk[k], col="red", linetype="dashed")
    }
    p1 <- p1 + geom_vline(xintercept = meanSN(alpha), col="green")

    p2 <- ggplot(data.frame(k = 1:K, prob=pk), aes(x = k)) +
      geom_point(aes(y=prob)) +
      geom_linerange(aes(ymax=prob), ymin=0) +
      ylab('Probability') + ylim(0, (max(pk)))
    p2 <- p2 + ggtitle(TeX(paste(disc_type, "$\\hat{X}$ with $\\xi =$", xi, "and $\\omega =$", omega)))
    grid.arrange(p1, p2)
  } else {
    p1 <- "None"
    p2 <- "None"
  }

  rk <- get_new_rk(xk, fX)
  return(list("pk"=pk, "rk"=rk, "xk"=xk, "p1"=p1, "p2"=p2))
}

newtons_method <- function(alpha, pk, xk,
                           x = matrix(c(0, 1)), # initial guess
                           n = 100, # number of iterations
                           trace=FALSE,
                           for_real=c(0,0)) {
  K <- length(pk)
  mu <- meanSN(alpha) # because:    X_0 ~ X - mu    ----->      f_X_0(x) = f_X(x + mu)

  cdf_X <- function(x) {
    if (x == Inf) {
      return(1)
    } else {
      return(psn(x + mu, 0, 1, alpha))
    }
  }
  pdf_X <- function(x) { dsn(x + mu, 0, 1, alpha) }

  hk <- function(x, k) {
    return(
      cdf_X(x[2]*xk[k + 1] - x[2]*x[1])
      - cdf_X(x[2]*xk[k] - x[2]*x[1]) - pk[k])
  }
  dhk_u <- function(x, k) { # partial derivative for u
    if (k == 1) { # we are on the edge p1
      return(pdf_X(x[2]*xk[k + 1] - x[2]*x[1])*(-x[2]))
    } else if (k == K) { # we are on the edge p_K
      return(-pdf_X(x[2]*xk[k] - x[2]*x[1])*(-x[2]))
    } else {
      return(
        pdf_X(x[2]*xk[k + 1] - x[2]*x[1])*(-x[2])
        - pdf_X(x[2]*xk[k] - x[2]*x[1])*(-x[2])
      )
    }
  }
  dhk_v <- function(x, k) { # partial derivative for v
    if (k == 1) { # we are on the edge p_1
      return(pdf_X(x[2]*xk[k + 1] - x[2]*x[1])*(xk[k + 1] - x[1]))
    } else if (k == K) { # we are on the edge p_K
      return(-pdf_X(x[2]*xk[k] - x[2]*x[1])*(xk[k] - x[1]))
    } else {
      return(
        pdf_X(x[2]*xk[k + 1] - x[2]*x[1])*(xk[k + 1] - x[1])
        - pdf_X(x[2]*xk[k] - x[2]*x[1])*(xk[k] - x[1])
      )
    }
  }
  k_range <- 1:5
  # funtion to find roots
  f <- function(x) {
    matrix(sapply(k_range, function(k) hk(x, k)))
  }

  # Jacobian matrix column wise!
  Df <- function(x) {
    matrix(c(sapply(k_range, function(k) dhk_u(x, k)),
             sapply(k_range, function(k) dhk_v(x, k))),
           ncol=2)
  }

  x_trace <- c(x[1])
  y_trace <- c(x[2])
  # if (trace) {
  #   print("0")
  #   print(x)
  # }
  for (i in 1:n) {
    b <- f(x) # eval f
    A <- Df(x) # eval Jacobian

    # can be unstable, we use SVD
    A.svd <- svd(A)
    A_diag <- diag(1 / A.svd$d)
    d <- A.svd$v %*% A_diag %*% t(A.svd$u) %*% (-b) # solve linear least squares norm(A*d + b) = min

    # iteration step adaptive method
    while ((x + d)[2] < 0) { # as long as step is so big, we jump to negative omega
      d <- d/2 # we adjust it (exponentially smaller)
      if (trace) {
        print("adjusted d")
      }
    }

    x <- x + 0.2*d # 0.2
    if (trace) {
      x_trace <- c(x_trace, x[1])
      y_trace <- c(y_trace, x[2])
      # print(i)
      # print(x)
      # cat("\n")
    }
    if (norm(d, "2") < 1e-15) {
      break
    }
  }
  if (trace) {
    # draw contour
    xlen <- 50
    ylen <- 50
    xgrid <- seq(-3, 3, length.out = xlen) # -5, 5
    ygrid <- seq(0.1, 3, length.out = ylen) # 0.1, 10
    zvals <- matrix(NA, ncol = xlen, nrow = ylen)
    for (i in 1:xlen) {
      for (j in 1:ylen) {
        zvals[i, j] <- norm(f(matrix(c(xgrid[i], ygrid[j]))) , "2")
      }
    }
    contour(x = xgrid, y = ygrid, z = zvals,
            col="gray42",
            xlab = TeX("$u$"), ylab = TeX("v"))
    grid(col = "lightgray", lty = "dotted")
    points(x_trace, y_trace, pch=20, col="blue")
  }
  return(x)
}

estimate_parameters <- function(pk,
                                K=5,
                                disc_type="OD",
                                testing=list("yes"=FALSE, "alpha"=NA),
                                trace=FALSE) {
  if (!testing$yes) {
    rng <- c(-10, - 8, - 6, seq(-4, 4, 0.2), 6, 8, 10) # range for alpha
    errors <- c()
    for (alpha_guess in rng) {
      sl <- simulate_Likert(disc_type, K, 0, 1, alpha_guess)
      xk_est <- c(-Inf, sl$xk, Inf)
      sol_xi_omega <- newtons_method(alpha_guess, pk, xk_est)
      omega_est <- 1/sol_xi_omega[2]
      xi_est <- sol_xi_omega[1]
      fX_0 <- function(x) { dsn(x + meanSN(alpha_guess), 0, 1, alpha_guess) }
      fX_alpha <- function(x) { fX_0((x - xi_est)/omega_est)*(1/omega_est) }
      pk_hat <- get_pk_of_xk(xk_est[2:5], fX_alpha)
      error <- norm(matrix(pk_hat - pk), "2")
      errors <- c(errors, error)
    }
    # on sample solve it, by finding peaks and valleys, because it can flip
    ers_pv <- sign(diff(errors)) + 1
    peaks_and_valleys <- c() # P's and V's
    pv_indexes <- c() # where they occur
    epv_previous <- ers_pv[1]
    for (i in 1:length(ers_pv)) {
      epv <- ers_pv[i]
      if (epv > epv_previous) {
        peaks_and_valleys <- c(peaks_and_valleys, "V")
        pv_indexes <- c(pv_indexes, i)
      } else if (epv < epv_previous) {
        peaks_and_valleys <- c(peaks_and_valleys, "P")
        pv_indexes <- c(pv_indexes, i)
      }
      epv_previous <- epv
    }
    # rng[1:4]: -10, ..., -4
    # rng[(length(rng) - 3):length(rng)]: 4, ..., 10
    new_peaks_and_valleys <- c()
    for (i in 1:length(pv_indexes)) {
      condition <- !(pv_indexes[i] %in% 1:4) & !(pv_indexes[i] %in% (length(rng) - 3):length(rng))
      if (condition) {
        new_peaks_and_valleys <- c(new_peaks_and_valleys, peaks_and_valleys[i])
      }
    }
    peaks_and_valleys <- new_peaks_and_valleys

    if (trace) {
      print("Detected peaks and valleys:")
      print(peaks_and_valleys)
    }
    if (length(peaks_and_valleys) == 3) {
      if (sum(peaks_and_valleys == c("V", "P", "V")) == 3) {
        # must take the first valley and leave the tail, that can go down
        min1 <- sort(errors[1:(length(rng) - 4)])[1]
        min2 <- sort(errors[1:(length(rng) - 4)])[2]
        min1 <- which(errors == min1)
        min2 <- which(errors == min2)
        if (rng[min1] < rng[min2]) {
          alpha_est <- rng[min1]
        } else {
          alpha_est <- rng[min2]
        }
      } else if (sum(peaks_and_valleys == c("P", "V", "P")) == 3) {
        # alpha value is very high
        alpha_est <- 10
      }
    } else if (length(peaks_and_valleys) == 2) {
      if (sum(peaks_and_valleys == c("V", "P")) == 2) {
        # take valley and leave the tail, that can go down
        alpha_est <- rng[which.min(errors[1:(length(rng) - 4)])]
      }
    } else if (length(peaks_and_valleys) == 1) {
      if (peaks_and_valleys == "V") {
        alpha_est <- rng[which.min(errors)]
      } else {
        alpha_est <- rng[which.max(errors)]
        alpha_est
      }
    } else {
      alpha_est <- rng[which.min(errors)]
    }
  } else {
    alpha_est <- testing$alpha
  }

  # calculate solution again for omega, xi (could be optimized)
  sl <- simulate_Likert(disc_type, K, 0, 1, alpha_est)
  xk_est <- c(-Inf, sl$xk, Inf)
  sol_xi_omega <- newtons_method(alpha_est, pk, xk_est,
                                     x=matrix(c(0, 1)),
                                     trace=FALSE)
  omega_est <- 1/sol_xi_omega[2]
  xi_est <- sol_xi_omega[1]

  if (trace) {
    par(mfrow=c(1,2))

    # show Gauss-Newton trace
    sol_xi_omega <- newtons_method(alpha_est, pk, xk_est,
                                       x=matrix(c(0, 1)),
                                       trace=TRUE)

    rng <- -10:10 # range for alpha just for testing

    errors <- c()
    for (alpha_guess in rng) {
      sl <- simulate_Likert(disc_type, K, 0, 1, alpha_guess)
      xk_est <- c(-Inf, sl$xk, Inf)
      sol_xi_omega <- newtons_method(alpha_guess, pk, xk_est)
      omega_est_0 <- 1/sol_xi_omega[2]
      xi_est_0 <- sol_xi_omega[1]
      fX_0 <- function(x) { dsn(x + meanSN(alpha_guess), 0, 1, alpha_guess) }
      fX_alpha <- function(x) { fX_0((x - xi_est_0)/omega_est_0)*(1/omega_est_0) }
      pk_hat <- get_pk_of_xk(xk_est[2:5], fX_alpha)
      error <- norm(matrix(pk_hat - pk), "2")
      errors <- c(errors, error)
    }

    # show alpha trace
    #pdf(paste(path, "alpha_trace_2_3_shift_0_2_scale_0_8", ".pdf", sep=""), height = 10, width = 7.5)
    par(mar=c(rep(8, 4)))
    plot(rng, errors, "l",
         xlab=TeX("$\\alpha"),
         ylab=TeX(paste("\\sum_{k} (\\hat{p}_{k} - p_{k})^{2}$", sep=" ")),
         col="blue",
         cex.axis=1.2,
         cex.lab = 1.5)
    axis(side = 1, at = round(alpha_est,2), cex.axis=1.2, col="red")
    grid(col = "lightgray", lty = "dotted")
    abline(h = 0, col="gray")
    par(mfrow=c(1,1))

    #dev.off()
  }

  return(list("alpha_est"=alpha_est, "xi_est"=xi_est, "omega_est"=omega_est))
}

RD_test <- function(x1, x2, # samples of Likert data
                    testing=list("yes"=FALSE, "alpha"=0, "pk1"=0, "pk2"=0), # if we want theoretic result (infinite population all responses)
                    K=5,
                    disc_type="OD",
                    correction=500 # we underestime for omega on the sample, so we do a simple correction
                    ) {
  # if for some level k, p_k is not present in Likert sample x_i (e.g. p_1), set p_k to 0
  fill_levels <- function(x, K) {
    pk <- prop.table(table(x))
    levels <- as.character(1:K)
    for (k in 1:K) {
      if (!(levels[k] %in% names(pk))) {
        pk[[levels[k]]] <- 0
      }
    }
    pk_new <- rep(0, K)
    for (k in 1:K) {
      pk_new[k] <- pk[[levels[k]]]
    }
    return(pk_new)
  }

  if (!testing$yes) {
    n1 <- length(x1)
    n2 <- length(x2)

    pk1 <- fill_levels(x1, K)
    pk2 <- fill_levels(x2, K)
  } else { # theoretic (only for testing)
    n1 <- 10^4
    n2 <- 10^4

    pk1 <- testing$pk1
    pk2 <- testing$pk2
  }

  # estimate first group parameters
  testing$yes <- TRUE # TODO: speed it up: precalculate OD thresholds
  theta1 <- estimate_parameters(pk = pk1, K = K, disc_type = disc_type, testing=testing)
  alpha1_est <- theta1$alpha_est
  mu1_est <- theta1$xi_est
  sd1_est <- theta1$omega_est * sqrt(varSN(alpha1_est)) * (n1 + correction)/n1

  # estimate second group parameters
  theta2 <- estimate_parameters(pk = pk2, K = K, disc_type = disc_type, testing=testing)
  alpha2_est <- theta2$alpha_est
  mu2_est <- theta2$xi_est
  sd2_est <- theta2$omega_est * sqrt(varSN(alpha2_est)) * (n2 + correction)/n2

  p_val <- get_p_val_t_test(mu1_est, sd1_est, mu2_est, sd2_est, n1, n2)
  return(p_val)
}

# test reconstruction idea
get_rejected_matrix <- function(fX_name, disc_type, nsim=10, theoretic=FALSE,
                                prior="unif", alpha=0.05, K=3, n=1000, noise=0.01) {



  p_reject_RD <- matrix(0, nrow=4, ncol=5)
  colnames(p_reject_RD) <- c("T-Test", "Mann-W", "Bootstr.", "RDT-EW", "RDT-OD")
  rownames(p_reject_RD) <- c("0.0, 1", "0.0, 2", "0.2, 1", "0.2, 2")

  for (xi_omega in rownames(p_reject_RD)) {
    # extract xi, omega from row name
    xi_omega_split <- strsplit(xi_omega, ", ")[[1]]
    xi <- as.numeric(xi_omega_split[1])
    omega <- as.numeric(xi_omega_split[2])

    if (theoretic) {
      # for testing we feed actual pk to calculate and estimate functions
      sl1 <- simulate_Likert(fX_name, disc_type, K, 0, omega, FALSE, prior)
      pk1 <- sl1$pk

      sl2 <- simulate_Likert(fX_name, disc_type, K, xi, 1, FALSE, prior)
      pk2 <- sl2$pk
    } else {
      # otherwise we feed zeros to calculate and estimate functions and they use estimates of pk from Likert sample
      pk1 <- rep(0, K)
      pk2 <- rep(0, K)
    }

    for (i in 1:nsim) {
      # create pair of rLikert samples x1 ~ X1 and x2 ~ X2
      # where: X1 ~ Disc(1/omega X)
      # and:   X2 ~ Disc(X - xi)
      # or equivalently: omega * xk_1 and xk_2 + xi
      x1 <- rLikert(n, fX_name, disc_type, K, 0,     omega, prior)
      x2 <- rLikert(n, fX_name, disc_type, K, xi, 1,      prior)

      # calculate p-values using all three approaches: on rLikert responses, EW reconstructed and OD reconstructed
      # t test on simulated Likert responses
      #p_val <- t.test(x1, x2)$p.value
      p_reject_RD[xi_omega, 1] <- p_reject_RD[xi_omega, 1] +
        is_rejected(x = x1, y = x2, test_type = "t", fX_name = fX_name)

      # Mann-Whitney test on simulated Likert responses
      p_reject_RD[xi_omega, 2] <- p_reject_RD[xi_omega, 2] +
        is_rejected(x = x1, y = x2, test_type = "Mann-Whitney-U", fX_name = fX_name)

      # Bootstrap test on simulated Likert responses
      p_reject_RD[xi_omega, 3] <- p_reject_RD[xi_omega, 3] +
        is_rejected(x = x1, y = x2, test_type = "bootstrap", fX_name = fX_name)

      # reverse EW discretization test: (x, y, test_type, fX_name, disc_type="OD", K=5, pk1=0, pk2=0, alpha=0.05)
      p_reject_RD[xi_omega, 4] <- p_reject_RD[xi_omega, 4] +
        is_rejected(x = x1, y = x2, test_type = "RDT", disc_type = "EW", fX_name = fX_name, K = K, pk1 = pk1, pk2 = pk2)

      # reverse OD discretization test
      p_reject_RD[xi_omega, 5] <- p_reject_RD[xi_omega, 5] +
        is_rejected(x = x1, y = x2, test_type = "RDT", disc_type = "OD", fX_name = fX_name, K = K, pk1 = pk1, pk2 = pk2)
    }
  }
  #p_reject_RD/nsim
  return(p_reject_RD/nsim)
}


##########################################################################################
### AUXILIARY FUNCTIONS

# implementation of Lloyd-Max algorithm for simulation of Likert variables from latent variables
# returns xk (breaks), pk (\hat(X) probabilities) and rk (representatives)
# returns also mean and MSE list of approximation for convergence plots
apply_Lloyd_Max <- function(fX, K) {
  rk <- seq(-5, 5, length.out = K) # start with arbitrary representatives in domain of fX

  niter <- 10 # TODO: niter = 10 is probably fine, convergence is fast!
  means_of_pk <- c()
  MSEs_of_fX_pk <- c()
  for (i in 1:niter) {
    xk <- get_new_xk(rk) # calculate new breaks
    rk <- get_new_rk(xk, fX) # calculate new representatives

    # calculate pk's, mean and distortion measure for pk
    pk <- get_pk_of_xk(xk, fX)
    means_of_pk <- c(means_of_pk, get_mean_of_pk(pk))
    MSEs_of_fX_pk <- c(MSEs_of_fX_pk, get_MSE_fX_pk(xk, rk, fX))
  }
  return(list("xk_est"=xk, "pk_est"=pk, "rk_est"=rk, "means_of_pk"=means_of_pk, "MSEs"=MSEs_of_fX_pk))
}

# calculate MSE
get_MSE_fX_pk <- function(xk, rk, fX) {
  K <- length(rk) # generalize it so we can test it
  xk <- c(-Inf, xk, Inf)

  MSE <- 0
  for (k in 1:K) {
    lower_bound <- xk[k]
    upper_bound <- xk[k+1]
    integrand <- function(x) { ((x - rk[k])^2)*fX(x)  }
    MSE <- MSE + integrate(integrand, lower=lower_bound, upper=upper_bound)[[1]]
  }
  return(MSE)
}

# calculate new breaks xk from representatives rk
get_new_xk <- function(rk) {
  K <- (length(rk) - 1) # generalized to arbitrary K
  xk <- rep(0, K)
  for (k in 1:K) {
    xk[k] <- (rk[k] + rk[k+1])/2
  }
  return(xk)
}

# calculate new representatives rk
get_new_rk <- function(xk, fX) {
  K <- length(xk) + 1
  xk <- c(-Inf, xk, Inf)
  rk <- rep(0, K)
  for (k in 1:K) {
    f_upper <- function(x) { x*fX(x) }
    f_lower <- function(x) { fX(x) }
    I_upper <- integrate(f_upper, lower=xk[k], upper=xk[k+1])[[1]]
    I_lower <- integrate(f_lower, lower=xk[k], upper=xk[k+1])[[1]]
    rk[k] <- I_upper/I_lower
  }
  return(rk)
}

# calculate probabilities pk from breaks xk
get_pk_of_xk <- function(xk, fX) {
  K <- length(xk) + 1
  xk <- c(-Inf, xk, Inf)
  pk <- rep(0, K)
  # TODO: maybe faster to just use psn instead of integrate
  for (k in 1:K) {
    lower_bound <- xk[k]
    upper_bound <- xk[k+1]
    pk[k] <- integrate(fX, lower=lower_bound, upper=upper_bound)[[1]]
  }
  return(pk)
}

get_mean_of_pk <- function(pk) {
  domain <- (1:length(pk)) # generalized it so we can test it
  return(sum(pk * domain))
}
get_var_of_pk <- function(pk) {
  domain <- (1:length(pk)) # generalized it so we can test it
  m <- get_mean_of_pk(pk)
  return(sum(pk * (domain - m)^2))
}

plot_Likert <- function(alpha) {
  plots <- list()
  K <- 5
  i <- 1
  for (omega in c(1, 2)) {
    for (xi in c(0, 0.2)) {
      for (disc_type in c("EW", "OD")) {
        sl <- simulate_Likert(disc_type, K, xi, omega, alpha, TRUE)
        plots[[i]] <- list(sl$p1, sl$p2)
        i <- i + 1
      }
    }
  }
  grid.arrange(plots[[1]][[1]], plots[[1]][[2]],
               plots[[2]][[1]], plots[[2]][[2]],
               plots[[3]][[1]], plots[[3]][[2]],
               plots[[4]][[1]], plots[[4]][[2]],
               plots[[5]][[1]], plots[[5]][[2]],
               plots[[6]][[1]], plots[[6]][[2]],
               ncol=4)
}

# for random sample
rLikert <- function(n, pk) {
  K <- length(pk)
  return(sample(x = 1:K, size = n, replace = TRUE, prob = pk))
}

# returns if test of difference of means rejects null hypothesis
is_rejected <- function(x1, x2, test_type, testing) {
  alpha <- 0.05
  if (test_type == "T-Test") {
    return(as.numeric(t.test(x1, x2)$p.value < alpha))
  } else if (test_type == "RD-Test") {
    return(as.numeric(RD_test(x1, x2, testing) < alpha))
  } else if (test_type == "MWW") {
    return(as.numeric(wilcox.test(x1, x2)$p.value < alpha))
  } else if (test_type == "Bootstr.") {
    return(as.numeric(bootstrap.diff.of.means(x1, x2) < alpha))
  } else if (test_type == "Permut.") {
    return(as.numeric(permutation.test(x1, x2) < alpha))
  } else {
    print("Wrong test type.")
    return(NA)
  }
}

# Welch's t-test
get_p_val_t_test <- function(mu1_est, sd1_est, mu2_est, sd2_est, n1, n2) {
  s <- sqrt( (sd1_est^2)/n1 + (sd2_est^2)/n2 )
  df1 <- (((sd1_est^2)/n1)^2)/(n1 - 1)
  df2 <- (((sd2_est^2)/n2)^2)/(n2 - 1)
  df <- ( ((sd1_est^2)/n1 + (sd2_est^2)/n2)^2 ) / (df1 + df2)
  t <- (mu1_est - mu2_est)/s; t
  if (t < 0) {
    return(pt(t, df, lower.tail = TRUE)*2)
  } else {
    return(pt(t, df, lower.tail = FALSE)*2)
  }
}

# based on chapter 16 of Efron's and Tibshirani's An Introduction to the bootstrap (page 220-224)
bootstrap.diff.of.means <- function(x, y, nsim=100) {
  pooled <- c(x,y)
  xt <- x - mean(x) + mean(pooled)
  yt <- y - mean(y) + mean(pooled)

  boot.t <- rep(0, nsim)
  for (i in 1:nsim){
    sample.x <- sample(xt, replace=TRUE)
    sample.y <- sample(yt, replace=TRUE)
    boot.t[i] <- t.test(sample.x, sample.y)$statistic
  }
  t.stat <- t.test(x,y)$statistic
  p.bootstrap <- (1 + sum(abs(boot.t) >= abs(t.stat))) / (nsim + 1)
  return(p.bootstrap)
}

permutation.test <- function(x, y, nsim=100) {
  pool <- c(x, y)
  len_x <- length(x)
  len_y <- length(y)
  obs_diff_p <- mean(x) - mean(y)
  sampl_dist_p <- NULL
  for (i in 1:nsim) {
    resample <- sample(1:length(pool), length(pool))
    x_perm = pool[resample][1:len_x]
    y_perm = pool[resample][(len_x + 1):length(pool)]
    sampl_dist_p[i] = mean(x_perm) - mean(y_perm)
  }
  p.permute <- (1 + sum(abs(sampl_dist_p) >= abs(obs_diff_p))) / (nsim + 1)
  return(p.permute)
}

# mean and variance for skew normal distribution for alpha
deltaSN <- function(alpha=alpha) {
  return(alpha / (sqrt(1 + alpha^2)))
}
meanSN <- function(alpha=alpha) {
  return(deltaSN(alpha) * sqrt(2/pi))
}
varSN <- function(alpha=alpha) {
  return(1 - 2*(deltaSN(alpha)^2)/pi)
}
scale_shift_SN_values <- function(x, omega, xi, alpha) {
  return( (x - meanSN(alpha))/omega + meanSN(alpha) - xi )
}
scale_shift_SN_breaks <- function(x, omega, xi, alpha) {
  return( (x - meanSN(alpha))*omega + meanSN(alpha) + xi )
}

# print out tables
color_alpha <- function(p_vals_mat, delta) {
  pp_vals_mat <- p_vals_mat
  p_vals_mat <- round(p_vals_mat, 2)
  p_vals_mat <- data.frame(p_vals_mat)
  colnames(p_vals_mat) <- c("T-Test", "RD-Test", "MWW", "Bootstr.", "Permut.")
  p_vals_mat$`T-Test` <- as.character(p_vals_mat$`T-Test`)
  p_vals_mat$`RD-Test` <- as.character(p_vals_mat$`RD-Test`)
  p_vals_mat$`MWW` <- as.character(p_vals_mat$`MWW`)
  p_vals_mat$`Bootstr.` <- as.character(p_vals_mat$`Bootstr.`)
  p_vals_mat$`Permut.` <- as.character(p_vals_mat$`Permut.`)
  p_vals_mat
  p_vals_mat[p_vals_mat == "-1"] <- "" # RD-test only works for Likert data

  if (delta == 0) {
    bad <- pp_vals_mat > 0.05
    p_vals_mat[bad] <- sapply(p_vals_mat[bad], function(x) paste("\\textcolor{red}{", x, "}", sep=""))
  } else {
    bad <- (pp_vals_mat < 0.8)
    p_vals_mat[bad] <- sapply(p_vals_mat[bad], function(x) paste("\\textcolor{red}{", x, "}", sep=""))
  }


  return(p_vals_mat)
}

output_TeX_table <- function(delta, n) {
  sink(paste(path, "compare_2_means_results.tex", sep=""))
  cat(
    "\\documentclass[5p,times]{elsarticle}
    \\usepackage[slovene]{babel}
    \\usepackage[utf8]{inputenc}
    \\geometry{textheight=25cm}
    \\usepackage[usenames, dvipsnames]{color}


    \\begin{document}
    \\pagenumbering{gobble}
    ")
  if ((delta == 0) & (n == 1000)) {
    cat(
      "\\setcounter{table}{2}
      "
    )
  } else if ((delta == 0.2) & (n == 1000)) {
    cat(
      "\\setcounter{table}{3}
      "
    )
  } else if ((delta == 0) & (n == 500)) {
    cat(
      "\\setcounter{table}{4}
      "
    )
  } else {
    cat(
      "\\setcounter{table}{5}
      "
    )
  }
  cat(
    "\\begin{center}
    \\begin{table*}
    \\small
    \\begin{tabular}{lcc}
    & $\\lambda = 1$ & $\\lambda = 0.5$ \\\\
    \\vspace{2mm}
    ")

  for (alpha in c(0, 3)) {
    cat("\\vspace{2mm}")
    cat(paste("$\\alpha = ", alpha, "$ &

              "))

    for (lambda in c(1, 0.5)) {
      cat("\\begin{tabular}{lrrrrr}  ")
      file_name <- paste("p_vals_mat", "_delta=", delta, "_alpha=", alpha, "_lambda=", lambda, "_n=", n, ".rds", sep="")
      p_vals_mat <- readRDS(file_name)
      p_vals_mat <- color_alpha(p_vals_mat, delta)
      print(xtable(p_vals_mat, type="latex", hline.after = c(-1, 0),
                   digits = rep(0, 6), align = "lrrrrr"),
            table.placement=NULL, only.contents = TRUE,
            sanitize.text.function = function(x) x)
      cat("\\end{tabular} ")
      if (lambda == 1) {
        cat("&
            ")
      } else {
        cat("
            \\\\
            ")
      }
    }
  }

  if (n == 1000) {
    if (delta == 0) {
      cat(
        "\\\\
\\end{tabular}
\\centering
\\caption{Rejection of null for $\\delta = 0$ and $n=1000$.}

\\end{table*}
\\end{center}"
      )
    } else if (delta == 0.2) {
      cat(
        "\\\\
\\end{tabular}
\\centering
\\caption{Rejection of null for $\\delta = 0.2$ and $n=1000$.}
\\end{table*}
\\end{center}"
      )
    }
  } else if (n == 500) {
    if (delta == 0) {
    cat(
      "\\\\
\\end{tabular}
\\centering
\\caption{Rejection of null for $\\delta = 0$ and $n=500$.}

\\end{table*}
\\end{center}"
    )
  } else if (delta == 0.2) {
    cat(
      "\\\\
\\end{tabular}
\\centering
\\caption{Rejection of null for $\\delta = 0.2$ and $n=500$.}
\\end{table*}
\\end{center}"
    )
   }
 }

  cat("

      \\end{document}
      ")
  sink()
}
