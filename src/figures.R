rm(list = ls()) 
path <- "/nfs/general/repos/latent-variable-reconstruction/src"
setwd(path)
set.seed(1)
source("functions.R")

#######################################################
## DEFINE THE PROBLEM
dp <- define_problem(alpha=0, 0, 0.5, 0, "OD", trace = TRUE) # alpha, delta, lambda, eta, disc_type, trace
print(dp$mu_X1_hat - dp$mu_X2_hat)

h0 <- 4
w0 <- 7
pdf("../data/X_0.pdf", height = h0, width = w0)
dp$p10
dev.off()

pdf("../data/X_0_hat.pdf", height = h0, width = w0)
dp$p20
dev.off()

pdf("../data/X_1.pdf", height = h0, width = w0)
dp$p11
dev.off()

pdf("../data/X_1_hat.pdf", height = h0, width = w0)
dp$p21
dev.off()

pdf("../data/X_2.pdf", height = h0, width = w0)
dp$p12
dev.off()

pdf("../data/X_2_hat.pdf", height = h0, width = w0)
dp$p22
dev.off()


#######################################################
## HOOPS EXAMPLE FOR MOTIVATION
d <- 3 # depth of water tank

# Naive approach
rk <- seq(0, d, length.out = 7) 
rk[2:6]
# 0.5 1.0 1.5 2.0 2.5

# EW approach
xk <- seq(2 - 3/2, 2 + 3/2, length.out = 6)
xk <- c(0, xk[2:5], 3)
d - xk
rk <- get_hoops(xk[2:5])
d - rk

# correct approach
# calculate centroids (hoops)
get_hoops <- function(xk) {
  fX <- function(x) { (2/9)*x }
  K <- length(xk) + 1
  xk <- c(0, xk, 3)
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
F_X <- function(x) { (x^2)/9 }
q_F_X <- function(q) { 3*sqrt(q) }
xk <- q_F_X(seq(0.2, 0.8, 0.2))
rk <- get_hoops(xk)

d - xk
# [1] 1.6583592 1.1026334 0.6762100 0.3167184
d - rk
# 2.1055728 1.3646051 0.8822421 0.4921625 0.1554175
round(abs(diff(d - rk)), 3)



##########################################################################################
### COMPARE 2 MEANS

## SMALL SIMULATION TEST
alpha <- 0 # normal case
n <- 1000 # reflects the sample size of many contemporary commercial, academic and omnibus surveys
nsim <- 1000
eta <- 0.1
dp <- define_problem(alpha, 0.2, 0.5, eta, "OD", trace = TRUE) # alpha, delta, lambda, eta, disc_type, trace

# no correction
get_p_vals_mat_parallel(alpha, eta, n, nsim, "OD") # alpha, eta, n, nsim, disc_type
# T-Test RD-Test
# 0.0, 1    0.047   0.091
# 0.0, 0.5  0.902   0.113
# 0.2, 1    0.990   0.988
# 0.2, 0.5  0.603   1.000

#get_p_vals_mat_sequential(alpha, eta, n, nsim, "OD") # non-parallel version, same results, just for testing
get_p_vals_mat_parallel(alpha, eta, n, nsim, "OD") # alpha, eta, n, nsim, disc_type
# T-Test RD-Test
# 0.0, 1    0.044   0.007
# 0.0, 0.5  0.899   0.024
# 0.2, 1    0.985   0.906
# 0.2, 0.5  0.622   0.995

get_p_vals_mat_parallel(alpha, eta, n, nsim, "EW") # alpha, eta, n, nsim, disc_type
# T-Test RD-Test
# 0.0, 1    0.056   0.016
# 0.0, 0.5  0.276   0.017
# 0.2, 1    0.984   0.910
# 0.2, 0.5  0.926   0.989

## BIG SIMULATION TEST

alpha <- 3
delta <- 0
lambda <- 0.5
eta <- 0.1
n <- 1000
nsim <- 1000

eta <- 0.1
nsim <- 1000
for (n in c(1000, 500)) {
  for (delta in c(0, 0.2)) {
    for (alpha in c(0, 3)) {
      for (lambda in c(1, 0.5)) {
        p_vals_mat <- get_p_vals_mat_parallel_big(alpha, delta, lambda, eta, n, nsim)
        file_name <- paste("p_vals_mat", "_delta=", delta, "_alpha=", alpha, "_lambda=", lambda, "_n=", n, ".rds", sep="")
        saveRDS(p_vals_mat, file_name)
      }
    }
  }
}

for (n in c(1000, 500)) {
  for (delta in c(0, 0.2)) {
    output_TeX_table(delta, n)
    system("pdflatex compare_2_means_results.tex")
    system(paste("pdfcrop compare_2_means_results.pdf", "out.pdf"))
    if (delta == 0) {
      system(paste("pdftk out.pdf cat 2-end output compare_means_delta=", delta, "_n=", n, ".pdf", sep="")) 
    } else {
      system(paste("pdftk out.pdf cat 2-end output compare_means_delta=", "0_2", "_n=", n, ".pdf", sep="")) 
    }
  }
}




