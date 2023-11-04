# Test package on intro example
log_lik <- function(psi0, psi1, Ymat, score = TRUE, finf = TRUE, expected = TRUE)
{
  n1 <- nrow(Ymat)
  n2 <- ncol(Ymat)
  
  YtY <- rowSums(Ymat^2)
  Yt1 <- rowSums(Ymat)
  sum_YtY <- sum(YtY)
  sum_Yt1sq <- sum(Yt1^2)
  
  ll <- -0.5 * n1 * log(psi0 + n2 * psi1) - 0.5 * n1 * (n2 - 1) * log(psi0) -
    0.5 * sum_YtY / psi0 + 0.5 * psi1 * sum_Yt1sq  / (psi0^2 + n2 * psi1 * psi0)
  s <- rep(NA, 2)
  if(score){
    s[1] <- -n1 * 0.5 / (psi0 + n2 * psi1) - 0.5 * n1 * (n2 - 1) / psi0 + 0.5 * sum_YtY / psi0^2 - 
      0.5 * psi1 *  sum_Yt1sq  * (2 * psi0 + n2 * psi1) / (psi0^2 + n2 * psi1 * psi0)^2
    s[2] <- -0.5 * n1 * n2 / (psi0 + n2 * psi1) + 0.5 * sum_Yt1sq  / (psi0 + n2 * psi1)^2
  }
  H <- matrix(NA, 2, 2)
  if(finf){
    if(expected){
      sum_Yt1sq <- n1 * (n2 * psi0 + n2^2 * psi1)
      sum_YtY <- n1 * n2 * (psi1 + psi0)
    }
    t1 <- -n1 / (psi0 + psi1 * n2)^2
    t2 <- - n1 * (n2 - 1) / psi0^2
    t3 <- 2 * sum_YtY / psi0^3
    t4 <- -2 * psi1 * sum_Yt1sq * (2 * psi0 + n2 * psi1)^2 / (psi0^2 + n2 * psi1 * psi0)^3
    t5 <- 2 * psi1 * sum_Yt1sq / (psi0^2 + n2 * psi1 * psi0)^2
    H[1, 1] <- -0.5 * (t1 + t2 + t3 + t4 + t5)
    
    t1 <- -n1 * n2 / (psi0 + n2 * psi1)^2
    t2 <- 2 * sum_Yt1sq / (psi0 + n2 * psi1)^3
    H[1, 2] <- -0.5 * (t1 + t2)
    H[2, 1] <- H[1, 2]
    
    t1 <- -n1 * n2^2 / (psi0 + n2 * psi1)^2
    t2 <- 2 * n2 * sum_Yt1sq /  (psi0 + n2 * psi1)^3
    H[2, 2] <- -0.5 * (t1 + t2)
  }
  return(list("ll" = ll, "score" = s, "finf" = -H))
}
gen_dat <- function(n1, n2, psi0, psi1){
  Ymat <- matrix(rnorm(n1 * n2, sd = sqrt(psi0)), n1, n2)
  Ymat + rnorm(n1, sd = sqrt(psi1)) # Recycling
}

set.seed(12658)
n1 <- 100
n2 <- 3
psi <- c(0.3, 0.02)
Y_test<- gen_dat(n1 = n1, n2 = n2, psi0 = psi[1], psi1 = psi[2])
D <- tibble(response = as.vector(t(Y_test)), clust = rep(1:n1, each = n2))
fit <- lme4::lmer(response ~ 0 + (1|clust), data = D)

Z <- getME(fit, "Z")
y <- getME(fit, "y")
H <- Matrix::Diagonal(n1)
Psi0 <- Matrix::Diagonal(n1, psi[2] / psi[1])
print(limestest::loglik_psi(Z = Z, ZtZXe = crossprod(Z, cbind(Z, y)), e = y, H = H, Psi0 = Psi0, psi0 = psi[1], expected = FALSE))
print(log_lik(psi[1], psi[2], Y_test, expected = FALSE))


