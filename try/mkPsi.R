# Use components of findbars and mkReTrms in lme4 to create covariance matrix
library(lme4)

################################################################################
## These are the lme4::findbars examples
################################################################################
findbars(f1 <- Reaction ~ Days + (Days | Subject))
## These two are equivalent:% tests in ../inst/tests/test-doubleVertNotation.R
findbars(y ~ Days + (1 | Subject) + (0 + Days | Subject))
findbars(y ~ Days + (Days || Subject))
findbars(~ 1 + (1 | batch / cask))
################################################################################



# This is a modified example from mkReTrms
data("Pixel", package="nlme")
Pixel$New <- rnorm(nrow(Pixel))
# mform <- pixel ~ day + I(day^2) + (day * New | Dog) + (1 | Side/Dog)
mform <- pixel ~ day + I(day^2) + (day * New | Dog) + (1 | Side/Dog)
bars <- findbars(mform)
fr <- model.frame(subbars(mform), data = Pixel)
fit <- lme4::lmer(mform, data = Pixel)

# Make Psi
make_Psi <- function(lmerfit){
  # Lambdat has the right structure, but not the same entries as Psi
  Psi_half <- lme4::getME(lmerfit, "Lambdat")

  # Extract variances and covariances of random effects ordered as in the covmat
  # lower.tri may seem wrong since we are creaing upper tringular Psi,
  # but appears to conform with the indexing in getME(, "Lind")
  psi <- as.data.frame(VarCorr(lmerfit), order = "lower.tri")$vcov

  # Separate error variance and covariance matrix for random effects
  psi0 <- psi[length(psi)]
  psi <- psi[-length(psi)]

  # Fill in the lower-triangular part of Psi with the extracted elements
  param_idx <- lme4::getME(lmerfit, "Lind")
  Psi_half@x <- psi[param_idx]
  Matrix::forceSymmetric(Psi_half, uplo = "U")
}

# Make list of H matrices
make_H_list <- function(lmerfit)
{
  # Psi, and hence H, has the same structure as Lambda
  H <- lme4::getME(lmerfit, "Lambdat")
  param_idx <- lme4::getME(lmerfit, "Lind")

  # Replace values by parameter index
  H@x <- param_idx
  # Return list of H matrices
  lapply(seq_len(getME(lmerfit, "m")),
         function(i){
           M <- H
           M@x <- as.numeric(i == M@x)
           Matrix::forceSymmetric(Matrix::drop0(M), uplo = "U")
         })
}

# THE BELOW CAN BE USED TO CREATE AND INDEXING FOR PSI IF LAMBDAT IS UNAVAILABLE
# Based on the code in mkReTrms
#
# get_Psi_idx <- function(bars, fr, drop.unused.levels = TRUE, reorder.terms = TRUE,
#                  reorder.vars = FALSE)
# {
# if (!length(bars))
#   stop("No random effects terms specified in formula",
#          call. = FALSE)
#   stopifnot(is.list(bars), vapply(bars, is.language, NA), inherits(fr, "data.frame"))
#   names(bars) <- lme4:::barnames(bars)
#   term.names <- vapply(bars, deparse1, "")
#
#   blist <- lapply(bars, lme4:::mkBlist, fr,
#                   drop.unused.levels = drop.unused.levels,
#                   reorder.vars = reorder.vars)
#   # Number of levels per factor per random effect term
#   nl <- vapply(blist, `[[`, 0L, "nl")
#
#   # Order random effect terms by increasing number of levels
#   if (reorder.terms) {
#     if (any(diff(nl) > 0)) {
#       ord <- rev(order(nl))
#       blist <- blist[ord]
#       nl <- nl[ord]
#       term.names <- term.names[ord]
#     }
#   }
#
#   # Number of columns in Z
#   q <- sum(sapply(blist, function(x)nrow(x$sm)))
#   # Ztlist <- lapply(blist, `[[`, "sm")
#   # Zt <- do.call(rbind, Ztlist)
#   # names(Ztlist) <- term.names
#   # q <- nrow(Zt)
#
#
#   cnms <- lapply(blist, `[[`, "cnms")
#   browser()
#   # How many random effects are there per random effect part
#   nc <- lengths(cnms)
#
#   # How many parameters are there per random effect part (number of parameters
#   # in a covariance matrix)
#   nth <- as.integer((nc * (nc + 1))/2)
#
#   # How many columns per random effect part
#   nb <- nc * nl
#   if (sum(nb) != q) {
#     stop(sprintf("total number of RE (%d) not equal to nrow(Zt) (%d)",
#                  sum(nb), q))
#   }
#   # Offsets used for indexing correct submatrices of Psi
#   boff <- cumsum(c(0L, nb))
#   thoff <- cumsum(c(0L, nth))
#
#   make_sparse_core <- function(i) {
#     mm <- matrix(seq_len(nb[i]), ncol = nc[i], byrow = TRUE)
#     dd <- diag(nc[i])
#     ltri <- lower.tri(dd, diag = TRUE)
#     ii <- row(dd)[ltri]
#     jj <- col(dd)[ltri]
#     data.frame(i = as.vector(mm[, ii]) + boff[i],
#                j = as.vector(mm[, jj]) + boff[i],
#                x = as.double(rep.int(seq_along(ii), rep.int(nl[i], length(ii)))
#                              + thoff[i]))
#   }
#
#   # This gives an indexing matrix for lower triangle of RE covmat Psi
#   do.call(sparseMatrix, do.call(rbind, lapply(seq_along(blist),
#                                                  make_sparse_core)))
# }
