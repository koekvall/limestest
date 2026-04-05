#' Confidence interval for a single covariance parameter
#'
#' Computes a score-based confidence interval for a single covariance parameter
#' in a linear mixed model fitted with lme4, by inverting the score test
#' statistic. The signed score profile is evaluated on a grid and the CI bounds
#' are located by linear interpolation at the +/- critical value crossings.
#'1
#' @param lmerfit An \code{lmerMod} object from fitting a linear mixed model
#'   using \code{lme4::lmer}.
#' @param test_idx Single positive integer specifying which covariance parameter
#'   to compute a CI for. Indexes into the vector returned by
#'   \code{as.data.frame(VarCorr(lmerfit), order = "lower.tri")$vcov}, with
#'   error variance last.
#' @param level Numeric confidence level in (0, 1). Default is \code{0.95}.
#' @param step_size Positive numeric. Step size for the outward search (on the
#'   parameter scale). If \code{NULL} (default), set automatically to
#'   \code{SE / 20} where SE is derived from the expected information, giving
#'   roughly 20 steps per Wald CI half-width.
#' @param num_points Positive integer. Maximum number of steps in each
#'   direction. Default is 500. Increase if the CI bound is not found.
#' @param REML Logical or \code{NULL}. If \code{NULL} (default), the estimation
#'   method is taken from \code{lmerfit}.
#' @param expected Logical. If \code{TRUE} (default), use expected Fisher
#'   information.
#' @param known_idx Integer vector or \code{NULL}. Covariance parameters to
#'   treat as fixed (known) when profiling over nuisance parameters. See
#'   \code{\link{score_profile}} for details.
#' @param return_profile Logical. Currently unused; reserved for future use.
#' @param ... Additional arguments passed to the trust-region optimizer.
#'
#' @return If \code{return_profile = FALSE}, a named numeric vector with
#'   elements \code{lower} and \code{upper}. If \code{return_profile = TRUE},
#'   a list with elements \code{ci} (the named vector), \code{null_values}
#'   (grid of null values evaluated), and \code{stat} (signed score profile).
#'
#' @details
#' The confidence interval is constructed as the inversion of the one-dimensional
#' signed score test statistic. Specifically, the CI at level \eqn{1 - \alpha}
#' is \deqn{\{\psi^{(1)} : |T_n(\psi^{(1)})| \leq z_{1-\alpha/2}\},}
#' where \eqn{T_n} is the signed score statistic and \eqn{z_{1-\alpha/2}} is the
#' corresponding standard normal quantile.
#'
#' The signed score profile is evaluated on a uniform grid around the MLE, with
#' nuisance parameters optimized at each grid point using warm starts (see
#' \code{\link{score_profile}}). CI bounds are then found by linear interpolation
#' between adjacent grid points that bracket the \eqn{\pm z_{1-\alpha/2}}
#' crossings. Increase \code{num_points} or \code{search_radius} if a warning is
#' issued about the profile not crossing the critical value.
#'
#' @seealso \code{\link{ci_all_lmer}}, \code{\link{score_profile}}
#'
#' @examples
#' \donttest{
#' library(lme4)
#' fit <- lmer(Reaction ~ Days + (Days | Subject), data = sleepstudy)
#'
#' # 95% CI for the random intercept variance (parameter 1)
#' ci_lmer(fit, test_idx = 1)
#'
#' # 95% CI for the random slope variance (parameter 3)
#' ci_lmer(fit, test_idx = 3)
#' }
#' @export
ci_lmer <- function(lmerfit, test_idx, level = 0.95, step_size = NULL,
                    num_points = 500L, REML = NULL, expected = TRUE,
                    known_idx = NULL, return_profile = FALSE, ...) {

  if (!inherits(lmerfit, "lmerMod"))
    stop("lmerfit must be an lmerMod object from lme4::lmer")

  assertthat::assert_that(
    is.numeric(test_idx), length(test_idx) == 1L,
    test_idx == floor(test_idx), test_idx >= 1L,
    msg = "test_idx must be a single positive integer"
  )
  assertthat::assert_that(
    is.numeric(level), length(level) == 1L, level > 0, level < 1,
    msg = "level must be a single number in (0, 1)"
  )
  assertthat::assert_that(
    is.numeric(num_points), length(num_points) == 1L, num_points >= 2L,
    msg = "num_points must be a single integer >= 2"
  )

  # Extract model components
  Y       <- lme4::getME(lmerfit, "y")
  X       <- lme4::getME(lmerfit, "X")
  Z       <- lme4::getME(lmerfit, "Z")
  Hlist   <- get_Hlist_lmer(lmerfit)
  if (is.null(REML)) REML <- lme4::getME(lmerfit, "REML") != 0
  # For ML, only pass geometry (XtX, XtZ, ZtZ) — residuals are recomputed from
  # beta inside loglikelihood, so they don't go stale during the outward search.
  if (REML) {
    precomp <- get_precomp_lmer(lmerfit, REML = TRUE)
  } else {
    precomp <- list(XtX = as.matrix(crossprod(X)),
                    XtZ = as.matrix(crossprod(X, Z)),
                    ZtZ = methods::as(crossprod(Z), "generalMatrix"))
  }
  psi_hat <- get_psi_hat_lmer(lmerfit)
  r       <- length(psi_hat)
  p       <- ncol(X)

  assertthat::assert_that(
    test_idx <= r,
    msg = paste0("test_idx must not exceed the number of covariance parameters (", r, ")")
  )

  # For ML fits with fixed effects, the full parameter vector is c(beta, psi)
  # and test_idx must be shifted by p to index into psi.
  if (!REML && p > 0) {
    b_hat      <- lme4::fixef(lmerfit)
    theta_hat  <- c(b_hat, psi_hat)
    test_idx_  <- p + test_idx
    known_idx_ <- if (is.null(known_idx)) NULL else p + known_idx
  } else {
    theta_hat  <- psi_hat
    test_idx_  <- test_idx
    known_idx_ <- known_idx
    b_hat      <- NULL
  }

  # Determine step size from expected information if not provided.
  # Use SE/20 so roughly 20 steps cover one Wald CI half-width on each side.
  if (is.null(step_size)) {
    ll <- loglikelihood(psi = psi_hat, b = b_hat, Y = Y, X = X, Z = Z,
                        Hlist = Hlist, REML = REML,
                        get_val = FALSE, get_score = FALSE, get_inf = TRUE,
                        get_beta = (!REML && p > 0),
                        expected = TRUE, precomp = precomp)
    se_approx <- tryCatch(
      sqrt(solve(ll$inf_mat)[test_idx_, test_idx_]),
      error = function(e) sqrt(1 / ll$inf_mat[test_idx_, test_idx_])
    )
    step_size <- se_approx / 40
  }

  z_crit <- stats::qnorm((1 + level) / 2)

  # Search outward in both directions from the MLE, warm-starting each nuisance
  # optimisation from the previous step's solution. Stop as soon as the signed
  # score profile crosses the critical value on that side.
  lower <- .outward_bound(
    theta_hat, test_idx_, z_crit, direction = -1L,
    step_size = step_size, max_steps = as.integer(num_points),
    Y = Y, X = X, Z = Z, Hlist = Hlist,
    REML = REML, expected = expected, known_idx = known_idx_,
    precomp = precomp, p = p, ...
  )
  upper <- .outward_bound(
    theta_hat, test_idx_, z_crit, direction =  1L,
    step_size = step_size, max_steps = as.integer(num_points),
    Y = Y, X = X, Z = Z, Hlist = Hlist,
    REML = REML, expected = expected, known_idx = known_idx_,
    precomp = precomp, p = p, ...
  )

  if (is.finite(lower) && is.finite(upper) && lower >= upper) {
    warning("Lower bound is not less than upper bound. ",
            "Consider decreasing step_size or increasing num_points.")
  }

  vc <- as.data.frame(lme4::VarCorr(lmerfit), order = "lower.tri")
  param_name <- paste0(
    vc$var1[test_idx],
    ifelse(is.na(vc$var2[test_idx]), "", paste0(":", vc$var2[test_idx])),
    " | ", vc$grp[test_idx]
  )

  ci <- matrix(c(lower, upper), nrow = 1L,
               dimnames = list(param_name, c("lower", "upper")))

  if (return_profile)
    warning("return_profile not supported with outward search; returning CI only.")

  ci
}


#' Confidence intervals for covariance parameters
#'
#' Computes score-based confidence intervals for each covariance parameter in a
#' linear mixed model fitted with lme4, by calling \code{\link{ci_lmer}} for
#' each parameter.
#'
#' @param lmerfit An \code{lmerMod} object from fitting a linear mixed model
#'   using \code{lme4::lmer}.
#' @param test_idx Integer vector specifying which parameters to compute CIs for.
#'   If \code{NULL} (default), CIs are computed for all covariance parameters
#'   except the error variance (i.e., all random effect variances and
#'   covariances).
#' @param level Numeric confidence level in (0, 1). Default is \code{0.95}.
#' @param ... Additional arguments passed to \code{\link{ci_lmer}}.
#'
#' @return A matrix with one row per parameter and columns \code{lower} and
#'   \code{upper}. Row names are of the form \code{"var | grouping_factor"} or
#'   \code{"var1:var2 | grouping_factor"} for covariance parameters.
#'
#' @seealso \code{\link{ci_lmer}}
#'
#' @examples
#' \donttest{
#' library(lme4)
#' fit <- lmer(Reaction ~ Days + (Days | Subject), data = sleepstudy)
#'
#' # 95% CIs for all random-effect variance/covariance parameters
#' ci_all_lmer(fit)
#' }
#' @export
ci_all_lmer <- function(lmerfit, test_idx = NULL, level = 0.95, ...) {
  if (!inherits(lmerfit, "lmerMod"))
    stop("lmerfit must be an lmerMod object from lme4::lmer")

  psi_hat <- get_psi_hat_lmer(lmerfit)
  r <- length(psi_hat)

  if (is.null(test_idx)) test_idx <- seq_len(r - 1L)  # Exclude error variance

  do.call(rbind, lapply(test_idx, function(i) {
    ci_lmer(lmerfit, test_idx = i, level = level, ...)
  }))
}


# Internal helper ---------------------------------------------------------

# Find where the numeric vector y (evaluated at x) crosses `target` from above
# (i.e., where y - target transitions from positive to negative).
# Returns the interpolated x value at the first such crossing if side = "lower",
# or the last such crossing if side = "upper".
# Issues a warning and returns +/-Inf if no crossing is found.
.find_crossing <- function(x, y, target) {
  f <- y - target
  # Indices where f is positive and next point is negative (downward crossing)
  down <- which(f[-length(f)] > 0 & f[-1] < 0)

  if (length(down) == 0L) {
    direction <- if (target > 0) "lower" else "upper"
    warning("Score profile does not cross the critical value on the ", direction,
            " side. Consider increasing search_radius or num_points.")
    return(if (target > 0) -Inf else Inf)
  }

  # Lower bound: use first crossing; upper bound: use last crossing
  idx <- if (target > 0) down[1L] else down[length(down)]

  # Linear interpolation
  x[idx] - f[idx] * (x[idx + 1L] - x[idx]) / (f[idx + 1L] - f[idx])
}


# Search outward from the MLE in one direction until the signed score profile
# crosses the critical value, then interpolate to find the CI bound.
#
# direction: -1L to search left (lower bound), +1L to search right (upper bound)
# target:  +z_crit for lower bound, -z_crit for upper bound
#
# At each step the nuisance optimizer warm-starts from the previous solution,
# which maintains feasibility automatically. If the starting combination
# (new test value, previous nuisance solution) is infeasible, the step is
# halved up to 20 times before giving up.
.outward_bound <- function(psi_hat, test_idx, z_crit, direction,
                           step_size, max_steps,
                           Y, X, Z, Hlist, REML, expected, known_idx,
                           precomp, p = 0L, ...) {

  target    <- if (direction == -1L) z_crit else -z_crit
  r         <- length(psi_hat)
  exclude   <- unique(c(test_idx, known_idx))
  opt_idx   <- seq_len(r)[-exclude]

  .ll_val <- function(theta) {
    # For ML with fixed effects, the first p elements are beta
    b_   <- if (!REML && p > 0L) theta[seq_len(p)] else NULL
    psi_ <- if (!REML && p > 0L) theta[-seq_len(p)] else theta
    tryCatch(
      loglikelihood(psi = psi_, b = b_, Y = Y, X = X, Z = Z,
                    Hlist = Hlist, REML = REML,
                    get_val = TRUE, get_score = FALSE, get_inf = FALSE,
                    expected = TRUE, precomp = precomp)$value,
      error = function(e) -Inf
    )
  }

  theta      <- psi_hat          # current parameter vector (warm-start state)
  prev_val   <- psi_hat[test_idx]
  prev_stat  <- 0                # stat at MLE is 0 by definition
  step       <- step_size

  for (i in seq_len(max_steps)) {

    # Propose next test value
    theta_prop <- theta
    theta_prop[test_idx] <- theta[test_idx] + direction * step

    # If starting combination is infeasible, halve the step up to 20 times
    n_halve <- 0L
    while (!is.finite(.ll_val(theta_prop)) && n_halve < 20L) {
      step <- step / 2
      theta_prop[test_idx] <- theta[test_idx] + direction * step
      n_halve <- n_halve + 1L
    }
    if (n_halve == 20L) {
      # True feasibility boundary reached before crossing z_crit
      side <- if (direction == -1L) "lower" else "upper"
      warning("Hit feasibility boundary before crossing the critical value on ",
              "the ", side, " side. CI bound may be at the boundary of the ",
              "parameter space.")
      return(if (direction == -1L) -Inf else Inf)
    }

    # Optimise nuisance parameters from the current warm-start
    if (length(opt_idx) > 0L) {
      opt <- tryCatch(
        maximize_loglik(start_val = theta_prop, opt_idx = opt_idx,
                        Y = Y, X = X, Z = Z, Hlist = Hlist,
                        expected = expected, REML = REML,
                        precomp = precomp, ...)$arg,
        error = function(e) NULL
      )
      if (is.null(opt)) {
        step <- step / 2
        next
      }
      theta_prop <- opt
    }

    # Compute signed score statistic
    stat <- tryCatch(
      as.numeric(score_stat(theta = theta_prop, test_idx = test_idx,
                            Y = Y, X = X, Z = Z, Hlist = Hlist,
                            REML = REML, expected = expected,
                            efficient = TRUE, signed = TRUE,
                            known_idx = known_idx, precomp = precomp)),
      error = function(e) NA_real_
    )

    if (is.na(stat)) {
      step <- step / 2
      next
    }

    # Reset step to full size after any halvings
    step <- step_size

    # Check for crossing of the target critical value
    # (prev_stat - target) and (stat - target) have opposite signs at a crossing
    if ((prev_stat - target) * (stat - target) < 0) {
      # Linear interpolation between previous and current test values
      x1 <- prev_val;  y1 <- prev_stat - target
      x2 <- theta_prop[test_idx]; y2 <- stat - target
      return(x1 - y1 * (x2 - x1) / (y2 - y1))
    }

    # Advance warm-start state
    theta     <- theta_prop
    prev_val  <- theta_prop[test_idx]
    prev_stat <- stat
  }

  side <- if (direction == -1L) "lower" else "upper"
  warning("CI ", side, " bound not found within ", max_steps, " steps. ",
          "Consider increasing num_points or decreasing step_size.")
  if (direction == -1L) -Inf else Inf
}
