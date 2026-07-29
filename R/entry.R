#' ZINB_GP
#' @description Fits the zero-inflated negative-binomial Gaussian-process model
#' described in <https://doi.org/10.1016/j.jspi.2023.106098>. The supplied
#' spatial and temporal design and distance matrices determine which Gaussian
#' processes are included.
#'
#' @param X Fixed-effect design matrix with one row per observation.
#' @param y Non-negative integer count response.
#' @param coords Spatial coordinates used by the legacy spNNGP implementation.
#' @param nsim Total number of MCMC iterations; must exceed `burn`.
#' @param burn Number of burn-in iterations.
#' @param use_count_gp Whether to include GP random effects in the count component.
#' @param use_inflation_gp Whether to include GP random effects in the zero-inflation component.
#' @param thin Store every `thin`-th iteration after burn-in.
#' @param kern Kernel function accepting a distance matrix and length scale.
#' @param save_ypred Whether to save posterior predictive draws.
#' @param print_iter Print progress every `print_iter` iterations.
#' @param print_progress Whether to print MCMC progress.
#' @param Vs Spatial random-effect design matrix with one row per observation.
#' @param Vt Temporal random-effect design matrix with one row per observation.
#' @param Ds Spatial distance matrix; its diagonal must be zero.
#' @param Dt Temporal distance matrix; its diagonal must be zero.
#' @param ltPrior List with `max`, `mh_sd`, `a`, and `b` for temporal length-scale prior and proposal controls.
#' @param lsPrior List with `max`, `mh_sd`, `a`, and `b` for spatial length-scale prior and proposal controls.
#' @param sigmaPrior List with `a` and `b` inverse-gamma prior parameters for GP scales.
#' @param noisePrior List with `a`, `b`, and `mh_sd` for the GP noise-ratio prior and proposal.
#' @param mh_sd_r Proposal standard deviation for the negative-binomial dispersion parameter.
#' @return A list containing posterior MCMC draws:
#' \describe{
#'   \item{Alpha}{Fixed-effect coefficients for the zero-inflation component.}
#'   \item{Beta}{Fixed-effect coefficients for the count component.}
#'   \item{A, B}{Spatial and temporal random effects for the zero-inflation component.}
#'   \item{C, D}{Spatial and temporal random effects for the count component.}
#'   \item{L1t, L2t}{Temporal GP length scales for the zero-inflation and count components.}
#'   \item{Sigma1t, Sigma2t}{Temporal GP scale parameters.}
#'   \item{Phi_bin, Phi_nb}{Spatial GP length-scale parameters.}
#'   \item{Sigma1s, Sigma2s}{Spatial GP scale parameters.}
#'   \item{R}{Negative-binomial dispersion parameter.}
#'   \item{at_risk}{At-risk indicator draws for each observation.}
#'   \item{Y_pred}{Posterior predictive draws, or `NULL` when `save_ypred` is `FALSE`.}
#' }
#' @export
ZINB_GP <- function(X, y, coords, nsim = 5000, burn = 1000, use_count_gp = TRUE, use_inflation_gp = FALSE, thin = 1, kern = NULL, save_ypred = FALSE, print_iter = 100, print_progress = FALSE, Vs = NULL, Vt = NULL, Ds = NULL, Dt = NULL, ltPrior = NULL, lsPrior = NULL, sigmaPrior = NULL, noisePrior = NULL, mh_sd_r = NULL) 
{
    errMsg <- "Error: must specify at least 1 GP to use. Use optimization GLM software like INLA, MASS, glmmTMB, pscl, etc. to fit the model instead.\n"
    if (Vs == NULL)
    {
        if (Vt == NULL)
        {
            # No GPS
            cat(errMsg)
            return(NULL)
        }
        else
        {
            # Only temporal GPs
            if (use_count_gp)
            {
                if (use_inflation_gp)
                {
                    # Both types of GP, temporal only
                    results <- ZINB_GP_spatial(X, y, Vt, Dt, nsim, burn, thin, save_ypred, print_iter, print_progress, ltPrior, sigmaPrior, noisePrior, mh_sd_r, kern)
                    toReturn <- list(Alpha = results$Alpha, Beta = results$Beta, B = results$A, D = results$C, L1t = results$L1s, Sigma1t = results$Sigma1s, Noise1t = results$Noise1s, L2t = results$L2s, Sigma2t = results$Sigma2s, Noise2t = results$Noise2s, R = results$R)
                    return(toReturn)
                }
                else
                {
                    # Only count GP, temporal only
                    results <- ZINB_GP_spatial_count(X, y, Vt, Dt, nsim, burn, thin, save_ypred, print_iter, print_progress, ltPrior, sigmaPrior, noisePrior, mh_sd_r, kern)
                    toReturn <- list(Alpha = results$Alpha, Beta = results$Beta, D = results$C, L2t = results$L2s, Sigma2t = results$Sigma2s, Noise2t = results$Noise2s, R = results$R)
                    return(toReturn)
                }
            }
            else
            {
                # No count GP
                if (use_inflation_gp)
                {
                    # Inflation GP only
                    results <- ZINB_GP_spatial_inflation(X, y, Vt, Dt, nsim, burn, thin, save_ypred, print_iter, print_progress, ltPrior, sigmaPrior, noisePrior, mh_sd_r, kern)
                    toReturn <- list(Alpha = results$Alpha, Beta = results$Beta, B = results$A, L1t = results$L1s, Sigma1t = results$Sigma1s, Noise1t = results$Noise1s, R = results$R)
                    return()
                }
                else
                {
                    # No GPS
                    cat(errMsg)
                    return(NULL)
                }
            }
        }
    }
    else
    {
        if (Vt == NULL)
        {
            # Only spatial GPs
            if (use_count_gp)
            {
                # Do use count gps
                if (use_inflation_gp)
                {
                    # Both GP types, spatial only
                    return(ZINB_GP_spatial(X, y, Vs, Ds, nsim, burn, thin, save_ypred, print_iter, print_progress, lsPrior, sigmaPrior, noisePrior, mh_sd_r, kern))
                }
                else
                {
                    # Only count GP, spatial only
                    return(ZINB_GP_spatial_count(X, y, Vs, Ds, nsim, burn, thin, save_ypred, print_iter, print_progress, lsPrior, sigmaPrior, noisePrior, mh_sd_r, kern))
                }
            }
            else
            {
                # No count GPS
                if (use_inflation_gp)
                {
                    # Inflation gp only, spatial only
                    return(ZINB_GP_spatial_inflation(X, y, Vs, Ds, nsim, burn, thin, save_ypred, print_iter, print_progress, lsPrior, sigmaPrior, noisePrior, mh_sd_r, kern))
                }
                else
                {
                    # No GPS
                    cat(errMsg)
                    return(NULL)
                }
            }
        }
        else
        {
            # All GPs
            if (use_count_gp)
            {
                # Do use count gps
                if (use_inflation_gp)
                {
                    # Full GPs
                    return(ZINB_GP_orig(X, y, coords, Vs, Vt, Ds, Dt, nsim, burn, thin, save_ypred, print_iter, print_progress, ltPrior, lsPrior, sigmaPrior, noisePrior, mh_sd_r, kern))
                }
                else
                {
                    # Count GP only
                    return(ZINB_GP_count(X, y, coords, Vs, Vt, Ds, Dt, nsim, burn, thin, save_ypred, print_iter, print_progress, ltPrior, lsPrior, sigmaPrior, noisePrior, mh_sd_r, kern))
                }
            }
            else
            {
                # No Count GPS
                if (use_inflation_gp)
                {
                    return(ZINB_GP_inflation(X, y, coords, Vs, Vt, Ds, Dt, nsim, burn, thin, save_ypred, print_iter, print_progress, ltPrior, lsPrior, sigmaPrior, noisePrior, mh_sd_r, kern))
                }
                else
                {
                    # No GPS
                    cat(errMsg)
                    return(NULL)
                }
            }
        }
    }
}
