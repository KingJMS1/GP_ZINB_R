#' ZINB_GP
#' @description Run the ZINB NNGP model described in https://doi.org/10.1016/j.jspi.2023.106098.
#'
#' @param X Other Predictor variables
#' @param y Zero inflated count response
#' @param coords Spatial coordinates for NNGP
#' @param nsim How long to run MCMC in total, must be greater than burn, defaults to 5000 iterations
#' @param burn How long to run MCMC before saving samples, defaults to 1000 iterations
#' @param use_count_gp Whether or not to include GP random effects in the count model, defaults to true, recommended
#' @param use_inflation_gp Whether or not to include GP random effects in the zero-inflation portion of the model, defaults to false, not recommended
#' @param thin How often to save MCMC samples, default is 1, saves every iteration. Increase this if running out of memory.
#' @param kern Kernel function, takes a distance matrix and length scale, returns evaluated kernel. e.g. function(dist, ls) {return(exp(-dist / (ls^2)))}
#' @param save_ypred Whether or not to output the predicted values of the response, defaults to false
#' @param print_iter If print_progress is set, how often to print the iteration number of the MCMC chain, defaults to every 100 iterations.
#' @param print_progress Whether or not to print the iteration number of the MCMC chain, defaults to True
#' @param Vs   Spatially varying predictor variables (e.g. one-hot indication of which location this is for varying intercept), wrapped in sparseMatrix from Matrix R package. Will be multiplied by the spatial random effects for prediction.
#' @param Vt   Temporal varying predictor variables, wrapped in sparseMatrix from Matrix R package. Will be multiplied by the temporal random effects for prediction.
#' @param Ds   Spatial distance matrix, diagonal should be 0, off diagonal is distance between elements i and j in space, inputs to the spatial NNGP kernel
#' @param Dt   Temporal distance matirx, diagonal should be 0, off diagonal is distance between elements i and j in time, inputs to the temporal GP kernel
#' @param ltPrior Parameters for a gamma prior and MH update controls for temporal lengthscale: e.g. list(max=50, mh_sd=3, a=1, b=0.001), must contain all listed values.
#' @param lsPrior Parameters for a gamma prior and MH update controls for spatial lengthscale: e.g. list(max=50, mh_sd=3, a=1, b=0.001), must contain all listed values.
#' @param sigmaPrior Parameters for inverse-gamma prior for sigma e.g. list(a=0.01, b=0.1)
#' @param noisePrior Parameters for beta prior for kernel signal to noise ratio, along with MH proposal controls, e.g. list(a=1.5, b=1.5, mh_sd=0.2)
#' @param mh_sd_r MH standard deviation for proposal distribution for r, change if r seems to be walking too slowly. Default is 0.4.
#' @return A List of the following sampled values:          
#' \itemize{
#'      \item {\strong{Alpha:} } {Model coefficients for logit model}
#'      \item {\strong{Beta:} } {Model coefficients for NB model}
#'      \item {\strong{A:} } {Portion of spatial random effect in the logit model explained by kernel}
#'      \item {\strong{B:} } {Portion of temporal random effect in the logit model explained by kernel}
#'      \item {\strong{C:} } {Portion of spatial random effect in the NB model explained by kernel}
#'      \item {\strong{D:} } {Portion of temporal random effect in the NB model explained by kernel}
#'      \item {\strong{L1t:} } {Length scale for temporal kernel in logit model, i.e. } \eqn{e^{-\frac{d^{2}}{2 l_{1t}^{2}}}} 
#'      \item {\strong{Sigma1t:} } {Kernel scale parameter for above kernel, i.e. } \eqn{\sigma_{1t}^{2}e^{.}}
#'      \item {\strong{L2t:} } {Length scale for temporal kernel in NB model, i.e. } \eqn{e^{-\frac{d^{2}}{2 l_{1t}^{2}}}}
#'      \item {\strong{Sigma2t:} } {Kernel scale parameter for above kernel, i.e. } \eqn{\sigma_{2t}^{2}e^{.}}
#'      \item {\strong{Phi_bin:} } {Length scale for spatial kernel in logit model, i.e. } \eqn{e^{-\Phi_{bin}d^{2}}}
#'      \item {\strong{Sigma1s:} } {Square root of multiplier for spatial kernel in logit model}
#'      \item {\strong{Phi_nb:} } {Length scale for spatial kernel in NB model, i.e. } \eqn{e^{-\Phi_{nb}d^{2}}}
#'      \item {\strong{Sigma2s:} } {Square root of multiplier for spatial kernel in NB model}
#'      \item {\strong{R:} } {Dispersion parameter for Negative Binomial distribution.}
#'      \item {\strong{at_risk:} } {At risk indicator for each observation}
#'      \item {\strong{Y_pred:} } {Predictions, sampled from the posterior distribution at each iteration, NULL if save_ypred is false}
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
                    toReturn <- list(Alpha = results$Alpha, Beta = results$Beta, D = results$C, L2t = results$L2s, Sigma2t = results$Sigma2s, Noise2t = results$Noise2s, R = R)
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
                    return(ZINB_GP(X, y, coords, Vs, Vt, Ds, Dt, nsim, burn, thin, save_ypred, print_iter, print_progress, ltPrior, lsPrior, sigmaPrior, noisePrior, mh_sd_r, kern))
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