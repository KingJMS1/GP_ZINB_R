#' Draw a Zero-Inflated Negative-Binomial Response
#'
#' Generates one posterior predictive count for each pair of at-risk and count
#' linear predictors using the parameterization fitted by `ZINB_GP()`.
#'
#' @param eta_at_risk Numeric vector of linear predictors for the Bernoulli
#'   at-risk component.
#' @param eta_count Numeric vector of linear predictors for the
#'   negative-binomial count component.
#' @param r Positive negative-binomial dispersion parameter.
#'
#' @return A non-negative integer vector with one predictive count per linear
#'   predictor pair.
#' @keywords internal
draw_zinb_response <- function(eta_at_risk, eta_count, r) {
    eta_at_risk <- c(eta_at_risk)
    eta_count <- c(eta_count)
    N <- length(eta_at_risk)

    if (length(eta_count) != N) {
        stop("`eta_at_risk` and `eta_count` must have the same length.")
    }

    at_risk <- stats::rbinom(N, 1, sigmoid(eta_at_risk))
    y <- integer(N)

    if (any(at_risk == 1)) {
        y[at_risk == 1] <- stats::rnbinom(
            sum(at_risk == 1),
            size = r,
            prob = sigmoid(-eta_count[at_risk == 1])
        )
    }

    y
}


#' @title Numerically stable sigmoid function
#' @description Computes the inverse-logit transformation while clipping its
#' output away from zero and one.
#' @param eta Numeric vector or matrix of linear predictors.
#' @return `eta` transformed to probabilities in `[1e-6, 1 - 1e-6]`.
sigmoid <- function(eta) {
    eta <- pmin(700, eta)
    return(pmax(1e-6, pmin(1 - 1e-6, exp(eta) / (1 + exp(eta))))) 
}

#' @title Predict at new locations and times
#' @description Predicts from fitted Gaussian-process model output at new
#' spatial locations and times.
#' @param X Fixed-effect design matrix for the prediction observations.
#' @param Ds_new Distances between new and observed spatial locations.
#' @param Dt_new Distances between new and observed temporal locations.
#' @param Vs_new Spatial random-effect design matrix for prediction observations.
#' @param Vt_new Temporal random-effect design matrix for prediction observations.
#' @param output Fitted model output containing posterior parameter draws.
#' @return A prediction object. This function is currently a placeholder and
#'   does not return calculated predictions.
predict <- function(X, Ds_new, Dt_new, Vs_new, Vt_new, output) {
    l1t <- output$L1t
    l2t <- output$L2t
    phi_bin <- output$L1s
    phi_nb <- output$L2s
}
