#' @title Draw from the full posterior predictive distribution
#' @description Draw a zero-inflated negative-binomial response conditional on
#' sampled fixed effects and spatial and temporal random effects.
#' @param X Fixed-effect design matrix.
#' @param alpha Fixed-effect coefficients for the zero-inflation component.
#' @param beta Fixed-effect coefficients for the count component.
#' @param Vs Spatial random-effect design matrix.
#' @param Vt Temporal random-effect design matrix.
#' @param a Spatial random-effect draw for the zero-inflation component.
#' @param b Temporal random-effect draw for the zero-inflation component.
#' @param c Spatial random-effect draw for the count component.
#' @param d Temporal random-effect draw for the count component.
#' @param r Negative-binomial dispersion parameter.
#' @return A vector of posterior predictive count draws.
#' @importFrom stats rbinom rnbinom
estimate <- function(X, alpha, beta, Vs, Vt, a, b, c, d, r) {
    N <- nrow(X)
    
    # Calculate etas
    eta1 <- X %*% alpha + Vs %*% a + Vt %*% b
    p_at_risk <- sigmoid(eta1) # at-risk probability
    u <- rbinom(N, 1, p_at_risk) # at-risk indicator
    if (ncol(X) == 1) {
        eta2 <- X[u == 1, ] * beta + Vs[u == 1, ] %*% c + Vt[u == 1, ] %*% d # Linear predictor for count part
    } else {
        eta2 <- X[u == 1, ] %*% beta + Vs[u == 1, ] %*% c + Vt[u == 1, ] %*% d # Linear predictor for count part
    }
    psi <- sigmoid(eta2) # Prob of success
    mu <- r * psi / (1 - psi) # NB mean
    y <- rep(0, N) # Response
    y[u == 1] <- rnbinom(n=length(psi), size=r, prob=(1 - psi)) # Draw from posterior
    return(y)
}

#' @title Draw posterior predictions without inflation random effects
#' @description Draw a zero-inflated negative-binomial response conditional on
#' fixed effects and count-component spatial and temporal random effects.
#' @param X Fixed-effect design matrix.
#' @param alpha Fixed-effect coefficients for the zero-inflation component.
#' @param beta Fixed-effect coefficients for the count component.
#' @param Vs Spatial random-effect design matrix.
#' @param Vt Temporal random-effect design matrix.
#' @param c Spatial random-effect draw for the count component.
#' @param d Temporal random-effect draw for the count component.
#' @param r Negative-binomial dispersion parameter.
#' @return A vector of posterior predictive count draws.
estimate_noab <- function(X, alpha, beta, Vs, Vt, c, d, r) {
    N <- nrow(X)
    
    # Calculate etas
    eta1 <- X %*% alpha
    p_at_risk <- sigmoid(eta1) # at-risk probability
    u <- rbinom(N, 1, p_at_risk) # at-risk indicator
    if (ncol(X) == 1) {
        eta2 <- X[u == 1, ] * beta + Vs[u == 1, ] %*% c + Vt[u == 1, ] %*% d # Linear predictor for count part
    } else {
        eta2 <- X[u == 1, ] %*% beta + Vs[u == 1, ] %*% c + Vt[u == 1, ] %*% d # Linear predictor for count part
    }
    psi <- sigmoid(eta2) # Prob of success
    mu <- r * psi / (1 - psi) # NB mean
    y <- rep(0, N) # Response
    y[u == 1] <- rnbinom(n=length(psi), size=r, prob=(1 - psi)) # Draw from posterior
    return(y)
}


#' @title Draw posterior predictions without count random effects
#' @description Draw a zero-inflated negative-binomial response conditional on
#' fixed effects and zero-inflation-component spatial and temporal random effects.
#' @param X Fixed-effect design matrix.
#' @param alpha Fixed-effect coefficients for the zero-inflation component.
#' @param beta Fixed-effect coefficients for the count component.
#' @param Vs Spatial random-effect design matrix.
#' @param Vt Temporal random-effect design matrix.
#' @param a Spatial random-effect draw for the zero-inflation component.
#' @param b Temporal random-effect draw for the zero-inflation component.
#' @param r Negative-binomial dispersion parameter.
#' @return A vector of posterior predictive count draws.
estimate_nocd <- function(X, alpha, beta, Vs, Vt, a, b, r) {
    N <- nrow(X)
    
    # Calculate etas
    eta1 <- X %*% alpha + Vs %*% a + Vt %*% b
    p_at_risk <- sigmoid(eta1) # at-risk probability
    u <- rbinom(N, 1, p_at_risk) # at-risk indicator
    if (ncol(X) == 1) {
        eta2 <- X[u == 1, ] * beta # Linear predictor for count part
    } else {
        eta2 <- X[u == 1, ] %*% beta # Linear predictor for count part
    }
    psi <- sigmoid(eta2) # Prob of success
    mu <- r * psi / (1 - psi) # NB mean
    y <- rep(0, N) # Response
    y[u == 1] <- rnbinom(n=length(psi), size=r, prob=(1 - psi)) # Draw from posterior
    return(y)
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
    phi_bin <- output$Phi_bin
    phi_nb <- output$Phi_nb
}
