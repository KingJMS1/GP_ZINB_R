#' estimate
#' @description Draw y from posterior distribution
estimate <- function(X, alpha, beta, Vs, Vt, a, b, c, d, sigma_eps1s, sigma_eps1t, sigma_eps2s, sigma_eps2t, r) {
    N <- nrow(X)
    
    # Draw new epsilon for new predictions of new observations
    eps1s <- rnorm(ncol(Vs), mean = 0, sd = sigma_eps1s)
    eps1t <- rnorm(ncol(Vt), mean = 0, sd = sigma_eps1t)
    eps2s <- rnorm(ncol(Vs), mean = 0, sd = sigma_eps2s)
    eps2t <- rnorm(ncol(Vt), mean = 0, sd = sigma_eps2t)
    
    eta1 <- X %*% alpha + Vs %*% (a + eps1s) + Vt %*% (b + eps1t)
    p_at_risk <- sigmoid(eta1) # at-risk probability
    u <- rbinom(N, 1, p_at_risk) # at-risk indicator
    if (ncol(X) == 1) {
        eta2 <- X[u == 1, ] * beta + Vs[u == 1, ] %*% c + Vt[u == 1, ] %*% d + Vs[u == 1, ] %*% eps2s + Vt[u == 1, ] %*% eps2t # Linear predictor for count part
    } else {
        eta2 <- X[u == 1, ] %*% beta + Vs[u == 1, ] %*% (c + eps2s) + Vt[u == 1, ] %*% (d + eps2t) # Linear predictor for count part
    }
    psi <- sigmoid(eta2) # Prob of success
    mu <- r * psi / (1 - psi) # NB mean
    y <- rep(0, N) # Response
    y[u == 1] <- rnbinom(n=length(psi), size=r, prob=(1 - psi)) # Draw from posterior
    return(y)
}

#' sigmoid
#' @description Compute sigmoid function, clip properly to prevent infinity/nan
sigmoid <- function(eta) {
    eta <- pmin(700, eta)
    return(pmax(1e-6, pmin(1 - 1e-6, exp(eta) / (1 + exp(eta))))) 
}

#' predict
#' @description Predict at new locations/times, provides predictions along with estimated variance
predict <- function(X, Ds_new, Dt_new, Vs_new, Vt_new, output) {
    l1t <- output$L1t
    l2t <- output$L2t
    phi_bin <- output$Phi_bin
    phi_nb <- output$Phi_nb
}
