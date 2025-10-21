#' noise_mix
#' @description Create the following matrix: ratio * A + (1-ratio) * I
#' 
#' @param A Matrix to normalize (square)
#' @param noise_ratio Noise mixing ratio
noise_mix <- function(A, noise_ratio) {
    return(noise_ratio * A + diag(1 - noise_ratio, nrow=nrow(A)))
}

#' kern
#' @description Create the kernel matrix e^(-dist / (ls^2))
#' 
#' @param dist Distance matrix
#' @param ls length scale
kern <- function(dist, ls) {
    return(exp(-dist / (ls^2)))
}


#' update_ls_sigma_noise
#' @description Update kernel parameters for a GP
#' 
#' @param ls Current length scale
#' @param sigma Current sigma
#' @param noise_ratio Current noise ratio
#' @param gpdraw Last draw from the gp with these parameters
#' @param K Current kernel matrix
#' @param D Distance matrix
#' @param lsprior prior information for length scale, needs mh_sd, max, a, b
#' @param sigmaPrior prior information for sigma, needs a, b
#' @param noisePrior prior information for noise_ratio, needs mh_sd, a, b
#' #' @return A List of the following sampled values:          
#' \itemize{
#'      \item {\strong{ls:} } {Length scale}
#'      \item {\strong{sigma:} } {sigma}
#'      \item {\strong{noise_ratio:} } {noise ratio}
#'      \item {\strong{K:} } {Kernel matrix}
#'      \item {\strong{K_inv:} } {Inverse of kernel matrix}
#' }
#' @importFrom stats dgamma
#' @improtFrom stats dbeta
#' @importFrom stats runif
#' @importFrom mvtnorm dmvnorm
update_ls_sigma_noise <- function(ls, sigma, noise_ratio, gpdraw, K, D, lsPrior, sigmaPrior, noisePrior) {
    # update ls
    proposal <- max(min(rnorm(1, ls, lsprior$mh_sd), lsprior$max - 1), 1e-6)
    if (TRUE) {
        K_star <- sigma^2 * noise_mix(kern(D, proposal), noise_ratio)
        
        # Calculate model likelihood
        likelihood_ls <- dmvnorm(gpdraw, mean = rep(0, length(gpdraw)), sigma = K_star, log = TRUE) -
            dmvnorm(gpdraw, mean = rep(0, length(gpdraw)), sigma = K, log = TRUE)
        
        # Calculate prior likelihood
        prior_ls <- dgamma(x = proposal, shape = lsprior$a, rate = lsprior$b, log = TRUE) -
            dgamma(x = ls, shape = lsprior$a, rate = lsprior$b, log = TRUE)
        
        posterior_ls <- likelihood_ls + prior_ls

        if (!is.na(posterior_ls)) {
            if (log(runif(1)) < posterior_ls) {
                ls <- proposal
                K <- K_star
            }
        }
    }
    
    # Update noise ratio
    proposal <- max(min(rnorm(1, noise_ratio, noisePrior$mh_sd), 1 - 1e-7), 0 + 1e-7)
    if (TRUE) {
        K_star <- sigma^2 * noise_mix(kern(D, ls), proposal)
        
        # Calculate model likelihood
        likelihood_nr <- dmvnorm(gpdraw, mean = rep(0, length(gpdraw)), sigma = K_star, log = TRUE) -
            dmvnorm(gpdraw, mean = rep(0, length(gpdraw)), sigma = K, log = TRUE)
        
        # Calculate prior likelihood
        prior_nr <- dbeta(x = proposal, shape1 = noisePrior$a, shape2 = noisePrior$b, log = TRUE) -
            dbeta(x = ls, shape1 = noisePrior$a, shape2 = noisePrior$b, log = TRUE)
        
        posterior_nr <- likelihood_nr + prior_nr

        if (!is.na(posterior_nr)) {
            if (log(runif(1)) < posterior_nr) {
                noise_ratio <- proposal
            }
        }
    }

    ## update sigma1t
    K_nosigma <- noise_mix(kern(D, ls^2), noise_ratio)
    K_nosigma_inv <- forceSymmetric(solve(Kt_bin_nosigma))
    a_new <- sigmaPrior$a + 0.5 * length(gpdraw)
    b_new <- sigmaPrior$b + 0.5 * (t(gpdraw) %*% K_nosigma %*% gpdraw)
    sigma.sq <- rinvgamma(n = 1, shape = a_new, scale = b_new[1,1])
    sigma <- sqrt(sigma.sq)

    K <- K_nosigma * sigma.sq
    K_inv <- K_nosigma_inv * (1 / sigma.sq)

    return(list(ls=ls, sigma=sigma, noise_ratio=noise_ratio, K=K, K_inv=K_inv))
}

#' ZINB_GP
#' @description Run the ZINB NNGP model described in https://doi.org/10.1016/j.jspi.2023.106098.
#'
#' @param X Other Predictor variables
#' @param y Zero inflated count response
#' @param coords Spatial coordinates for NNGP
#' @param Vs   Spatially varying predictor variables (e.g. one-hot indication of which location this is for varying intercept), wrapped in sparseMatrix from Matrix R package. Will be multiplied by the spatial random effects for prediction.
#' @param Vt   Temporal varying predictor variables, wrapped in sparseMatrix from Matrix R package. Will be multiplied by the temporal random effects for prediction.
#' @param Ds   Spatial distance matrix, diagonal should be 0, off diagonal is distance between elements i and j in space, inputs to the spatial NNGP kernel
#' @param Dt   Temporal distance matirx, diagonal should be 0, off diagonal is distance between elements i and j in time, inputs to the temporal GP kernel
#' @param M    How many neighbors to allow in the spatial NNGP algorithm, defaults to 10.
#' @param nsim How long to run MCMC in total, must be greater than burn.
#' @param burn How long to run MCMC before saving samples.
#' @param thin How often to save MCMC samples, default is 1, saves every iteration.
#' @param save_ypred Whether or not to output the predicted values at every iteration
#' @param print_iter How often to print the iteration number of the MCMC chain.
#' @param print_progress Whether or not to print the iteration number of the MCMC chain.
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
#' @importFrom MASS glm.nb
#' @importFrom mvtnorm rmvnorm
#' @importFrom BayesLogit rpg
#' @importFrom Matrix bdiag
#' @importFrom Matrix sparseMatrix
#' @importFrom msm rtnorm
#' @importFrom msm dtnorm
#' @importFrom stats rnorm
#' @importFrom LaplacesDemon rinvgamma
#' @importFrom stats runif
#' @importFrom Matrix forceSymmetric
ZINB_GP <- function(X, y, coords, Vs, Vt, Ds, Dt, M = 10, nsim, burn, thin = 1, save_ypred = FALSE, print_iter = 100, print_progress = FALSE) {
    # TODO: Break down the Gibbs sampling and test all steps independently
    # TODO: Remove the need to compute Ds, Dt manually, take in coords for both instead so you can NNGP with large datasets

    # X is the design matrix with dimension N*p
    # x is the vector with length N
    # y is the count response with length N
    n <- nrow(coords) # number of clusters
    N <- nrow(X) # number of observations
    p <- ncol(X) # dimension of alpha and beta
    n_time_points <- ncol(Vt)

    # Use squared exponential kernel
    Ds <- Ds * Ds
    Dt <- Dt * Dt

    # Find reasonable bounds for GP length scales
    param_bounds <- gp_param_bounds(Ds, Dt)
    lsmax <- param_bounds$lsmax
    ltmax <- param_bounds$ltmax

    ##########
    # Priors #
    ##########

    ####### priors for alpha and beta ######
    T0a <- T0b <- diag(100, p)
    sd_r <- 0.02 # TODO: Why so low.

    ####### kernel hyperparameters  ######
    ltPrior <- list(max=ltmax, mh_sd=3, a=1, b=0.001)
    lsPrior <- list(max=lsmax, mh_sd=3, a=1, b=0.001)
    sigmaPrior <- list(a=0.01, b=0.1)
    noisePrior <- list(a=1.5, b=1.5, mh_sd=0.05)

    # Model init
    r <- 1
    y1 <- rep(0, N) # At risk indicator (this is W in paper)
    y1[y > 0] <- 1 # If y>0, then at risk w.p. 1
    q <- rep(.5, N) # 1-p=1/(1+exp(X%*%alpha)), used for updating y1


    #########
    # Inits #
    #########

    #################
    # Fixed Effects #
    #################
    y_ind <- rep(0, N) # convert y to a two class indicator and use logistic regression
    y_ind[y != 0] <- 1

    m1 <- glm(y_ind ~ 0 + X, family = "binomial")
    alpha1 <- m1$coefficients # initial for alpha in the binary component
    m2 <- glm.nb(y[y != 0] ~ 0 + X[y != 0, ])
    beta <- m2$coefficients # initial for beta in the binary component

    eta1 <- X %*% alpha1 + 0
    eta2 <- X %*% beta + 0 # Use all n observations
    p_at_risk <- sigmoid(eta1) # at-risk probability

    q <- 1 / (1 + exp(eta2))
    theta <- p_at_risk * (q^r) / (p_at_risk * (q^r) + 1 - p_at_risk) # Conditional prob that y1=1 given y=0 -- i.e. Pr(chance zero|observed zero)
    y1[y == 0] <- rbinom(sum(y == 0), 1, theta[y == 0]) # If y=0, then draw a "chance zero" w.p. theta, otherwise y1=1

    m1 <- glm(y1 ~ 0 + X, family = "binomial")
    alpha <- as.matrix(m1$coefficients) # initial for alpha in the binary component

    noise_ratio_t1 <- 0.5
    noise_ratio_t2 <- 0.5
    noise_ratio_s1 <- 0.5
    noise_ratio_s2 <- 0.5

    ##########################
    # Spatial Random Effects #
    ##########################
    l1s <- l2s <- 1
    sigma1s <- sigma2s <- sqrt(b_sigmas / (a_sigmas - 1))
    Ks_bin <- sigma1s^2 * noise_mix(kern(Ds, l1s), noise_ratio_s1)
    Ks_bin_inv <- forceSymmetric(solve(Ks_bin))
    Ks_nb <- sigma2s^2 * noise_mix(kern(Ds, l2s), noise_ratio_s2)
    Ks_nb_inv <- forceSymmetric(solve(Ks_nb))
    a <- t(rmvnorm(n = 1, sigma = Ks_bin))
    c <- t(rmvnorm(n = 1, sigma = Ks_nb))

    #################
    # Temporal Random Effects #
    #################
    sigma1t <- sigma2t <- sqrt(b_sigmat / (a_sigmat - 1))
    l1t <- l2t <- 1
    Kt_bin <- sigma1t^2 * noise_mix(kern(Dt, l1t), noise_ratio_t1)
    Kt_bin_inv <- forceSymmetric(solve(Kt_bin))
    Kt_nb <- sigma2t^2 * noise_mix(kern(Dt, l2t), noise_ratio_t2)
    Kt_nb_inv <- forceSymmetric(solve(Kt_nb))
    b <- t(rmvnorm(n = 1, sigma = Kt_bin))
    d <- t(rmvnorm(n = 1, sigma = Kt_nb))

    ############
    # Num Sims #
    ############
    lastit <- (nsim - burn) / thin # Last stored value

    #########
    # Store #
    #########
    Beta <- Alpha <- matrix(0, lastit, p)
    R <- rep(0, lastit)
    A <- C <- matrix(0, lastit, n)
    B <- D <- matrix(0, lastit, n_time_points)
    L1t <- Sigma1t <- Noise1t <- rep(0, lastit)
    L2t <- Sigma2t <- Noise2t <- rep(0, lastit)
    L1s <- Sigma1s <- Noise1s <- rep(0, lastit)
    L2s <- Sigma2s <- Noise2s <- rep(0, lastit)
    if (save_ypred == TRUE) {
        Y_pred <- matrix(NA, lastit, N)
        y1s <- matrix(NA, lastit, N)
    }


    ########
    # MCMC #
    ########
    XV <- cbind(X, Vs, Vt)
    for (i in 1:nsim)   {
        # Ensure these are all updated from last iteration
        Sigma0_bin.inv <- as.matrix(bdiag(
            Ks_bin_inv,
            Kt_bin_inv
        ))
        Sigma0_nb.inv <- as.matrix(bdiag(
            Ks_nb_inv,
            Kt_nb_inv
        ))
        T0_bin <- as.matrix(bdiag(T0a, Sigma0_bin.inv))
        T0_nb <- as.matrix(bdiag(T0b, Sigma0_nb.inv))

        # Update latent variable z
        mu <- X %*% alpha + Vs %*% a + Vt %*% b
        w <- rpg(N, 1, mu[, 1])
        z <- (y1 - 1 / 2) / w

        # Update alpha, a, b
        svd_vinv <- svd(crossprod(sqrt(w) * XV) + T0_bin)
        m <- solve_svd(svd_vinv, (t(sqrt(w) * XV) %*% (sqrt(w) * z)))
        alphaab <- c(mvn_sample_svd(svd_vinv, m))
        alpha <- alphaab[1:p]
        a <- alphaab[(p + 1):(p + n)]
        b <- alphaab[-(1:(p + n))]

        # Update at-risk indicator y1 (W in paper)
        eta1 <- X %*% alpha + Vs %*% a + Vt %*% b
        eta2 <- X %*% beta + Vs %*% c + Vt %*% d # Use all n observations
        pi <- sigmoid(eta1) # at-risk probability
        q <- 1 / (1 + exp(eta2)) # Pr(y=0|y1=1)
        theta <- pi * (q^r) / (pi * (q^r) + 1 - pi) # Conditional prob that y1=1 given y=0 -- i.e. Pr(chance zero|observed zero)
        y1[y == 0] <- rbinom(sum(y == 0), 1, theta[y == 0]) # If y=0, then draw a "chance zero" w.p. theta, otherwise y1=1
        N1 <- sum(y1)

        # Update r, TODO: Check this.
        rnew <- rtnorm(1, r, sd_r, lower = 0) # Treat r as continuous
        ratio <- sum(dnbinom(y[y1 == 1], rnew, q[y1 == 1], log = TRUE)) - sum(dnbinom(y[y1 == 1], r, q[y1 == 1], log = TRUE)) +
            dtnorm(r, rnew, sd_r, 0, log = TRUE) - dtnorm(rnew, r, sd_r, 0, log = TRUE) # Uniform Prior for R
        # Proposal not symmetric
        if (log(runif(1)) < ratio) {
            r <- rnew
        }

        # update l1t, sigma1t, noise_ratio_t1
        out <- update_ls_sigma_noise(l1t, sigma1t, noise_ratio_t1, b, Kt_bin, Dt, ltPrior, sigmaPrior)
        l1t <- out$ls
        sigma1t <- out$sigma
        noise_ratio_t1 <- out$noise_ratio
        Kt_bin <- out$K
        Kt_bin_inv <- out$K_inv

        # update l1s, sigma1s, noise_ratio_s1
        out <- update_ls_sigma_noise(l1s, sigma1s, noise_ratio_s1, a, Ks_bin, Ds, lsPrior, sigmaPrior)
        l1s <- out$ls
        sigma1s <- out$sigma
        noise_ratio_s1 <- out$noise_ratio
        Ks_bin <- out_K
        Ks_bin_inv <- out$K_inv

        # Update beta
        eta <- X[y1 == 1, ] %*% beta + Vs[y1 == 1, ] %*% c + Vt[y1 == 1, ] %*% d
        w <- rpg(N1, y[y1 == 1] + r, eta) # Polya weights
        z <- (y[y1 == 1] - r) / (2 * w)

        # Update beta, c, d
        svd_vinv <- svd(crossprod(sqrt(w) * XV[y1 == 1, ]) + T0_nb)
        m <- solve_svd(svd_vinv, (t(sqrt(w) * XV[y1 == 1, ]) %*% (sqrt(w) * z)))
        betacd <- c(mvn_sample_svd(svd_vinv, m))
        beta <- betacd[1:p]
        c <- betacd[(p + 1):(p + n)]
        d <- betacd[-(1:(p + n))]

        # update l2t, sigma2t, noise_ratio_t2
        out <- update_ls_sigma_noise(l2t, sigam2t, noise_ratio_t2, d, Kt_nb, Dt, ltPrior, sigmaPrior)
        l2t <- out$ls
        sigma2t <- out$sigma
        noise_ratio_t2 <- out$noise_ratio
        Kt_nb <- out$K
        Kt_nb_inv <- out$K_inv

        # update l2s, sigma2s, noise_ratio_s2
        out <- update_ls_sigma_noise(l2s, sigma2s, noise_ratio_s2, c, Ks_nb, Ds, lsPrior, sigmaPrior)
        l2s <- out$ls
        sigma2s <- out$sigma
        noise_ratio_s2 <- out$noise_ratio
        Ks_nb <- out$K
        Ks_nb_inv <- out$K_inv

        # Store
        if ((i > burn) && (i %% thin == 0)) {
            j <- (i - burn) / thin
            
            Alpha[j, ] <- alpha
            Beta[j, ] <- beta # fixed effects

            A[j, ] <- a
            B[j, ] <- b
            C[j, ] <- c
            D[j, ] <- d # random effects
            
            L1t[j] <- l1t
            Noise1t[j] <- noise_ratio_t1
            Sigma1t[j] <- sigma1t # temporal hyperparameters
            
            L2t[j] <- l2t
            Noise2t[j] <- noise_ratio_t2
            Sigma2t[j] <- sigma2t # temporal hyperparameters
            
            L1s[j] <- l1s
            Noise1s[j] <- noise_ratio_s1
            Sigma1s[j] <- sigma1s # spatial hyperparameters
            
            L2s[j] <- l2s
            Noise2s[j] <- noise_ratio_s2
            Sigma2s[j] <- sigma2s # spatial hyperparameters
            
            R[j] <- r
            if (save_ypred) {
                Y_pred[j, ] <- estimate(X, alpha, beta, Vs, Vt, a, b, c, d, r)
                y1s[j, ] <- y1
            }
        }
        if ((i %% print_iter == 0) && (print_progress)) print(i)
    }
    # Put the results into a list
    results <- list(
        Alpha = Alpha, Beta = Beta, A = A, B = B, C = C, D = D,
        L1t = L1t, Sigma1t = Sigma1t, Noise1t=Noise1t, L2t = L2t, Sigma2t = Sigma2t, Noise2t=Noise2t,
        L1s = L1s, Sigma1s = Sigma1s, Noise1s=Noise1s, L2s = L2s, Sigma2s = Sigma2s, Noise2s=Noise2s,
        R = R
    )
    if (save_ypred) {
        temp <- list(Y_pred = Y_pred, at_risk = y1s)
        results <- append(results, temp)
    }
    return(results)
}
