library(mvtnorm)
library(Matrix)
library(ZINB.GP)

# Generate a number of samples at a number of spatial locations, return design matrices and total number of observations, avg_obs describes how many observations there area at each spatiotemporal point
make_Vs_Vt <- function(num_spatial, num_temporal, avg_obs) {
    n_time_points <- num_temporal # Number of temporal units
    n_unit_mat <- matrix(rpois(num_spatial * n_time_points, avg_obs), nrow = num_spatial, byrow = TRUE) # sample around avg_obs observations per sampling unit (both space and time)

    N <- sum(n_unit_mat)
    id <- c()
    for (i in seq_len(nrow(n_unit_mat))) {
        id <- c(id, rep(i, sum(n_unit_mat[i, ])))
    }

    tp_seq <- c()
    for (j in seq_len(ncol(n_unit_mat))) {
        tp_seq <- c(tp_seq, rep(j, sum(n_unit_mat[, j])))
    }

    # spatial design matrix
    Vs <- as.matrix(sparseMatrix(i = 1:N, j = id, x = rep(1, N)))

    # temporal design matrix
    Vt <- as.matrix(sparseMatrix(i = 1:N, j = tp_seq, x = rep(1, N)))

    return(list(Vs = Vs, Vt = Vt, N = N))
}

make_spatial_effects <- function(phi_nb, phi_bin, sigma_bin_s, sigma_nb_s, coords, noise_ratio) {
    ##########################
    # Spatial Random Effects #
    ##########################
    Ds <- as.matrix(dist(coords))

    # noise_mix(A,ratio) computes ratio*A + (1 - ratio)*I
    Ks_bin <- sigma_bin_s^2 * noise_mix(exp(-(Ds^2) / (phi_bin^2)), noise_ratio)
    Ks_nb <- sigma_nb_s^2 * noise_mix(exp(-(Ds^2) / (phi_nb^2)), noise_ratio)
    a <- t(rmvnorm(n = 1, sigma = Ks_bin))
    c <- t(rmvnorm(n = 1, sigma = Ks_nb))

    return(list(a = a, c = c, Ds = Ds))
}

make_temporal_effects <- function(l1t, l2t, sigma1t, sigma2t, n_time_points, noise_ratio) {
    ###########################
    # Temporal Random Effects #
    ###########################
    w <- matrix(1:n_time_points, ncol = 1) * 5
    Dt <- as.matrix(dist(w))
    
    # noise_mix(A,ratio) computes ratio*A + (1 - ratio)*I
    Kt_bin <- sigma1t^2 * noise_mix(exp(-(Dt^2) / (l1t^2)), noise_ratio)
    Kt_nb <- sigma2t^2 * noise_mix(exp(-(Dt^2) / (l2t^2)), noise_ratio)
    
    b <- t(rmvnorm(n = 1, sigma = Kt_bin))
    d <- t(rmvnorm(n = 1, sigma = Kt_nb))
    
    return(list(b = b, d = d, Dt = Dt))
}

# This file can be run via `Rscript example_usage.R`
dir.create("simulations", showWarnings = FALSE)
obs_nums <- c(2, 4, 8, 16, 24, 32, 48, 64)
for (obs_num in obs_nums)
{
    cat(paste("\nStarting: ", obs_num, "\n"))
    for (kk in 1:100)
    {
        #################
        # Generate Data #
        #################
        num_spatial <- 30
        num_temporal <- 10

        # Get Spatial and temporal design matrices, and total number of observations
        out <- make_Vs_Vt(num_spatial, num_temporal, obs_num)
        Vs <- out$Vs
        Vt <- out$Vt
        N  <- out$N
        spatial_noise <- 0.5
        temporal_noise <- 0.2

        coords <- cbind(runif(num_spatial), runif(num_spatial)) * 100
        x <- rnorm(N, 0, 1)
        X <- as.matrix(x) # Design matrix, can add additional covariates (e.g., race, age, gender)
        X <- cbind(1, X)
        p <- ncol(X)


        phi_nb <- 1
        phi_bin <- 2
        sigma_bin_s <- 1
        sigma_nb_s <- 1
        out <- make_spatial_effects(phi_nb, phi_bin, sigma_bin_s, sigma_nb_s, coords, spatial_noise)
        a <- out$a
        c <- out$c
        Ds <- out$Ds

        l1t <- 2
        l2t <- 3
        sigma1t <- 0.5
        sigma2t <- 0.5
        out <- make_temporal_effects(l1t, l2t, sigma1t, sigma2t, num_temporal, temporal_noise)
        b <- out$b
        d <- out$d
        Dt <- out$Dt

        #################
        # Fixed Effects #
        #################
        # Binomial Part
        alpha <- c(-0.25, 0.25)

        # Count Part
        beta <- c(.5, -.25)

        #######################
        # Binomial Simulation #
        #######################
        eta1 <- X %*% alpha + Vs %*% a + Vt %*% b

        p_at_risk <- exp(eta1) / (1 + exp(eta1)) # 1-pr("structural zero")
        u <- rbinom(N, 1, p_at_risk[, 1]) # at-risk indicator

        #################
        # NB Simulation #
        #################
        res <- Vs %*% c + Vt %*% d
        eta2 <- X[u == 1, ] %*% beta + res[u == 1, ] # Linear predictor for count part
        N1 <- sum(u == 1)

        r <- 1 # NB dispersion
        psi <- exp(eta2) / (1 + exp(eta2)) # Prob of success

        mu <- r * psi / (1 - psi) # NB mean
        y <- rep(0, N) # Response
        y[u == 1] <- rnbinom(N1, r, mu = mu[, 1]) # If at risk, draw from NB
        save(list=c("X", "y", "coords", "Vs", "Vt", "Ds", "Dt", "u", "a", "b", "c", "d"), file=paste("simulations/ob_", obs_num, "_rep_", kk, ".rda", sep=""))
    }
}
