simulate_entry_fixture <- function() {
    set.seed(271828)

    n_space <- 4
    n_time <- 4
    replicates <- 3

    cells <- expand.grid(
        spatial = seq_len(n_space),
        temporal = seq_len(n_time),
        replicate = seq_len(replicates)
    )

    N <- nrow(cells)
    Vs_full <- diag(n_space)[cells$spatial, , drop = FALSE]
    Vt_full <- diag(n_time)[cells$temporal, , drop = FALSE]

    # Baseline contrasts match the model fitted by ZINB_GP().
    Vs <- Vs_full[, -1, drop = FALSE]
    Vt <- Vt_full[, -1, drop = FALSE]

    coords <- rbind(
        c(0, 0),
        c(1000, 0),
        c(0, 1000),
        c(1000, 1000)
    )
    time_coords <- matrix(seq(0, 3000, length.out = n_time), ncol = 1)

    Ds <- as.matrix(stats::dist(coords))
    Dt <- as.matrix(stats::dist(time_coords))

    noisy_covariance <- function(distance, length_scale, sigma, kappa) {
        correlation <- exp(-(distance^2) / length_scale^2)
        sigma^2 * (
            kappa * correlation +
                (1 - kappa) * diag(nrow(distance))
        )
    }

    spatial_distance <- Ds[-1, -1, drop = FALSE]
    temporal_distance <- Dt[-1, -1, drop = FALSE]

    a <- drop(mvtnorm::rmvnorm(
        1,
        sigma = noisy_covariance(spatial_distance, 900, 0.4, 0.6)
    ))
    c_effect <- drop(mvtnorm::rmvnorm(
        1,
        sigma = noisy_covariance(spatial_distance, 700, 0.4, 0.6)
    ))
    b <- drop(mvtnorm::rmvnorm(
        1,
        sigma = noisy_covariance(temporal_distance, 1200, 0.3, 0.6)
    ))
    d <- drop(mvtnorm::rmvnorm(
        1,
        sigma = noisy_covariance(temporal_distance, 1000, 0.3, 0.6)
    ))

    x <- as.numeric(scale(seq_len(N)))
    X <- cbind("(Intercept)" = 1, x = x)

    alpha <- c(0.4, -0.2)
    beta <- c(3, 0.15)
    dispersion <- 0.25

    eta_at_risk <- drop(X %*% alpha + Vs %*% a + Vt %*% b)
    at_risk <- stats::rbinom(N, 1, stats::plogis(eta_at_risk))

    eta_count <- drop(
        X %*% beta + Vs %*% c_effect + Vt %*% d
    )
    y <- integer(N)
    y[at_risk == 1] <- stats::rnbinom(
        sum(at_risk == 1),
        size = dispersion,
        mu = dispersion * exp(eta_count[at_risk == 1])
    )

    # The initialization requires both observed zeros and positive counts.
    if (!any(y > 0)) {
        y[1] <- 1L
    }
    if (!any(y == 0)) {
        y[1] <- 0L
    }

    list(
        X = X,
        y = y,
        coords = coords,
        Vs = Vs,
        Vt = Vt,
        Ds = Ds,
        Dt = Dt
    )
}
