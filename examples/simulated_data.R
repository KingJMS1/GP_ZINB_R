library(mvtnorm)
library(Matrix)
library(ZINB.GP)

# This file can be run via `Rscript examples/simulated_data.R`.
set.seed(20260729)

# Generate observations at space-time cells and return the full indicator
# matrices. The model uses the first spatial and temporal levels as baselines,
# so their columns are removed below after the distance matrices are formed.
make_Vs_Vt <- function(num_spatial, num_temporal, avg_obs) {
    n_unit_mat <- matrix(
        rpois(num_spatial * num_temporal, avg_obs),
        nrow = num_spatial
    )

    N <- sum(n_unit_mat)
    cells <- expand.grid(
        spatial = seq_len(num_spatial),
        temporal = seq_len(num_temporal)
    )
    id <- rep(cells$spatial, times = as.vector(n_unit_mat))
    tp_seq <- rep(cells$temporal, times = as.vector(n_unit_mat))

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
    Ds_reduced <- Ds[-1, -1, drop = FALSE]

    # noise_mix(A,ratio) computes ratio*A + (1 - ratio)*I
    Ks_bin <- sigma_bin_s^2 * noise_mix(
        exp(-(Ds_reduced^2) / (phi_bin^2)),
        noise_ratio
    )
    Ks_nb <- sigma_nb_s^2 * noise_mix(
        exp(-(Ds_reduced^2) / (phi_nb^2)),
        noise_ratio
    )
    a <- t(rmvnorm(n = 1, sigma = Ks_bin))
    c <- t(rmvnorm(n = 1, sigma = Ks_nb))

    return(list(a = a, c = c, Ds = Ds))
}

make_temporal_effects <- function(l1t, l2t, sigma1t, sigma2t, n_time_points, noise_ratio) {
    ###########################
    # Temporal Random Effects #
    ###########################
    w <- matrix(0:(n_time_points - 1), ncol = 1) * 50
    Dt <- as.matrix(dist(w))
    Dt_reduced <- Dt[-1, -1, drop = FALSE]
    
    # noise_mix(A,ratio) computes ratio*A + (1 - ratio)*I
    Kt_bin <- sigma1t^2 * noise_mix(
        exp(-(Dt_reduced^2) / (l1t^2)),
        noise_ratio
    )
    Kt_nb <- sigma2t^2 * noise_mix(
        exp(-(Dt_reduced^2) / (l2t^2)),
        noise_ratio
    )
    
    b <- t(rmvnorm(n = 1, sigma = Kt_bin))
    d <- t(rmvnorm(n = 1, sigma = Kt_nb))
    
    return(list(b = b, d = d, Dt = Dt))
}

#################
# Generate Data #
#################
num_spatial <- 30
num_temporal <- 10

# Get Spatial and temporal design matrices, and total number of observations
out <- make_Vs_Vt(num_spatial, num_temporal, 20)
Vs <- out$Vs
Vt <- out$Vt
N  <- out$N
spatial_noise <- 0.5
temporal_noise <- 0.2

coords <- cbind(runif(num_spatial), runif(num_spatial)) * 1000
x <- rnorm(N, 0, 1)
X <- as.matrix(x) # Design matrix, can add additional covariates (e.g., race, age, gender)
X <- cbind(1, X)
p <- ncol(X)


phi_nb <- 250
phi_bin <- 350
sigma_bin_s <- 1
sigma_nb_s <- 1
out <- make_spatial_effects(phi_nb, phi_bin, sigma_bin_s, sigma_nb_s, coords, spatial_noise)
a <- out$a
c <- out$c
Ds <- out$Ds

l1t <- 100
l2t <- 150
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
beta <- c(0.5, -0.25)

# ZINB_GP_orig() uses baseline contrasts: the first spatial and temporal
# effects are fixed at zero, and the GP priors are defined on the remaining
# effects. The distance matrices remain full because ZINB_GP_orig() removes
# their first rows and columns internally.
Vs <- Vs[, -1, drop = FALSE]
Vt <- Vt[, -1, drop = FALSE]

#######################
# Binomial Simulation #
#######################
eta1 <- X %*% alpha + Vs %*% a + Vt %*% b

p_at_risk <- plogis(eta1) # 1 - Pr("structural zero")
u <- rbinom(N, 1, p_at_risk[, 1]) # at-risk indicator

#################
# NB Simulation #
#################
res <- Vs %*% c + Vt %*% d
eta2 <- X[u == 1, ] %*% beta + res[u == 1, ] # Linear predictor for count part
N1 <- sum(u == 1)

r <- 1 # NB dispersion
psi <- plogis(eta2) # NB success parameter under the package convention

mu <- r * psi / (1 - psi) # NB mean
y <- rep(0, N) # Response
y[u == 1] <- rnbinom(N1, r, mu = mu[, 1]) # If at risk, draw from NB


#################
# Run the Model #
#################
# Run for a short time for demo purposes
cat("Running Model\n")
output <- ZINB_GP(
    X = X,
    y = y,
    coords = coords,
    Vs = Vs,
    Vt = Vt,
    Ds = Ds,
    Dt = Dt,
    nsim = 1000,
    burn = 200,
    thin = 1,
    save_ypred = TRUE,
    print_progress = TRUE,
    use_count_gp = TRUE,
    use_inflation_gp = TRUE
)
predictions <- output$Y_pred
sim_alpha <- output$Alpha
sim_beta <- output$Beta
sim_spatial_noise1 <- output$Noise1s

# Examine coefficients for regressions
cat("\nLogistic Regression Coefficients:\n")
alpha
apply(sim_alpha, 2, function(x) quantile(x, probs=c(0.025, 0.5, 0.975)))

cat("\nNegative Binomial Coefficients:\n")
beta
apply(sim_beta, 2, function(x) quantile(x, probs=c(0.025, 0.5, 0.975)))

# Examine how often various samples are at risk
cat("\nAt risk 'probabilities':\n")
at_risk <- output$at_risk
sim_p_at_risk <- apply(at_risk, 2, mean)
sim_p_at_risk[1:20]
cat("\nActual at risk:\n")
u[1:20]

cat("\nEstimated Spatial Noise:\n")
mean(sim_spatial_noise1)
cat("\nActual Spatial Noise:\n")
spatial_noise
