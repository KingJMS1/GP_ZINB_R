source("../R/ZINB_NNGP_final.R")

# Number of locations and time points
n_locs <- 30    # spatial locations
n_times <- 1000   # time points
n <- n_locs * n_times  # total observations

# Generate spatial coordinates for 10 locations
coords <- cbind(runif(n_locs, 0, 10), runif(n_locs, 0, 10))

# Create a time vector (assuming evenly spaced time points)
time_points <- seq(1, n_times)

# Compute the full spatial distance matrix for all observations
Ds <- as.matrix(dist(coords))
Dt <- as.matrix(dist(time_points)) * 100

l1t <- 1000
l2t <- 1
l1s <- 1
l2s <- 1
sigma1t <- 10000
sigma2t <- 1
sigma1s <- 1
sigma2s <- 1
Kt_bin <- sigma1t^2 * normalize(exp(-Dt / (l1t^2)))
Kt_nb <- sigma2t^2 * normalize(exp(-Dt / (l2t^2)))
Ks_bin <- sigma1s^2 * normalize(exp(-l1s * Ds))
Ks_nb <- sigma2s^2 * normalize(exp(-l2s * Ds))

print(gp_param_bounds(Ds, Dt))