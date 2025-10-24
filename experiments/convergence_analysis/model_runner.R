library(ZINB.GP)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1)
{
    stop("Must provide rep number I am responsible for.", call. = FALSE)
}
repNum <- as.integer(args[1]) + 1

# Load: X, y, coords, Vs, Vt, Ds, Dt, u, a, b, c, d
obs_nums <- c(2, 4, 8, 16, 24, 32, 48, 64)
for (obs_num in obs_nums)
{
    cat(paste("\nRunning", obs_num, "on rep", repNum, "\n"))
    load(paste("simulations/", "ob_", obs_num, "_rep_", repNum, ".rda", sep=""))
    output <- ZINB_GP(X, y, coords, Vs, Vt, Ds, Dt, M = 10, 1000, 200, 1, TRUE, print_progress = TRUE)
    save(list=c("output"), file=paste("simulations/out_", obs_num, "_rep_", repNum, ".rda", sep=""))
    rm(list=c("X", "y", "coords", "Vs", "Vt", "Ds", "Dt", "u", "a", "b", "c", "d", "output"))
}


# # #################
# # # Run the Model #
# # #################
# # # Run for a short time for demo purposes
# cat("Running Model\n")
# output <- ZINB_GP(X, y, coords, Vs, Vt, Ds, Dt, M = 10, 1000, 200, 1, TRUE, print_progress = TRUE)
# predictions <- output$Y_pred
# sim_alpha <- output$Alpha
# sim_beta <- output$Beta
# sim_spatial_noise1 <- output$Noise1s

# # # Examine coefficients for regressions
# cat("\nLogistic Regression Coefficients:\n")
# alpha
# apply(sim_alpha, 2, function(x) quantile(x, probs=c(0.025, 0.5, 0.975)))

# cat("\nNegative Binomial Coefficients:\n")
# beta
# apply(sim_beta, 2, function(x) quantile(x, probs=c(0.025, 0.5, 0.975)))

# # Examine how often various samples are at risk
# cat("\nAt risk 'probabilities':\n")
# at_risk <- output$at_risk
# sim_p_at_risk <- apply(at_risk, 2, mean)
# sim_p_at_risk[1:20]
# cat("\nActual at risk:\n")
# u[1:20]

# cat("\nEstimated Spatial Noise:\n")
# mean(sim_spatial_noise1)
# cat("\nActual Spatial Noise:\n")
# spatial_noise