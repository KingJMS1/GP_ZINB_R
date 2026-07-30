test_that("prediction inputs preserve unique new space-time levels", {
    observed_space <- rbind(
        c(0, 0),
        c(1000, 0),
        c(0, 1000)
    )
    observed_time <- c(0, 100, 200)
    new_space <- rbind(
        c(250, 250),
        c(250, 250),
        c(750, 750)
    )
    new_time <- c(300, 400, 300)

    inputs <- make_prediction_inputs(
        coords = observed_space,
        time_coords = observed_time,
        coords_new = new_space,
        time_coords_new = new_time
    )

    expect_equal(
        inputs$Vs_new,
        rbind(c(1, 0), c(1, 0), c(0, 1))
    )
    expect_equal(
        inputs$Vt_new,
        rbind(c(1, 0), c(0, 1), c(1, 0))
    )
    expect_equal(dim(inputs$Ds_new), c(4L, 4L))
    expect_equal(dim(inputs$Dt_new), c(4L, 4L))
    expect_equal(inputs$Ds_new, t(inputs$Ds_new))
    expect_equal(inputs$Dt_new, t(inputs$Dt_new))
})

test_that("GP conditioning reproduces a noiseless effect at the same point", {
    set.seed(1)
    draw <- ZINB.GP:::draw_conditional_gp(
        effect = 1.75,
        distance = matrix(0, nrow = 2, ncol = 2),
        length_scale = 1,
        sigma = 2,
        noise_ratio = 1,
        kern = kernel
    )

    expect_equal(draw, 1.75)
})

test_that("predict conditions every active GP and draws new responses", {
    fixture <- simulate_entry_fixture()
    time_coords <- matrix(seq(0, 3000, length.out = 4), ncol = 1)
    new_space <- rbind(
        c(250, 250),
        c(250, 250),
        c(750, 750)
    )
    new_time <- c(3500, 4000, 3500)
    inputs <- make_prediction_inputs(
        coords = fixture$coords,
        time_coords = time_coords,
        coords_new = new_space,
        time_coords_new = new_time
    )
    X_new <- cbind(
        "(Intercept)" = 1,
        x = c(-0.5, 0, 0.5)
    )

    set.seed(101)
    fit <- ZINB_GP(
        X = fixture$X,
        y = fixture$y,
        coords = fixture$coords,
        Vs = fixture$Vs,
        Vt = fixture$Vt,
        Ds = fixture$Ds,
        Dt = fixture$Dt,
        nsim = 5,
        burn = 1,
        thin = 2,
        use_count_gp = TRUE,
        use_inflation_gp = TRUE
    )
    set.seed(202)
    prediction <- do.call(
        stats::predict,
        c(list(object = fit, X = X_new), inputs)
    )

    expect_s3_class(fit, "zinb_gp_fit")
    expect_s3_class(prediction, "zinb_gp_prediction")
    expect_setequal(
        names(prediction),
        c("Y_pred", "eta_at_risk", "eta_count", "A", "B", "C", "D")
    )
    expect_equal(dim(prediction$Y_pred), c(2L, 3L))
    expect_equal(dim(prediction$eta_at_risk), c(2L, 3L))
    expect_equal(dim(prediction$eta_count), c(2L, 3L))
    for (effect in c("A", "B", "C", "D")) {
        expect_equal(dim(prediction[[effect]]), c(2L, 2L))
        expect_true(all(is.finite(prediction[[effect]])))
    }
    expect_true(all(is.finite(prediction$eta_at_risk)))
    expect_true(all(is.finite(prediction$eta_count)))
    expect_true(all(prediction$Y_pred >= 0))
    expect_equal(prediction$Y_pred, floor(prediction$Y_pred))
})
