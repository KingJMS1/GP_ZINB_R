entry_variations <- list(
    list(
        name = "spatial and temporal GPs in both components",
        use_space = TRUE,
        use_time = TRUE,
        use_count_gp = TRUE,
        use_inflation_gp = TRUE,
        effects = c("A", "B", "C", "D")
    ),
    list(
        name = "spatial and temporal GPs in the count component",
        use_space = TRUE,
        use_time = TRUE,
        use_count_gp = TRUE,
        use_inflation_gp = FALSE,
        effects = c("C", "D")
    ),
    list(
        name = "spatial and temporal GPs in the inflation component",
        use_space = TRUE,
        use_time = TRUE,
        use_count_gp = FALSE,
        use_inflation_gp = TRUE,
        effects = c("A", "B")
    ),
    list(
        name = "spatial GPs in both components",
        use_space = TRUE,
        use_time = FALSE,
        use_count_gp = TRUE,
        use_inflation_gp = TRUE,
        effects = c("A", "C")
    ),
    list(
        name = "a spatial GP in the count component",
        use_space = TRUE,
        use_time = FALSE,
        use_count_gp = TRUE,
        use_inflation_gp = FALSE,
        effects = "C"
    ),
    list(
        name = "a spatial GP in the inflation component",
        use_space = TRUE,
        use_time = FALSE,
        use_count_gp = FALSE,
        use_inflation_gp = TRUE,
        effects = "A"
    ),
    list(
        name = "temporal GPs in both components",
        use_space = FALSE,
        use_time = TRUE,
        use_count_gp = TRUE,
        use_inflation_gp = TRUE,
        effects = c("B", "D")
    ),
    list(
        name = "a temporal GP in the count component",
        use_space = FALSE,
        use_time = TRUE,
        use_count_gp = TRUE,
        use_inflation_gp = FALSE,
        effects = "D"
    ),
    list(
        name = "a temporal GP in the inflation component",
        use_space = FALSE,
        use_time = TRUE,
        use_count_gp = FALSE,
        use_inflation_gp = TRUE,
        effects = "B"
    )
)

for (variation in entry_variations) {
    test_that(variation$name, {
        fixture <- simulate_entry_fixture()

        args <- list(
            X = fixture$X,
            y = fixture$y,
            coords = fixture$coords,
            nsim = 5,
            burn = 1,
            thin = 2,
            save_ypred = TRUE,
            print_progress = FALSE,
            use_count_gp = variation$use_count_gp,
            use_inflation_gp = variation$use_inflation_gp
        )

        if (variation$use_space) {
            args$Vs <- fixture$Vs
            args$Ds <- fixture$Ds
        }
        if (variation$use_time) {
            args$Vt <- fixture$Vt
            args$Dt <- fixture$Dt
        }

        set.seed(variation$use_space * 100 + variation$use_time * 10 +
            variation$use_count_gp * 2 + variation$use_inflation_gp)
        fit <- do.call(ZINB_GP, args)

        expected_names <- c("Alpha", "Beta")
        if (variation$use_inflation_gp && variation$use_space) {
            expected_names <- c(
                expected_names, "A", "L1s", "Sigma1s", "Noise1s"
            )
        }
        if (variation$use_inflation_gp && variation$use_time) {
            expected_names <- c(
                expected_names, "B", "L1t", "Sigma1t", "Noise1t"
            )
        }
        if (variation$use_count_gp && variation$use_space) {
            expected_names <- c(
                expected_names, "C", "L2s", "Sigma2s", "Noise2s"
            )
        }
        if (variation$use_count_gp && variation$use_time) {
            expected_names <- c(
                expected_names, "D", "L2t", "Sigma2t", "Noise2t"
            )
        }
        expected_names <- c(
            expected_names, "R", "Y_pred", "at_risk"
        )

        expect_type(fit, "list")
        expect_setequal(names(fit), expected_names)
        expect_true(all(variation$effects %in% names(fit)))

        expect_equal(nrow(fit$Alpha), 2)
        expect_equal(nrow(fit$Beta), 2)
        expect_equal(length(fit$R), 2)
        expect_equal(dim(fit$Y_pred), c(2L, length(fixture$y)))
        expect_equal(dim(fit$at_risk), c(2L, length(fixture$y)))
        expect_true(all(fit$Y_pred >= 0))
        expect_equal(fit$Y_pred, floor(fit$Y_pred))
        expect_true(all(is.finite(fit$Alpha)))
        expect_true(all(is.finite(fit$Beta)))
        expect_true(all(is.finite(fit$R)))

        spatial_effects <- intersect(variation$effects, c("A", "C"))
        for (effect in spatial_effects) {
            expect_equal(dim(fit[[effect]]), c(2L, ncol(fixture$Vs)))
            expect_true(all(is.finite(fit[[effect]])))
        }

        temporal_effects <- intersect(variation$effects, c("B", "D"))
        for (effect in temporal_effects) {
            expect_equal(dim(fit[[effect]]), c(2L, ncol(fixture$Vt)))
            expect_true(all(is.finite(fit[[effect]])))
        }
    })
}

test_that("combinations with no active GP are explicitly unsupported", {
    fixture <- simulate_entry_fixture()
    base_args <- list(
        X = fixture$X,
        y = fixture$y,
        coords = fixture$coords,
        nsim = 4,
        burn = 2,
        thin = 1
    )

    unsupported <- list(
        list(use_count_gp = TRUE, use_inflation_gp = TRUE),
        list(use_count_gp = TRUE, use_inflation_gp = FALSE),
        list(use_count_gp = FALSE, use_inflation_gp = TRUE),
        list(use_count_gp = FALSE, use_inflation_gp = FALSE),
        list(
            Vs = fixture$Vs,
            Ds = fixture$Ds,
            use_count_gp = FALSE,
            use_inflation_gp = FALSE
        ),
        list(
            Vt = fixture$Vt,
            Dt = fixture$Dt,
            use_count_gp = FALSE,
            use_inflation_gp = FALSE
        ),
        list(
            Vs = fixture$Vs,
            Vt = fixture$Vt,
            Ds = fixture$Ds,
            Dt = fixture$Dt,
            use_count_gp = FALSE,
            use_inflation_gp = FALSE
        )
    )

    for (args in unsupported) {
        expect_output(
            fit <- do.call(ZINB_GP, c(base_args, args)),
            "must specify at least 1 GP"
        )
        expect_null(fit)
    }
})

test_that("spatial coordinates, designs, and distances must align", {
    fixture <- simulate_entry_fixture()
    args <- list(
        X = fixture$X,
        y = fixture$y,
        coords = fixture$coords[-1, , drop = FALSE],
        Vs = fixture$Vs,
        Vt = fixture$Vt,
        Ds = fixture$Ds,
        Dt = fixture$Dt,
        nsim = 2,
        burn = 1,
        use_count_gp = TRUE,
        use_inflation_gp = TRUE
    )

    expect_error(
        do.call(ZINB_GP, args),
        "nrow\\(coords\\) must equal ncol\\(Vs\\) \\+ 1"
    )

    args$coords <- fixture$coords
    args$Ds <- fixture$Ds[, -1, drop = FALSE]
    expect_error(
        do.call(ZINB_GP, args),
        "Ds must be a square matrix"
    )

    args$Ds <- fixture$Ds
    args$Ds[1, 2] <- args$Ds[1, 2] + 1
    expect_error(
        do.call(ZINB_GP, args),
        "Ds must be symmetric"
    )
})

test_that("zero-mixture probabilities remain finite for extreme predictors", {
    pi <- ZINB.GP:::sigmoid(c(-1000, 1000))
    q <- ZINB.GP:::sigmoid(-c(-1000, 1000))
    theta <- ZINB.GP:::zero_at_risk_probability(pi, q, r = 2)

    expect_true(all(is.finite(theta)))
    expect_true(all(theta >= 0 & theta <= 1))
})
