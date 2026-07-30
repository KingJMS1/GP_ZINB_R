# Helper functions that are used before the model is run

#' make_y_Vs_Vt
#' @description Create y, along with spatial and temporal design matrices from an observation matrix.
#' @param obs_matrix s by t matrix, where s is the number of locations, t is the number of times, each entry of the matrix is a nonnegative integer.
#' @return A list with the following elements:
#' \describe{
#'   \item{y}{The observation matrix flattened in column-major order.}
#'   \item{Vs}{Spatial design matrix mapping observations to locations.}
#'   \item{Vt}{Temporal design matrix mapping observations to times.}
#' }
#' @examples
#' counts <- matrix(c(0, 1, 2, 0, 1, 3), nrow = 2)
#' design <- make_y_Vs_Vt(counts)
#' design$y
#' design$Vs
#' design$Vt
#' @export
#' @importFrom Matrix sparseMatrix
make_y_Vs_Vt <- function(obs_matrix) {
    # Assumes rows are space, columns are time
    obs_matrix <- as.matrix(unname(obs_matrix))
    n_temporal <- ncol(obs_matrix)
    n_spatial <- nrow(obs_matrix)
    N <- n_spatial * n_temporal

    # Create y, Vs, Vt
    y <- as.vector(obs_matrix)
    Vt <- as.matrix(sparseMatrix(i = 1:N, j = rep(1:n_temporal, each=n_spatial), x=rep(1, N)))
    Vs <- as.matrix(sparseMatrix(i = 1:N, j = rep(1:n_spatial, n_temporal), x=rep(1,N)))

    # Sacrifice to the intercept gods
    Vt <- Vt[,2:ncol(Vt)]
    Vs <- Vs[,2:ncol(Vs)]

    return(list(Vs=Vs, Vt=Vt, y=y))
}

#' Build Inputs for Prediction at New GP Coordinates
#'
#' Creates the augmented distance matrices needed for Gaussian-process
#' conditioning and indicator matrices mapping prediction observations to
#' unique new spatial and temporal levels. The first observed level in each
#' dimension is the baseline omitted by `ZINB_GP()` and is therefore excluded
#' from the augmented distance matrix.
#'
#' @param coords Observed spatial coordinate matrix, including the baseline as
#'   its first row. Use `NULL` for a temporal-only model.
#' @param time_coords Observed temporal coordinate matrix, including the
#'   baseline as its first row. A numeric vector is treated as a one-column
#'   matrix. Use `NULL` for a spatial-only model.
#' @param coords_new Spatial coordinates for the prediction observations, with
#'   one row per row of the new fixed-effect design matrix. Repeated rows share
#'   one predicted spatial random effect.
#' @param time_coords_new Temporal coordinates for the prediction observations,
#'   with one row per row of the new fixed-effect design matrix. A numeric
#'   vector is treated as a one-column matrix. Repeated rows share one predicted
#'   temporal random effect.
#'
#' @return A list containing `Ds_new`, `Dt_new`, `Vs_new`, and `Vt_new`.
#'   Each non-`NULL` distance matrix covers the observed nonbaseline levels
#'   followed by the unique new levels. Each design matrix has one row per
#'   prediction observation and one column per unique new level.
#' @examples
#' observed_space <- rbind(c(0, 0), c(1000, 0), c(0, 1000))
#' observed_time <- c(0, 100)
#' new_space <- rbind(c(500, 500), c(500, 500), c(750, 250))
#' new_time <- c(150, 200, 150)
#' prediction_inputs <- make_prediction_inputs(
#'   coords = observed_space,
#'   time_coords = observed_time,
#'   coords_new = new_space,
#'   time_coords_new = new_time
#' )
#' prediction_inputs$Vs_new
#' prediction_inputs$Vt_new
#' @export
make_prediction_inputs <- function(
    coords = NULL,
    time_coords = NULL,
    coords_new = NULL,
    time_coords_new = NULL
) {
    # Normalize one coordinate dimension to the matrix form used by dist().
    as_coordinate_matrix <- function(x, name) {
        if (is.null(x)) {
            return(NULL)
        }

        if (is.vector(x) && is.numeric(x)) {
            x <- matrix(x, ncol = 1L)
        } else {
            x <- as.matrix(x)
        }

        if (nrow(x) < 1L || ncol(x) < 1L || any(!is.finite(x))) {
            stop("`", name, "` must be a non-empty matrix of finite coordinates.")
        }
        storage.mode(x) <- "double"

        x
    }

    # Full-precision keys let repeated coordinate rows share one GP level.
    coordinate_key <- function(x) {
        apply(x, 1L, function(row) {
            paste(sprintf("%.17g", row), collapse = "\r")
        })
    }

    # Build the conditioning distance matrix and prediction-level indicators
    # for either the spatial or temporal dimension.
    build_dimension <- function(observed, new, observed_name, new_name) {
        observed <- as_coordinate_matrix(observed, observed_name)
        new <- as_coordinate_matrix(new, new_name)

        if (is.null(observed) && is.null(new)) {
            return(list(distance = NULL, design = NULL, n = NULL))
        }

        if (is.null(observed) || is.null(new)) {
            stop("`", observed_name, "` and `", new_name, "` must be supplied together.")
        }

        if (ncol(observed) != ncol(new)) {
            stop("`", observed_name, "` and `", new_name, "` must have the same number of columns.")
        }

        if (nrow(observed) < 2L) {
            stop("`", observed_name, "` must include a baseline and at least one nonbaseline level.")
        }

        # Retain each new level once, then map observation rows back to it.
        keys <- coordinate_key(new)
        keep <- !duplicated(keys)
        unique_new <- new[keep, , drop = FALSE]
        level <- match(keys, keys[keep])
        design <- diag(nrow(unique_new))[level, , drop = FALSE]

        # Fitted random effects omit the first observed (baseline) level.
        augmented <- rbind(observed[-1L, , drop = FALSE], unique_new)

        list(
            distance = as.matrix(stats::dist(augmented)),
            design = design,
            n = nrow(new)
        )
    }

    # Construct each active GP dimension independently.
    spatial <- build_dimension(coords, coords_new, "coords", "coords_new")
    temporal <- build_dimension(
        time_coords,
        time_coords_new,
        "time_coords",
        "time_coords_new"
    )

    # Spatial and temporal coordinates must describe the same prediction rows.
    prediction_rows <- c(spatial$n, temporal$n)
    prediction_rows <- prediction_rows[!vapply(prediction_rows, is.null, logical(1))]

    if (length(prediction_rows) == 0L) {
        stop("At least one pair of observed and new coordinate inputs is required.")
    }

    if (length(unique(prediction_rows)) != 1L) {
        stop("`coords_new` and `time_coords_new` must have the same number of rows.")
    }

    # Keep names aligned with the arguments accepted by predict().
    list(
        Ds_new = spatial$distance,
        Dt_new = temporal$distance,
        Vs_new = spatial$design,
        Vt_new = temporal$design
    )
}


#' gp_param_bounds
#' @description Finds reasonable upper/lower bounds on gp parameters ensuring matrices remain pd invertible
#' 
#' @param Ds Spatial distance matrix
#' @param Dt Temporal distance matrix
#' @param kernel Kernel function accepting a distance matrix and a length scale.
#' @param tolerance Numerical tolerance used to assess matrix invertibility.
#' @return A list containing the upper bounds `ltmax` and `lsmax` for the
#'   temporal and spatial length scales, respectively.
gp_param_bounds <- function(Ds, Dt, kernel, tolerance = 1e-10) {
    smin <- 1
    Ks <- kernel(Ds, 1 / smin)
    # TODO: Add try catch to solve for better error messages
    err <- sqrt(sum(((solve(Ks) %*% Ks) - diag(1, nrow=nrow(Ks)))^2))
    while(err < tolerance) {
        smin <- smin / 2
        Ks <- kernel(Ds, 1 / smin)
        err <- sqrt(sum(((solve(Ks) %*% Ks) - diag(1, nrow=nrow(Ks)))^2))
    }
    smin <- smin * 2
    if (smin > 0.01) {
        stop("Ds causes ill-conditioned kernel matrix, try increasing distances between spatial coordinates, e.g. Ds <- 100 * Ds")
    }
    
    tmin <- 1
    Kt <- kernel(Dt, 1 / tmin)
    err <- sqrt(sum(((solve(Kt) %*% Kt) - diag(1, nrow=nrow(Kt)))^2))
    while(err < tolerance)
    {
        tmin <- tmin / 2
        Kt <- kernel(Dt, 1 / tmin)
        err <- sqrt(sum(((solve(Kt) %*% Kt) - diag(1, nrow=nrow(Kt)))^2))
    }
    tmin <- tmin * 2
    if (tmin > 0.01) {
        stop("Dt casuses ill-conditioned kernel matrix, try increasing distances between temporal coordinates, e.g. Dt <- Dt * 100")
    }
    
    ltmax <- (1 / tmin)
    lsmax <- (1 / smin)
    return(list(ltmax=ltmax, lsmax=lsmax))
}
