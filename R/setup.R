# Helper functions that are used before the model is run

#' make_y_Vs_Vt
#' @description Create y, along with spatial and temporal design matrices from an observation matrix.
#' @param obs_matrix s by t matrix, where s is the number of locations, t is the number of times, each entry of the matrix is a nonnegative integer.
#' @return A List of the following values:          
#' \itemize{
#'      \item {\strong{y:} } {Flattened version of the observation matrix, flattened in column-major order.}
#'      \item {\strong{Vs:} } {Spatial design matrix, indicates which elements in y correspond with which positions in space}
#'      \item {\strong{Vt:} } {Temporal design matrix, indicates which elements in y correspond with which positions in time}
#' }
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


#' gp_param_bounds
#' @description Finds reasonable upper/lower bounds on gp parameters ensuring matrices remain pd invertible
#' 
#' @param Ds Spatial distance matrix
#' @param Dt Temporal distance matrix
#' @return Minimum values for l*s, maximum values for l*t 
gp_param_bounds <- function(Ds, Dt, tolerance = 1e-10) {
    smin <- 1
    Ks <- exp(-Ds * smin)
    # TODO: Add try catch to solve for better error messages
    err <- sqrt(sum(((solve(Ks) %*% Ks) - diag(1, nrow=nrow(Ks)))^2))
    while(err < tolerance) {
        smin <- smin / 2
        Ks <- exp(-Ds * smin)
        err <- sqrt(sum(((solve(Ks) %*% Ks) - diag(1, nrow=nrow(Ks)))^2))
    }
    smin <- smin * 2
    if (smin > 0.01) {
        stop("Ds causes ill-conditioned kernel matrix, try increasing distances between spatial coordinates, e.g. Ds <- 100 * Ds")
    }
    
    tmin <- 1
    Kt <- exp(-Dt * tmin)
    err <- sqrt(sum(((solve(Kt) %*% Kt) - diag(1, nrow=nrow(Kt)))^2))
    while(err < tolerance)
    {
        tmin <- tmin / 2
        Kt <- exp(-Dt * tmin)
        err <- sqrt(sum(((solve(Kt) %*% Kt) - diag(1, nrow=nrow(Kt)))^2))
    }
    tmin <- tmin * 2
    if (tmin > 0.01) {
        stop("Dt casuses ill-conditioned kernel matrix, try increasing distances between temporal coordinates, e.g. Dt <- Dt * 100")
    }
    
    ltmax <- sqrt(1 / tmin)
    lsmax <- sqrt(1 / smin)
    return(list(ltmax=ltmax, lsmax=lsmax))
}