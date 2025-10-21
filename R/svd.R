# SVD helper functions

#' solve_svd
#' @description Use SVD to solve a linear system Ax=b
#' 
#' @param A_svd SVD of A
#' @param b b
#' @return x
solve_svd <- function(A_svd, b, threshold = 1e-12) {
    # Get parts that have nonzero eigenvalues
    rank <- sum(A_svd$d > threshold)
    Upart <- A_svd$u[,1:rank]
    Vpart <- A_svd$v[,1:rank]
    dpart <- A_svd$d[1:rank]
    
    # Compute solution
    result <- Vpart %*% ((1 / dpart) * crossprod(Upart, b))
    return(result)
}

#' mvn_sample_svd
#' @description Use SVD for precision matrix to sample from multivariate normal distribution
#' 
#' @param P_svd SVD of precision matrix
#' @param mu Mean of MVN to draw from
#' @param entropy Draw from MVN of correct size (can be used to draw all mvns at once for efficiency)
mvn_sample_svd <- function(P_svd, mu, entropy = NULL, threshold = 1e-12) {
    if (is.null(entropy)) {
        entropy <- as.matrix(rnorm(length(mu)))
    }

    P_half_svd <- list(u=P_svd$u, d=sqrt(P_svd$d), v=P_svd$u)
    varPart <- solve_svd(P_half_svd, entropy, threshold = sqrt(threshold))
    return(mu + varPart)
}