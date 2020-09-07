#' Gradient for weighted Unifrac
#'
#' If d_{x,T}(y) denotes the weighted Unifrac distance between x and y
#' with tree T, computes the gradient of d_{x,T} evaluated at y.
#' @param x A p-vector.
#' @param y A p-vector.
#' @param A A p x B matrix, with A_{ij} = 1 indicating that taxon i
#' descends from branch j.
#'
#' @return A p-vector giving the gradient.
#' @export
wuf_dist_deriv <- function(x, y, A, L) {
    p = length(y)
    sy = sum(y)
    sx = sum(x)
    py = y / sy
    px = x / sx
    by = py %*% A
    bx = px %*% A
    Jh = sy^(-1) * (diag(rep(1, p)) + outer(py, rep(1, p)))
    Jg = t(A)
    Jf = matrix(ifelse(by >= bx, 1, -1) * L, nrow = 1)
    grad = Jf %*% Jg %*% Jh
    return(as.vector(grad))
}


#' Weighted Unifrac distance
#'
#' @param X A n x p matrix, n samples, p taxa.
#' @param A A p x B marix, A_{ij} = 1 indicating that taxon i descends
#' from branch j.
#'
#' @return An n x n distance matrix, with element i,j giving the
#' weighted Unifrac distance between X[i,] and X[j,].
#' @export
wuf_dist <- function(X, A, L) {
    P = sweep(X, MARGIN = 1, STATS = rowSums(X), FUN = "/")
    bp_matrix = P %*% sweep(A, MARGIN = 2, STATS = L, FUN = "*")
    return(dist(bp_matrix, method = "manhattan"))
}

#' Unweighted Unifrac distance
#'
#' @param X A n x p matrix, n samples, p taxa.
#' @param A A p x B marix, A_{ij} = 1 indicating that taxon i descends
#' from branch j.
#' @param L A vector giving branch lengths.
#'
#' @return An n x n distance matrix, with element i,j giving the
#' weighted Unifrac distance between X[i,] and X[j,].
#' @export
uf_dist <- function(X, A, L) {
    bp_matrix = X %*% A
    br_ind_matrix = bp_matrix > 0
    br_ind_matrix_scaled = sweep(br_ind_matrix, MARGIN = 2, STATS = L, FUN = "*")
    return(dist(br_ind_matrix_scaled, method = "manhattan") / sum(L))
}

#' Generalized gradient for unweighted Unifrac
#'
#' If d_{x,T}(y) denotes the unweighted Unifrac distance between x and y
#' with tree T, computes the generalized gradient of d_{x,T} evaluated at y.
#' @param x A p-vector.
#' @param y A p-vector.
#' @param A A p x B matrix, with A_{ij} = 1 indicating that taxon i
#' descends from branch j.
#' @param L A vector of branch lengths.
#' @param delta_min Smallest perturbation.
#' @param positive Logical, positive perturbation?
#'
#' @return A p-vector giving the gradient.
#' @export
uf_dist_deriv <- function(x, y, A, L, delta_min = 1, positive = TRUE,
                          return_dist_difference = FALSE) {
    y = as.vector(y)
    x = as.vector(x)
    if(positive) {
        rho = ifelse(y == 0, 1, Inf)
    } else {
        rho = ifelse(y == 0, Inf, y)
    }
    dxy = as.numeric(uf_dist(rbind(x, y), A, L))
    deriv = sapply(seq_along(y), function(i) {
        if(rho[i] == Inf) return(0)
        if(positive) {
            delta = uf_delta(x, y, A, L, i, rho = rho[i])
        } else {
            delta = -uf_delta(x, y, A, L, i, rho = -rho[i])
        }
        if(return_dist_difference) {
            return((-1) * delta * (delta + 2 * dxy))
        } else {
            return(delta_min * rho[i]^(-1) * (-1) * delta * (delta + 2 * dxy))
        }
    })
    return(deriv)
}


#' Pre-computations for unifrac biplots
#'
#' Saves matrices that we need to access many times when computing the
#' unifrac biplot axes.
#'
#' @param X A samples by species matrix
#' @param A The branches by species ancestry matrix
#' @param L Branch-length vector
#' 
#' @export
uf_lb_computations <- function(X, A, L) {
    A_lgc = A != 0
    X_lgc_sparse = as(X != 0, "lgCMatrix")
    X_lgc = X != 0
    sample_indicator_list = list()
    for(i in 1:nrow(X)) {
        sample_indicator_list[[i]] = Matrix::t(X_lgc_sparse[i,,drop=FALSE])
    }
    sample_species_combinations = expand.grid(1:nrow(X), 1:ncol(X))
    ## this is a list where sample_sensitivity_indicators[[j]][[k]]
    ## has A^{(jk)}, the matrix indicating whether the contribution of
    ## the bth branch to the unifrac distance between sample j and any
    ## other sample will change if taxon k is increased
    sample_sensitivity_indicators = list()
    for(i in 1:nrow(X)) {
        sample_sensitivity_indicators[[i]] = list()
    }
    for(i in 1:nrow(sample_species_combinations)) {
        sample = as.numeric(sample_species_combinations[i,1])
        species = as.numeric(sample_species_combinations[i,2])
        d = as.logical(A_lgc[,species]) & !as.logical(A_lgc %*% X_lgc[sample,])
        diagonal = Diagonal(n = nrow(A), x = d)
        sample_sensitivity_indicators[[sample]][[species]] = drop0(diagonal %*% A != 0)
    }
    uf_distances = dist(sweep((X_lgc %*% Matrix::t(A_lgc) > 0), MARGIN = 2, STATS = L, FUN = "*"), "manhattan") / sum(L)
    return(list(
        sample_indicator_list = sample_indicator_list,
        sample_sensitivity_indicators = sample_sensitivity_indicators,
        uf_distances = uf_distances,
        L = L,
        T = sum(L)
    ))

}

#' Computes biplot axes around a point
#'
#' @param X A samples by species data matrix.
#' @param sample_idx The point around which to compute biplot axes.
#' @param species_indices The species for which to compute biplot axes.
#' @param delta_min Minimum value for perturbation.
#' @param positive Logical, perturbation in the positive direction?
#' @param uf_cached The output from uf_lb_computations
#' @param mds_matrices The output from make_mds_matrices
#'
#' @return A species by n_axes matrix with biplot axes.
#' @export
uf_lb <- function(X, sample_idx, species_indices = 1:ncol(X),
                           delta_min = 1, positive = TRUE, uf_cached,
                           mds_matrices, n_axes = 2) {
    if(class(X) != "Matrix") {
        X = as(X, "Matrix")
    }
    generalized_gradients =
        make_uf_generalized_gradients(X, sample_idx, species_indices,
                                      delta_min, positive = positive,
                                      uf_cached = uf_cached)
    biplot_axes =
        t(.5 * diag(mds_matrices$Lambda[1:n_axes]^(-1)) %*%
            t(mds_matrices$Y[,1:n_axes,drop=FALSE]) %*% generalized_gradients)
    return(biplot_axes)
}

#' Computes change in uf distance after perturbation
#'
#' Computes d(xi,xj+rho e_k) - d(xi,xj)
#'
#' @param uf_cache The output from uf_lb_computations.
delta_from_A_lb <- function(uf_cache, i, j, k) {
    x_j = x_j_pert = uf_cache$sample_indicator_list[[j]]
    ## next five lines should be equivalent to x_j_pert[k,1] = TRUE
    if(all(x_j_pert@i != (k-1))) {
        x_j_pert@x = rep(TRUE, length(x_j_pert@x) + 1)
        x_j_pert@p[2] = x_j_pert@p[2] + 1L
        x_j_pert@i = sort(c(x_j_pert@i, as.integer(k-1)))
    }
    A_sens = uf_cache$sample_sensitivity_indicators[[j]][[k]]
    A_sens_x_i = A_sens %*% uf_cache$sample_indicator_list[[i]]
    A_sens_x_j = A_sens %*% x_j
    A_sens_x_j_pert = A_sens %*% x_j_pert
    d1 = uf_num_from_indicators(A_sens_x_i, A_sens_x_j, uf_cache$L)
    d2 = uf_num_from_indicators(A_sens_x_i, A_sens_x_j_pert, uf_cache$L)
    (d2 - d1) / sum(uf_cache$T)
}

sym_diff <- function(a,b) unique(c(setdiff(a,b), setdiff(b,a)))

#' Unifrac numerator
#'
#' More efficient, less transparent version of numerator computation.
#'
#' @param sm1 A p x 1 sparse matrix of class "dgCMatrix". Assumed to
#' have zeros dropped.
#' @param sm2 A p x 1 sparse matrix of class "dgCMatrix" Assumed to
#' have zeros dropped.
#' @param L A numeric p-vector.
#' @return Sum of unique branch lengths.
uf_num_from_indicators <- function(sm1, sm2, L) {
    #sm1 = as(sm1, "dgTMatrix")
    #sm2 = as(sm2, "dgTMatrix")
    #sm1 = drop0(sm1)
    #sm2 = drop0(sm2)
    sum_indices = sym_diff(sm1@i, sm2@i)
    ## +1 because Matrix uses zero indexing
    sum(L[sum_indices + 1])
}


make_uf_generalized_gradients <- function(X, sample_idx, species_indices,
                                          delta_min, positive, uf_cached) {
    X_lgc = as(X != 0, "lgTMatrix")
    if(positive) {
        rho = ifelse(!X_lgc[sample_idx,], delta_min, Inf)
    } else {
        rho = -ifelse(!X_lgc[sample_idx,], Inf, X[sample_idx,])
    }
    uf_distances = as.matrix(uf_cached$uf_distances)
    Xt_lgc = Matrix::t(X_lgc)
    generalized_gradients = matrix(0, nrow(X), length(species_indices))
    for(i in 1:nrow(X)) {
        for(k in 1:length(species_indices)) {
            species_idx = species_indices[k]
            if(rho[species_idx] == Inf) {
                generalized_gradients[i, k] = 0
            } else {
                dxy = uf_distances[i, sample_idx]
                delta = delta_from_A_lb(uf_cache, i, sample_idx, species_idx)
                generalized_gradients[i, species_idx] =
                    delta_min / rho[species_idx] * (-1) * delta * (delta + 2 * dxy)
            }
        }
    }
    return(generalized_gradients)
}
#' Computes delta for unifrac
#'
#' Computes d(x,y) - d(x,y + rho e_i) for d the unweighted Unifrac
#' distance.
#'
#' @param A A species by branches ancestry matrix.
#' @param L A branch length vector.
uf_delta <- function(x, y, A, L, i, rho) {
    i_s_ancestors = which(A[i,] > 0)
    A_sub = A[,i_s_ancestors,drop=FALSE]
    L_sub = L[i_s_ancestors]
    y_pert = y
    y_pert[i] = y_pert[i] + rho
    original = sum(abs(((t(x) %*% A_sub) > 0) - ((t(y) %*% A_sub) > 0)) * L_sub)
    ##new = sum(abs(((t(x) %*% A_sub) > 0) - (t(y_pert) %*% A_sub) > 0) * L_sub)
    new = sum(abs(((t(x) %*% A_sub) > 0) - ((t(y_pert) %*% A_sub) > 0)) * L_sub)
    delta = (original - new) / sum(L)
}
