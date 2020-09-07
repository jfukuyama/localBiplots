## library(Matrix)
## library(mdessence)
## library(phyloseq)
## library(ape)
## library(gridExtra)
## library(ggtree)
## library(viridis)
## library(rmutil)
#' Make descendant matrix
#'
#' Make a matrix describing the ancestry structure of a tree. Element
#' (i,j) indicates whether leaf i is a descendant of node j.
#'
#' @param tree A tree object of class phylo
#' @return A matrix describing the ancestry structure of a tree. 
#'
#' @importFrom Matrix Matrix
#' @keywords internal
makeDescendantMatrix <- function(tree) {
    etc = edgesToChildren(tree$edge)
    n = length(tree$tip.label)
    m = tree$Nnode
    A = Matrix::Matrix(data = 0, nrow = n, ncol = m + n)
    fillA = function(node, n) {
        if(node <= n) {
            A[node,node] <<- 1
        }
        else {
            children = etc[[node]]
            for(c in children) {
                fillA(c, n)
            }
            A[,node] <<- Matrix::rowSums(A[,children,drop=FALSE])
        }
    }
    fillA(n+1, n)
    return(A)
}

#' Makes a hash table with nodes and their children
#'
#' Takes the edge matrix from a phylo-class object and turns it into a
#' list where the entries are nodes and the elements are vectors with
#' the children of the nodes.
#' @keywords internal
edgesToChildren <- function(edges) {
    ## ape requires that the tips are numbered 1:n, nodes are numbered
    ## n+1:n+m (n tips, m internal nodes), root is numbered n+1. The
    ## first column of the edge matrix contains only the internal nodes
    etc = vector("list", max(edges))
    for(i in 1:nrow(edges)) {
        parent = edges[i,1]
        child = edges[i,2]
        if(is.null(etc[[parent]]))
            etc[[parent]] = child
        else
            etc[[parent]] = c(etc[[parent]], child)
    }
    return(etc)
}

make_example_data = function(n_per_group, log_2_p, deep_split_scaling = 1,
                             overall_scaling = 1, s = 1,
                             sample_partition_type = "confounded",
                             species_partition_type = "alternating_tips") {
    tree = make_tree(log_2_p)
    data_matrices = make_data_matrix(n_per_group = n_per_group, p = 2^log_2_p,
                                     deep_split_scaling = deep_split_scaling,
                                     overall_scaling = overall_scaling,
                                     s = s, 
                                     sample_partition_type = sample_partition_type,
                                     species_partition_type = species_partition_type)
    variable_descriptions = make_variable_descriptions(log_2_p, varnames = colnames(data_matrices$X), species_partition_type = species_partition_type)
    return(list(tree = tree,
                X = data_matrices$X,
                C = data_matrices$C,
                mean_matrix = data_matrices$mean_matrix,
                species_partitions = data_matrices$species_partitions,
                sample_partitions = data_matrices$sample_partitions,
                variable_descriptions = variable_descriptions))
}

make_variable_descriptions = function(log_2_p, varnames, species_partition_type) {
    p = 2^log_2_p
    partitions = make_partitions(p, species_partition_type = species_partition_type)
    df = data.frame(variable = varnames, shallow = partitions$shallow, deep = partitions$deep,
                    var_description = paste0(partitions$deep, partitions$shallow))
    return(df)
}

make_tree = function(log_2_p) {
    tree = stree(n = 2^log_2_p, type = "balanced")
    tree$edge.length = rep(1, nrow(tree$edge))
    return(tree)
}


make_partitions = function(p, species_partition_type = "alternating_tips") {
    if(species_partition_type == "alternating_tips") {
        uf_partition = seq(1, p, by = 2)
        uf_partition_vec = 1:p %in% uf_partition
    } else if(species_partition_type == "one_clade") {
        uf_partition_vec = c(rep(1, p/8), rep(0, p/8), rep(.5, p * .75))
    }
    wuf_partition = 1:(p/2)
    wuf_partition_vec = 1:p %in% wuf_partition
    return(list(shallow = uf_partition_vec, deep = wuf_partition_vec))
}

make_sample_partitions = function(n_per_group, type) {
    a = rep(0:1, each = n_per_group)
    if(type == "confounded") {
        b = a
    } else if(type == "orthogonal") {
        b = rep(0:1, n_per_group)
    }
    return(list(a = a, b = b))
}


#' Make simulated data
#'
#' @param s The dispersion factor for the double poisson. Smaller
#'     values of s correspond to larger variance.
make_data_matrix = function(n_per_group, p, deep_split_scaling, overall_scaling, s,
                            sample_partition_type,
                            species_partition_type) {
    species_partitions = make_partitions(p, species_partition_type = species_partition_type)
    sample_partitions = make_sample_partitions(n_per_group, type = sample_partition_type)
    m1 = overall_scaling * (outer(sample_partitions$a, species_partitions$shallow) +
                           outer(1 - sample_partitions$a, 1 - species_partitions$shallow))
    m2 = deep_split_scaling * (outer(sample_partitions$b, species_partitions$deep) +
                            outer(1 - sample_partitions$b, 1 - species_partitions$deep))
    mean_matrix = m1 * exp(m2)
    C = matrix(rdoublepois(n = length(mean_matrix), m = mean_matrix, s = s), nrow = nrow(mean_matrix))
    colnames(C) = colnames(mean_matrix) = paste0("t", 1:p)
    X = sweep(C, 1, STATS = rowSums(C), FUN = "/")
    return(list(X = X, C = C, mean_matrix = mean_matrix, species_partitions = species_partitions, sample_partitions = sample_partitions))
}

getBranchLengths = function(tree) {
    branchLengths = rep(0, length(tree$edge.length))
    branchLengths[tree$edge[,2]] = tree$edge.length
    return(branchLengths)
}

make_phyloseq_from_example_data = function(example_data, mean_matrix = FALSE) {
    colnames(example_data$mean_matrix) = example_data$tree$tip.label
    tax_data = data.frame(lapply(example_data$species_partitions, as.character),
                          row.names = example_data$tree$tip.label)
    tax_data$var_description = paste0(example_data$species_partitions$deep,
                                      example_data$species_partitions$shallow)
    if(mean_matrix) {
        ps = phyloseq(phy_tree(example_data$tree),
                      otu_table = otu_table(example_data$mean_matrix, taxa_are_rows = FALSE),
                      tax_table = tax_table(as.matrix(tax_data)))

    } else {
        ps = phyloseq(phy_tree(example_data$tree),
                      otu_table = otu_table(example_data$X, taxa_are_rows = FALSE),
                      tax_table = tax_table(as.matrix(tax_data)))
    }
    return(ps)
}
