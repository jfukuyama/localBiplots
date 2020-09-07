#' Interactive biplot
#'
#' @param X A data matrix, samples as rows.
#' @param dist Either a string describing one of the supported
#'     distances or a function that takes a matrix and returns the
#'     distances between the rows of the matrix.
#' @param dist_deriv Either NULL (if dist is a string describing one
#'     of the supported distances) or a function that takes two
#'     vectors and computes \eqn{\frac{\partial}{\partial y_j}d(x,y)}.
#' @param k The number of embedding dimensions.
#' @param sample_data A data frame, containing extra information about
#'     the samples to be appended to the embeddings for plotting.
#' @param sample_mapping The output from aes_string. Defaut is
#'     aes_string(x = "Axis1", y = "Axis2"); add in additional
#'     aesthetics to map variables in sample_data.  in sample_data.
#' @param sample_facet Output from facet_grid or facet_wrap. Can facet
#'     on variables in sample_data.
#' @param var_data A data frame, containing extra information about
#'     the variables to be appended to the local biplot axes for
#'     plotting.
#' @param var_mapping The output from aes_string. Default is
#'     aes_string(x = "Axis1", y = "Axis2"); add in additional
#'     aesthetics to map variabes in var_data.
#' @param layout The widths of the sample plot window and the variable
#'     plot window. Widths are out of 12.
#' @param width The width of the plots, in pixels.
#' @param res The resolution for the plots.
#'
#' @import shiny
#' @import ggplot2
#' @export
interactive_biplot <- function(X, dist, dist_deriv = NULL, k = 2,
                               sample_data = NULL,
                               sample_mapping = aes_string(x = "Axis1", y = "Axis2"),
                               sample_facet = NULL,
                               var_data = NULL,
                               var_mapping = aes_string(x = "Axis1", y = "Axis2"),
                               layout = c(6,6),
                               width = 600,
                               res = 90) {
    samples = NULL
    ui <- fluidPage(
        headerPanel("Interactive Local Biplot"),
        mainPanel(width = 11,
                  fluidPage(
                      fluidRow(
                          column(layout[1], plotOutput("plot_samples", click = "plot_click")),
                          column(layout[2], plotOutput("plot_variables"))),
                      fluidRow(column(1), actionButton("done", "Done"))
                  )
                  )
    )

    server <- function(input, output, session) {
        dist_fns = make_dist_fns(dist, dist_deriv)
        mds_matrices = make_mds_matrices(X, dist_fns$dist_fn)
        embedding_and_sample_name = data.frame(mds_matrices$Y, sample = paste0("Original", 1:nrow(mds_matrices$Y)))
        lb_df = compute_lb_samples(mds_matrices, dist_fns, k = k, samples = samples)
        varnames = get_variable_names(mds_matrices$X)
        output$plot_samples = renderPlot({
            if(!is.null(sample_data))
                p = ggplot(data.frame(embedding_and_sample_name, sample_data), sample_mapping) +
                    geom_point()
                else
                    p = ggplot(embedding_and_sample_name, sample_mapping) +
                        geom_point()
                if(!is.null(sample_facet))
                    p = p + sample_facet
                p + coord_fixed() + get_square_limits(p)
        }, width = 600, res = 90)
        output$plot_variables = renderPlot({
            biplot_center = nearPoints(embedding_and_sample_name, input$plot_click, maxpoints = 1)
            sample_name = biplot_center$sample
            if(length(sample_name) > 0) {
                if(!(sample_name %in% lb_df$sample)) {
                    sample_num = as.numeric(gsub("Original", "", sample_name))
                    lb_df = rbind(lb_df,
                                  compute_lb_samples(mds_matrices, dist_fns, k, samples = sample_num))
                }
                if(!is.null(var_data)) {
                    p = ggplot(data.frame(subset(lb_df, sample == sample_name), var_data), var_mapping) +
                        geom_point()
                } else {
                    p = ggplot(subset(lb_df, sample == sample_name), var_mapping) +
                        geom_point()
                }
                p + coord_fixed() + get_square_limits(p)
            }
        }, width = 600, res = 90)

        observeEvent(input$done, {
            stopApp(lb_df)
        })
    }
    
    runGadget(ui, server)
}

get_square_limits <- function(p) {
    build = ggplot_build(p)
    y_limits = build$layout$panel_scales_y[[1]]$range$range
    x_limits = build$layout$panel_scales_x[[1]]$range$range
    y_range = y_limits[2] - y_limits[1]
    x_range = x_limits[2] - x_limits[1]
    if(y_range > x_range) {
        x_limits = c(mean(x_limits) - y_range / 2, mean(x_limits) + y_range / 2)
        return(xlim(x_limits))
    } else {
        y_limits = c(mean(y_limits) - x_range / 2, mean(y_limits) + x_range / 2)
        return(ylim(x_limits))
    }
}

get_variable_names <- function(m) {
    if(is.null(colnames(m))) {
        return(paste("Var", 1:ncol(m), sep = ""))
    }
    return(colnames(m))
}
