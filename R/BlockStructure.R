#' @title Ecological Block Struture plot with ggplot2 package
#'
#' @description Make a Structure Plot split into ecological blocks and ordered
#' in each block by an ordering metadata
#'
#' @param omega Cluster membership probabilities of each sample. Usually a
#' sample by cluster matrix in the Topic model output. The cluster weights
#' sum to 1 for each sample.
#' @param blocker_metadata a factor metadata used for creating blocks of Structure plots.
#' @param order_metadata A quantitative metadata used for ordering sites in each block.
#' @param palette A vector of colors assigned to the clusters. First color in
#' the vector is assigned to the cluster labeled as 1, and second color in the
#' vector is assigned to the cluster labeled as 2, etc. The number of colors
#' must be the same or greater than the number of clusters. The clusters not
#' assigned a color are filled with white in the figure. In addition, the
#' recommended choice of color palette here is RColorBrewer, for instance
#' RColorBrewer::brewer.pal(8, "Accent") or RColorBrewwer::brewer.pal(9, "Set1").
#' @param yaxis_label Axis label for the samples.
#' @param split_line Control parameters for line splitting the batches in the
#' plot.
#' @param axis_tick Control parameters for x-axis and y-axis tick sizes.
#' @param plot_labels A logical parameter, if TRUE the function plots the axis labels.
#' @param levels_decreasing if TRUE, the sites in the block are ordered in decreasing order
#'                        of magnitude of ordering metadata. Default is TRUE.
#' @param round_off The ordering metadata is rounded off to these many digits.
#' @param layout The graph layout for plotting the block Structure output.
#' @param plot_order The order of blocks within the specified layout for plotting the block Structure output.
#' @param panel_title_size The size of the title for layout panel
#' @param panel_title_font The font of the title for the layout panel
#' @param main_title the title of the Block Structure model plot.
#'
#'
#' @return Plots the Block Structure plot visualization of the GoM model
#'
#' @examples
#'
#' library(HimalayanBirdsAbundance)
#' data("HimalayanBirdsAbundance")
#' library(Biobase)
#' new_counts <- t(exprs(HimalayanBirdsAbundance));
#' metadata <- pData(HimalayanBirdsAbundance);
#' elevation_metadata=metadata$Elevation;
#' east_west_dir = metadata$WorE;
#' topic_clus <- maptpx::topics(new_counts, K=2, tol=0.1)
#' omega <- topic_clus$omega
#' BlockStructure(omega, blocker_metadata = east_west_dir,
#'               order_metadata = elevation_metadata,
#'               yaxis_label = "Elevation",
#'               levels_decreasing = FALSE)
#'
#' @import CountClust
#' @import grid
#' @import gridExtra
#' @import maptpx
#' @import slam
#'
#' @export



BlockStructure = function( omega,
                           blocker_metadata,
                           order_metadata,
                           palette = c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7"),
                           yaxis_label = "Y",
                           split_line=list(split_lwd = 1,
                                           split_col = "white"),
                           axis_tick = list(axis_ticks_length = .1,
                                            axis_ticks_lwd_y = .1,
                                            axis_ticks_lwd_x = .1,
                                            axis_label_size = 7,
                                            axis_label_face = "bold"),
                           plot_labels = TRUE,
                           levels_decreasing=TRUE,
                           round_off=1,
                           layout=c(1,2),
                           plot_order=c(1,2),
                           panel_title_size=10,
                           panel_title_font=4,
                           main_title="Block Structure Plot"){

  split_indices <- split(1:dim(omega)[1], as.factor(blocker_metadata))
  split_struct <- list()
  for(l in 1:length(split_indices)){

    order_split <- round(order_metadata[split_indices[[l]]], round_off);
    omega_split <- omega[split_indices[[l]],]
    if(levels_decreasing){
      order_split_ordered <- order_split[order(order_split, decreasing=TRUE)]
      omega_split_ordered <- omega_split[order(order_split, decreasing=TRUE),]
    }else{
      order_split_ordered <- order_split[order(order_split, decreasing=FALSE)]
      omega_split_ordered <- omega_split[order(order_split, decreasing=FALSE),]
    }
    annotation <- data.frame(
      sample_id = paste0("X", c(1:NROW(omega_split_ordered))),
      tissue_label = factor(order_split_ordered,
                            levels = unique(order_split_ordered) ) );

    split_struct[[l]] <- CountClust::StructureGGplot(omega = omega_split_ordered,
                                                     annotation = annotation,
                                                     figure_title = names(split_indices)[l],
                                                     palette = palette,
                                                     yaxis_label = yaxis_label,
                                                     split_line=split_line,
                                                     order_sample = FALSE,
                                                     axis_tick = axis_tick,
                                                     plot_labels=TRUE)

  }

  if((layout[1]*layout[2]) != length(split_struct)){
    stop("The layout size does not match with number of blocks")
  }
  if(length(plot_order) != length(split_indices)){
    stop("The plot order length does not match with number of blocks")
  }

  do.call("grid.arrange",
          args = list(grobs=split_struct[plot_order],
                      ncol = layout[2],
                      nrow = layout[1],
                      top=textGrob(main_title,
                                   gp=gpar(fontsize=panel_title_size,
                                           font=panel_title_font))))

}
