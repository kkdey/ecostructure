#' @title Ecological Block Struture plot with ggplot2 package
#'
#' @description Make a Structure Plot split into ecological blocks by a blocker
#' metadata and ordered in each block by an ordering metadata
#'
#' @param omega Matrix of Cluster membership probabilities of each sample
#'              obtained from \code{ecos_fit()}. The row sums 
#'              sum to 1 for each sample.
#' @param blocker_metadata a factor metadata used for creating blocks of 
#'                         Structure plots.
#' @param order_metadata A quantitative metadata used for ordering sites 
#'                       in each block.
#' @param palette A vector of colors assigned to the clusters. First color in
#'                the vector is assigned to the cluster labeled as 1, and second
#'                color in the vector is assigned to the cluster labeled as 2, 
#'                etc. The number of colors must be the same or greater than the
#'                number of clusters. The recommended choice of color palette 
#'                here is RColorBrewer, for instance 
#'                RColorBrewer::brewer.pal(8, "Accent") or 
#'                RColorBrewwer::brewer.pal(9, "Set1").
#' @param structure_control Control parameters for the Block Structure plot.
#'                          Fixes the title, axis labels, tick sizes, splitting
#'                          line characteristics and panel orientation.
#' @param layout The graph layout for plotting the Block Structure output.
#'               Is missing, automatically determined as a square 
#'               configuration.
#' 
#' @return Plots the Block Structure plot visualization of the GoM model
#'
#' @examples
#'
#' data("himalayan_birds")
#' species_abundance_counts <- t(exprs(himalayan_birds));
#' site_metadata <- pData(himalayan_birds);
#' elevation_metadata=site_metadata$Elevation;
#' east_west_dir = site_metadata$WorE;
#' topic_clus <- ecos_fit(species_abundance_counts, K = 2, tol = 0.1)
#' ecos_blocks(topic_clus$omega,
#'             blocker_metadata = east_west_dir,
#'             order_metadata = elevation_metadata)
#'
#' @importFrom CountClust StructureGGplot
#' @import grid
#' @importFrom gridExtra grid.arrange
#'
#' @export



ecos_blocks = function(omega,
                       blocker_metadata,
                       order_metadata,
                       palette = c("#E69F00", "#56B4E9", "#009E73", "#F0E442", 
                                       "#0072B2", "#D55E00", "#CC79A7"),
                       structure_control = list(),
                       layout){
  
  if(!is.factor(blocker_metadata)){
    stop("the blocker_metadata must be a factor variable")
  }
  
  if(!is.numeric(order_metadata)){
    stop("the order_metadata must be a numeric variable")
  }
  
  num_levels <- length(levels(blocker_metadata))
  
  if(missing(layout)){
    layout <- c(floor(sqrt(num_levels)), ceiling(sqrt(num_levels)))
  }
  
  structure_control_default <- list(split_line=list(split_lwd = 1,
                                                    split_col = "white"),
                                    axis_tick = list(axis_ticks_length = .1,
                                                     axis_ticks_lwd_y = .1,
                                                     axis_ticks_lwd_x = .1,
                                                     axis_label_size = 7,
                                                     axis_label_face = "bold"),
                                    plot_labels = TRUE,
                                    levels_decreasing=FALSE,
                                    order_sample=TRUE,
                                    round_off=1,
                                    panel_title_size=10,
                                    panel_title_font=4,
                                    main_title="Block Structure Plot",
                                    yaxis_label = "Y")
  
  structure_control <- modifyList(structure_control_default, structure_control)
  
  split_indices <- split(1:dim(omega)[1], as.factor(blocker_metadata))
  split_struct <- list()
  
  for(l in 1:length(split_indices)){

    order_split <- round(order_metadata[split_indices[[l]]], structure_control$round_off);
    omega_split <- omega[split_indices[[l]],]
    if(structure_control$levels_decreasing){
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
                          yaxis_label = structure_control$yaxis_label,
                          split_line=structure_control$split_line,
                          order_sample = structure_control$order_sample,
                          axis_tick = structure_control$axis_tick,
                          plot_labels=structure_control$plot_labels)

  }

  do.call(gridExtra::grid.arrange,
          args = list(grobs=split_struct,
                      ncol = layout[2],
                      nrow = layout[1],
                      top=grid::textGrob(structure_control$main_title,
                      gp=gpar(fontsize=structure_control$panel_title_size,
                      font=structure_control$panel_title_font))))

}
