#' @title A bi-Y plot of topic proportions and a metadata across sites, ordered by another metadata
#'
#' @description a bi-Y axis plot of topic proportions (or any site metadata) across sites arranged
#' preserving relative ordering of some metadata of interest (say elevation) along one Y axis
#' and another metadata plotted along the other Y axis. The goal is to figure out
#' drop/jump in topic proportions across sites is related to the two metadata.
#'
#'  @param annotation A data frame that contains four objects
#'          \itemize{
#'               \item{\code{x_names}} {is a N vector of site names}
#'               \item{\code{x}} {a N vector of metadata that we use to order sites along x axis}
#'               \item{\code{y}} {a N vector of topic grades of membership or metadata plotted along one of the Y axis}
#'               \item{\code{y2}}{ a N vector of metadata we want to compare with grades of membership plotted along second Y axis}
#'           }
#'
#'  @param  margin  The values of the margin on the bottom, left, top and right.
#'  @param  col.pts1 The color of the points for \code{annotation$y}.
#'  @param  pch.pts1 The size of the points for \code{annotation$y}. Defaults to 20.
#'  @param  lwd.pts1 The lwd of the points plot for \code{annotation$y}. Defaults to 3.
#'  @param  lty.pts1 The lty of the points plot for \code{annotation$y}. Defaults to 1
#'  @param  cex.pts1 The lty of the points plot for \code{annotation$y}. Defaults to 1
#'  @param  cex.axis The size of the X axis labels
#'  @param  col.pts2 The color of the segments generated from distance metric
#'  @param  pch.pts2 The size of the points for \code{annotation$y2}. Defaults to 18.
#'  @param  lwd.pts2 The lwd of the points plot for \code{annotation$y2}. Defaults to 3
#'  @param  lty.pts2 The lty of the points plot for \code{annotation$y2}. Defaults to 1
#'  @param  cex.pts2 The lty of the points plot for \code{annotation$y2}. Defaults to 1
#'  @param  las Orientation of the X axis labels. Defaults to 2.
#'  @param  text.x The X-axis label
#'  @param  text.y1 The label on the first Y axis
#'  @param  text.y2 The label on the second Y axis.
#'  @param  ylim   The range of values for the first Y axis. Defaults to NULL in which case
#'                 the range is determined from data
#'  @param  y2lim The range of values for the second Y axis. Defaults to NULL in which case
#'                 the range is determined from data
#'  @param  line.x The gap between the X-axis label and the axis.
#'  @param  line.y The gap between the first Y-axis label and the first Y-axis
#'  @param  line.y2 The gap between the second Y-axis label and the second Y axis
#'  @param  round_off the rounding factor used to report the X-axis metadata with names
#'  @param  legend.pos the position of the legend
#'
#'  @return Produces a bi-Y plot of \code{annotation$y} and  \code{annotation$y2}
#'          against \code{annotation$x} in the X-axis.
#'
#'  @examples
#'  annotation = data.frame(x_names = c(paste0("A",1:5), paste0("B",1:5)),
#'  x = c(0.5,2.0, 3.2, 4.6, 6.3,  23.5, 26.4, 28.5, 29.6, 31.8),
#'  y1 = c(0.4, 0.3, 0.4, 0.35, 0.4, 0.8, 0.85, 0.9, 0.8, 0.75),
#'  y2 =c(5, 6.6, 4, 5.2, 20, 3.4, 5.6, 4.5, 8, 10))
#'   TopicMetaMeta(annotation)
#'
#'  @export

#############  Plot with 2 Y -axes in R (ecology + evolution) ###################


TopicMetaMeta = function(     annotation,
                              margin=c(5,5,2,5),
                              col.pts1 = "red3",
                              pch.pts1=20,
                              lwd.pts1=3,
                              lty.pts1=1,
                              cex.pts1=1,
                              cex.axis=0.5,
                              col.pts2 = "blue",
                              pch.pts2 = 18,
                              lwd.pts2 = 2,
                              lty.pts2 = 3,
                              cex.pts2 = 1.5,
                              las=2,
                              text.x = "X",
                              text.y1 = "Y",
                              text.y2 = "Y2",
                              ylim = NULL,
                              y2lim = NULL,
                              line.x=3,
                              line.y=3,
                              line.y2=3,
                              round_off=0,
                              legend.pos="topleft" ){

  if(is.null(ylim)){
    par(mar = margin)
    with(annotation, plot(x = x, y = y1, col=col.pts,
                          ylab=NA, xaxt="n", xlab=NA,
                          pch=pch.pts, lwd=lwd.pts, lty=lty.pts))
  }
  if(!is.null(ylim)){
    par(mar = margin)
    with(annotation, plot(x = x, y = y1, col=col.pts, ylim=ylim,
                          ylab=NA, xaxt="n", xlab=NA,
                          pch=pch.pts, lwd=lwd.pts, lty=lty.pts))
  }

  axis(1, at=annotation$x, paste0(annotation$x_names, "(",
          round(annotation$x,round_off), ")"),
          cex.axis=cex.axis, las=las)
  mtext(side=2, line=line.y, text.y1)
  mtext(side=1, line=line.x, text.x)

  if(is.null(y2lim)){
    y2lim = c(0, (max(annotation$y2)+0.5));
  }
  par(new = T)
  with(annotation, plot(x = x, y = rep(0, length(x)),
                        ylim=y2lim, type="n", axes=F, xlab="", ylab=""))
  with(annotation, points(x = x, y = y2, col=col.pts2,
                          ylab=NA, xaxt="n", xlab=NA,
                          pch=pch.pts2, lwd=lwd.pts2, lty=lty.pts2,
                          cex=cex.pts2))

  axis(side=4, ylim=y2lim)
  mtext(side = 4, line = line.y2, text.y2)

  legend(legend.pos,
         legend=c(text.y2, text.y1),
         pch=c(pch.pts2, pch.pts1),
         col=c(col.pts2, col.pts1))
}


#################   Example  #############################


