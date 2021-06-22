#' Tree Clustering Plot
#'
#'\code{plotTreeClust} plots a tree clustering from a cluster model
#'
#'@description Based on a function written by dmontaner at cipf.es
#'@param cluster cluster model
#'@param title plot title
#'@param palette color palette
#'@importFrom ggdendro dendro_data segment
#'@import ggplot2
#'@export

plotTreeClust <- function(cluster, title,
                          palette=c("#808080", "#FFA500", "#009E73", "#56B4E9", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")) {
  ##
  dendr <- ggdendro::dendro_data(cluster, type = "rectangle")
  ##
  clases <- as.character(cluster$clase)
  clust.df <- data.frame(label = cluster$labels, Clusters = factor(clases))
  ##
  dendr[["labels"]] <- merge(dendr[["labels"]], clust.df, by = "label")
  ##
  ggplot() +
    geom_segment(data = segment(dendr), aes(x=x, y=y, xend=xend, yend=yend)) +
    geom_text(data = label(dendr), aes(x, y, label = label, hjust = 0, color = Clusters), size = 3) +
    coord_flip() + scale_y_reverse(expand = c(0.2, 0)) +
    ggtitle(title) +
    theme(axis.line.y=element_blank(),
          axis.ticks.y=element_blank(),
          axis.text.y=element_blank(),
          axis.title.y=element_blank(),
          axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          panel.background=element_rect(fill="white"),
          panel.grid=element_blank(),
          plot.title = element_text(hjust = 0.5)) +
    scale_colour_manual(values=palette) +
    scale_fill_manual(values=palette)
}
