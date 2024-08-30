#' Creates a PCA from a counts matrix
#'
#' @param counts Numeric matrix with the count values to plot.
#' @param targets Dataframe including the metadata annotation of each sample and the variables of interest. Used to access the condition information to color.
#' @param PCs Principal components to plot. By default plots PC1 and PC2, but non-default PCs can be plotted if specified; e.g. c(2,3) will plot PC2 and PC3. Default = NULL
#' @param label Whether to show the labels of the points of the PCA. Default = FALSE
#' @param title Title for each plot. Default = NULL
#' @param color_condition Name of the variable in `targets` that includes the conditions to be used for coloring the PCA. Default = NULL
#' @param shape_condition Name of the variable in `targets` that includes the conditions to be used for the shape of the PCA. Passing shape = FALSE makes plot without points. In this case, label is turned on unless otherwise specified. Default = NULL
#' @param colors Vector with the colors assigned to each condition (in order of the unique(conditions)). Default = NULL
#' @param fmtPlot Format for the image file, "png" or "pdf" (static plot); or "html" (interactive). If none specified image will pop up in R session. Default = ""
#' @param resultsDir Output directory. Default = current working directory (getwd())
#' @param fileName Name of the output file, without extension. Default = "PCAPlot"
#' @param interactive Creates a PCA plot that shows the label of each point as you hover over it. Must be TRUE to be saved as html. Can also show plot in session. Default = FALSE
#'
#' @import ggplot2
#' @import ggfortify
#' @import RColorBrewer
#'
#' @return PCA image
#' @export makePCA

makePCA <-function(counts, targets, PCs=NULL, label=FALSE, title=NULL, color_condition=NULL, shape_condition=NULL, colors=NULL, resultsDir = getwd(), fileName = "PCAPlot", fmtPlot = "", interactive=FALSE, ...) {

  # principal components to plot
  if (is.null(PCs)) {
    x = 1
    y = 2
  } else {
    x = PCs[1]
    y = PCs[2]
  }

  if(is.null(colors)) {
    colors <- brewer.pal(8, "Dark2")
  }

  pca.filt <- prcomp(t(counts), scale=T)

  p <- autoplot(pca.filt, x=x, y=y, data=targets, label=label, colour=color_condition, shape=shape_condition, main=title, ...) + theme_classic() + scale_color_manual(values=colors)

  if (interactive) {
    pca_data <- as.data.frame(pca.filt$x)
    pca_data$label <- rownames(pca_data)  # Assuming rownames are your labels
    pca_data$color_cond <- targets[,color_condition]
    pca_data$shape_cond <- targets[,shape_condition]

    explained_variance <- pca.filt$sdev^2 / sum(pca.filt$sdev^2) * 100
    pc_labels <- paste0("PC", 1:length(explained_variance), " (", round(explained_variance, 2), "%)")
    # Create a ggplot object
    p <- ggplot(pca_data, aes(x = pca_data[,x], y = pca_data[,y], colour = color_cond, text = label)) +
      geom_point(size=1) +
      theme_classic() +
      scale_color_manual(values = colors) +
      labs(colour="Condition") +
      xlab(pc_labels[1]) +
      ylab(pc_labels[2])

    if(!is.null(shape_condition)) {
      p <- ggplot(pca_data, aes(x = pca_data[,x], y = pca_data[,y], colour = color_cond, shape = shape_cond, text = label)) +
        geom_point(size=1) +
        theme_classic() +
        scale_color_manual(values = colors) +
        labs(colour="Condition") +
        xlab(pc_labels[1]) +
        ylab(pc_labels[2])
    }

    # Convert ggplot to an interactive plotly object
    p_interactive <- ggplotly(p, tooltip = "text")
  }

  if(fmtPlot == "") {
    if (interactive) {
      print(p_interactive)
      return(p_interactive)
    } else {
      print(p)
      return(p)
    }

  } else if (fmtPlot %in% c("png","pdf")) {
    ggsave(file.path(resultsDir,paste0(fileName,".",fmtPlot)), plot=p, width=8, height=6)
  } else if (fmtPlot == "html") {
    if (interactive==FALSE) {
      stop("To save the plot in html, you must set parameter interactive=TRUE")
    } else {
      htmlwidgets::saveWidget(as_widget(p_interactive), file.path(resultsDir, paste(fileName, fmtPlot, sep = ".")))

    }

  }


}





#' Creates a PCA 3D from a numerical matrix
#'
#' @param estimates Numeric matrix with the intensity values to plot.
#' @param resultsDir Output directory. Default = NULL
#' @param fileName Name of the output file, without extension. Default = NULL
#' @param fmtPlot Format for the image file, "png" or "pdf" (default). If none specified image will pop up in R session
#' @param title Title for each plot. Default = NULL
#' @param conditions Vector with the different conditions. Default = NULL
#' @param colors Vector with the colors assigned to each condition (in order of the unique(conditions)). Default = NULL
#' @param dim Dimensions of the PCA plot. Possible values are 2 or 3. Default = 3
#' @param dist Distance from labels to dots. Default = 2
#' @param legend.pos position of the legend. Default "bottomright"
#'
#' @import scatterplot3d
#' @return PCA image
#' @export makePCA3D

makePCA3D <-function(estimates, resultsDir = NULL, fileName = NULL, fmtPlot = "pdf",
                   title = NULL, conditions = NULL, colors = NULL, dist = 2, dim = 3,
                   legend.pos = "bottomright")
{

    labels <- colnames(estimates)
    parameters <- setParameters(labels)
    summary(pca.filt <- prcomp(t(estimates), center = TRUE, scale. = TRUE ))
    vars <- apply(pca.filt$x, 2, var)
    props <- vars / sum(vars)
    PCAvec <- round(props*100, 2)

    if (fmtPlot == "pdf") {

      pdf(file = file.path(resultsDir, paste(fileName,"pdf", sep=".")))

    } else if (fmtPlot == "png") {

      png(file = file.path(resultsDir, paste(fileName, "png", sep = ".")),
          width = parameters$wid, height = parameters$hei, res = parameters$res)

    }

    if (is.null(conditions)){

      colors = "black"

    } else if (is.null(colors)) {

      palette = c(brewer.pal(8, "Dark2"), brewer.pal(12, "Paired"))
      list1 <- unique(as.character(sort(conditions)))
      ColVect <- palette #20 colors by default
      list2 <- ColVect[1:length(unique(conditions))]
      map = setNames(list2, list1)
      colors <- map[conditions]

    } else {

      list1 <- unique(as.character(conditions))
      list2 <- unique(colors)
      map = setNames(list2, list1)
      colors <- map[conditions]

    }

    if (dim==3)  {

      pca3d <- scatterplot3d(x = pca.filt$x[, 1], y = pca.filt$x[, 2], z = pca.filt$x[, 3],
                             xlab = sprintf("PC1 %.0f%%", PCAvec[1]),
                             ylab = sprintf("PC2 %.0f%%", PCAvec[2]),
                             zlab = sprintf("PC3 %.0f%%", PCAvec[3]),
                             main = 'PCA', col.grid = "lightblue",
                             cex.symbols = parameters$ce+0.2, color = colors, pch = 16)
      text(pca3d$xyz.convert(pca.filt$x+dist), labels = rownames(pca.filt$x),
           cex = parameters$ce+0.2, col = colors)
      if (!is.null(conditions)) legend(legend.pos, legend = list1, col = list2, pch = 16, ncol = 1, cex = 0.9)

    } else {

      pca2d <- plot(x = pca.filt$x[,1],y = pca.filt$x[,2],
                    xlab = sprintf("PC1 %.0f%%", PCAvec[1]),
                    ylab = sprintf("PC2 %.0f%%", PCAvec[2]),
                    main = 'PCA', col = colors, pch=20)

      text(pca.filt$x[, 1]+dist, pca.filt$x[, 2]+dist, labels = rownames(pca.filt$x),
           cex = parameters$ce+0.2, col = colors)
      if (!is.null(conditions)) legend(legend.pos, legend = list1, col = list2, pch = 16, ncol = 1, cex = 0.9)

    }

   if (fmtPlot %in% c("pdf", "png")) dev.off()
}
