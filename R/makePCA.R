#' Creates a PCA 3D from a numerical matrix
#' 
#' @param estimates: Numeric matrix with the intensity values to plot
#' @param resultsDir: Output directory
#' @param fileName: name of the output file, without extension
#' @param fmtPlot: Format for the image file, "png" or "pdf" (default). If none specified image will pop up in R session
#' @param title: Title for each plot
#' @param conditions: Vector with the different conditions
#' @param colors: Vector with the colors assigned to each condition (in order of the unique(conditions)). Default = NULL
#' @param dim: Dimensions of the PCA plot. Default = 3. Possible values are 2 or 3
#' @param dist: Distance from labels to dots. Default = 2
#' 
#' @param legend.pos position of the legend. Default "bottomright"
#' 
#' @import scatterplot3d
#' @return PCA image
#' @export makePCA

makePCA <-function(estimates, resultsDir = NULL, fileName = NULL, fmtPlot = "pdf",
                   title, conditions = NULL, colors = NULL, dist = 2, dim = 3,
                   legend.pos = "bottomright")
{
    
    labels <- colnames(estimates)
    parameters <- setParameters(labels)
    summary(pca.filt <- prcomp(t(estimates), center = TRUE, scale. = TRUE )) 
    vars <- apply(pca.filt$x, 2, var)  
    props <- vars / sum(vars)
    PCAvec <- round(props*100, 2)
    
    if (fmtPlot == "pdf") {
      
      pdf(file=file.path(resultsDir, paste(fileName,"pdf", sep=".")))
      
    } else if (fmtPlot == "png") {
      
      png(file.path(resultsDir, paste(title, "PCA.png", sep = "_")),
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
