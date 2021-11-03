#' Creates a pdf with different dendrograms (hierarchical clusters) combining a distance and a linkage method
#'
#' @param estimates: Numeric matrix with the intensity values to plot
#' @param resultsDir: Output directory. Default = NULL
#' @param fileName: name of the output file, without extension. Default = NULL
#' @param fmtPlot: Format for the image file, pdf" (default). If none specified images will be generated in R session
#' @param title: Title for each plot. Default = NULL
#' @param method: linkage method to cluster samples. Values given by parameter method of hclust function. Default = "ward.D2"
#' @param conditions: Vector with the different conditions. Default = NULL
#' @param colors: Vector with the colors assigned to each condition (in order of the unique(conditions)). Default = NULL
#'
#' @import RColorBrewer
#' @return a pdf and/or png with the hierarchical cluster
#' @export oneCluster

oneCluster <-function(estimates, resultsDir = NULL, fileName = NULL, fmtPlot = "pdf",
                      title = NULL, method = "ward.D2", conditions = NULL, colors = NULL)
{

    labels <- colnames(estimates)
    parameters <- setParameters(labels)

    use.cor = "pairwise.complete.obs"

    clust.cor.ward <- hclust(as.dist(1-cor(estimates,use=use.cor)), method = method)
    xlab <- paste("Correlation", method, sep = "-")

    #plots
    if (fmtPlot == "pdf") {

      pdf(file = file.path(resultsDir, paste(fileName, "pdf", sep = ".")))

    } else if (fmtPlot == "png") {

      png(file=file.path(resultsDir, paste(fileName,"png", sep=".")),
          width=parameters$wid,height=parameters$hei,res=parameters$res)

    }

    if (is.null(conditions)) {

        opt<-par(cex.main=1, cex = parameters$ce, cex.axis = parameters$ce)
          plot(clust.cor.ward, main=title, hang=-1, xlab=xlab, ylab="", sub = "")
        par(opt)

    } else {

      if(is.null(colors)){

        palette=c(brewer.pal(8, "Dark2"), brewer.pal(12, "Paired"))
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

      opt<-par(cex.main=1,cex.axis=0.8, cex=0.8)
      clust.cor.ward <- colorCluster(clust.cor.ward, colors, parameters$ce)
      plot(clust.cor.ward, main=title, xlab=xlab)
      legend("topright",legend=list1, cex=parameters$ce+0.2,fill=list2)

    }

    if (fmtPlot %in% c("pdf", "png")) dev.off()

    return(clust.cor.ward)
}
