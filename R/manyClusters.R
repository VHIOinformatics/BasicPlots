#' Creates a pdf with different dendrograms (hierarchical clusters) combining a distance and a linkage method
#'
#' @param estimates Numeric matrix with the intensity values to plot. The matrix should have samples as columns and genes/transcripts/exons as rows
#' @param resultsDir Output directory. Default = NULL
#' @param fileName name of the output file, without extension. Default = "manyClusters"
#' @param fmtPlot Format for the image file. If none specified images will be generated in R session. Default = "pdf"
#' @param title Title for each plot. Default = NULL
#' @param conditions Vector with the different conditions. Default = NULL
#' @param colors Vector with the colors assigned to each condition. Default = NULL
#'
#' @import RColorBrewer
#' @return a pdf with the hierarchical clusters
#' @export manyClusters

manyClusters <-function(estimates, resultsDir = NULL, fileName = "manyClusters", fmtPlot = "pdf",
                        title = NULL, conditions = NULL, colors = NULL)
{

    labels <- colnames(estimates)
    parameters <- setParameters(labels)

    use.cor="pairwise.complete.obs"

    #all clusters
    clust.cor.ward<- hclust(as.dist(1-cor(estimates,use=use.cor)), method = "ward.D2")
    clust.cor.average <- hclust(as.dist(1-cor(estimates,use=use.cor)), method = "average")
    clust.cor.complete <- hclust(as.dist(1-cor(estimates,use=use.cor)), method = "complete")
    clust.euclid.ward <- hclust(dist(t(estimates)), method = "ward.D2")
    clust.euclid.average <- hclust(dist(t(estimates)), method = "average")
    clust.euclid.complete <- hclust(dist(t(estimates)), method = "complete")

    #plots
    if (fmtPlot == "pdf") pdf(file=file.path(resultsDir, paste(fileName, "pdf", sep = ".")))

    if (is.null(conditions)) {

        opt<-par(cex.main=1, cex = parameters$ce, cex.axis = parameters$ce)
          plot(clust.cor.ward, main = title, hang =-1, xlab = "Correlation-Ward", ylab="", sub = "")
          plot(clust.cor.average, main = title, hang =-1, xlab = "Correlation-Average", ylab="", sub = "")
          plot(clust.cor.complete, main = title, hang =-1, xlab = "Correlation-Complete", ylab="", sub = "")
          plot(clust.euclid.ward, main = title, hang =-1, xlab = "Euclidean-Ward", ylab="", sub = "")
          plot(clust.euclid.average, main = title, hang =-1, xlab = "Euclidean-Average", ylab="", sub = "")
          plot(clust.euclid.complete, main = title, hang =-1, xlab = "Euclidean-Complete", ylab="", sub = "")
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

        list1 <- unique(as.character(sort(conditions)))
        list2 <- unique(colors)
        map = setNames(list2, list1)
        colors <- map[conditions]

      }

      opt<-par(cex.main = 1, cex.axis = 0.8, cex = 0.8)
      clust.cor.ward <- colorCluster(clust.cor.ward, colors, parameters$ce)
      plot(clust.cor.ward, main = title, xlab = "ward.D2")
      legend("topright", legend = list1, cex = parameters$ce+0.2, fill = list2)

      clust.cor.average <- colorCluster(clust.cor.average, colors, parameters$ce)
      plot(clust.cor.average, main = title, xlab = "Average")
      legend("topright", legend = list1, cex = parameters$ce+0.2, fill = list2)

      clust.cor.complete <- colorCluster(clust.cor.complete, colors, parameters$ce)
      plot(clust.cor.complete, main = title, xlab = "Complete")
      legend("topright", legend = list1, cex = parameters$ce+0.2, fill = list2)

      clust.euclid.ward <- colorCluster(clust.euclid.ward, colors, parameters$ce)
      plot(clust.euclid.ward, main = title, xlab = "Euclidean Ward")
      legend("topright",legend = list1, cex = parameters$ce+0.2,fill = list2)

      clust.euclid.average <- colorCluster(clust.euclid.average, colors, parameters$ce)
      plot(clust.euclid.average, main = title, xlab = "Euclidean Average")
      legend("topright",legend = list1, cex = parameters$ce+0.2, fill = list2)

      clust.euclid.complete <- colorCluster(clust.euclid.complete, colors, parameters$ce)
      plot(clust.euclid.complete, main = title, xlab = "Euclidean Complete")
      legend("topright",legend = list1, cex = parameters$ce+0.2, fill = list2)
      par(opt)

    }

    if (fmtPlot == "pdf") dev.off()

    return(list(Corr.ward = clust.cor.ward, Corr.avg = clust.cor.average,
                Corr.compl = clust.cor.complete, Euclid.ward = clust.euclid.ward,
                Euclid.avg = clust.euclid.average, Euclid.compl = clust.euclid.complete))

}
