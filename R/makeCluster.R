#' Creates a pdf or png with a dendrogram (hierarchical clusters) combining a distance and a linkage method
#'
#' @param estimates Numeric matrix with the intensity values to plot
#' @param resultsDir Output directory. Default = NULL
#' @param fileName name of the output file, without extension. Default = NULL
#' @param fmtPlot Format for the image file, pdf or png. If none specified images will be generated in R session. Default = "pdf"
#' @param title Title for each plot. Default = NULL
#' @param distance Distance measure for the hierarchical clustering. Either correlation or euclidean. Default = "correlation"
#' @param use.cor Used for correlation distance: method for computing covariances in the presence of missing values. Default = "pairwise.complete.obs"
#' @param cor.method Correlation coefficient (or covariance) to be computed. Default = "pearson"
#' @param clust.method Linkage method to cluster samples. Values given by parameter method of hclust function. Default = "ward.D2"
#' @param conditions Vector with the different conditions. Default = NULL
#' @param colors Vector with the colors assigned to each condition (in order of the unique(conditions)). Default = NULL
#' @param wid Plot parameter: width of the device (in pixels). Internal default (3500) will be used if not specified. Default = NULL
#' @param hei Plot parameter: height of the device (in pixels). Internal default (3500) will be used if not specified. Default = NULL
#' @param res Plot parameter: nominal resolution in ppi. Internal default (400) will be used if not specified. Default = NULL
#' @param ce Plot parameter: character expansion is the numerical value giving the amount by which plotting text and symbols should be magnified relative to the default. Will be defined internally according to character length and number of samples if not specified. Default = NULL
#' @param ... Other elements of plots
#'
#' @import RColorBrewer
#' @return a pdf and/or png with the hierarchical cluster
#' @export oneCluster

oneCluster <- function(estimates, resultsDir = NULL, fileName = NULL, fmtPlot = "pdf", title = NULL, distance="correlation", use.cor = "pairwise.complete.obs", cor.method = "pearson", clust.method = "ward.D2", conditions = NULL, colors = NULL, wid = NULL, hei = NULL, res = NULL, ce = NULL, ...)
{

  # parameters
  labels <- colnames(estimates)
  parameters <- setParameters(labels)

  if (is.null(wid)) wid <- parameters$wid
  if (is.null(hei)) hei <- parameters$hei
  if (is.null(res)) res <- parameters$res
  if (is.null(ce))  ce  <- parameters$ce

  # make clustering
  if (distance == "correlation") {
    clust <- hclust(as.dist(1 - cor(estimates, use = use.cor, method = cor.method)), method = clust.method)
    xlab <- paste("Correlation", clust.method, sep = "-")

  } else if(distance == "euclidean") {
    clust <- hclust(dist(t(estimates)), method = clust.method)
    xlab <- paste("Euclidean", clust.method, sep="-")
  }


  # initiate plots
  if (fmtPlot == "pdf") {

    pdf(file = file.path(resultsDir, paste(fileName, "pdf", sep = ".")), width = wid/res,
        height = hei/res)

  } else if (fmtPlot == "png") {

    png(file = file.path(resultsDir, paste(fileName,"png", sep=".")),
        width = wid, height = hei, res = res)

  }


  if (is.null(conditions)) { # plot with no conditions defined

    opt <- par(cex.main = 1, cex = ce, cex.axis = ce)
    plot(clust, main = title, xlab = xlab, ylab = "", sub = "", ...)
    par(opt)

  } else { # plot with defined conditions

    unique_cond <- unique(as.character(sort(conditions)))


    if (is.null(colors)) { # if user didn't add colors, we set our custom defaults

      custom_palette = c("#288A6C", "#CF7937", "#8682BB", "#BB4F86", "#E3BC51","#8CB369", "#A6CEE3","#F1AEAD", "#7C7B7B","#1F78B4", "#D8BFD8", "#854849" ,"#5D7344","#10438B", "#FDBF6F", "#61A175", "#BB4042", "#F18684", "#5C3F7C", "#D0A376", "#EDED90","#B2DF8A")


      color_list <- custom_palette[seq_along(unique_cond)]

    } else {

      color_list <- unique(colors)

    }

    # associate labels to colors
    map = setNames(color_list, unique_cond)
    colors <- map[conditions]


    opt <- par(cex.main = 1, cex.axis = ce, cex = ce)
    clust <- colorCluster(clust, colors, ce)
    plot(clust, main = title, xlab = xlab, ...)
    legend("topright",legend = unique_cond, cex = ce + 0.2, fill = color_list)

  }

  if (fmtPlot %in% c("pdf", "png")) dev.off()

  return(clust)
}


#' Creates a pdf with different dendrograms (hierarchical clusters) combining a distance and a linkage method
#'
#' @param estimates Numeric matrix with the intensity values to plot. The matrix should have samples as columns and genes/transcripts/exons as rows
#' @param use.cor Used for correlation distance: method for computing covariances in the presence of missing values. Default = "pairwise.complete.obs"
#' @param cor.method Correlation coefficient (or covariance) to be computed. Default = "pearson"
#' @param resultsDir Output directory. Default = NULL
#' @param fileName name of the output file, without extension. Default = "manyClusters"
#' @param fmtPlot Format for the image file. If none specified images will be generated in R session. Default = "pdf"
#' @param title Title for each plot. Default = NULL
#' @param conditions Vector with the different conditions. Default = NULL
#' @param colors Vector with the colors assigned to each condition. Default = NULL
#' @param wid Plot parameter: width of the device (in pixels). Internal default (3500) will be used if not specified. Default = NULL
#' @param hei Plot parameter: height of the device (in pixels). Internal default (3500) will be used if not specified. Default = NULL
#' @param res Plot parameter: nominal resolution in ppi. Internal default (400) will be used if not specified. Default = NULL
#' @param ce Plot parameter: character expansion is the numerical value giving the amount by which plotting text and symbols should be magnified relative to the default. Will be defined internally according to character length and number of samples if not specified. Default = NULL
#'
#' @import RColorBrewer
#' @return a pdf with the hierarchical clusters
#' @export manyClusters

manyClusters <- function(estimates, use.cor = "pairwise.complete.obs", cor.method = "pearson", resultsDir = NULL, fileName = "manyClusters", fmtPlot = "pdf", title = NULL, conditions = NULL, colors = NULL, wid = NULL, hei = NULL, res = NULL, ce = NULL)
{

  labels <- colnames(estimates)
  parameters <- setParameters(labels)

  if (is.null(wid)) wid <- parameters$wid
  if (is.null(hei)) hei <- parameters$hei
  if (is.null(res)) res <- parameters$res
  if (is.null(ce))  ce  <- parameters$ce

  # compute all clusters
  clust_list <- list(
    "Correlation Ward.D2" = hclust(as.dist(1 - cor(estimates, use = use.cor, method = cor.method)), method = "ward.D2"),
    "Correlation Average" = hclust(as.dist(1 - cor(estimates, use = use.cor, method = cor.method)), method = "average"),
    "Correlation Complete" = hclust(as.dist(1 - cor(estimates, use = use.cor, method = cor.method)), method = "complete"),
    "Euclidean Ward.D2" = hclust(dist(t(estimates)), method = "ward.D2"),
    "Euclidean Average" = hclust(dist(t(estimates)), method = "average"),
    "Euclidean Complete" = hclust(dist(t(estimates)), method = "complete")
  )

  # initiate plots
  if (fmtPlot == "pdf") {
    pdf(file = file.path(resultsDir, paste(fileName, "pdf", sep = ".")),
        width = wid/res,
        height = hei/res)
  }

  if (!is.null(conditions)) {

    unique_cond <- unique(as.character(sort(conditions)))

    if (is.null(colors)) { # if user didn't add colors, we set our custom defaults

      custom_palette = c("#288A6C", "#CF7937", "#8682BB", "#BB4F86", "#E3BC51","#8CB369", "#A6CEE3","#F1AEAD", "#7C7B7B","#1F78B4", "#D8BFD8", "#854849" ,"#5D7344","#10438B", "#FDBF6F", "#61A175", "#BB4042", "#F18684", "#5C3F7C", "#D0A376", "#EDED90","#B2DF8A")

      color_list <- custom_palette[seq_along(unique_cond)]

    } else {

      color_list <- unique(colors)

    }

    # map conditions to colors
    map = setNames(color_list, unique_cond)
    sample_colors <- map[conditions]

  }

  # plot all clusters
  opt <- par(cex.main = 1, cex.axis = ce, cex = ce)

  for (clustplot in names(clust_list)) {
    dend <- clust_list[[clustplot]]

    if (!is.null(conditions)) { # colors only added if conditions are provided
      dend <- colorCluster(dend, sample_colors, ce)
    }

    plot(dend, main = title, xlab = clustplot, ylab = "", sub = "")

    if (!is.null(conditions)) { # legend is added if conditions are provided
      legend("topright", legend = unique_cond, cex = ce + 0.2, fill = color_list)
    }

    # update the cluster object in the list in case it was colored
    clust_list[[clustplot]] <- dend
  }

  par(opt)


  if (fmtPlot == "pdf") dev.off()

  return(clust_list)

}


#' Colors an hclust object according to a vector of conditions

#' @param hclus an hclust object obtained with the hclust function
#' @param condition a numeric vector of the length=numb of samples, each number is a condition(ex: condition=c(1,1,1,1,2,2,2,2))
#' @param ce cex parameter that is set in function of the length of the characters of the vector
#'
#' @return returns the cluster with colored leafs
#' @export colorCluster

colorCluster <- function(hclus, condition, ce)
{

  sampleDendrogram <- as.dendrogram(hclus)
  names(condition) <- hclus$labels

  sampleDendrogram <- dendrapply(sampleDendrogram, function(estimates, batch) {
    if (is.leaf(estimates)) {
      label <- attr(estimates, "label")
      attr(estimates, "nodePar") <- list(lab.col = as.vector(batch[label]), pch = "", cex=0.8, lab.cex = ce)

    }
    estimates

  }, condition)

}
