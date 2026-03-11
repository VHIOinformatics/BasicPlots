#' Function to make heat map from a numeric matrix using pheatmap package

#' @param matrix Numeric matrix with the intensity values to plot, with row.names (usually genes) and col.names (samples)
#' @param resultsDir Output directory of the plot. Default = NULL
#' @param fileName Name of the output file, without extension. Default = NULL
#' @param fmtPlot Format for the image file, "png" or "pdf". If none specified image will pop up in R session. Default = ""
#' @param scaleColors Vector with three colors to compose the color scale. Default = c("blue","white","red")
#' @param title Title of the plot. Default = ""
#' @param clusterRow Whether to cluster rows, or hclust object. Equivalent to parameter cluster_row of pheatmap. Default = TRUE
#' @param clusterCol Whether to cluster columns, or hclust object. Equivalent to parameter cluster_col of pheatmap. Default = TRUE
#' @param clustMethod Clustering method. Default "ward.D2". Values given by parameter clustering_method of pheatmap function. Default = "ward.D2"
#' @param clustDistCols Distance for clustering columns. Values given by parameter clustering_distance_cols of pheatmap function. Default = "correlation"
#' @param scale Character indicating if the values should be centered and scaled. Values given by parameter scale of pheatmap function. Default = "row"
#' @param annot Dataframe that specifies the annotations shown for columns. rownames of the data.frame must be the colnames of the matrix object. Default = NA
#' @param annotColors List with as many elements as conditions. Each element must be a vector with the color assignment for each category of the condition. Default = NA (assigned internally)
#' @param showRowNames Indicates whether you want to show rownames or not. Each row defines the features for a specific row. A maximum of 150 genes is possible, if this number is exceeded the parameter will automatically be turned off. Default = TRUE
#' @param fontSizeRows Size of row labels. Computed depending on the number of rows of matrix. Default = 10
#' @param showRowCluster Whether to show the tree of the row clustering (genes). Default = FALSE
#' @param showColCluster Whether to show the tree of the column clustering (samples). Default = TRUE
#' @param wid Plot parameter: width of the device (in pixels). Internal default (3500) will be used if not specified. Default = NULL
#' @param hei Plot parameter: height of the device (in pixels). Internal default (3500) will be used if not specified. Default = NULL
#' @param res Plot parameter: nominal resolution in ppi. Internal default (400) will be used if not specified. Default = NULL
#' @param ... Additional parameters of pheatmap

#' @return The heatmap plot is created in the "resultsDir" with the name "fileName"

#' @export
#' @import pheatmap
#' @import gplots

makeHeatmap <- function(matrix, resultsDir = NULL, fileName = NULL, fmtPlot = "", title = "", scaleColors = c("blue","white","red"), clusterRow = TRUE, clusterCol = TRUE, clustMethod = "ward.D2", clustDistCols = "correlation", scale = "row", annot = NA, annotColors = NA, showRowNames = TRUE, fontSizeRows = 10, showRowCluster = FALSE, showColCluster = TRUE, wid = NULL, hei = NULL, res = NULL, ...)
{

  heatcol <- colorRampPalette(scaleColors, space = "rgb")

  # check how to plot the gene names
  genes.l <- nrow(matrix) # number of genes to plot

  if (showRowNames & fontSizeRows == 10){ # if the fontSizeRows is specified (diff to default) this will be ignored

    if (genes.l < 40) {
      fontSizeRows = 10
    } else if (genes.l< 70) {
      fontSizeRows = 8
    } else if (genes.l < 100) {
      fontSizeRows = 6
    } else if (genes.l < 150) {
      fontSizeRows = 4
    } else {
      showRowNames = FALSE
      warning(paste(genes.l, "rows are too many rows in the heatmap for cex.row to be defined. Maximum is 150 rows"), call. = FALSE)
    }
  }

  # plot parameters
  parameters <- setParameters(colnames(matrix))
  if (is.null(wid)) wid <- parameters$wid
  if (is.null(hei)) hei <- parameters$hei
  if (is.null(res)) res <- parameters$res

  if (fmtPlot == "pdf") {

    pdf(file = file.path(resultsDir, paste(fileName, "pdf", sep = ".")), width = wid/res, height = hei/res)

  } else if (fmtPlot == "png") {


    png(file = file.path(resultsDir, paste(fileName,"png", sep = ".")),
        width = wid, height = hei, res = res)

  }

  if (showRowCluster) {
    treeheightRow = 50
  } else {
    treeheightRow = 0
  }

  if (showColCluster) {
    treeheightCol = 50
  } else {
    treeheightCol = 0
  }

  if (is.na(annotColors)) { # if colors not provided use our palette

    custom_palette = c("#288A6C", "#CF7937", "#8682BB", "#BB4F86", "#E3BC51","#8CB369", "#A6CEE3","#F1AEAD", "#7C7B7B","#1F78B4", "#D8BFD8", "#854849" ,"#5D7344","#10438B", "#FDBF6F", "#61A175", "#BB4042", "#F18684", "#5C3F7C", "#D0A376", "#EDED90","#B2DF8A")

    conditions <- levels(annot[,1]) # must have only 1 condition which must be a factor
    myCols <- custom_palette[1:length(conditions)]
    names(myCols) <- conditions # named list with all conditions of the variable
    condCols <- myCols[names(myCols) %in% annot[,1]] # make sure we get only the comparison of interest
    annotColors <- list(condCols)
    names(annotColors) <- colnames(annot)[1] # make a colored list with same column name as the given annot df

  }

  pheatm <- pheatmap(matrix, color = heatcol(256), main = title,
                     scale = scale, cluster_row = clusterRow, cluster_cols = clusterCol,
                     clustering_method = clustMethod, clustering_distance_cols = clustDistCols,
                     show_rownames = showRowNames, treeheight_row = treeheightRow,
                     treeheight_col = treeheightCol,
                     fontsize_row = fontSizeRows,
                     annotation_col = annot, annotation_colors = annotColors, ...)

  if (fmtPlot %in% c("pdf", "png")) dev.off()

  if (fmtPlot == "") print(pheatm)



}
