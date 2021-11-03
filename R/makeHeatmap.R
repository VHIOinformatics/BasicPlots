#' Function to make heat map from a numeric matrix using pheatmap package

#' @param matrix: Numeric matrix with the intensity values to plot, with row.names and col.names
#' @param resultsDir: Output directory. Default = NULL
#' @param fileName: name of the output file, without extension. Default = NULL
#' @param fmtPlot: Format for the image file, "png" or "pdf" (default). If none specified image will pop up in R session
#' @param scaleColors: Vector with three colors to compose the color scale. Default = c("blue","white","red")
#' @param title: Title for each plot. Default = NULL
#' @param clustMethod: Clustering method. Default "ward.D2". Values given by parameter clustering_method of pheatmap function. Default = "ward.D2"
#' @param clustDistCols: Distance for clustering columns. Default "correlation". Values given by parameter clustering_distance_cols of pheatmap function
#' @param scale: Character indicating if the values should be centered and scaled. Values given by parameter scale of pheatmap function. efault = "row"
#' @param annot: data.frame that specifies the annotations shown for columns. rownames of the data.frame must be the colnames of the matrix object. Default = NA
#' @param annotColors: list with as many elements as conditions. Each element must be a vector with the color assignment for each category of the condition. Default = NA (assigned internally)
#' @param showRownames: Indicates whether you want to show rownames or not. Each row defines the features for a specific row. Default = TRUE
#' @param fontSizeRows: Size of row labels. Computed depending on the number of rows of matrix. Default = 10

#' @return The plot is created in the "resultsDir" with the name "fileName"

#' @export
#' @import pheatmap
#' @import gplots

makeHeatmap <- function(matrix, resultsDir = NULL, fileName = NULL, fmtPlot = "pdf",
                        title = NULL, scaleColors = c("blue","white","red"),
                        clustMethod = "ward.D2", clustDistCols = "correlation", scale = "row",
                        annot = NA, annotColors = NA, showRowNames = TRUE, fontSizeRows = 10, ...)
{

  heatcol<-colorRampPalette(scaleColors, space = "rgb")

  genes.l <-nrow(matrix)

  if (showRowNames & fontSizeRows==10){

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

  if (fmtPlot == "pdf") {

    pdf(file = file.path(resultsDir, paste(fileName, "pdf", sep = ".")))

  } else if (fmtPlot == "png") {

    parameters <- setParameters(colnames(matrix))
    png(file = file.path(resultsDir, paste(fileName,"png", sep = ".")),
        width = parameters$wid, height = parameters$hei,res = parameters$res)

  }

  pheatm<-pheatmap(matrix, color = heatcol(256), main = title,
         scale = scale, cluster_row = T, cluster_cols = T,
         clustering_method = clustMethod, clustering_distance_cols = clustDistCols,
         show_rownames = showRowNames, treeheight_row = 0,  fontsize_row = fontSizeRows,
         annotation_col = annot, annotation_colors = annotColors)

  if (fmtPlot %in% c("pdf", "png")) dev.off()

}
