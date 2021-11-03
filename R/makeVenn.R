#' Creates a Venn diagram based on several gene lists
#'
#' @param listGenes: named list with gene vectors, with 2 to 4 elements. Names should be short strings of up to 10 characters
#' @param resultsDir: Output directory. Default = NULL
#' @param fileName: name of the output file, without extension. Default = NULL
#' @param fmtPlot: Format for the image file, "png" or "pdf" (default). If none specified image will pop up in R session
#' @param title: Title for each plot. Default = NULL
#' @param mkExcel: whether to generate an excel with gene intersection and specific sections. Default = TRUE
#' @param colors: vector of colors to paint the different sections of the Venn diagram
#'
#' @import Vennerable
#' @import colorfulVennPlot
#' @import openxlsx
#' @import RColorBrewer
#' @return a file
#' @export makeVenn

makeVenn <-  function (listGenes, resultsDir = NULL, fileName = NULL, fmtPlot = "pdf",
                       title = NULL, mkExcel = TRUE, colors)
{

  colors <- c(brewer.pal(12, "Set3"), brewer.pal(4, "Pastel2"))
  nList <- length(listGenes)
  vtest <- Venn(listGenes)
  vennData <- sapply(vtest@IntersectionSets, function(x) length(unlist(x)))

  # plot and save file
  if (fmtPlot == "png") {

    png(file.path(resultsDir, paste("VennDiagram", fileName, "png", sep = ".")))

  } else if (fmtPlot == "pdf") {

    pdf(file.path(resultsDir, paste("VennDiagram", fileName, "pdf", sep = ".")))

  }

  if (nList==2){

  plotVenn2d(vennData, rev(names(listGenes)), Colors = colors, Title = title, shrink = 1)

  } else if (nList==3){

    plotVenn3d(vennData, names(listGenes), Colors = colors, Title = title, shrink = 1)

  } else if (nList==4){

    plotVenn4d(vennData[-1], names(listGenes), Colors = colors, Title = title, shrink = 1) #notice we need to remove first element in this case

  }

  if (fmtPlot %in% c("pdf", "png")) dev.off()

if (mkExcel) {
    hs1 <- createStyle(fgFill = "#737373", halign = "CENTER", textDecoration = "Bold",
                       border = "Bottom", fontColour = "white")
    wb <- createWorkbook()

    nSheets = length(vtest@IntersectionSets)

    for (i in 2:(nSheets-1)){
    nameVector <- names(listGenes)[as.logical(as.numeric(unlist(strsplit(names(vtest@IntersectionSets[i]),
                                                                         split = ""))))]
    nameSheet <- paste(nameVector, collapse = ".")

    addWorksheet(wb, sheetName =  nameSheet)

    writeData(wb, unlist(vtest@IntersectionSets[[i]]),
              sheet =  nameSheet, startRow = 1, startCol = 1,
              headerStyle = hs1)
    }
    #intersection
    addWorksheet(wb, sheetName = "Common")

    writeData(wb, unlist(vtest@IntersectionSets[[nSheets]]),
              sheet = "Common", startRow = 1, startCol = 1,
              headerStyle = hs1)

    saveWorkbook(wb, file = file.path(resultsDir,
                                      paste("VennGenes", fileName, "xlsx", sep = ".")),
                 overwrite = TRUE)
  }

}

