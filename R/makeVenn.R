#' Creates a Venn diagram based on several gene lists
#'
#' @param listGenes Named list with gene vectors, with 2 to 4 elements. Names should be short strings of up to 10 characters
#' @param resultsDir Output directory. Default = NULL
#' @param fileName Name of the output file, without extension. Default = NULL
#' @param fmtPlot Format for the image file, "png" or "pdf". If none specified image will pop up in R session. Default = ""
#' @param title Title for each plot. Default = NULL
#' @param mkExcel Whether to generate an excel with gene intersection and specific sections. Default = TRUE
#' @param colors Vector of colors to paint the different main conditions of the Venn diagram. Intersecting sections will be colored with color overlap. If none specified, pre-defined default colors will be used.
#' @param percentage Whether to show the percentage of genes for each set. Default = FALSE
#'
#' @import Vennerable
#' @import ggplot2
#' @import ggvenn
#' @import openxlsx

#' @return Venn diagram plot and list of genes in each set (in excel format) if specified
#' @export makeVenn

makeVenn <-  function (listGenes, resultsDir = NULL, fileName = NULL, fmtPlot = "",
                       title = NULL, mkExcel = TRUE, colors = NULL, percentage = FALSE)
{

  if (is.null(colors)) { #default colors
    colors=c("#c33c54","#254e70","#43aa8b","#ffc857")
  }

  nList <- length(listGenes)
  vtest <- Venn(listGenes)
  vennData <- sapply(vtest@IntersectionSets, function(x) length(unlist(x)))

  # plot and save file
  venn_plot <- ggvenn(
    data = listGenes,
    show_elements = FALSE,    # Optional: Set to TRUE if you want to display the elements
    fill_color = colors,  # Custom colors if desired
    stroke_size = 0.5,
    set_name_size = 6,
    text_size = 5,
    show_percentage = percentage   # Optional: Set to TRUE to show percentages in intersections
  )

  if (fmtPlot %in% c("png", "pdf")) {
    ggsave(venn_plot, filename= file.path(resultsDir, paste("VennDiagram",fileName, fmtPlot, sep = ".")), device = fmtPlot, bg="white", width=8, height=8, units="in")
  } else {
    print(venn_plot)
  }


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

