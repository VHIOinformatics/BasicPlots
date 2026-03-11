#' Creates a Venn diagram based on several gene lists
#'
#' @param listGenes Named list with gene vectors, with 2 to 4 elements. Names should be short strings of up to 10 characters
#' @param resultsDir Output directory. Default = NULL
#' @param fileName Name of the output file, without extension. Default = NULL
#' @param fmtPlot Format for the image file, "png" or "pdf". If none specified image will pop up in R session. Default = ""
#' @param title Title for each plot. Default = NULL
#' @param mkExcel Whether to generate an excel with gene intersection and specific sections. Default = TRUE
#' @param colors Vector of colors to paint the different main conditions of the Venn diagram. Intersecting sections will be colored with color overlap. If none specified, pre-defined default colors will be used.
#' @param labelSize Text size of the name label of each set. Equivalent to set_name_size from ggvenn. Default = 6
#' @param setSize Text size for the set number contents. Equivalent to text_size from ggvenn. Default = 5
#' @param percentage Whether to show the percentage of genes for each set. Default = FALSE
#' @param wid Plot parameter: width of the device (in inches). Default = 8
#' @param hei Plot parameter: height of the device (in inches). Default = 8
#' @param res Plot parameter: nominal resolution in ppi. Default = 400
#' @param ... Additional parameters to add to ggvenn function
#'
#' @import ggplot2
#' @import ggvenn
#' @import openxlsx

#' @return Venn diagram plot and list of genes in each set (in excel format) if specified
#' @export makeVenn

makeVenn <-  function (listGenes, resultsDir = NULL, fileName = NULL, fmtPlot = "",
                       title = NULL, mkExcel = TRUE, colors = NULL, labelSize = 6, setSize = 5, percentage = FALSE, wid = 8, hei = 8, res = 400, ...)
{

  if (is.null(colors)) { # default colors (only 4 because more than 4 comparisons is not advised)
    colors=c("#c33c54","#254e70","#43aa8b","#ffc857")
  }

  nList <- length(listGenes)
  if (nList > 4) {
    stop("Cannot make a Venn diagram with more than 4 comparisons as it would be too complex to interpret.")
  } else if (nList == 1) {
    stop("Cannot make a Venn diagram with only one list as there is nothing to compare.")
  }

  # create and plot ggvenn
  venn_plot <- ggvenn(
    data = listGenes,
    show_elements = FALSE,    # Optional: Set to TRUE if you want to display the elements
    fill_color = colors,  # Custom colors if desired
    stroke_size = 0.5,
    set_name_size = labelSize,
    text_size = setSize,
    show_percentage = percentage,   # Optional: Set to TRUE to show percentages in intersections
    ...
    )

  # extract data from ggvenn object
  venn_df <- as.data.frame(venn_plot@data)
  rownames(venn_df) <- venn_df$`_key`
  venn_df <- venn_df[,-1]
  # we have a table of whether each gene is in one comparison or the other (TRUE/FALSE)
  venn_df$pattern <- apply(venn_df, 1, function(x) {
    paste(as.numeric(x), collapse = "")
  }) # we assign a value of the pattern (00 in neither, 01 in second but not first, 10 in first but not second, 11 in both)
  #table(venn_df$pattern)
  intersectionSets <- split(rownames(venn_df), venn_df$pattern) # we get a list of the genes for each pattern

  # save file
  if (fmtPlot %in% c("png", "pdf")) {
    ggsave(filename = file.path(resultsDir, paste("VennDiagram", fileName, fmtPlot, sep = ".")),
           plot = venn_plot,
           device = fmtPlot,
           bg = "white",
           width = wid,
           height = hei,
           units = "in")
  } else {
    print(venn_plot)
  }


  # make excel with genes
  if (mkExcel) {
    hs1 <- createStyle(fgFill = "#737373", halign = "CENTER", textDecoration = "Bold",
                       border = "Bottom", fontColour = "white")
    wb <- createWorkbook()

    nSheets = length(intersectionSets)

    for (i in 1:(nSheets-1)){ # started from 2 before bc 00 was included, but not now. nSheets-1 keeps from adding the common 11, which is added later
      nameVector <- names(listGenes)[as.logical(as.numeric(unlist(strsplit(names(intersectionSets[i]), split = ""))))]
      nameSheet <- paste(nameVector, collapse = ".")

      addWorksheet(wb, sheetName = nameSheet)

      writeData(wb, unlist(intersectionSets[[i]]),
                sheet = nameSheet, startRow = 1, startCol = 1,
                headerStyle = hs1)
    }
    # intersection
    addWorksheet(wb, sheetName = "Common")

    writeData(wb, unlist(intersectionSets[[nSheets]]),
              sheet = "Common", startRow = 1, startCol = 1,
              headerStyle = hs1)

    saveWorkbook(wb, file = file.path(resultsDir,
                                      paste("VennGenes", fileName, "xlsx", sep = ".")),
                 overwrite = TRUE)
  }

}

