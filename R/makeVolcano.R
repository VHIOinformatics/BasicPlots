#' Function to make volcano plot from a topTable extracted from limma

#' @param dataV data.frame as the output obtained through function topTable() in limma package
#' @param resultsDir Output directory. Default = volcanoDir
#' @param fileName Name of the output file, without extension. Default = NULL
#' @param fmtPlot Format for the image file, "png" or "pdf" (static plot); or "html" (interactive). If none specified ("") image will pop up in R session. Default = ""
#' @param scaleColors Vector with three colors to compose the color scale. Default = c("blue","grey","red")
#' @param title Title for each plot. Default = NULL
#' @param yVal Value to plot on the y axis. Possible values are p.value or p.adjust. Values will be transformed to minus(log2(yVal)). Default = p.adjust
#' @param p.val P-value cutoff for topTable. Default = NULL
#' @param p.adj Adjusted p-value cutoff for topTable. Default = 0.05
#' @param thres.logFC LogFC cutoff for topTable. Default = 1
#' @param coefs Coefficients of fit.main results. Will be used for the limits of x in the volcano plot. Default = fit.main$coefficients
#' @param topGenes Number of genes to label. Default = 20
#' @param geneLabel Variable with the label of genes to be shown in plot. Default = Geneid
#' @param sortLabel Sorting criteria for the label highlight. Specify "logFC" to plot the labels of the top most different genes in absolute fold change. Otherwise, plot will show the labels of the top most significant (according to previously defined cutoff). Default = p.value
#' @param annotFile Dataframe that includes the equivalence of Geneid and geneLabel, if geneLabel is specified. Default = NULL
#' @param interactive Prints an interactive plot showing the geneLabel, FC and (adj)p.value when hovering the points. Shows plot in session. For saving interactive plot, use fmtPlot="html". Default = FALSE
#'
#' @return The plot is created in the "resultsDir" with the name "fileName" or in session if the format is not specified.

#' @export

#' @import ggplot2
#' @import ggrepel
#' @import plotly
#' @import htmlwidgets

makeVolcano <- function(topTable, resultsDir = volcanoDir, fileName = NULL, fmtPlot = "", title = NULL, scaleColors = c("blue", "grey", "red"), yVal="p.adjust", p.val = NULL, p.adj = 0.05, thres.logFC=1, coefs=fit.main$coefficients ,topGenes = 20, maxOverlaps = 25, geneLabel = "Geneid", sortLabel="p.value", annotFile = NULL, interactive = FALSE, ...)
{

  colorS <- scaleColors

  maxxlim <- max(abs(coefs))
  xlim <- xlim(-ceiling(maxxlim),ceiling(maxxlim))
   ## ceiling rounds up to nearest integer

  #topTable <- topTable(fit.main.CAF, n = Inf, coef = i, adjust = "fdr")
  dataV <- topTable
  dataV <- dataV %>% mutate(Geneid = rownames(dataV), minlogp = -(log10(P.Value)), minlogadjp = -(log10(adj.P.Val)), FC = ifelse(logFC>0, 2^logFC, -(2^abs(logFC))))

  if(is.null(p.val)) { # if p.value is null we use p.adjust as significance threshold
    dataV <- dataV %>%
    mutate(sig = ifelse(adj.P.Val<p.adj & logFC > thres.logFC, "UP", ifelse(adj.P.Val<p.adj & logFC < (-thres.logFC), "DN","n.s")))

    if (sortLabel == "logFC") {
      dataV$abslogFC <- abs(dataV$logFC)
      dataV <- dataV[order(dataV$abslogFC, decreasing = TRUE),] # order sorted by logFC
    } else {
      dataV <- dataV[order(dataV$adj.P.Val),] # order sorted by p.adj
    }

  } else { # if p.value is not null, we use p.value as threshold
    dataV <- dataV %>%
      mutate(sig = ifelse(P.Value<p.val & logFC > thres.logFC, "UP", ifelse(P.Value<p.val & logFC < (-thres.logFC), "DN","n.s")))


    if (sortLabel == "logFC") {
      dataV$abslogFC <- abs(dataV$logFC)
      dataV <- dataV[order(dataV$abslogFC, decreasing = TRUE),] # sort by logFC
    } else {
      dataV <- dataV[order(dataV$P.Value),] # order sorted by p.val
    }
  }

  if (geneLabel != "Geneid") { # if gene label is not Geneid then we will merge the annotation so it can be used
    if (!is.null(annotFile)) {
      dataV <- merge(annotFile[,c("Geneid",geneLabel)],dataV,by="Geneid")
      if (is.null(p.val)) {
        if (sortLabel == "logFC") {
          dataV <- dataV[order(dataV$abslogFC, decreasing = TRUE),]
        } else { # per pvalue, default
          dataV <- dataV[order(dataV$adj.P.Val),]
        }

      } else { # pvalue instead of padj
        if (sortLabel == "logFC") {
          dataV <- dataV[order(dataV$abslogFC, decreasing = TRUE),]
        } else { # per pvalue, default
          dataV <- dataV[order(dataV$P.Value),]
        }
      }
    } else {
      stop("Please, provide an annotation file with the specified 'geneLabel' and Geneid, using 'annotFile' parameter.\n")
    }
  }

  if (yVal == "p.adjust") {
    ylim <- ylim(0,ceiling(max(dataV$minlogadjp)))
    topData = head(dataV[dataV$sig != "n.s",],topGenes) # the genes that will be labeled in plot
    p <- ggplot(data=dataV, aes(x=logFC, y=minlogadjp )) +
      geom_point(alpha = 1, size= 1, aes(col = sig,
                                         text= paste(dataV[,c(geneLabel)],
                                                     "\nFC =", round(FC,3),
                                                     "\nadj.pval =", formatC(adj.P.Val,format="e",digits=3)))) +
      scale_color_manual(values = colorS) + xlim + ylim + labs(col=" ") +
      geom_vline(xintercept = thres.logFC, linetype= "dotted") + geom_vline(xintercept = -thres.logFC, linetype= "dotted") +
      geom_hline(yintercept = -log10(p.adj), linetype= "dotted")  +  theme_bw() + ggtitle(title)
    # for interactive, expressions don't work
    intp <- p + xlab("log2FC") + ylab("-log10(adj.pval)")
    intp <- ggplotly(intp, tooltip = "text")
    # but it does work with normal ggplot so we maintain the labs after
    p <- p + xlab(expression("log"[2]*"FC")) + ylab(expression("-log"[10]*"(adj.pval)")) +
      geom_text_repel(data = topData, aes(label = topData[,c(geneLabel)]), max.overlaps = maxOverlaps) # we hard label for only non-interactive plot

  } else {
    ylim <- ylim(0,ceiling(max(dataV$minlogp)))
    topData = head(dataV[dataV$sig != "n.s",],topGenes)
    p <- ggplot(data=dataV, aes(x=logFC, y=minlogp )) +
      geom_point(alpha = 1, size= 1, aes(col = sig,
                                         text= paste(dataV[,c(geneLabel)],
                                                     "\nFC =", round(FC,3),
                                                     "\np.val =", formatC(P.Value,format="e",digits=3)))) +
      scale_color_manual(values = colorS) + xlim + ylim + labs(col=" ") +
      geom_vline(xintercept = thres.logFC, linetype= "dotted") + geom_vline(xintercept = -thres.logFC, linetype= "dotted") +
      geom_hline(yintercept = -log10(p.val), linetype= "dotted")  +  theme_bw() + ggtitle(title)
    # for interactive, expressions don't work
    intp <- p + xlab("log2FC") + ylab("-log10(p.val)")
    intp <- ggplotly(intp, tooltip = "text")
    # but it does work with normal ggplot so we maintain the labs after
    p <- p + xlab(expression("log"[2]*"FC")) + ylab(expression("-log"[10]*"(p.val)")) +
      geom_text_repel(data = topData, aes(label = topData[,c(geneLabel)]), max.overlaps = maxOverlaps)

  }


  if(fmtPlot %in% c("png","pdf")) {
    ggsave(p, filename= file.path(resultsDir, paste(fileName, fmtPlot, sep = ".")), device = fmtPlot, ...)
  } else if (fmtPlot == "html") {
    htmlwidgets::saveWidget(as_widget(intp), file.path(resultsDir, paste(fileName, fmtPlot, sep = ".")))

  } else {
    if(interactive == TRUE) {
      print(intp)
    } else {
      print(p)
    }
  }


}
