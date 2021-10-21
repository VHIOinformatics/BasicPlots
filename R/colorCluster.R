#' Colors an hclust object according to a vector of conditions

#' @param hclus: an hclust object obtained with the hclust function
#' @param condition: a numeric vector of the length=numb of samples, each number is a condition(ex: condition=c(1,1,1,1,2,2,2,2))
#' @param Ce: cex parameter that is set in function of the length of the characters of the vector
#' 
#' @return returns the cluster with colored leafs

colorCluster <- function(hclus, condition, ce) 
{
  
  sampleDendrogram <- as.dendrogram(hclus)
  names(condition) <- hclus$labels
  
  sampleDendrogram <- dendrapply(sampleDendrogram, function(estimates, batch) {
    if(is.leaf(estimates)){
      label <- attr(estimates, "label")
      attr(estimates, "nodePar") <- list(lab.col = as.vector(batch[label]), pch = "", cex=0.8, lab.cex = ce)
      
    }
    estimates
    
  }, condition)
  
} 