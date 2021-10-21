#' Estimates the parameters for the different graphics
#' 
#' @param labels: Vector with the sample names
#' @return a list of the parameters needed to perform graphics

setParameters <- function (labels) 
{
    
    long <- length(labels)
    
    wid <- 3500
    hei <- 3500
    res <- 400
    
    if (max(nchar(labels)) < 10) {
        ce <- 0.8
    }else if (max(nchar(labels)) < 15){
        ce <- 0.6
    }else {
        ce <- 0.4
    }

    if (long < 35) {
        if (max(nchar(labels)) < 10) {
            ce <- 0.8
        }else if (max(nchar(labels)) < 15){
            ce <- 0.6
        }else {
            ce <- 0.4
        }
    } else if (long < 90) {
        ce <- 0.4
    } else if (long < 140){
        ce <- 0.3
    } else {
        ce <- 0.2
    }
    return(list(wid = wid, hei = hei, res = res, ce = ce))
}
