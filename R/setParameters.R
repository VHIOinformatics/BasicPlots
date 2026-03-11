#' Estimates the parameters for the different graphics
#'
#' @param labels Vector with the sample names
#' @return a list of the parameters needed to perform graphics
#' @export

setParameters <- function(labels)
{

    n_samples <- length(labels) # number of samples
    max_label <- max(nchar(labels)) # maximum length of the label

    wid <- 3500
    hei <- 3500
    res <- 400

    # define text size according to the total number of samples
    if (n_samples < 35) {
        if (max_label < 10) {
            ce <- 0.8
        } else if (max_label < 15) {
            ce <- 0.6
        } else {
            ce <- 0.4
        }
    } else if (n_samples < 90) {
        ce <- 0.4
    } else if (n_samples < 140) {
        ce <- 0.3
    } else {
        ce <- 0.2
    }
    return(list(wid = wid, hei = hei, res = res, ce = ce))
}
