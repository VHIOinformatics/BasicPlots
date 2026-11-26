#' Estimates the parameters for the different graphics
#'
#' The font size is kept readable and the page grows with the number of samples
#' (and with the length of the labels), instead of shrinking the text on a fixed
#' size canvas. Only if the page reaches max_in the font size is reduced.
#'
#' @param labels Vector with the sample names
#' @param base_ce Preferred character expansion for the labels. Default = 0.8
#' @param tree_in Height (inches) reserved for the tree itself. Default = 5
#' @param min_in Minimum page width (inches). Default = 7
#' @param max_in Maximum page width (inches); above this the font is reduced. Default = 30
#' @return a list of the parameters needed to perform graphics
#' @export

setParameters <- function(
    labels,
    base_ce = 0.8,
    tree_in = 5,
    min_in = 7,
    max_in = 30
) {
    n_samples <- length(labels) # number of samples
    max_label <- max(nchar(labels)) # maximum length of the label

    res <- 400
    pointsize <- 12 # default pointsize of pdf() and png()

    ce <- base_ce

    # space (inches) that each label needs along the x axis (text height + 25% gap)
    per_leaf_in <- function(ce) ce * pointsize * 1.25 / 72
    side_in <- 2 # room for the y axis and the right margin

    # page width grows with the number of samples; if it reaches max_in, shrink the font
    wid_in <- n_samples * per_leaf_in(ce) + side_in
    if (wid_in > max_in) {
        ce <- (max_in - side_in) / (n_samples * pointsize * 1.25 / 72)
        wid_in <- max_in
    }
    wid_in <- max(wid_in, min_in)

    # page height = tree + room for the (vertical) labels + title/xlab.
    # 0.6 * font size is a conservative average character width; exact margins are set in fitMargins()
    label_in <- max_label * 0.6 * ce * pointsize / 72
    hei_in <- tree_in + label_in + 1.5

    return(list(
        wid = round(wid_in * res),
        hei = round(hei_in * res),
        res = res,
        ce = ce
    ))
}
