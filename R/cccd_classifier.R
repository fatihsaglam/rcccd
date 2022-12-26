#' @title  Class Cover Catch Digraph Classifier
#'
#' @description Class Cover Catch Digraph Classification.
#'
#' @param x asd.
#' @param y asd
#' @param proportion asd
#'
#' @details
#' asd
#'
#' @return asd.
#'  \item{asd}{asd}
#'
#' @author Fatih Saglam, saglamf89@gmail.com
#'
#' @importFrom  RANN nn2
#' @importFrom  Rfast Dist
#'
#' @references
#' asd
#'
#' @examples
#'
#' rnorm(100)
#'
#' @rdname cccd_classifier
#' @export

cccd_classifier <- function(x, y, proportion = 1) {

  p <- ncol(x)
  class_names <- levels(y)
  k_class <- length(class_names)

  i_dominant_list <- vector(mode = "list", length = k_class)
  x_dominant_list <- vector(mode = "list", length = k_class)
  radii_dominant_list <- vector(mode = "list", length = k_class)

  i = 1
  proportion = 1
  for (i in 1:k_class) {
    i_main <- which(y == class_names[i])
    x_main <- x[i_main,]
    x_other <- x[-i_main,]

    n_main <- nrow(x_main)

    dist_main2other <- RANN::nn2(data = x_other, query = x_main, k = 1)$nn.dist
    dist_main2main <- Rfast::Dist(x = x_main)
    M <- dist_main2main < c(dist_main2other)
    M <- matrix(as.numeric(M), n_main)

    cover <- rep(0, n_main)
    thresh <- n_main*proportion

    i_dominant <- f_cover_rcpp(cover = cover,
                 thresh = thresh,
                 dist_main2main = dist_main2main,
                 dist_main2other = dist_main2other,
                 M = M)$i_dominant

    i_dominant_list[[i]] <- i_dominant
    x_dominant_list[[i]] <- x_main[i_dominant,,drop = FALSE]
    radii_dominant_list[[i]] <- dist_main2other[i_dominant,]
  }

  results <- list(
    i_dominant_list = i_dominant_list,
    x_dominant_list = x_dominant_list,
    radii_dominant_list = radii_dominant_list,
    class_names = class_names,
    k_class = k_class
  )
  class(results) <- "cccd_classifier"
  return(results)
}








