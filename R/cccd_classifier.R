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
  proportions <- c()

  for (i in 1:k_class) {
    i_main <- which(y == class_names[i])
    x_main <- x[i_main,]
    x_other <- x[-i_main,]

    n_main <- nrow(x_main)

    dist_main2other <- RANN::nn2(data = x_other, query = x_main, k = 1)$nn.dist # nearest opposite class sample distance
    dist_main2main <- Rfast::Dist(x = x_main) # main class distance matrix
    M <- dist_main2main < c(dist_main2other) # main class observers which are nearer than opposite class nearest sample
    M <- matrix(as.numeric(M), n_main)

    cover <- rep(0, n_main) # cover vector
    thresh <- n_main*proportion # threshold

    m_dominant <- f_cover_rcpp(cover = cover,
                               thresh = thresh,
                               dist_main2main = dist_main2main,
                               dist_main2other = dist_main2other,
                               M = M)

    i_dominant_list[[i]] <- m_dominant$i_dominant # dominant main class indexes
    x_dominant_list[[i]] <- x_main[m_dominant$i_dominant,,drop = FALSE] # dominant main class samples
    radii_dominant_list[[i]] <- dist_main2other[m_dominant$i_dominant,] # radius of dominant main class samples
    proportions[i] <- m_dominant$cover_proportion
  }

  results <- list(
    i_dominant_list = i_dominant_list,
    x_dominant_list = x_dominant_list,
    radii_dominant_list = radii_dominant_list,
    class_names = class_names,
    k_class = k_class,
    proportions = proportions
  )
  class(results) <- "cccd_classifier"
  return(results)
}








