#' @title  Random Walk Class Cover Catch Digraph Classifier
#'
#' @description Random Walk Class Cover Catch Digraph Classification.
#'
#' @param x asd.
#' @param y asd
#' @param method asd
#' @param m asd
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
#' @importFrom  Rfast dista
#'
#' @references
#' asd
#'
#' @examples
#'
#' rnorm(100)
#'
#' @rdname rwcccd_classifier
#' @export

rwcccd_classifier <- function(x, y, method = NULL, m = 1, proportion = 1) {
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
    n_other <- nrow(x_other)

    n_main <- nrow(x_main)

    dist_main2other <- Rfast::dista(xnew = x_main, x = x_other)
    dist_main2main <- Rfast::Dist(x = x_main)
    dist_main2all <- cbind(dist_main2main, dist_main2other)

    # w <- c(rep(n_other/n_main, n_main), rep(-1, n_other))
    w <- c(rep(1/n_main, n_main), rep(-1/n_other, n_other))
    # w <- c(rep(1/n, n_main), rep(-1/n, n_other))

    radii <- radii_rw(dist_main2all = dist_main2all,
                      dist_main2main = dist_main2main,
                      w = w,
                      m = m,
                      p = p)

    M <- dist_main2main < c(radii)
    M <- matrix(as.numeric(M), n_main)

    cover <- rep(0, n_main)
    thresh <- n_main*proportion

    m_dominant <- f_cover_rcpp(cover = cover,
                               thresh = thresh,
                               M = M,
                               dist_main2main = dist_main2main,
                               dist_main2other = radii)

    i_dominant_list[[i]] <- m_dominant$i_dominant
    radii_dominant_list[[i]] <- radii[m_dominant$i_dominant]
    x_dominant_list[[i]] <- x_main[m_dominant$i_dominant,]
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
