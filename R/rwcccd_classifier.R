#' @title  Random Walk Class Cover Catch Digraph Classifier
#'
#' @description Random Walk Class Cover Catch Digraph Classification.
#'
#' @param x feature matrix or dataframe.
#' @param y class factor variable.
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

rwcccd_classifier <- function(x, y, method = "a", m = 1, proportion = 0.99) {
  p <- ncol(x)
  n <- nrow(x)
  class_names <- levels(y)
  k_class <- length(class_names)

  i_dominant_list <- vector(mode = "list", length = k_class)
  x_dominant_list <- vector(mode = "list", length = k_class)
  radii_dominant_list <- vector(mode = "list", length = k_class)
  proportions <- c()
  T_score_list <- vector(mode = "list", length = k_class)

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

    w <- switch(method,
                a = c(rep(1, n_main), rep(-1, n_other)),
                b = c(rep(n_other/n_main, n_main), rep(-1, n_other))
    )

    cover_main <- rep(0, n_main)
    cover_all <- rep(0, n)
    code_dominant <- rep(0, n_main)
    radii <- rep(NA, n_main)
    T_score <- rep(NA, n_main)

    while (Inf) {
      T_scores <- rep(NA, n_main)
      rs <- rep(NA, n_main)
      i_uncovered <- which(cover_all == 0)
      i_potentialdom <- which(code_dominant == 0 & cover_main == 0)

      if (length(i_potentialdom) < 1) {
        break
      }

      for (j in (1:n_main)[i_potentialdom]) {
        dd <- dist_main2all[j,i_uncovered]
        ii_sort <- sort.int(dist_main2all[j,i_uncovered], index.return = TRUE)$ix
        dd_sorted <- dd[ii_sort]
        R_all <- cumsum(w[ii_sort])
        P_all <- m*dd_sorted^p
        R_penalized_all <- R_all - P_all
        i_R_penalized_max <- which.max(R_penalized_all)
        rs[j] <- dd_sorted[i_R_penalized_max]
        n_uncovered <- sum(dd > rs[j])
        T_scores[j] <- R_all[i_R_penalized_max] - rs[j]/max(dist_main2main[j,]) * (n_uncovered/2)
      }

      i_selected <- which.max(T_scores)
      T_score_selected <- T_scores[i_selected]
      r_selected <- rs[i_selected]
      i_cover_main <- which(dist_main2main[i_selected,cover_main == 0] < r_selected)
      if (length(i_cover_main) > 0) {
        cover_main[i_cover_main] <- 1
      }
      i_cover_all <- which(dist_main2all[i_selected,i_uncovered] < r_selected)
      if (length(i_cover_all) > 0) {
        cover_all[i_cover_all] <- 1
      }
      code_dominant[i_selected] <- 1
      radii[i_selected] <- r_selected
      T_score[i_selected] <- T_score_selected

      if (sum(cover_main)/n_main >= proportion) {
        break
      }

      if (sum(code_dominant) == n_main) {
        break
      }

      if (sum(cover_main) == n_main) {
        break
      }
      cat("\r", "class:", i, "cover:", sum(cover_main)/n_main, "n_main:", n_main, "n_dom:", sum(code_dominant), "          ")
    }
    cat("\n")

    i_dominant_list[[i]] <- which(code_dominant == 1)
    radii_dominant_list[[i]] <- radii[i_dominant_list[[i]]]
    x_dominant_list[[i]] <- x_main[i_dominant_list[[i]],,drop = FALSE]
    T_score_list[[i]] <- T_score[i_dominant_list[[i]]]
    proportions[i] <- sum(cover_main)/n_main
  }

  results <- list(
    i_dominant_list = i_dominant_list,
    x_dominant_list = x_dominant_list,
    radii_dominant_list = radii_dominant_list,
    T_score_list = T_score_list,
    class_names = class_names,
    k_class = k_class,
    proportions = proportions
  )
  class(results) <- "rwcccd_classifier"
  return(results)
}
