#' @title Pure and Proper Class Cover Catch Digraph Ensemble Classifier
#'
#' @description \code{pcccd_ensemble_classifier} fits an Ensemble Pure and Proper Class Cover Catch
#' Digraph (PCCCD) classification model.
#'
#' @param x feature matrix or dataframe.
#' @param y class factor variable.
#' @param n_model an integer. Number of weak classifiers.
#' @param n_var an integer. number of variables in weak classifiers.
#' @param replace a bool. Should replacement be used in data sampling
#' @param prop_sample a value between 0 and 1. Proportion the number of resampled
#' samples to the number of samples in x.
#' @param min_proportion Minimum proportion of cover proportion in weak classifiers.
#' @param max_proportion Maximum proportion of cover proportion in weak classifiers.
#' @param verbose prints information. Default is FALSE.
#'
#' @details
#' Bagging framework for PCCCD.
#'
#' @return an object of "cccd_classifier" which includes:
#'  \item{i_dominant_list}{dominant sample indexes.}
#'  \item{x_dominant_list}{dominant samples from feature matrix, x}
#'  \item{radii_dominant_list}{Radiuses of the circle for dominant samples}
#'  \item{class_names}{class names}
#'  \item{k_class}{number of classes}
#'  \item{proportions}{proportions each class covered}
#'
#' @author Fatih Saglam, saglamf89@gmail.com
#'
#' @importFrom  stats runif
#'
#' @examples
#' n <- 1000
#'
#' @rdname pcccd_ensemble_classifier
#' @export

pcccd_ensemble_classifier <-
  function(x,
           y,
           n_model = 30,
           n_var = ncol(x),
           replace = FALSE,
           prop_sample = ifelse(replace, 1, 0.67),
           min_proportion = 0.7,
           max_proportion = 1,
           verbose = FALSE) {
    if (min_proportion < 0 | min_proportion > 1) {
      stop("min_proportion must be in range [0,1]")
    }

    if (max_proportion < 0 | max_proportion > 1) {
      stop("max_proportion must be in range [0,1]")
    }

    if (!is.matrix(x) & !is.data.frame(x)) {
      stop("x must be a matrix or data.frame")
    }

    if (is.data.frame(x)) {
      x <- as.matrix(x)
    }

    class_names <- levels(y)
    k_class <- length(class_names)
    n <- nrow(x)
    p <- ncol(x)

    m_list <- list()
    var_list <- list()
    for (i in 1:n_model) {
      i_sample <- sample(1:n, size = n * prop_sample, replace = replace)
      i_var <- sample(1:p, size = n_var)

      x_boot <- x[i_sample, i_var, drop = FALSE]
      y_boot <- y[i_sample]

      m_list[[i]] <-
        pcccd_classifier(
          x = x_boot,
          y = y_boot,
          proportion = runif(n = 1, min_proportion, max_proportion)
        )
      var_list[[i]] <- i_var
      if (verbose) {
        cat("\r", round(i/n_model*100, 2), "%")
      }
    }
    results <- list(m_list = m_list,
                    var_list = var_list,
                    k_class = k_class,
                    n_model = n_model,
                    class_names = class_names)
    class(results) <- "pcccd_ensemble_classifier"
    return(results)
  }
