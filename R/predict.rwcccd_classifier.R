#' @title  Random Walk Class Cover Catch Digraph Predictor
#'
#' @description Random Walk Class Cover Catch Digraph Prediction
#'
#' @param object asd.
#' @param newdata asd
#' @param type asd
#' @param e asd
#' @param ... asd
#'
#' @details
#' asd
#'
#' @return asd.
#'  \item{asd}{asd}
#'
#' @author Fatih Saglam, saglamf89@gmail.com
#'
#'
#' @references
#' asd
#'
#' @importFrom proxy dist
#' @importFrom Rfast colMins
#'
#'
#' @examples
#'
#' rnorm(100)
#'
#' @rdname predict.rwcccd_classifier
#' @export

predict.rwcccd_classifier <- function(object, newdata, type = "pred", e = 0, ...) {
  x_dominant_list <- object$x_dominant_list
  radii_dominant_list <- object$radii_dominant_list
  T_score_list <- object$T_score_list
  class_names <- object$class_names
  k_class <- object$k_class

  x <- newdata
  n <- nrow(x)

  dist_prop <- matrix(data = NA, nrow = n, ncol = k_class)

  for (i in 1:k_class) {
    dist_x2dom <- as.matrix(proxy::dist(x_dominant_list[[i]], x))
    prop_x2dom <- (dist_x2dom/radii_dominant_list[[i]])^(T_score_list[[i]]^e)
    dist_prop[,i] <- Rfast::colMins(prop_x2dom, value = TRUE)
  }
  prob <- 1 - t(apply(dist_prop, 1, function(m) m/sum(m)))

  if (type == "prob") {
    colnames(prob) <- class_names
    return(prob)
  }
  if (type == "pred") {
    pred <- factor(class_names[max.col(prob)], levels = class_names, labels = class_names)
    return(pred)
  }
}










