#' @title  Class Cover Catch Digraph Predictor
#'
#' @description Class Cover Catch Digraph Prediction
#'
#' @param object asd.
#' @param newdata asd
#' @param type asd
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
#' @importFrom  proxy dist
#'
#'
#' @examples
#'
#' rnorm(100)
#'
#' @rdname predict.cccd_classifier
#' @export

predict.cccd_classifier <- function(object, newdata, type = "pred", ...) {
  x_dominant_list <- object$x_dominant_list
  radii_dominant_list <- object$radii_dominant_list
  class_names <- object$class_names
  k_class <- object$k_class

  x <- newdata
  n <- nrow(x)

  dist_prop <- matrix(data = NA, nrow = n, ncol = k_class)

  for (i in 1:k_class) {
    dist_x2dom <- as.matrix(proxy::dist(x_dominant_list[[i]], x))
    prop_x2dom <- dist_x2dom/radii_dominant_list[[i]]
    dist_prop[,i] <- apply(prop_x2dom, 2, min)
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










