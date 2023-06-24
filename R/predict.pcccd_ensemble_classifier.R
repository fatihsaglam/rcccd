#' @title  Pure and Proper Class Cover Catch Digraph Ensemble Prediction
#'
#' @description \code{predict.pcccd_ensemble_classifier} makes prediction using
#' \code{pcccd_ensemble_classifier} object.
#'
#' @param object a \code{rwcccd_classifier} object
#' @param newdata newdata as matrix or dataframe.
#' @param type "pred" or "prob". Default is "pred". "pred" is class estimations,
#'  "prob" is \eqn{n\times k} matrix of class probabilities.
#' @param ... not used.
#'
#' @return a vector of class predictions (if type is "pred") or a \eqn{n\times p}
#' matrix of class probabilities (if type is "prob").
#'
#' @author Fatih Saglam, saglamf89@gmail.com
#'
#' @importFrom  stats predict
#'
#' @examples
#' n <- 1000
#'
#' @rdname predict.pcccd_ensemble_classifier
#' @export

predict.pcccd_ensemble_classifier <- function(object, newdata, type = "pred", ...) {
  m_list <- object$m_list
  var_list <- object$var_list
  k_class <- object$k_class
  n_model <- object$n_model
  class_names <- object$class_names

  x <- newdata
  n <- nrow(x)

  M_votes <- matrix(data = 0, nrow = n, ncol = k_class)

  for (i in 1:n_model) {
    x_selected <- x[, var_list[[i]], drop = FALSE]
    prob <-
      predict(object = m_list[[i]],
              newdata = x_selected,
              type = "prob")

    for (j in 1:n) {
      M_votes[j, which.max(prob[j,])] <-
        M_votes[j, which.max(prob[j,])] + 1
    }
  }

  prob <- t(apply(M_votes, 1, function(m)
    m / sum(m)))

  if (type == "prob") {
    colnames(prob) <- class_names
    return(prob)
  }
  if (type == "pred") {
    pred <-
      factor(class_names[max.col(prob)], levels = class_names, labels = class_names)
    return(pred)
  }
}










