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
  if (!(type %in% c("pred", "prob"))) {
    stop("type must be 'pred' or 'prob'")
  }
  x <- check_newdata(newdata)

  n <- nrow(x)
  class_names <- object$class_names
  M_votes <- matrix(0, nrow = n, ncol = object$k_class)
  i_row <- seq_len(n)

  for (i in seq_len(object$n_model)) {
    prob <- predict(object = object$m_list[[i]],
                    newdata = x[, object$var_list[[i]], drop = FALSE],
                    type = "prob")
    i_vote <- cbind(i_row, max.col(prob, ties.method = "first"))
    M_votes[i_vote] <- M_votes[i_vote] + 1
  }

  prob <- M_votes/rowSums(M_votes)
  colnames(prob) <- class_names

  if (type == "prob") {
    return(prob)
  }
  prob_to_pred(prob, class_names)
}
