#' @title  Random Walk Class Cover Catch Digraph Prediction
#'
#' @description \code{predict.rwcccd_classifier} makes prediction using
#' \code{rwcccd_classifier} object.
#'
#' @param object a \code{rwcccd_classifier} object
#' @param newdata newdata as matrix or dataframe.
#' @param type "pred" or "prob". Default is "pred". "pred" is class estimations,
#'  "prob" is \eqn{n\times k} matrix of class probabilities.
#' @param e 0 or 1. Default is 0. Penalty based on \eqn{T} scores in
#' \code{rwcccd_classifier} object.
#' @param ... not used.
#'
#' @details
#' Estimations are based on nearest dominant neighbor in radius unit.
#' \code{e} argument is used to penalize estimations based on \eqn{T} scores in
#' \code{rwcccd_classifier} object.
#'
#' For detail, please refer to Priebe et al. (2001), Priebe et al. (2003),
#' and Manukyan and Ceyhan (2016).
#'
#' @return a vector of class predictions (if type is "pred") or a \eqn{n\times p}
#' matrix of class probabilities (if type is "prob").
#'
#' @author Fatih Saglam, saglamf89@gmail.com
#'
#' @references
#' Priebe, C. E., DeVinney, J., & Marchette, D. J. (2001). On the distribution
#' of the domination number for random class cover catch digraphs. Statistics &
#' Probability Letters, 55(3), 239–246. https://doi.org/10.1016/s0167-7152(01)00129-8
#'
#' Priebe, C. E., Marchette, D. J., DeVinney, J., & Socolinsky, D. A. (2003).
#' Classification Using Class Cover Catch Digraphs. Journal of Classification,
#' 20(1), 3–23. https://doi.org/10.1007/s00357-003-0003-7
#'
#' Manukyan, A., & Ceyhan, E. (2016). Classification of imbalanced data with a
#' geometric digraph family. Journal of Machine Learning Research, 17(1),
#' 6504–6543. https://jmlr.org/papers/volume17/15-604/15-604.pdf
#'
#'
#' @examples
#' n <- 1000
#' x1 <- runif(n, 1, 10)
#' x2 <- runif(n, 1, 10)
#' x <- cbind(x1, x2)
#' y <- as.factor(ifelse(3 < x1 & x1 < 7 & 3 < x2 & x2 < 7, "A", "B"))
#'
#' # testing the performance
#' i_train <- sample(1:n, round(n*0.8))
#'
#' x_train <- x[i_train,]
#' y_train <- y[i_train]
#'
#' x_test <- x[-i_train,]
#' y_test <- y[-i_train]
#'
#' m_rwcccd <- rwcccd_classifier(x = x_train, y = y_train)
#' pred <- predict(object = m_rwcccd, newdata = x_test, e = 0)
#'
#' # confusion matrix
#' table(y_test, pred)
#'
#' # test accuracy
#' sum(y_test == pred)/nrow(x_test)
#'
#' @rdname predict.rwcccd_classifier
#' @export

predict.rwcccd_classifier <- function(object, newdata, type = "pred", e = 0, ...) {
  if (!is.numeric(e) || length(e) != 1 || is.na(e) || e < 0 || e > 1) {
    stop("e must be between 0 and 1")
  }
  if (!(type %in% c("pred", "prob"))) {
    stop("type must be 'pred' or 'prob'")
  }
  x <- check_newdata(newdata)

  expo_list <- lapply(object$T_score_list, function(t) t^e)
  dist_prop <- class_dist_prop(x, object$x_dominant_list, object$radii_dominant_list, expo_list)
  prob <- dist_to_prob(dist_prop, object$class_names)

  if (type == "prob") {
    return(prob)
  }
  prob_to_pred(prob, object$class_names)
}
