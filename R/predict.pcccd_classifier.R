#' @title  Pure and Proper Class Cover Catch Digraph Prediction
#'
#' @description \code{predict.pcccd_classifier} makes prediction using \code{pcccd_classifier} object.
#'
#' @param object a \code{pcccd_classifier} object
#' @param newdata newdata as matrix or dataframe.
#' @param type "pred" or "prob". Default is "pred". "pred" is class estimations,
#'  "prob" is \eqn{n\times k} matrix of class probabilities.
#' @param ... not used.
#' @details
#' Estimations are based on nearest dominant neighbor in radius unit.
#'
#' For detail,
#' please refer to Priebe et al. (2001), Marchette and Socolinsky (2003), and
#' Manukyan and Ceyhan (2016).
#'
#' @return a vector of class predictions (if type is "pred") or a \eqn{n\times p}
#' matrix of class probabilities (if type is "prob").
#'
#' @author Fatih Saglam, saglamf89@gmail.com
#'
#' @references
#' Priebe, C. E., DeVinney, J. G., & Marchette, D. J. (2001). On the
#' distribution of the domination number for random class cover catch digraphs.
#' Statistics & Probability Letters, 55(3), 239-246.
#'
#' Marchette, C. E. P. D. J., & Socolinsky, J. G. D. D. A. (2003).
#' Classiﬁcation Using Class Cover Catch Digraphs. Journal of Classiﬁcation,
#' 20, 3-23.
#'
#' Manukyan, A., & Ceyhan, E. (2016). Classification of imbalanced
#' data with a geometric digraph family. The Journal of Machine Learning
#' Research, 17(1), 6504-6543.
#'
#' @importFrom proxy dist
#' @importFrom Rfast colMins
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
#' m_pcccd <- pcccd_classifier(x = x_train, y = y_train)
#' pred <- predict(object = m_pcccd, newdata = x_test)
#'
#' # confusion matrix
#' table(y_test, pred)
#'
#' # test accuracy
#' sum(y_test == pred)/nrow(x_test)
#'
#' @rdname predict.pcccd_classifier
#' @export

predict.pcccd_classifier <- function(object, newdata, type = "pred", ...) {
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










