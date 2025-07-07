#' @title Pure and Proper Class Cover Catch Digraph Classifier
#'
#' @description \code{pcccd_classifier} fits a Pure and Proper Class Cover Catch
#' Digraph (PCCCD) classification model.
#'
#' @param x feature matrix or dataframe.
#' @param y class factor variable.
#' @param proportion proportion of covered samples. A real number between \eqn{(0,1]}.
#' 1 by default. Smaller numbers results in less dominant samples.
#'
#' @details
#' Multiclass framework for PCCCD. PCCCD determines target class dominant points
#' set \eqn{S} and their circular cover area by determining balls
#' \eqn{B(x^{\text{target}}, r_{i})} with radii r using minimum amount of
#' dominant point which satisfies \eqn{X^{\text{non-target}}\cap \bigcup_{i}
#' B_{i} = \varnothing} (pure) and \eqn{X^{\text{target}}\subset \bigcup_{i}
#' B_{i}} (proper).
#'
#' This guarantees that balls of target class never covers any non-target
#' samples (pure) and balls cover all target samples (proper).
#'
#' For detail, please refer to Priebe et al. (2001), Priebe et al. (2003),
#' and Manukyan and Ceyhan (2016).
#'
#' Note: Much faster than \code{cccd} package.
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
#' @importFrom  RANN nn2
#' @importFrom  Rfast Dist
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
#' @examples
#' n <- 1000
#' x1 <- runif(n, 1, 10)
#' x2 <- runif(n, 1, 10)
#' x <- cbind(x1, x2)
#' y <- as.factor(ifelse(3 < x1 & x1 < 7 & 3 < x2 & x2 < 7, "A", "B"))
#'
#' m_pcccd <- pcccd_classifier(x = x, y = y)
#'
#' # dataset
#' plot(x, col = y, asp = 1)
#'
#' # dominant samples of first class
#' x_center <- m_pcccd$x_dominant_list[[1]]
#'
#' # radii of balls for first class
#' radii <- m_pcccd$radii_dominant_list[[1]]
#'
#' # balls
#' for (i in 1:nrow(x_center)) {
#' xx <- x_center[i, 1]
#' yy <- x_center[i, 2]
#' r <- radii[i]
#' theta <- seq(0, 2*pi, length.out = 100)
#' xx <- xx + r*cos(theta)
#' yy <- yy + r*sin(theta)
#' lines(xx, yy, type = "l", col = "green")
#' }
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
#' @rdname pcccd_classifier
#' @export

pcccd_classifier <- function(x, y, proportion = 1) {

  if (proportion < 0 | proportion > 1) {
    stop("proportion must be in range [0,1]")
  }

  if (!is.matrix(x) & !is.data.frame(x)) {
    stop("x must be a matrix or data.frame")
  }

  if (is.data.frame(x)) {
    x <- as.matrix(x)
  }

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

    dist_main2other <- RANN::nn2(data = x_other, query = x_main, k = 1)$nn.dist
    dist_main2main <- Rfast::Dist(x = x_main)
    M <- dist_main2main < c(dist_main2other)
    M <- matrix(as.numeric(M), n_main)

    thresh <- n_main*proportion

    m_dominant <- f_cover_pcccd(thresh = thresh, M = M)

    i_dominant_list[[i]] <- m_dominant$i_dominant
    x_dominant_list[[i]] <- x_main[m_dominant$i_dominant,,drop = FALSE]
    radii_dominant_list[[i]] <- dist_main2other[m_dominant$i_dominant,]
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
  class(results) <- "pcccd_classifier"
  return(results)
}








