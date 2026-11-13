# Internal helpers ------------------------------------------------------------

# Validates and standardizes the training data. Returns a list with a numeric
# matrix `x` and a factor `y`.
check_xy <- function(x, y) {
  if (!is.matrix(x) && !is.data.frame(x)) {
    stop("x must be a matrix or data.frame")
  }
  x <- as.matrix(x)
  if (!is.numeric(x)) {
    stop("x must be numeric")
  }
  if (anyNA(x)) {
    stop("x must not contain missing values")
  }
  if (!is.factor(y)) {
    stop("y must be a factor")
  }
  if (length(y) != nrow(x)) {
    stop("length of y must be equal to the number of rows of x")
  }
  if (anyNA(y)) {
    stop("y must not contain missing values")
  }
  if (sum(table(y) > 0) < 2) {
    stop("at least two classes must have observations")
  }
  list(x = x, y = y)
}

check_proportion <- function(proportion, name = "proportion") {
  if (!is.numeric(proportion) || length(proportion) != 1 || is.na(proportion) ||
      proportion < 0 || proportion > 1) {
    stop(name, " must be a single number in range [0,1]")
  }
}

check_newdata <- function(newdata) {
  if (!is.matrix(newdata) && !is.data.frame(newdata)) {
    stop("newdata must be a matrix or data.frame")
  }
  newdata <- as.matrix(newdata)
  if (!is.numeric(newdata)) {
    stop("newdata must be numeric")
  }
  newdata
}

# Scales class "distances" to class probabilities.
dist_to_prob <- function(dist_prop, class_names) {
  prob <- 1/(dist_prop + 1)
  prob <- prob/rowSums(prob)
  colnames(prob) <- class_names
  prob
}

prob_to_pred <- function(prob, class_names) {
  factor(class_names[max.col(prob, ties.method = "first")], levels = class_names)
}

# Distance (in radius unit) of each new sample to the nearest dominant ball of
# every class. Classes without any dominant sample get Inf.
class_dist_prop <- function(x, x_dominant_list, radii_dominant_list, expo_list = NULL) {
  k_class <- length(x_dominant_list)
  dist_prop <- matrix(Inf, nrow = nrow(x), ncol = k_class)
  for (i in seq_len(k_class)) {
    x_dom <- x_dominant_list[[i]]
    if (nrow(x_dom) == 0) next
    expo <- if (is.null(expo_list)) rep(1, nrow(x_dom)) else expo_list[[i]]
    dist_prop[, i] <- min_ratio_cpp(x_dom, x, radii_dominant_list[[i]], expo)
  }
  dist_prop
}
