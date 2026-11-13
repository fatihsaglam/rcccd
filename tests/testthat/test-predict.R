test_that("predict returns valid classes and probabilities", {
  set.seed(1)
  n <- 200
  x <- matrix(runif(n * 2, 1, 10), n)
  y <- factor(ifelse(3 < x[, 1] & x[, 1] < 7 & 3 < x[, 2] & x[, 2] < 7, "A", "B"))

  for (m in list(pcccd_classifier(x, y), rwcccd_classifier(x, y))) {
    prob <- predict(m, x, type = "prob")
    expect_equal(dim(prob), c(n, 2))
    expect_equal(colnames(prob), levels(y))
    expect_equal(unname(rowSums(prob)), rep(1, n))
    pred <- predict(m, x)
    expect_s3_class(pred, "factor")
    expect_equal(levels(pred), levels(y))
    expect_gt(mean(pred == y), 0.7)
  }
})

test_that("predict agrees with a direct computation", {
  set.seed(2)
  x <- matrix(rnorm(100 * 2), 100)
  y <- factor(rep(c("a", "b"), 50))
  x[y == "b", ] <- x[y == "b", ] + 2
  m <- pcccd_classifier(x, y)
  xt <- matrix(rnorm(20 * 2), 20) + 1

  dp <- sapply(1:2, function(i) {
    d <- rcccd:::cross_dist_cpp(m$x_dominant_list[[i]], xt)
    apply(d / m$radii_dominant_list[[i]], 2, min)
  })
  prob <- 1 / (dp + 1)
  prob <- prob / rowSums(prob)
  expect_equal(unname(predict(m, xt, type = "prob")), unname(prob))
})

test_that("ensemble classifier works", {
  set.seed(3)
  x <- matrix(rnorm(150 * 3), 150)
  y <- factor(sample(c("a", "b"), 150, TRUE, prob = c(.8, .2)))
  x[y == "b", 1] <- x[y == "b", 1] + 2
  m <- pcccd_ensemble_classifier(x, y, n_model = 7, n_var = 2)
  prob <- predict(m, x, type = "prob")
  expect_equal(unname(rowSums(prob)), rep(1, 150))
  expect_equal(levels(predict(m, x)), levels(y))
})

test_that("a class absent from the training data is handled", {
  set.seed(4)
  x <- matrix(rnorm(60 * 2), 60)
  y <- factor(rep(c("a", "b"), 30), levels = c("a", "b", "c"))
  m <- pcccd_classifier(x, y)
  expect_equal(unname(predict(m, x, type = "prob")[, "c"]), rep(0, 60))
})

test_that("input is validated", {
  x <- matrix(rnorm(40), 20)
  y <- factor(rep(c("a", "b"), 10))
  expect_error(pcccd_classifier(x, y, proportion = 2))
  expect_error(pcccd_classifier(x, as.character(y)), "factor")
  expect_error(pcccd_classifier(x, y[-1]), "length")
  expect_error(pcccd_classifier(x, factor(rep("a", 20))), "two classes")
  expect_error(rwcccd_classifier(x, y, method = "x"))
  expect_error(rwcccd_classifier(x, y, m = -1))
  m <- pcccd_classifier(x, y)
  expect_error(predict(m, x, type = "z"))
  expect_error(predict(m, "a"))
  expect_error(predict(rwcccd_classifier(x, y), x, e = 2))
})

test_that("legacy cover function keeps its old signature", {
  expect_true(is.function(rcccd:::f_cover_pcccd))
  expect_named(formals(rcccd:::f_cover_pcccd),
               c("cover", "thresh", "M", "dist_main2main", "dist_main2other"))
})
