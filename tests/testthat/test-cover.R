make_data <- function(seed, n = 120, p = 2, round_digits = NULL) {
  set.seed(seed)
  y <- factor(sample(c("a", "b", "c"), n, TRUE, prob = c(.6, .3, .1)))
  x <- matrix(rnorm(n * p), n)
  x[, 1] <- x[, 1] + as.integer(y)
  if (!is.null(round_digits)) x <- round(x, round_digits)
  list(x = x, y = y)
}

# Plain R reference of the pure and proper cover (greedy, static counts)
pcccd_ref <- function(x_main, radii, thresh) {
  n <- nrow(x_main)
  D <- as.matrix(stats::dist(x_main))
  M <- D < matrix(radii, n, n)           # row i: ball i
  n_covers <- rowSums(M)
  covered <- rep(FALSE, n)
  dom <- integer(0)
  while (sum(covered) < thresh) {
    n_covers[covered] <- -1
    i <- which.max(n_covers)
    dom <- c(dom, unname(i))
    covered[M[i, ]] <- TRUE
    covered[i] <- TRUE
  }
  list(i_dominant = dom, proportion = sum(covered) / n)
}

test_that("pcccd cover matches the R reference", {
  for (seed in 1:6) {
    d <- make_data(seed, round_digits = if (seed %% 2 == 0) 1 else NULL)
    for (prop in c(1, 0.8)) {
      m <- pcccd_classifier(d$x, d$y, proportion = prop)
      for (i in seq_len(m$k_class)) {
        i_main <- which(d$y == m$class_names[i])
        x_main <- d$x[i_main, , drop = FALSE]
        radii <- c(RANN::nn2(d$x[-i_main, , drop = FALSE], x_main, k = 1)$nn.dist)
        ref <- pcccd_ref(x_main, radii, length(i_main) * prop)
        expect_equal(m$i_dominant_list[[i]], ref$i_dominant)
        expect_equal(m$proportions[i], ref$proportion)
        expect_equal(m$radii_dominant_list[[i]], radii[ref$i_dominant])
      }
    }
  }
})

test_that("pcccd is pure and proper", {
  d <- make_data(10)
  m <- pcccd_classifier(d$x, d$y)
  expect_true(all(m$proportions == 1))
  for (i in seq_len(m$k_class)) {
    other <- d$x[d$y != m$class_names[i], , drop = FALSE]
    dist_to_other <- rcccd:::cross_dist_cpp(m$x_dominant_list[[i]], other)
    expect_true(all(dist_to_other >= m$radii_dominant_list[[i]]))
  }
})

test_that("rwcccd C++ matches the R implementation", {
  for (seed in 1:6) {
    d <- make_data(seed, round_digits = if (seed %% 2 == 0) 1 else NULL)
    for (method in c("default", "balanced")) {
      for (m in c(0, 1)) {
        a <- rwcccd_classifier(d$x, d$y, method = method, m = m, proportion = 0.95)
        b <- rwcccd_classifier_2(d$x, d$y, method = method, m = m, proportion = 0.95)
        expect_equal(a$i_dominant_list, b$i_dominant_list)
        expect_equal(a$radii_dominant_list, b$radii_dominant_list)
        expect_equal(a$T_score_list, b$T_score_list)
        expect_equal(a$proportions, b$proportions)
      }
    }
  }
})

test_that("a single predictor and a data.frame work", {
  d <- make_data(3, p = 1)
  m <- pcccd_classifier(d$x, d$y)
  expect_s3_class(m, "pcccd_classifier")
  m2 <- rwcccd_classifier(as.data.frame(d$x), d$y)
  expect_length(predict(m2, as.data.frame(d$x)), nrow(d$x))
})
