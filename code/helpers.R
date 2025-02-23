#' Perform test for the comparison of two observed proportions
#'
#' @description
#' When the sample size is large enough (nAp, nAq, nBp, and nBq should be ≥5),
#' the function performs a two-proportions z-test. A Fisher Exact probability
#' test is performed when the two independent samples are small in size.
#'
#' @param test_data A <data.frame>-like object with at least three columns: a
#'    with the names of the two groups to be compared, a column with the total
#'    count in the corresponding groups, a column with the count of observed.
#' @param count_col A <character> with the name of the column that contains the
#'    total counts
#' @param observed_col A <character> with the name of the column that contains
#'    the observed counts
#' @param group_col A <character> with the name of the column that contains
#'    the names of the two groups to be compared
#'
#' @return A <tibble> of one row with the test statistics
#' @export
#'
#' @examples
#' data <- readRDS(
#'   file.path(getwd(), "data", "example_data", "test_prop_comparison.RDS")
#' )
#' test_result <- compare_proportions(
#'   test_data = data,
#'   group_col = "statut",
#'   count_col = "nbre_femelles,
#'   observed_col = "inseminees"
#' )
compare_proportions <- function(test_data, group_col, count_col, observed_col) {
  # calculate the overall proportion of inseminated (p) and non-inseminated (q)
  # females
  p <- sum(test_data[[observed_col]]) / sum(test_data[[count_col]])
  q <- 1 - p

  # use the two-proportions z-test if the sample size (n) is large enough.
  # nAp, nAq, nBp and nBq should be ≥ 5
  # if not, use the Fisher Exact probability test
  if (all(c(test_data[[count_col]] * p, test_data[[count_col]] * q) >= 5)) {
    test_result <- rstatix::prop_test(
      x = test_data[[observed_col]],
      n = test_data[[count_col]],
      p = NULL,
      alternative = "two.sided",
      correct = TRUE,
      conf.level = 0.95,
      detailed = TRUE
    )
  } else {
    xtab <- as.table(rbind(
      c(test_data[[observed_col]][1],
        test_data[[count_col]][1] - test_data[[observed_col]][1]),
      c(test_data[[observed_col]][2],
        test_data[[count_col]][2] - test_data[[observed_col]][2]))
    )
    dimnames(xtab) <- list(
      groups = test_data[[group_col]],
      observed = c("yes", "no")
    )
    test_result <- rstatix::fisher_test(xtab, detailed = TRUE)
    attr(test_result, "test") <- "fisher_exact"
  }

  return(test_result)
}
