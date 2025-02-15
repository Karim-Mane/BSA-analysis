
compare_proportions <- function(test_data) {
  # calculate the overall proportion of inseminated (p) and non-inseminated (q)
  # females
  p <- sum(test_data$inseminees) / sum(test_data$nbre_femelles)
  q <- 1 - p

  # use the two-proportions z-test if the sample size (n) is large enough.
  # nAp, nAq, nBp and nBq should be ≥ 5
  # if not, use the Fisher Exact probability test
  if (all(c(test_data$nbre_femelles * p, test_data$nbre_femelles * q) >= 5)) {
    test_result <- rstatix::prop_test(
      x = test_data$inseminees,
      n = test_data$nbre_femelles,
      p = NULL,
      alternative = "two.sided",
      correct = TRUE,
      conf.level = 0.95,
      detailed = TRUE
    )
  } else {
    xtab <- as.table(rbind(
      c(test_data$inseminees[1], test_data$nbre_femelles[1] - test_data$inseminees[1]),
      c(test_data$inseminees[2], test_data$nbre_femelles[2] - test_data$inseminees[2]))
    )
    dimnames(xtab) <- list(
      statut = c("ctrl", "test"),
      insinee = c("yes", "no")
    )
    test_result <- rstatix::fisher_test(xtab, detailed = TRUE)
    attr(test_result, "test") <- "fisher_exact"
  }

  return(test_result)
}
