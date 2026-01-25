# Formula Interface Tests for synthdid package
# Addresses SUGGESTIONS.md Issue #8: Missing formula interface tests

test_that("basic formula interface works", {
  data(california_prop99)

  expect_no_error({
    result <- synthdid(PacksPerCapita ~ treated,
                       data = california_prop99,
                       index = c("State", "Year"))
  })

  # Should return synthdid class
  expect_s3_class(result, "synthdid")
  expect_s3_class(result, "synthdid_estimate")

  # Should have required attributes
  expect_true(!is.null(attr(result, "call")))
  expect_true(!is.null(attr(result, "formula")))
  expect_true(!is.null(attr(result, "index")))
  expect_true(!is.null(attr(result, "data_info")))
})

test_that("formula interface produces finite estimate", {
  data(california_prop99)

  result <- synthdid(PacksPerCapita ~ treated,
                     data = california_prop99,
                     index = c("State", "Year"))

  expect_true(is.finite(coef(result)))
})

test_that("formula interface print method works", {
  data(california_prop99)

  result <- synthdid(PacksPerCapita ~ treated,
                     data = california_prop99,
                     index = c("State", "Year"))

  expect_output(print(result), "Synthetic Difference-in-Differences Estimate")
  expect_output(print(result), "Treatment Effect")
  expect_output(print(result), "Units:")
  expect_output(print(result), "Time Periods:")
})

test_that("formula interface summary method works", {
  data(california_prop99)

  result <- synthdid(PacksPerCapita ~ treated,
                     data = california_prop99,
                     index = c("State", "Year"))

  expect_output(summary(result), "Treatment Effect Estimate")
  expect_output(summary(result), "Dimensions")
})

test_that("formula interface coef method works", {
  data(california_prop99)

  result <- synthdid(PacksPerCapita ~ treated,
                     data = california_prop99,
                     index = c("State", "Year"))

  coef_val <- coef(result)
  expect_true(is.numeric(coef_val))
  expect_equal(length(coef_val), 1)
  expect_equal(names(coef_val), "treated")
})

test_that("formula interface with different methods", {
  data(california_prop99)

  # synthdid
  result_sdid <- synthdid(PacksPerCapita ~ treated,
                          data = california_prop99,
                          index = c("State", "Year"),
                          method = "synthdid")

  # SC
  result_sc <- synthdid(PacksPerCapita ~ treated,
                        data = california_prop99,
                        index = c("State", "Year"),
                        method = "sc")

  # DID
  result_did <- synthdid(PacksPerCapita ~ treated,
                         data = california_prop99,
                         index = c("State", "Year"),
                         method = "did")

  # All should produce estimates
  expect_true(is.finite(coef(result_sdid)))
  expect_true(is.finite(coef(result_sc)))
  expect_true(is.finite(coef(result_did)))

  # Estimates should differ (different methods)
  expect_false(isTRUE(all.equal(coef(result_sdid), coef(result_sc))))
  expect_false(isTRUE(all.equal(coef(result_sdid), coef(result_did))))
})

test_that("formula interface with standard errors", {
  data(california_prop99)

  # This may be slow, so we use few replications
  result <- synthdid(PacksPerCapita ~ treated,
                     data = california_prop99,
                     index = c("State", "Year"),
                     se = TRUE,
                     se_method = "jackknife")

  # Should have SE attribute
  expect_true(!is.null(attr(result, "se")))
  expect_true(is.finite(attr(result, "se")))
})

test_that("formula interface confint method works", {
  data(california_prop99)

  result <- synthdid(PacksPerCapita ~ treated,
                     data = california_prop99,
                     index = c("State", "Year"),
                     se = TRUE,
                     se_method = "jackknife")

  ci <- confint(result)
  expect_true(is.matrix(ci))
  expect_equal(nrow(ci), 1)
  expect_equal(ncol(ci), 2)
  expect_equal(colnames(ci), c("Lower", "Upper"))
  expect_equal(rownames(ci), "treated")
})

test_that("formula interface without index uses first two columns", {
  data(california_prop99)

  # Reorder so unit and time are first
  data_reordered <- california_prop99[, c("State", "Year", "PacksPerCapita", "treated")]

  expect_message({
    result <- synthdid(PacksPerCapita ~ treated,
                       data = data_reordered)
  }, "Using first two columns")

  expect_true(is.finite(coef(result)))
})

test_that("formula interface errors on missing variables", {
  data(california_prop99)

  expect_error({
    synthdid(NonExistent ~ treated,
             data = california_prop99,
             index = c("State", "Year"))
  }, "not found in data")
})

test_that("formula interface errors on bad index", {
  data(california_prop99)

  expect_error({
    synthdid(PacksPerCapita ~ treated,
             data = california_prop99,
             index = "BadIndex")
  }, "must be a character vector of length 2")
})

test_that("formula interface errors on non-data.frame", {
  expect_error({
    synthdid(PacksPerCapita ~ treated,
             data = matrix(1:10, 2, 5),
             index = c("unit", "time"))
  }, "'data' must be a data.frame")
})

test_that("formula interface errors on non-formula", {
  data(california_prop99)

  expect_error({
    synthdid("not a formula",
             data = california_prop99,
             index = c("State", "Year"))
  }, "'formula' must be a formula")
})

test_that("formula interface update method works", {
  data(california_prop99)

  result_original <- synthdid(PacksPerCapita ~ treated,
                              data = california_prop99,
                              index = c("State", "Year"),
                              method = "synthdid")

  result_updated <- update(result_original, method = "did")

  # Should have different estimates
  expect_false(isTRUE(all.equal(coef(result_original), coef(result_updated))))

  # Updated result should use DID
  expect_equal(attr(result_updated, "estimator"), "did")
})

test_that("formula interface predict method works", {
  data(california_prop99)

  result <- synthdid(PacksPerCapita ~ treated,
                     data = california_prop99,
                     index = c("State", "Year"))

  # Counterfactual prediction
  cf <- predict(result, type = "counterfactual")
  expect_true(is.matrix(cf))

  # Treated outcomes
  treated_outcomes <- predict(result, type = "treated")
  expect_true(is.matrix(treated_outcomes))

  # Effect curve
  effect_curve <- predict(result, type = "effect")
  expect_true(is.numeric(effect_curve))
})

test_that("formula interface fitted method works", {
  data(california_prop99)

  result <- synthdid(PacksPerCapita ~ treated,
                     data = california_prop99,
                     index = c("State", "Year"))

  fitted_vals <- fitted(result)
  expect_true(is.matrix(fitted_vals))

  # Should have same dimensions as Y
  setup <- attr(result, "setup")
  expect_equal(dim(fitted_vals), dim(setup$Y))
})

test_that("formula interface residuals method works", {
  data(california_prop99)

  result <- synthdid(PacksPerCapita ~ treated,
                     data = california_prop99,
                     index = c("State", "Year"))

  # Control residuals
  resid_control <- residuals(result, type = "control")
  expect_true(is.matrix(resid_control))

  # Pretreatment residuals
  resid_pre <- residuals(result, type = "pretreatment")
  expect_true(is.matrix(resid_pre))

  # All residuals
  resid_all <- residuals(result, type = "all")
  expect_true(is.matrix(resid_all))
})

test_that("formula interface with small dataset", {
  # Create minimal panel data
  panel_data <- data.frame(
    unit = rep(c("A", "B", "C"), each = 4),
    time = rep(1:4, times = 3),
    outcome = rnorm(12, mean = 10),
    treatment = c(rep(0, 8), rep(0, 2), rep(1, 2))  # Unit C treated in periods 3-4
  )

  expect_no_error({
    result <- synthdid(outcome ~ treatment,
                       data = panel_data,
                       index = c("unit", "time"))
  })

  expect_true(is.finite(coef(result)))
})

test_that("formula interface preserves convergence information", {
  data(california_prop99)

  result <- synthdid(PacksPerCapita ~ treated,
                     data = california_prop99,
                     index = c("State", "Year"))

  # Should have convergence attribute
  expect_true(!is.null(attr(result, "convergence")))

  # Convergence info should be accessible
  conv_info <- synthdid_convergence_info(result)
  expect_true(!is.null(conv_info))
})

test_that("formula interface model.frame method works", {
  data(california_prop99)

  result <- synthdid(PacksPerCapita ~ treated,
                     data = california_prop99,
                     index = c("State", "Year"))

  mf <- model.frame(result)
  expect_true(is.list(mf))
  expect_true(!is.null(mf$call))
  expect_true(!is.null(mf$terms))
})
