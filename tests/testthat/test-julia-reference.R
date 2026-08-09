test_that("pure-R evaluators reproduce frozen Julia reference fixtures", {
  input_path <- system.file(
    "extdata", "julia_reference_input.csv", package = "RobustSandwichRV"
  )
  expected_path <- system.file(
    "extdata", "julia_reference_expected.csv", package = "RobustSandwichRV"
  )
  expect_true(nzchar(input_path))
  expect_true(nzchar(expected_path))
  input <- read.csv(input_path)
  expected <- read.csv(expected_path)
  cache <- prepare_model(
    input$y, input$d, cbind(1, input$x1, input$x2)
  )
  for (row in seq_len(nrow(expected))) {
    vcov <- expected$vcov[[row]]
    cluster <- if (grepl("CR", vcov)) input$cluster else NULL
    fast <- evaluate_confounder(
      cache, input$z, vcov = vcov, cluster = cluster
    )
    reference <- evaluate_confounder_reference(
      cache, input$z, vcov = vcov, cluster = cluster
    )
    for (field in c("tau", "se_classical", "se_robust", "omega",
                    "t_classical", "t_robust", "r2_d", "r2_y")) {
      expect_equal(fast[[field]], expected[[field]][[row]],
                   tolerance = 2e-11,
                   info = paste(vcov, field))
    }
    expect_equal(reference$tau, expected$reference_tau[[row]],
                 tolerance = 2e-11)
    expect_equal(reference$se_robust,
                 expected$reference_se_robust[[row]], tolerance = 2e-11)
    expect_equal(reference$omega, expected$reference_omega[[row]],
                 tolerance = 2e-11)
  }
})
