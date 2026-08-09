# Package index

## Main workflow

- [`robust_rv()`](https://JanAcedanski.github.io/RobustSandwichRV/reference/robust_rv.md)
  [`observed_control_rv()`](https://JanAcedanski.github.io/RobustSandwichRV/reference/robust_rv.md)
  : Main model-facing robustness-value workflow
- [`calibrated_rv()`](https://JanAcedanski.github.io/RobustSandwichRV/reference/calibrated_rv.md)
  : Calibrated robustness value
- [`rv_cal_bounded()`](https://JanAcedanski.github.io/RobustSandwichRV/reference/rv_cal_bounded.md)
  : Benchmark-bounded calibrated robustness value
- [`rv_cal_bounded_sequence()`](https://JanAcedanski.github.io/RobustSandwichRV/reference/rv_cal_bounded_sequence.md)
  : Evaluate nested benchmark bands with witness reuse
- [`conditional_omega_rv()`](https://JanAcedanski.github.io/RobustSandwichRV/reference/conditional_omega_rv.md)
  [`rv_att_point()`](https://JanAcedanski.github.io/RobustSandwichRV/reference/conditional_omega_rv.md)
  : Point-conditioned attainable robustness value

## Model preparation and evaluation

- [`prepare_model()`](https://JanAcedanski.github.io/RobustSandwichRV/reference/prepare_model.md)
  : Prepare a linear-regression sensitivity cache

- [`cache_from_lm()`](https://JanAcedanski.github.io/RobustSandwichRV/reference/cache_from_lm.md)
  :

  Create a sensitivity cache from an `lm` model

- [`evaluate_confounder()`](https://JanAcedanski.github.io/RobustSandwichRV/reference/evaluate_confounder.md)
  : Fast FWL evaluation of an omitted confounder

- [`evaluate_confounder_reference()`](https://JanAcedanski.github.io/RobustSandwichRV/reference/evaluate_confounder_reference.md)
  : Independent explicit full-regression evaluation

- [`evaluate_subspace()`](https://JanAcedanski.github.io/RobustSandwichRV/reference/evaluate_subspace.md)
  : Evaluate a residualized confounder subspace

- [`validate_sandwich_reference()`](https://JanAcedanski.github.io/RobustSandwichRV/reference/validate_sandwich_reference.md)
  :

  Validate a short-model covariance against `sandwich`

## Calibration and Cinelli-Hazlett geometry

- [`benchmark_calibration()`](https://JanAcedanski.github.io/RobustSandwichRV/reference/benchmark_calibration.md)
  : Construct a rank-normalized benchmark calibration
- [`benchmark_from_model()`](https://JanAcedanski.github.io/RobustSandwichRV/reference/benchmark_from_model.md)
  : Calibrate an observed scalar control or control group
- [`bounded_benchmark_band()`](https://JanAcedanski.github.io/RobustSandwichRV/reference/bounded_benchmark_band.md)
  : Construct a benchmark sandwich-ratio band
- [`alignment_directions()`](https://JanAcedanski.github.io/RobustSandwichRV/reference/alignment_directions.md)
  : Exact normalized score directions
- [`single_parameterization()`](https://JanAcedanski.github.io/RobustSandwichRV/reference/single_parameterization.md)
  : Construct a scalar confounder from partial correlations
- [`classical_ch_t()`](https://JanAcedanski.github.io/RobustSandwichRV/reference/classical_ch_t.md)
  : Conventional full-model t statistic from CH geometry
- [`conditional_omega_evaluator()`](https://JanAcedanski.github.io/RobustSandwichRV/reference/conditional_omega_evaluator.md)
  : Construct the estimator-specific attainable omega evaluator
- [`conditional_equal_strength_rv()`](https://JanAcedanski.github.io/RobustSandwichRV/reference/conditional_equal_strength_rv.md)
  : Equal-strength conditional robustness value
- [`conditional_protected_rv()`](https://JanAcedanski.github.io/RobustSandwichRV/reference/conditional_protected_rv.md)
  : Analytical protected conditional robustness value

## Conventional benchmarks

- [`point_estimate_rv()`](https://JanAcedanski.github.io/RobustSandwichRV/reference/point_estimate_rv.md)
  : Point-estimate robustness value

- [`conventional_rv()`](https://JanAcedanski.github.io/RobustSandwichRV/reference/conventional_rv.md)
  : Conventional Cinelli-Hazlett robustness value

- [`naive_robust_t_rv()`](https://JanAcedanski.github.io/RobustSandwichRV/reference/naive_robust_t_rv.md)
  : Incorrect robust-t plug-in benchmark

- [`compare_sensemakr()`](https://JanAcedanski.github.io/RobustSandwichRV/reference/compare_sensemakr.md)
  :

  Compare the classical implementation with `sensemakr`

## Multiple-confounder geometry

- [`reduce_score_subspace()`](https://JanAcedanski.github.io/RobustSandwichRV/reference/reduce_score_subspace.md)
  : Reduce an omitted-confounder block to at most two active directions
- [`verify_score_reduction()`](https://JanAcedanski.github.io/RobustSandwichRV/reference/verify_score_reduction.md)
  : Verify the two-direction reduction numerically

## Result methods

- [`print(`*`<robust_rv_result>`*`)`](https://JanAcedanski.github.io/RobustSandwichRV/reference/print.robust_rv_result.md)
  : Print a robustness-value result
- [`summary(`*`<robust_rv_result>`*`)`](https://JanAcedanski.github.io/RobustSandwichRV/reference/summary.robust_rv_result.md)
  : Summarize a robustness-value result
- [`plot(`*`<robust_rv_result>`*`)`](https://JanAcedanski.github.io/RobustSandwichRV/reference/plot.robust_rv_result.md)
  : Plot the analytical decision geometry
- [`tidy_rv()`](https://JanAcedanski.github.io/RobustSandwichRV/reference/tidy_rv.md)
  : Tidy one robustness-value result
