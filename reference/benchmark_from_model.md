# Calibrate an observed scalar control or control group

The observed group is deleted from the control matrix, its numerical
rank is measured, and the score-alignment change is divided by that rank
before it is transported to the one-index sensitivity problem.

## Usage

``` r
benchmark_from_model(model, treatment, benchmark, vcov = "HC1", cluster = NULL)
```

## Arguments

- model:

  An unweighted `lm` model.

- treatment:

  Name of the treatment model-matrix column.

- benchmark:

  One or more model term or model-matrix column names.

- vcov:

  HC0, HC1, CR0, or CR1.

- cluster:

  Optional cluster formula or labels.

## Value

An `rv_benchmark` with full calibration metadata.
