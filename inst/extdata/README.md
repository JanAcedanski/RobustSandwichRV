# Frozen Julia reference fixtures

`julia_reference_input.csv` is a deterministic linear-regression design. The
expected fast and explicit-refit results in `julia_reference_expected.csv`
were generated once with the independently audited Julia
`RobustSandwichRV` implementation under Julia 1.12.4.

The R package never invokes Julia. These files are ordinary immutable test
fixtures used to detect numerical drift in coefficients, HC0/HC1/CR0/CR1
standard errors, sandwich ratios, t statistics, and partial R-squared values.

SHA-256:

- input: `2B2BD4AA2E9ECAE228244B338EEDDF3FFF1BA4FB65405ADB817650ABBB521A01`
- expected: `52ADEAFBEB6F69E294218BAC83AF62918694FA6181305160DE79062AA300E2E6`

The generation scripts are retained in `data-raw/` for provenance, but are
excluded from the built R package.

