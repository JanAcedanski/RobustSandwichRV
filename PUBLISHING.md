# Publishing on GitHub

Repository: <https://github.com/JanAcedanski/RVwithSandwichSE>

1.  Update `DESCRIPTION`, `NEWS.md`, and the version-pinned installation
    command in both `README.md` and `README.Rmd`. Keep their main
    content identical.

2.  Run the quick-start demo and the package checks:

    ``` sh
    Rscript tools/check-local.R .
    R CMD build .
    R CMD check --no-manual RobustSandwichRV_0.5.0.tar.gz
    ```

3.  Commit the reviewed source and push `main`. Inspect the R-CMD-check
    and pkgdown workflows. The check matrix covers Windows, macOS, and
    Linux.

4.  Create tag `v0.5.0` at the tested source commit and publish a GitHub
    release with release notes, `RobustSandwichRV_0.5.0.tar.gz`, and
    `SHA256SUMS`.

5.  Verify installation from the tag in a clean R library and run
    `demo("quick-start", package = "RobustSandwichRV", ask = FALSE)`.

The pkgdown workflow deploys the documentation to the `gh-pages` branch.
The source archive is an R source package, not a platform-specific
binary. For the next release, replace `0.5.0` in these steps with the
new version.
