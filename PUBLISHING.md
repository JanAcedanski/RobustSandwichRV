# Publishing on GitHub

The repository is ready to publish after choosing the GitHub account or
organization that will own it.

1. Configure repository-specific links:

   ```sh
   Rscript tools/configure-github.R GITHUB_OWNER
   ```

2. Run the local checks:

   ```sh
   R CMD build .
   R CMD check --no-manual RobustSandwichRV_0.1.0.tar.gz
   ```

3. Create an empty public repository named `RobustSandwichRV`, then push:

   ```sh
   git init
   git add .
   git commit -m "Initial pure-R RobustSandwichRV release"
   git branch -M main
   git remote add origin https://github.com/GITHUB_OWNER/RobustSandwichRV.git
   git push -u origin main
   ```

The included workflows run `R CMD check` on Windows, macOS, and Linux and
deploy the pkgdown site to the `gh-pages` branch. Enable GitHub Pages for that
branch after the first successful site build.
