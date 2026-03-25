# Contributing to PolyVAR

Bug reports, feature requests, and pull requests are welcome.

## Reporting bugs
Open an issue at https://github.com/Deylab999/PolyVAR/issues
Include: R version, OS, package version, minimal reproducible example.

## Pull requests
1. Fork the repo
2. Create a branch: `git checkout -b my-fix`
3. Make changes and add tests
4. Run `R CMD check .` -- no errors or warnings
5. Submit a pull request

## Code style
- Follow existing Rcpp/R conventions in the package
- Comment all C++ functions with input/output types
- Keep per-SNP functions self-contained (no R callbacks from C++)
