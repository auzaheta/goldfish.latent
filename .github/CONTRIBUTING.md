# Contributing

Contributions to `goldfish.latent`, whether in the form of issue identification,
bug fixes, new code or documentation are encouraged and welcome:

* [Submit an issue](#issues)
* [Fix a bug or implement new features](#adding-new-code)
* [Document existing code](#documentation)

`goldfish.latent` extends [`goldfish`](https://github.com/stocnet/goldfish) with
Bayesian relational event models featuring latent variables (random effects and
Hidden Markov Models), using [Stan](https://mc-stan.org/) for inference via
`cmdstanr`.

Please note that this project is released with a
[Contributor Code of Conduct](../CODE_OF_CONDUCT.md). By contributing to this
project, you agree to abide by its terms.

## Issues

Please use the [issue tracker on GitHub](https://github.com/auzaheta/goldfish.latent/issues)
to identify problems or suggest new functionality, before submitting changes to
the code. We use issues to identify bugs and tasks, discuss feature requests, and
to track implementation of changes.

The most useful issues are ones that precisely identify a bug, or propose a test
that should pass but instead fails. When reporting a Stan-related problem, please
include your CmdStan version and a minimal reproducible data set.

## Adding new code

Independent or assigned code contributions are most welcome. When writing new R
code, please follow the [tidyverse style guide](https://style.tidyverse.org/).
The `lintr` and `goodpractice` packages help check adherence, and `styler` fixes
most formatting issues non-invasively (it also ships an RStudio Addin):

```r
# lint the package
lintr::lint_package(path = ".")

# good-practice checks (exclude the 80-character line-length check)
goodpractice::gp(path = ".", checks = all_checks()[-c(8)])

# fix styling for a single file
styler::style_file("filePath")
```

When editing Stan models in `inst/stan/`, keep the existing formatting (two-space
indentation, explicit blocks) and prefer the shared templates in
`inst/stan/templates/` over duplicating model code.

### Branches

We use two **main branches**:

1. `main` is reserved for fully functional releases. When `develop` reaches a
   stable point, a maintainer merges it into `main` and tags it with a release
   number.
2. `develop` reflects the latest development stage. Minor changes that enhance
   existing functionality can go here directly; new features that may break
   existing functionality belong on a supporting branch.

We use two types of **supporting branches**:

3. *Feature branches* develop new functionality. They branch off `develop` and
   are merged back into it (or deleted if abandoned).
4. *Hotfix branches* fix severe bugs in a release without pulling in unstable
   `develop` changes. Branch off `main` and prefix the name with `hotfix-`.

This branching model is based on
<https://nvie.com/posts/a-successful-git-branching-model/>.

#### Release (maintainer only)

```
git checkout main
git merge --no-ff develop
git tag -a vX.Y.Z -m "vX.Y.Z"
git push origin main --tags
```

#### Feature branch

```
git checkout -b myfeature develop
# ...make changes and commit...
git checkout develop
git merge --no-ff myfeature
git branch -d myfeature
git push origin develop
```

#### Hotfix branch

Increment the PATCH digit of the version: a hotfix for `v1.3.0` is named
`hotfix-v1.3.1` and released as `v1.3.1`.

```
git checkout -b hotfix-vX.Y.Z main
# ...fix the bug, bump DESCRIPTION + NEWS.md...
git checkout main
git merge --no-ff hotfix-vX.Y.Z
git tag -a vX.Y.Z -m "vX.Y.Z"
git push origin main --tags
git checkout develop
git merge --no-ff hotfix-vX.Y.Z
git branch -d hotfix-vX.Y.Z
```

### Commit messages

Commits that relate to an existing issue should mention the issue number
(preceded by `#`) in the message. When the issue reference is preceded by
`resolve`, `resolves`, `resolved`, `close`, `closes`, `closed`, `fix`, `fixes`,
or `fixed` (capitalised or not), GitHub updates the issue status automatically.

Mention the issue first, then describe what the change does; keep it to a single
line, with any ancillary changes after a comma:

```
fix #31 correct random-effects design matrix, update documentation
```

## Testing

We use the [testthat](https://testthat.r-lib.org/) package; tests live in
`tests/testthat/`. Verify that all tests pass before committing changes to
existing code:

```r
devtools::test()
```

Some tests fit Stan models and therefore require a working
[CmdStan](https://mc-stan.org/users/interfaces/cmdstan) installation
(`>= 2.29.2`) reachable by `cmdstanr`. When writing a new function, consider
adding a unit test for it. Conventions:

- A test file should cover one or more aspects of a single function.
- Name test files after the R file they cover, prefixed with `test-`, optionally
  postfixed with the function name (e.g. `test-write-json-methods.R`).
- Shared setup (sample networks, fixtures) belongs in a `helper-*.R` file; reuse
  existing test data rather than creating new data for every test.

## Documentation

Documentation is written with `roxygen2`; run `devtools::document()` after
changing any roxygen comment, export, or function signature. Another valuable
contribution is developing the vignettes/articles that illustrate what the
package adds. Please open an issue with proposals, or if existing documentation
is unclear.

## Versioning

The package follows [semantic versioning](https://semver.org/): versions use the
`Major.Minor.Patch` format. Bump the version in `DESCRIPTION` and add a `NEWS.md`
entry at each release milestone.

## CmdStan setup

Model fitting and several tests need CmdStan. The quickest route is via
`cmdstanr`:

```r
# install.packages("cmdstanr", repos = c("https://stan-dev.r-universe.dev"))
cmdstanr::check_cmdstan_toolchain()
cmdstanr::install_cmdstan()
```

If CmdStan is installed in a non-default location, point `cmdstanr` at it with
`cmdstanr::set_cmdstan_path("/path/to/cmdstan")`.
