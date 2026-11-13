# Authors' exact environments (macOS arm64)

Record of the conda environments in which the figures were originally produced (one per figure folder; all have the same
core packages: python 3.13.16, numpy 2.5.1, pandas 3.0.5, matplotlib 3.11.1, scipy 1.18.0, seaborn 0.13.2,
openpyxl 3.1.5, scikit-learn 1.9.0 where used).

- `<env>.yml` = `conda env export --from-history`
- `<env>.lock.txt` = `conda list --export` (platform osx-arm64); recreate exactly with
  `conda create -n <name> --file <env>.lock.txt` on an Apple-silicon Mac.

For normal use the single top-level `environment.yml` (env `sc-folding-figures`) is sufficient.
