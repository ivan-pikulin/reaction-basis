# Contributing

Keep the runtime package focused on exact stoichiometry. Energies, datasets,
quantum chemistry execution and unit systems belong to external consumers.

## Checks

```sh
uv sync --locked --group docs
uv run pytest
uv run ruff check .
uv run ruff format --check .
uv run sphinx-build -n -W --keep-going -b html docs docs/_build/html
uv run sphinx-build -n -W --keep-going -b doctest docs docs/_build/doctest
uv build
```

The test workflow covers Python 3.11–3.13. Documentation CI checks strict HTML
and doctests and deploys GitHub Pages from main after a successful build.
README examples run in pytest. Derive the documentation version from installed
metadata. Commit source and `uv.lock`, not generated HTML or doctrees.

## Preview

```sh
uv run python -m http.server 8765 --bind 127.0.0.1 --directory docs/_build/html
```

Open http://127.0.0.1:8765/. Rebuild after edits; this is a static server.

## Releases

Update the version and changelog, then push a matching `vX.Y.Z` tag. The release
workflow runs tests, strict documentation builds and an installed-wheel smoke
test before publishing the same distributions with PyPI Trusted Publishing.
It creates a GitHub Release after upload. To recover a failed upload, inspect
the PyPI files before rerunning; released versions cannot be overwritten.
