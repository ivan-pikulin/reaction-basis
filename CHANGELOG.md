# Changelog

## 0.1.0 — 2026-10-22

- Renamed distribution to `reaction-basis` and import package to `reaction_basis`.
- Made catalog ordering properties read-only and tightened input validation.
- Rejected leading formula coefficients and ambiguous isotope-like prefixes.

- Exact element/charge balancing and independent reaction spaces.
- Distinct species IDs for isomers, geometries and electronic states.
- Fixed-side balancing, kernel bases, bounded circuits and independent subsets.
- Rational JSON interchange; energetics remains external.
- Strict Sphinx builds, executable examples and randomized tests.

Pre-release API changes: `Species.from_formula` takes formula first with
keyword-only options. `Reaction.energy` was removed before publication.
