# Reaction Basis

Exact balancing and independent reaction sets for Python.

Build reactions from a species collection, balance specified sides, or extract
an independent subset of existing reactions. Element counts and charge are
conserved exactly; isomers and electronic states remain distinct by ID.

Python 3.11+ · MIT

[Documentation](https://ivan-pikulin.github.io/reaction-basis/) ·
[PyPI](https://pypi.org/project/reaction-basis/)

## Install

```sh
pip install reaction-basis
```

## Balance specified sides

```python
from reaction_basis import ReactionSpace, Species

space = ReactionSpace(
    [
        Species.from_formula("H2"),
        Species.from_formula("O2"),
        Species.from_formula("H2O"),
    ]
)
reaction = space.balance(reactants=["H2", "O2"], products=["H2O"])
print(reaction.format())  # 2 H2 + O2 -> 2 H2O
```

Ambiguous, impossible or incorrectly oriented balances raise an error.

## Generate an independent set

```python
formulas = ["H2", "O2", "H2O", "H2O2", "C", "CO", "CO2", "CH4", "CH3OH"]
space = ReactionSpace(Species.from_formula(f) for f in formulas)
reactions = space.basis()

assert len(reactions) == space.dimension == 6
for reaction in reactions:
    print(reaction.format())
```

`basis()` returns a full kernel basis: **N − rank(A)** independent reactions.
It does not optimize chemical simplicity or predict a mechanism.

## Choose the operation

| Method | Purpose |
|---|---|
| `balance(reactants=..., products=...)` | Uniquely balance specified sides |
| `basis()` | Generate a complete independent reaction basis |
| `circuits(max_species=4)` | Find support-minimal reactions within a size bound |
| `independent(reactions)` | Select an independent subset spanning the candidates |

Circuits can be dependent; support-limited candidates may not span the full
space. Enumeration has an explicit subset budget and fails rather than
silently truncating results.

## Keep species distinct

```python
ethanol = Species.from_formula("C2H6O", id="ethanol")
ether = Species.from_formula("C2H6O", id="dimethyl_ether")
ion = Species.from_formula("U", charge=4, multiplicity=3)
```

Default IDs encode composition, charge and optional multiplicity. Supply
explicit IDs for separate structures or source records; duplicate IDs fail.
The library does not identify structures from formulas.

## Use the results

`reaction.coefficients` is a read-only mapping from species IDs to exact SymPy
rationals: reactants negative, products positive. Use it with your own property tables.
`reaction.to_dict()` exports rational strings for lossless JSON interchange.

Coefficient normalization matters: doubling a reaction doubles its assigned
extensive quantities. The consumer owns units, provenance and uncertainty.

See the [guide](https://ivan-pikulin.github.io/reaction-basis/guide.html)
for algorithm guarantees and integration examples, or
[contributing](CONTRIBUTING.md) for development instructions.
