"""External consumer example: energetics is NOT part of the library API.

Run from the project: uv run python examples/external_energy.py
Energies are illustrative, all in the same unit and reference channel.
"""

from math import fsum, isfinite

from reaction_basis import ReactionSpace, Species


def main() -> float:
    space = ReactionSpace([Species.from_formula(f) for f in ("H2", "O2", "H2O")])
    reaction = space.balance(reactants=["H2", "O2"], products=["H2O"])
    energies = {"H2": -1.0, "O2": -2.0, "H2O": -4.0}

    # The consumer owns units, reference methods, convergence and uncertainty.
    # Missing species raise KeyError; there is no zero-fill or ignored record.
    values = [energies[key] for key in reaction.coefficients]
    if not all(isfinite(value) for value in values):
        raise ValueError("nonfinite species energy")
    delta_energy = fsum(
        float(coefficient) * energies[key] for key, coefficient in reaction.coefficients.items()
    )
    print(reaction.format())
    print(f"ΔE = {delta_energy:g} (same units as energies)")
    return delta_energy


if __name__ == "__main__":
    main()
