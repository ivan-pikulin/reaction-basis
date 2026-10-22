"""Run directly against an installed distribution, without development dependencies."""

from fractions import Fraction
from importlib.metadata import version

from reaction_basis import Reaction, ReactionSpace, Species


def main():
    assert version("reaction-basis")
    space = ReactionSpace(Species.from_formula(f) for f in ("H2", "O2", "H2O"))
    reaction = space.balance(reactants=["H2", "O2"], products=["H2O"])
    assert reaction.format() == "2 H2 + O2 -> 2 H2O"
    assert len(space.basis()) == space.dimension == 1
    assert Reaction({"a": -1, "b": Fraction(1, 2)}).primitive().to_dict() == {"a": "-2", "b": "1"}
    assert len(space.independent(space.circuits())) == 1
    space.validate(reaction)
    print(f"reaction-basis {version('reaction-basis')}: installed-distribution smoke test passed")


if __name__ == "__main__":
    main()
