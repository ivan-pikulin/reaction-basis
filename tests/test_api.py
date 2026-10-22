import inspect
import runpy
from pathlib import Path

import pytest

from reaction_basis import Reaction, ReactionSpace, Species


def test_single_required_argument_and_keyword_only_options():
    for factory in (Species, Species.from_formula):
        params = list(inspect.signature(factory).parameters.values())
        required = [p for p in params if p.default is inspect.Parameter.empty]
        assert len(required) == 1
        assert all(p.kind is inspect.Parameter.KEYWORD_ONLY for p in params[1:])
    with pytest.raises(TypeError):
        Species.from_formula("ethanol", "C2H6O")
    with pytest.raises(TypeError):
        Species({"H": 1}, "hydrogen")


def test_default_ids_are_stable_but_not_implicit_isomer_identity():
    assert Species.from_formula("OH2").id == Species({"H": 2, "O": 1}).id == "H2O"
    states = [Species.from_formula("U", charge=q) for q in (4, 5, 6)]
    assert len({s.id for s in states}) == 3
    assert (
        Species.from_formula("O2", multiplicity=1).id
        != Species.from_formula("O2", multiplicity=3).id
    )
    with pytest.raises(ValueError, match="duplicate"):
        ReactionSpace([Species.from_formula("OH2"), Species.from_formula("H2O")])
    space = ReactionSpace(
        [
            Species.from_formula("C2H6O", id="ethanol"),
            Species.from_formula("C2H6O", id="ether"),
        ]
    )
    assert space.dimension == 1
    assert Species.from_formula("H2O", label="water").label == "water"


def test_no_energy_in_core():
    assert not hasattr(Reaction, "energy")


def test_invalid_shapes_fail_explicitly():
    with pytest.raises(TypeError, match="composition must be a mapping"):
        Species([("H", 2)])
    with pytest.raises(TypeError, match="coefficients must be a mapping"):
        Reaction([("H2", -1), ("H", 2)])
    with pytest.raises(TypeError, match="Reaction object"):
        ReactionSpace([]).validate({"H2": -1, "H": 2})
    with pytest.raises(ValueError, match="label"):
        Species.from_formula("H2", label=" ")
    with pytest.raises(TypeError, match="canonical"):
        Reaction({"a": -1, "b": 1}).primitive(canonical="false")


def test_space_catalog_is_read_only():
    space = ReactionSpace([Species.from_formula("H2")])
    for name in ("ids", "species", "elements", "constraints"):
        with pytest.raises(AttributeError):
            setattr(space, name, ())


@pytest.mark.parametrize("formula", ["13C", "2H2O", "13CH4"])
def test_ambiguous_prefix_is_not_an_isotope_parser(formula):
    with pytest.raises(ValueError, match="prefixes"):
        Species.from_formula(formula)


def test_balance_preserves_requested_direction():
    space = ReactionSpace([Species.from_formula(f) for f in ("H2", "O2", "H2O")])
    reaction = space.balance(reactants=["H2", "O2"], products=["H2O"])
    assert reaction.to_dict() == {"H2": "-2", "H2O": "2", "O2": "-1"}
    assert (
        space.balance(reactants=["H2O"], products=["H2", "O2"]).to_dict()
        == reaction.reversed().to_dict()
    )
    with pytest.raises(ValueError, match="incompatible"):
        space.balance(reactants=["H2O", "H2"], products=["O2"])
    with pytest.raises(ValueError, match="unknown"):
        space.balance(reactants=["unknown"], products=["H2"])
    with pytest.raises(ValueError, match="distinct"):
        space.balance(reactants=["H2"], products=["H2"])
    with pytest.raises(TypeError):
        space.balance(reactants="H2", products=["H2O"])


def test_balance_impossible_ambiguous_and_zero_participant():
    space = ReactionSpace([Species.from_formula(f) for f in ("C", "O2", "CO", "CO2", "H2")])
    with pytest.raises(ValueError, match="ambiguous"):
        space.balance(reactants=["C", "O2"], products=["CO", "CO2"])
    with pytest.raises(ValueError, match="no balanced"):
        space.balance(reactants=["C"], products=["O2"])
    with pytest.raises(ValueError, match="not all"):
        space.balance(reactants=["C", "O2", "H2"], products=["CO2"])


def test_external_example():
    path = Path(__file__).parents[1] / "examples" / "external_energy.py"
    example = runpy.run_path(str(path))
    assert example["main"]() == -4.0
