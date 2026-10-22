import json
from fractions import Fraction

import pytest
from hypothesis import given
from hypothesis import strategies as st
from sympy import Matrix

from reaction_basis import Reaction, ReactionSpace, Species, parse_formula


def water():
    return ReactionSpace(
        [
            Species.from_formula("H2", id="h2"),
            Species.from_formula("O2", id="o2"),
            Species.from_formula("H2O", id="water"),
        ]
    )


def test_water():
    space = water()
    assert (space.rank, space.dimension) == (2, 1)
    (r,) = space.basis()
    assert r.to_dict() == {"h2": "-2", "o2": "-1", "water": "2"}
    assert Reaction(json.loads(json.dumps(r.to_dict()))).to_dict() == r.to_dict()
    assert len(space.independent([r, r.reversed(), r])) == 1


def test_duplicate_compositions_preserve_identities():
    space = ReactionSpace(
        [Species.from_formula("C2H6O", id=k) for k in ["ethanol", "ether", "state3"]]
    )
    assert space.dimension == 2
    basis = space.basis()
    assert len(space.independent(basis)) == 2
    assert set().union(*(r.coefficients.keys() for r in basis)) == set(space.ids)
    assert len(space.circuits()) == 3


def test_charge_and_electron():
    space = ReactionSpace(
        [
            Species.from_formula("U", id="u4", charge=4, multiplicity=3),
            Species.from_formula("U", id="u5", charge=5, multiplicity=2),
            Species.from_formula("U", id="u6", charge=6, multiplicity=1),
        ]
    )
    assert space.dimension == 1
    space.validate(Reaction({"u4": 1, "u5": -2, "u6": 1}))
    with pytest.raises(ValueError, match="charge"):
        space.validate(Reaction({"u4": -1, "u5": 1}))
    redox = ReactionSpace(
        [
            Species({}, id="e", charge=-1),
            Species.from_formula("H", id="h"),
            Species.from_formula("H", id="p", charge=1),
        ]
    )
    assert redox.dimension == 1
    redox.validate(redox.basis()[0])


@pytest.mark.parametrize(
    "formula,expected",
    [
        ("UO2(NO3)2", {"U": 1, "O": 8, "N": 2}),
        ("K4[Fe(CN)6]", {"K": 4, "Fe": 1, "C": 6, "N": 6}),
        ("CuSO4·5H2O", {"Cu": 1, "S": 1, "O": 9, "H": 10}),
    ],
)
def test_formulas(formula, expected):
    assert parse_formula(formula) == expected


@pytest.mark.parametrize(
    "formula",
    ["", "Xx2", "H0", "H2O+", "(H2]", "()", "H2O.", "H2O·", "[H", "2", "C1.5", "H 2", "H02"],
)
def test_bad_formulas(formula):
    with pytest.raises(ValueError):
        parse_formula(formula)


def test_exact_rational_and_orientation():
    r = Reaction({"a": Fraction(-1, 3), "b": "2/3"})
    assert r.primitive().to_dict() == {"a": "-1", "b": "2"}
    assert r.reversed().primitive().to_dict() == r.primitive().to_dict()
    assert r.reversed().primitive(canonical=False).to_dict() == {"a": "1", "b": "-2"}
    assert r.format() == "1/3 a -> 2/3 b"
    with pytest.raises(TypeError):
        Reaction({"a": -1.0, "b": 1})


def test_invalid_inputs():
    with pytest.raises(ValueError):
        ReactionSpace([Species.from_formula("H", id="same"), Species.from_formula("O", id="same")])
    with pytest.raises(ValueError):
        Species({"H": 0}, id="bad")
    with pytest.raises(TypeError):
        Species({"H": True}, id="bad")
    with pytest.raises(ValueError):
        Species({}, id="bad")
    with pytest.raises(ValueError):
        Reaction({"a": 1})
    with pytest.raises(ValueError, match="unknown"):
        water().validate(Reaction({"unknown": -1, "water": 1}))
    with pytest.raises(ValueError, match="unbalanced"):
        water().validate(Reaction({"h2": -1, "water": 1}))


def test_empty_and_immutable():
    empty = ReactionSpace([])
    assert (
        empty.rank,
        empty.dimension,
        empty.basis(),
        empty.circuits(),
        empty.independent([]),
    ) == (0, 0, (), (), ())
    counts = {"H": 2}
    s = Species(counts, id="h2")
    counts["H"] = 99
    assert s.composition["H"] == 2
    with pytest.raises(TypeError):
        s.composition["H"] = 3
    space = water()
    m = space.matrix
    m[0, 0] = 99
    assert space.matrix[0, 0] != 99


def test_budget():
    space = ReactionSpace([Species.from_formula("H2", id=str(i)) for i in range(50)])
    with pytest.raises(ValueError, match="budget"):
        space.circuits(max_subsets=1)


def test_circuits_are_minimal_and_span_small_example():
    space = ReactionSpace(
        [
            Species.from_formula("CO", id="co"),
            Species.from_formula("CO2", id="co2"),
            Species.from_formula("O2", id="o2"),
            Species.from_formula("C", id="c"),
        ]
    )
    circuits = space.circuits()
    assert len(space.independent(circuits)) == space.dimension
    for reaction in circuits:
        space.validate(reaction)
        support = [space.ids.index(k) for k in reaction.coefficients]
        submatrix = space.matrix[:, support]
        assert len(support) - submatrix.rank() == 1
        for i in range(len(support)):
            smaller = support[:i] + support[i + 1 :]
            assert space.matrix[:, smaller].rank() == len(smaller)


def test_normalization_and_electron_boundary():
    half = Reaction({"h2": -1, "o2": Fraction(-1, 2), "water": 1})
    assert half.primitive().coefficients["water"] == 2 * half.coefficients["water"]
    assert half.to_dict() != half.primitive().to_dict()
    with pytest.raises(ValueError, match="electron"):
        Species({}, id="positive_void", charge=1)


@given(st.lists(st.tuples(st.integers(1, 6), st.integers(1, 6)), min_size=1, max_size=9))
def test_random_exact_kernel(compositions):
    space = ReactionSpace(
        [Species({"C": c, "H": h}, id=str(i)) for i, (c, h) in enumerate(compositions)]
    )
    basis = space.basis()
    assert len(basis) == space.dimension
    assert len(space.independent(basis)) == len(basis)
    for reaction in basis:
        assert all(v == 0 for v in space.residual(reaction).values())
    reversed_space = ReactionSpace(reversed(space.species))
    assert [r.to_dict() for r in reversed_space.basis()] == [r.to_dict() for r in basis]
    if basis:
        vectors = Matrix.hstack(
            *(Matrix([r.coefficients.get(k, 0) for k in space.ids]) for r in basis)
        )
        assert vectors.rank() == space.dimension
