"""Linear algebra over Q; no floating-point balancing or identity-by-composition."""

from __future__ import annotations

import re
from collections.abc import Iterable, Mapping
from dataclasses import dataclass, field
from fractions import Fraction
from functools import reduce
from itertools import combinations
from math import comb, gcd, lcm
from types import MappingProxyType

from sympy import Matrix, Rational, zeros

_ELEMENTS = frozenset(
    "H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Cs Ba La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn Fr Ra Ac Th Pa U Np Pu Am Cm Bk Cf Es Fm Md No Lr Rf Db Sg Bh Hs Mt Ds Rg Cn Nh Fl Mc Lv Ts Og".split()
)


def _integer(value: object, name: str, *, positive: bool = False) -> int:
    if isinstance(value, bool) or not isinstance(value, int):
        raise TypeError(f"{name} must be an integer")
    if positive and value <= 0:
        raise ValueError(f"{name} must be positive")
    return value


def parse_formula(formula: str) -> dict[str, int]:
    """Parse strict neutral formulas, nested ()/[], and middle-dot hydrates.

    Charges, isotopes, fractional occupancies and ASCII-dot hydrates are not
    inferred. Specify charge separately; use explicit composition for other cases.
    """
    if not isinstance(formula, str) or not formula:
        raise ValueError("formula must be a nonempty string")
    result: dict[str, int] = {}
    for part_index, part in enumerate(formula.split("·")):
        tokens = re.findall(r"[A-Z][a-z]?|[1-9][0-9]*|[()\[\]]", part)
        if not tokens or "".join(tokens) != part:
            raise ValueError(f"unsupported formula syntax: {formula!r}")
        pos = 0
        factor = 1
        if tokens[0].isdigit():
            if part_index == 0:
                raise ValueError("leading coefficients and isotope prefixes are unsupported")
            factor = int(tokens[0])
            pos += 1

        def group(close: str | None = None, tokens: list[str] = tokens) -> dict[str, int]:
            nonlocal pos
            counts: dict[str, int] = {}
            while pos < len(tokens):
                token = tokens[pos]
                if token in (")", "]"):
                    if token != close:
                        raise ValueError("mismatched brackets")
                    pos += 1
                    if not counts:
                        raise ValueError("empty group")
                    return counts
                pos += 1
                if token in ("(", "["):
                    unit = group(")" if token == "(" else "]")
                elif token in _ELEMENTS:
                    unit = {token: 1}
                else:
                    raise ValueError(f"unknown element or misplaced count: {token}")
                n = 1
                if pos < len(tokens) and tokens[pos].isdigit():
                    n = int(tokens[pos])
                    pos += 1
                for element, count in unit.items():
                    counts[element] = counts.get(element, 0) + count * n
            if close is not None:
                raise ValueError("unclosed group")
            if not counts:
                raise ValueError("empty formula component")
            return counts

        for element, count in group().items():
            result[element] = result.get(element, 0) + factor * count
    return result


@dataclass(frozen=True)
class Species:
    """Composition with optional identity and state; duplicate IDs fail explicitly.

    Default ID encodes composition, charge and multiplicity. Isomers with the
    same state need explicit IDs; they are never silently merged.

    Args:
        composition: Known element symbols mapped to positive integer counts.
            An empty mapping is allowed only for an electron (charge -1).
        id: Unique key. Defaults to sorted composition with charge/state suffixes.
        charge: Integer total charge, in elementary-charge units.
        label: Optional display label; does not determine identity.
        multiplicity: Optional positive integer spin multiplicity; not conserved.

    Raises:
        ValueError: Invalid ID, element, count, state or empty composition.
        TypeError: Charge, count or multiplicity is not an integer.
    """

    composition: Mapping[str, int]
    id: str | None = field(default=None, kw_only=True)
    charge: int = field(default=0, kw_only=True)
    label: str | None = field(default=None, kw_only=True)
    multiplicity: int | None = field(default=None, kw_only=True)

    def __post_init__(self) -> None:
        if self.id is not None and (not isinstance(self.id, str) or not self.id.strip()):
            raise ValueError("species ID must be a nonempty string")
        _integer(self.charge, "charge")
        if self.multiplicity is not None:
            _integer(self.multiplicity, "multiplicity", positive=True)
        if self.label is not None and (not isinstance(self.label, str) or not self.label.strip()):
            raise ValueError("label must be a nonempty string")
        if not isinstance(self.composition, Mapping):
            raise TypeError("composition must be a mapping of elements to counts")
        counts = dict(self.composition)
        for element, count in counts.items():
            if element not in _ELEMENTS:
                raise ValueError(f"unknown element: {element}")
            _integer(count, f"count of {element}", positive=True)
        if not counts and self.charge != -1:
            raise ValueError("empty composition is supported only for an electron (charge -1)")
        object.__setattr__(self, "composition", MappingProxyType(counts))
        if self.id is None:
            identifier = (
                "".join(e + (str(counts[e]) if counts[e] != 1 else "") for e in sorted(counts))
                or "e"
            )
            if self.charge:
                identifier += f"[q={self.charge:+d}]"
            if self.multiplicity is not None:
                identifier += f"[m={self.multiplicity}]"
            object.__setattr__(self, "id", identifier)

    @classmethod
    def from_formula(
        cls,
        formula: str,
        *,
        id: str | None = None,
        charge: int = 0,
        label: str | None = None,
        multiplicity: int | None = None,
    ) -> Species:
        """Create a species from a strict formula and keyword-only options.

        Args:
            formula: Neutral formula syntax; pass charge separately.
            id: Optional unique key; required to distinguish same-state isomers.
            charge: Integer total charge.
            label: Display label; defaults to the supplied formula.
            multiplicity: Optional positive spin multiplicity.

        Returns:
            A species with a read-only composition mapping.

        Raises:
            ValueError: Unsupported formula or invalid species values.
            TypeError: Invalid integer fields.
        """
        return cls(
            parse_formula(formula),
            id=id,
            charge=charge,
            label=formula if label is None else label,
            multiplicity=multiplicity,
        )


def _rational(value: object) -> Rational:
    if isinstance(value, (bool, float)):
        raise TypeError("use exact integers, Fraction, Rational or rational strings, not floats")
    if not isinstance(value, (int, Fraction, Rational, str)):
        raise TypeError("unsupported coefficient type")
    try:
        return Rational(value)
    except (ValueError, TypeError, ZeroDivisionError) as exc:
        raise ValueError(f"invalid rational coefficient: {value!r}") from exc


@dataclass(frozen=True)
class Reaction:
    """Signed exact stoichiometry: products positive, reactants negative.

    Args:
        coefficients: Species IDs mapped to integers, rational strings,
            Fraction or SymPy Rational values. Zero entries are removed.

    Raises:
        ValueError: Invalid coefficients or a reaction without both sides.
        TypeError: Float, bool or unsupported coefficient types.

    Construction does not establish balance; use a reaction space to validate.
    """

    coefficients: Mapping[str, int | Fraction | Rational | str]

    def __post_init__(self) -> None:
        if not isinstance(self.coefficients, Mapping):
            raise TypeError("coefficients must be a mapping of species IDs to exact values")
        values = {}
        for key, value in self.coefficients.items():
            if not isinstance(key, str) or not key.strip():
                raise ValueError("coefficient keys must be species IDs")
            rational = _rational(value)
            if rational:
                values[key] = rational
        if (
            not values
            or not any(v > 0 for v in values.values())
            or not any(v < 0 for v in values.values())
        ):
            raise ValueError("reaction needs nonzero reactants and products")
        object.__setattr__(self, "coefficients", MappingProxyType(values))

    def primitive(self, *, canonical: bool = True) -> Reaction:
        """Smallest integer coefficients; optional orientation independent of input order."""
        if not isinstance(canonical, bool):
            raise TypeError("canonical must be a bool")
        denominator = lcm(*(int(v.q) for v in self.coefficients.values()))
        integers = {k: int(v * denominator) for k, v in self.coefficients.items()}
        divisor = reduce(gcd, (abs(v) for v in integers.values()))
        sign = 1
        if canonical and integers[min(integers)] > 0:
            sign = -1
        return Reaction({k: sign * v // divisor for k, v in integers.items()})

    def reversed(self) -> Reaction:
        """Return a new reaction with all signs reversed."""
        return Reaction({k: -v for k, v in self.coefficients.items()})

    def to_dict(self) -> dict[str, str]:
        """Export sorted IDs and rational strings for JSON; restore with Reaction(data)."""
        return {k: str(v) for k, v in sorted(self.coefficients.items())}

    def format(self, *, labels: Mapping[str, str] | None = None) -> str:
        """Format both sides, using optional ID-to-label replacements.

        Labels affect display only. Terms remain sorted by species ID.
        """
        labels = labels or {}

        def side(positive: bool) -> str:
            terms = []
            for key, value in sorted(self.coefficients.items()):
                if (value > 0) == positive:
                    n = abs(value)
                    terms.append(("" if n == 1 else f"{n} ") + labels.get(key, key))
            return " + ".join(terms)

        return f"{side(False)} -> {side(True)}"


class ReactionSpace:
    """Conservation matrix and deterministic exact kernel basis.

    Species sorted by ID; elements sorted by symbol; charge always conserved.
    Spin multiplicity is metadata, not an additive conserved quantity.

    Args:
        species: Finite iterable of species with unique IDs. Empty input is valid.

    Raises:
        ValueError: Duplicate species IDs.
        TypeError: Entries are not Species objects.

    Catalog attributes (species, ids, elements, constraints) are for read access.
    """

    def __init__(self, species: Iterable[Species]):
        entries = tuple(species)
        if any(not isinstance(s, Species) for s in entries):
            raise TypeError("expected Species objects")
        if len({s.id for s in entries}) != len(entries):
            raise ValueError("duplicate species ID")
        self._species = tuple(sorted(entries, key=lambda s: s.id))
        self._ids = tuple(s.id for s in self.species)
        self._elements = tuple(sorted({e for s in self.species for e in s.composition}))
        self._constraints = self.elements + ("charge",)
        self._matrix = (
            Matrix(
                [[s.composition.get(e, 0) for s in self.species] for e in self.elements]
                + [[s.charge for s in self.species]]
            )
            if entries
            else zeros(1, 0)
        )

    @property
    def species(self) -> tuple[Species, ...]:
        """Read-only catalog sorted by species ID."""
        return self._species

    @property
    def ids(self) -> tuple[str, ...]:
        """Read-only column order of the conservation matrix."""
        return self._ids

    @property
    def elements(self) -> tuple[str, ...]:
        """Read-only element symbols in alphabetical order."""
        return self._elements

    @property
    def constraints(self) -> tuple[str, ...]:
        """Read-only row order (elements followed by charge)."""
        return self._constraints

    @property
    def matrix(self) -> Matrix:
        """Copy of the conservation matrix; rows follow constraints, columns follow ids."""
        return self._matrix.copy()

    @property
    def rank(self) -> int:
        """Exact rank of the conservation matrix."""
        return self._matrix.rank()

    @property
    def dimension(self) -> int:
        """Dimension of the balanced space (species count minus rank)."""
        return len(self.species) - self.rank

    def _vector(self, reaction: Reaction) -> Matrix:
        if not isinstance(reaction, Reaction):
            raise TypeError("expected a Reaction object")
        unknown = set(reaction.coefficients) - set(self.ids)
        if unknown:
            raise ValueError(f"unknown species IDs: {sorted(unknown)}")
        return Matrix([reaction.coefficients.get(key, 0) for key in self.ids])

    def residual(self, reaction: Reaction) -> dict[str, Rational]:
        """Return exact element/charge residuals; all are zero for a balanced reaction.

        Raises:
            ValueError: Unknown species IDs or an empty species space.
        """
        vector = self._vector(reaction)
        if not self.ids:
            raise ValueError("empty species space")
        return dict(zip(self.constraints, self._matrix * vector, strict=True))

    def validate(self, reaction: Reaction) -> None:
        """Require known IDs and exact balance; return None on success.

        Raises:
            ValueError: Unknown species, empty space or nonzero residuals.
        """
        nonzero = {k: v for k, v in self.residual(reaction).items() if v}
        if nonzero:
            raise ValueError(f"unbalanced reaction: {nonzero}")

    def basis(self) -> tuple[Reaction, ...]:
        """Return a complete independent kernel basis as primitive integer reactions.

        The smallest participating ID has a negative coefficient. Size equals
        dimension; no claim of chemical simplicity or preferred direction.
        """
        reactions = tuple(
            Reaction(dict(zip(self.ids, vector, strict=True))).primitive()
            for vector in self._matrix.nullspace()
        )
        for reaction in reactions:
            self.validate(reaction)
        return reactions

    def independent(self, reactions: Iterable[Reaction]) -> tuple[Reaction, ...]:
        """Select first independent input reactions, preserving their normalization."""
        reactions = tuple(reactions)
        for reaction in reactions:
            self.validate(reaction)
        if not reactions:
            return ()
        vectors = Matrix.hstack(*(self._vector(r) for r in reactions))
        return tuple(reactions[i] for i in vectors.rref()[1])

    def balance(self, *, reactants: Iterable[str], products: Iterable[str]) -> Reaction:
        """Balance specified sides when the solution is unique up to scale.

        All listed species must participate. Ambiguous, impossible or incorrectly
        oriented reactions fail instead of selecting an arbitrary solution.
        """
        if isinstance(reactants, str) or isinstance(products, str):
            raise TypeError("sides must be iterables of species IDs, not strings")
        left, right = tuple(reactants), tuple(products)
        if not left or not right:
            raise ValueError("both sides must be nonempty")
        ids = left + right
        if any(not isinstance(key, str) for key in ids):
            raise TypeError("sides must contain species IDs")
        if len(set(ids)) != len(ids):
            raise ValueError("species IDs must be distinct within and across sides")
        unknown = set(ids) - set(self.ids)
        if unknown:
            raise ValueError(f"unknown species IDs: {sorted(unknown)}")
        columns = [self.ids.index(key) for key in ids]
        kernel = self._matrix[:, columns].nullspace()
        if not kernel:
            raise ValueError("no balanced reaction exists for these species")
        if len(kernel) != 1:
            raise ValueError("ambiguous balance: more than one independent solution")
        vector = kernel[0]
        if not all(vector):
            raise ValueError("not all listed species can participate")
        if vector[0] > 0:
            vector = -vector
        if any(v >= 0 for v in vector[: len(left)]) or any(v <= 0 for v in vector[len(left) :]):
            raise ValueError("specified sides are incompatible with balance")
        result = Reaction(dict(zip(ids, vector, strict=True))).primitive(canonical=False)
        self.validate(result)
        return result

    def circuits(self, *, max_species: int = 4, max_subsets: int = 100_000) -> tuple[Reaction, ...]:
        """Enumerate support-minimal balanced reactions up to a support bound.

        These are NOT necessarily linearly independent. Enumeration is
        combinatorial; a finite budget fails explicitly, never returns a silently
        truncated result. Basis generation is the efficient default.
        """
        _integer(max_species, "max_species", positive=True)
        _integer(max_subsets, "max_subsets", positive=True)
        limit = min(max_species, len(self.ids), self.rank + 1)
        required = sum(comb(len(self.ids), k) for k in range(2, limit + 1))
        if required > max_subsets:
            raise ValueError(f"enumeration needs {required} subsets, budget is {max_subsets}")
        result = []
        for size in range(2, limit + 1):
            for subset in combinations(range(len(self.ids)), size):
                kernel = self._matrix[:, subset].nullspace()
                if len(kernel) == 1 and all(kernel[0]):
                    values = kernel[0]
                    if any(v > 0 for v in values) and any(v < 0 for v in values):
                        reaction = Reaction(
                            {self.ids[i]: v for i, v in zip(subset, values, strict=True)}
                        )
                        result.append(reaction.primitive())
        return tuple(result)
