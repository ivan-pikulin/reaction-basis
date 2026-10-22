Guide
=====

Choose an operation
-------------------

.. list-table::
   :header-rows: 1

   * - Method
     - Result
     - Boundary
   * - ``balance(reactants=..., products=...)``
     - Smallest integer balance in the requested direction
     - Requires a unique solution using every participant
   * - ``basis()``
     - Complete independent basis of the reaction space
     - Not optimized for chemical simplicity
   * - ``circuits(max_species=4)``
     - Support-minimal reactions within the size bound
     - Can be dependent; enumeration is combinatorial
   * - ``independent(reactions)``
     - First independent subset spanning the candidates
     - Completeness depends on the supplied candidates

``balance`` rejects ambiguous or impossible balances, zero participants,
unknown IDs and species repeated within or across sides. Cancel shared
participants beforehand. It does not choose an arbitrary solution in a
multidimensional kernel.

``independent`` validates every candidate and preserves its direction and
normalization. The result has minimum cardinality for that candidate span,
not necessarily the simplest chemistry. It spans the full species space only
when its length equals ``space.dimension``.

``circuits`` finds dependencies with no smaller dependent support.
``max_species`` counts distinct participants, not coefficient sums. The search
tests subsets up to size ``min(max_species, N, rank + 1)`` and fails before
enumeration if their count exceeds ``max_subsets`` (default 100000). This is
neither all integer reactions nor a Hilbert basis. A subset budget is not a
wall-clock or RAM guarantee; use ``basis`` for large collections.

Species and IDs
---------------

``Species`` holds composition, charge, ID, label and optional multiplicity.
Geometry and measured/calculated properties belong in external tables keyed
by ID. Defaults encode alphabetically sorted composition and optional charge
and multiplicity; labels do not affect identity.

.. doctest::

   >>> from reaction_basis import Species, parse_formula
   >>> Species.from_formula("OH2").id
   'H2O'
   >>> Species({"H": 2, "O": 1}).id
   'H2O'
   >>> Species.from_formula("U", charge=4, multiplicity=3).id
   'U[q=+4][m=3]'
   >>> parse_formula("UO2(NO3)2")
   {'U': 1, 'O': 8, 'N': 2}

A default ID is not a structure fingerprint or Hill formula. Use explicit IDs
for isomers, geometries and source records. Collisions raise an error rather
than merging species. Unspecified and explicit multiplicities produce distinct
default IDs; multiplicity is not an additive conservation constraint.

Formula syntax supports known elements, positive integer counts, nested
parentheses/brackets and middle-dot hydrates such as ``CuSO4·5H2O``.
Pass charges separately. Isotopes, fractional occupancies and ASCII-dot
hydrates and leading coefficients are unsupported. Explicit composition still undergoes element/count
validation. ``Species({}, id="electron", charge=-1)`` represents an electron;
other empty compositions are rejected.

Exact balance and normalization
----------------------------------------

The conservation matrix :math:`A` has species as columns and elements plus
charge as rows. Balanced reaction vectors satisfy:

.. math::

   A\nu = 0, \qquad \nu\in\mathbb{Q}^{N}, \qquad
   \dim\ker A = N-\operatorname{rank}(A).

Reactants are negative and products positive. Arithmetic is exact, without
float tolerances. ``Reaction`` accepts integers, rational strings, ``Fraction``
and SymPy rationals, but rejects floats and requires both sides. Construction
alone does not check balance: call ``space.validate(r)`` or inspect
``space.residual(r)``.

``basis`` returns primitive integer vectors, oriented with the smallest
participating ID negative. Ordering is deterministic, but renaming IDs can
change the basis or direction. A basis is not mathematically unique.

.. doctest::

   >>> from reaction_basis import Reaction
   >>> r = Reaction({"H2": -1, "O2": "-1/2", "water": 1})
   >>> r.primitive(canonical=False).to_dict()
   {'H2': '-2', 'O2': '-1', 'water': '2'}

``primitive`` clears denominators and divides by the coefficient GCD. Default
``canonical=True`` also chooses direction; use ``False`` to preserve it.
Changing normalization scales all assigned extensive quantities.

External data and serialization
----------------------------------------

Use the read-only ID-to-coefficient mapping with any additive species property.
Energetics is not part of the package:

.. doctest::

   >>> from math import fsum
   >>> from reaction_basis import ReactionSpace
   >>> space = ReactionSpace(Species.from_formula(f) for f in ["H2", "O2", "H2O"])
   >>> r = space.balance(reactants=["H2", "O2"], products=["H2O"])
   >>> energies = {"H2": -1.0, "O2": -2.0, "H2O": -4.0}  # Illustrative only.
   >>> fsum(float(nu) * energies[key] for key, nu in r.coefficients.items())
   -4.0
   >>> import json
   >>> Reaction(json.loads(json.dumps(r.to_dict()))).to_dict() == r.to_dict()
   True

The consumer checks finite values, units, methods, convergence and states.
Missing required properties must not be replaced with zero. Uncertainties
require covariance propagation, not summing standard deviations.
``to_dict`` exports exact rational strings, not a species catalog; preserve
the ID mapping and provenance separately. A complete consumer example is in
``examples/external_energy.py`` in the source checkout.

``space.matrix`` returns a copy. Columns follow ``space.ids`` and rows follow
``space.constraints``. Align a reaction using
``[r.coefficients.get(key, 0) for key in space.ids]`` after validation;
zero correctly means that a species does not participate.

Limits and diagnostics
----------------------

Balance is not proof of stability, spontaneity, a mechanism or a spin-allowed
transition. Isotope-, bond-type- and reservoir-specific constraints are not
implemented. Composition/coefficient mappings and the space's catalog and
ordering properties are read-only. The conservation matrix is copied on access.

``ValueError`` reports invalid values, duplicate/unknown IDs, imbalance,
ambiguous balancing and exceeded budgets. ``TypeError`` reports unsupported
types. There are no silent skips, approximate fallbacks or automatic retries.
