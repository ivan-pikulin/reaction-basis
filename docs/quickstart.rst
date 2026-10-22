Quick start
===========

Install
-------

Requires Python 3.11 or later.

.. code-block:: console

   pip install reaction-basis

Balance a reaction
------------------

Define species once, then pass their IDs as the two sides. All listed species
must participate; the result must be unique up to scale.

.. doctest::

   >>> from reaction_basis import ReactionSpace, Species
   >>> space = ReactionSpace([Species.from_formula(f) for f in ("H2", "O2", "H2O")])
   >>> reaction = space.balance(reactants=["H2", "O2"], products=["H2O"])
   >>> reaction.format()
   '2 H2 + O2 -> 2 H2O'
   >>> space.validate(reaction)  # Raises on imbalance; returns None on success.

Generate reactions from a collection
----------------------------------------

No reactant/product assignment is needed. A kernel basis spans all balanced
stoichiometric vectors over the supplied species.

.. doctest::

   >>> formulas = ["H2", "O2", "H2O", "H2O2", "C", "CO", "CO2", "CH4", "CH3OH"]
   >>> space = ReactionSpace(Species.from_formula(f) for f in formulas)
   >>> space.rank, space.dimension
   (3, 6)
   >>> reactions = space.basis()
   >>> len(reactions)
   6
   >>> all(not any(space.residual(r).values()) for r in reactions)
   True

Print the results with ``for r in reactions: print(r.format())``.
Canonical direction is a convention, not a prediction of spontaneity.

Prefer small reactions
----------------------

Enumerate minimal supports, then select an independent subset:

.. doctest::

   >>> candidates = space.circuits(max_species=4)
   >>> selected = space.independent(candidates)
   >>> len(selected) == space.dimension
   True

That equality holds for this example, not for every support-limited search.
See :doc:`guide` for guarantees and cost.

Keep isomers distinct
---------------------

Provide explicit IDs for separate structures with identical composition.

.. doctest::

   >>> ethanol = Species.from_formula("C2H6O", id="ethanol")
   >>> ether = Species.from_formula("C2H6O", id="dimethyl_ether")
   >>> ReactionSpace([ethanol, ether]).dimension
   1

Use ``charge=`` and ``multiplicity=`` for electronic states. See :doc:`guide`
for default IDs and the supported formula syntax.
