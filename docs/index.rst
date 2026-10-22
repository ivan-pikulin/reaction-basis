Reaction Basis
==============

Build balanced reactions and independent reaction sets with exact arithmetic.
Keep isomers and electronic states distinct through stable species IDs.

.. doctest::

   >>> from reaction_basis import ReactionSpace, Species
   >>> space = ReactionSpace(Species.from_formula(f) for f in ["H2", "O2", "H2O"])
   >>> space.balance(reactants=["H2", "O2"], products=["H2O"]).format()
   '2 H2 + O2 -> 2 H2O'

.. toctree::
   :maxdepth: 1

   quickstart
   guide
   api
