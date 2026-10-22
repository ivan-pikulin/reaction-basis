API reference
=============

Signatures and descriptions below are imported from the installed package.

Species
-------

.. autoclass:: reaction_basis.Species
   :members: from_formula

Reactions
---------

.. autoclass:: reaction_basis.Reaction
   :members: primitive, reversed, to_dict, format

Reaction spaces
---------------

.. autoclass:: reaction_basis.ReactionSpace
   :members: species, ids, elements, constraints, matrix, rank, dimension, residual, validate, basis, independent, balance, circuits

Formula parsing
---------------

.. autofunction:: reaction_basis.parse_formula
