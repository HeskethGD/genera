Genera documentation
====================

Genera provides arbitrary-precision numerical algorithms for algebraic curves,
Riemann theta functions, Abelian functions, and integrable systems. It uses
mpmath for numerical types, arithmetic, matrices, and precision contexts.

The public API is independent of the mpmath namespace. Functions use
``mpmath.mp`` by default and accept a keyword-only ``ctx`` argument when a
separate numerical context is required.

For an overview of the mathematical conventions and available methods, start
with :doc:`algebraic_curves` and :doc:`abelian`.

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   algebraic_curves
   abelian
   examples
   references
