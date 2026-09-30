Genera documentation
====================

Genera provides arbitrary-precision numerical algorithms for algebraic curves,
Abelian functions, and integrable systems. It uses mpmath for numerical types,
arithmetic, matrices, and precision contexts.

The mathematical API and its detailed documentation will be added during the
migration from the mpmath algebraic-curve development branch. Genera's
implementations receive an mpmath numerical context internally, allowing them
to respect the caller's working precision without modifying mpmath's public
namespace.

Riemann theta functions
-----------------------

.. autofunction:: genera.rtheta

.. autofunction:: genera.rtheta_jet

Kleinian functions
------------------

.. autofunction:: genera.kleinian_sigma

.. autofunction:: genera.kleinian_sigma_jet

.. autofunction:: genera.kleinian_zeta

.. autofunction:: genera.kleinian_p

.. autofunction:: genera.kleinian_baker_akhiezer

.. toctree::
   :maxdepth: 1

   references
