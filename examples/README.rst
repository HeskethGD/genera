Genera examples
===============

These examples demonstrate the algebraic-curve and Abelian-function APIs on
problems from the mathematical literature. They are version-controlled with
Genera but are not part of the installed ``genera`` package or wheel.

Run an example from the repository root using its module name, for example::

    python -m examples.manakov.manakov_kleinian_demo --steps 800 --samples 41

The dynamics examples compare their analytic Abelian-function solutions with
the shared arbitrary-precision fourth-order Runge--Kutta implementation in
``examples._rk4``. No plotting dependency is required.

Examples
--------

``bernatska``
    Trigonal genus-three periods and Kleinian identities following Bernatska.

``bobenko_reyman_semenov_tian_shansky``
    The working genus-three Kowalewski-top formula of Bobenko, Reyman and
    Semenov-Tian-Shansky. The experimental genus-two route is intentionally
    not included.

``kovalevskaya_original``
    The original theta-functional Kovalevskaya construction.

``manakov`` and ``neumann_moser``
    Finite-dimensional integrable systems checked against RK4 trajectories.

``onishi`` and ``onishi_trigonal_genus_3``
    Determinant and sigma-function identities.

``rn_desitter_9d``
    Genus-two inversion for nine-dimensional Reissner--Nordström--de Sitter
    geodesics, checked against radial RK4 integration.
