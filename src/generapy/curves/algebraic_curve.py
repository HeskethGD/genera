"""Object-oriented interface to numerical algebraic curves."""

import warnings

from .._context import resolve_context
from . import _operations, _records, charts


class Curve:
    """Numerical computations on a plane algebraic curve.

    Parameters
    ----------
    polynomial : ``Mapping[tuple[int, int], scalar]``
        Sparse coefficients ``{(i, j): coefficient}`` of ``F(x, y) = 0``,
        where ``i`` and ``j`` are nonnegative integer powers of ``x`` and
        ``y``. ``curve.polynomial`` returns a copy of this mapping.
    differentials_kind_1 : ``Iterable[Callable[[scalar, scalar], scalar]] | None``
        Ordered first-kind basis of ``g`` holomorphic forms, where ``g`` is
        the genus. Each callable ``f(x, y) -> scalar`` returns the numeric
        coefficient of ``dx`` in ``f(x, y) dx``. ``None`` selects an automatic
        basis where supported.
    differentials_kind_2 : ``Iterable[Callable[[scalar, scalar], scalar]] | None``
        Ordered second-kind basis of ``g`` meromorphic forms with zero
        residues. Each callable ``f(x, y) -> scalar`` returns the numeric
        coefficient of ``dx``. ``None`` selects a compatible automatic basis
        where supported; otherwise, methods needing this basis raise an error.
    ctx : mpmath context or ``None``
        Numerical context, defaulting to ``mpmath.mp``. A clone allows
        independent working precision.

    Notes
    -----
    Here ``scalar`` denotes a real or complex numeric value. Supplied bases
    are stored as tuples and used by all relevant methods.
    Definitions and examples are in :ref:`custom-differential-bases`.

    Examples
    --------
    A genus-one curve:

    >>> from generapy import Curve
    >>> # F(x, y) = y**2 + x - x**3
    >>> polynomial = {(0, 2): 1, (1, 0): 1, (3, 0): -1}
    >>> curve = Curve(polynomial=polynomial)
    >>> curve.genus
    1

    A genus-two curve with a mixed ``x*y`` term:

    >>> # F(x, y) = y**2 + 2*x*y + x**2 - x**5 + x - 1
    >>> polynomial = {(0, 2): 1, (1, 1): 2, (2, 0): 1,
    ...               (5, 0): -1, (1, 0): 1, (0, 0): -1}
    >>> Curve(polynomial=polynomial).genus
    2

    An independent numerical context:

    >>> from mpmath import mp
    >>> ctx = mp.clone()
    >>> ctx.dps = 50
    >>> curve = Curve(polynomial=polynomial, ctx=ctx)
    >>> curve.ctx.dps
    50
    """

    def __init__(self, polynomial, *, differentials_kind_1=None,
                 differentials_kind_2=None, ctx=None):
        ctx = resolve_context(ctx)
        self.ctx = ctx
        self._creation_state = _operations._curve_cache_state(ctx)
        self._warned_states = set()
        if not hasattr(polynomial, "items"):
            raise ValueError(
                "polynomial must be a sparse mapping from (x_power, y_power) to coefficients")
        self._polynomial = dict(polynomial.items())
        self._prepared, self._hyperelliptic_model = (
            _operations._normalize_algebraic_curve_input(ctx, self._polynomial))
        self._classified_states = {
            self._creation_state: _records._ClassifiedCurve(
                self._prepared, self._hyperelliptic_model)
        }
        self._differentials_kind_1 = (
            None if differentials_kind_1 is None else
            _operations._curve_differential_sequence(
                differentials_kind_1, "differentials_kind_1"))
        self._differentials_kind_2 = (
            None if differentials_kind_2 is None else
            _operations._curve_differential_sequence(
                differentials_kind_2, "differentials_kind_2"))
        self._custom_basis = (self._differentials_kind_1 is not None
                              or self._differentials_kind_2 is not None)
        self._first_kind_periods = {}

    @staticmethod
    def _copy_first_kind_periods(result):
        """Copy mutable matrices in a first-kind result."""
        return result._replace(
            omega=+result.omega,
            omega_prime=+result.omega_prime,
            tau=+result.tau,
        )

    def _check_precision(self):
        """Warn once for each numerical state different from construction."""
        state = _operations._curve_cache_state(self.ctx)
        if state != self._creation_state and state not in self._warned_states:
            warnings.warn(
                "the Curve context changed after construction; "
                "numerical stages will be recomputed for the current state, "
                "but inexact input coefficients retain their construction "
                "precision",
                UserWarning,
                stacklevel=3,
            )
            self._warned_states.add(state)

    def _call(self, function, *args, _basis=False, **kwargs):
        self._check_precision()
        state = _operations._curve_cache_state(self.ctx)
        classified = self._classified_states.get(state)
        if classified is None:
            prepared, hyperelliptic = (
                _operations._normalize_algebraic_curve_input(
                    self.ctx, self._polynomial))
            classified = _records._ClassifiedCurve(
                prepared, hyperelliptic)
            self._classified_states[state] = classified
        if _basis and self._custom_basis:
            classified = _records._ClassifiedCurve(classified.curve, None)
        return function(self.ctx, classified, *args, **kwargs)

    def _require_second_kind_basis(self, method):
        if self._differentials_kind_2 is None and (
                self._custom_basis or self._hyperelliptic_model is None):
            raise ValueError(
                f"{method} requires a second-kind basis; no compatible automatic "
                "basis is available for this Curve. Construct Curve(..., "
                "differentials_kind_2=...) with a compatible sequence of forms.")

    @property
    def differentials_kind_1(self):
        """The supplied first-kind basis, or ``None`` for automatic construction."""
        return self._differentials_kind_1

    @property
    def differentials_kind_2(self):
        """The supplied second-kind basis, or ``None`` if none was supplied."""
        return self._differentials_kind_2

    @property
    def polynomial(self):
        """Return a copy of the defining sparse polynomial mapping."""
        return dict(self._polynomial)

    @property
    def x_degree(self):
        """Degree of the defining polynomial in ``x``."""
        return self._prepared.x_degree

    @property
    def y_degree(self):
        """Degree of the defining polynomial in ``y``."""
        return self._prepared.y_degree

    @property
    def branch_locus(self):
        r"""Return the finite branch locus of a plane algebraic curve.

        The returned ``CurveBranchLocus`` record contains the degree of the
        ``x`` projection, the distinct finite branch values above which the
        projection ramifies, and the ascending coefficients of the
        y-derivative resultant whose roots they are.  Ramification above
        infinity is reported by :attr:`Curve.monodromy` instead,
        because it requires monodromy rather than the resultant alone.

        The lemniscatic curve :math:`y^2 = x^3 - x` has a two-sheeted
        projection with three finite branch values::

            >>> from generapy import Curve
            >>> locus = Curve({(0, 2): 1, (1, 0): 1, (3, 0): -1}).branch_locus
            >>> locus.degree
            2
            >>> locus.branch_values
            (mpf('-1.0'), mpf('0.0'), mpf('1.0'))
        """
        return self._call(_operations.branch_locus)

    @property
    def monodromy(self):
        r"""Return the monodromy of a plane algebraic curve over the x-line.

        The curve is continued numerically along guarded radial loops around
        the finite branch values, from an exterior base point chosen
        automatically.  A large outer loop supplies the monodromy at infinity
        geometrically.  The returned ``CurveMonodromy`` record contains the
        computational base point, its ordered fibre, the branch values, the
        counter-clockwise product-ordered local permutations, the permutation
        at infinity, the total ramification, the genus from
        Riemann--Hurwitz, the transitivity and product identities of the
        permutation system, and the minimum geometric clearance of the
        continuation paths.

        The sheet labels refer to the internally selected base fibre.  The
        routing continuations used to compute them are private.

        Each finite branch value of the lemniscatic curve
        :math:`y^2 = x^3 - x` exchanges its two sheets, as does infinity::

            >>> from generapy import Curve
            >>> monodromy = Curve({(0, 2): 1, (1, 0): 1, (3, 0): -1}).monodromy
            >>> monodromy.genus
            1
            >>> monodromy.permutations
            ((1, 0), (1, 0), (1, 0))
            >>> monodromy.infinity_permutation
            (1, 0)
        """
        return self._call(_operations.monodromy)

    @property
    def genus(self):
        r"""Genus of the compact curve.

        See :attr:`genus_data` for projection degree and ramification.
        """
        return self._call(_operations.genus_data).genus

    @property
    def genus_data(self):
        r"""Return the genus of a plane algebraic curve.

        Returns a ``CurveGenus`` record with the genus, projection degree, and
        total ramification. The Riemann-Hurwitz balance :math:`2g-2 = -2d+r`
        can be verified directly::

            >>> from generapy import Curve
            >>> Curve({(0, 2): 1, (1, 0): 1, (3, 0): -1}).genus_data
            CurveGenus(genus=1, degree=2, ramification=4)
        """
        return self._call(_operations.genus_data)

    @property
    def homology(self):
        r"""Return the homology marking selected by the curve's bound bases.

        Returns a ``CurveHomology`` record containing ``2*genus`` cycles in a
        canonical symplectic basis with the standard intersection form. The
        ``marking`` field identifies the computational engine: ``'baker'`` for
        automatic hyperelliptic bases, ``'geometric-polygon'`` for custom bases
        and other supported curves.

        The lemniscatic curve :math:`y^2 = x^3 - x` uses Baker marking::

            >>> from generapy import Curve
            >>> homology = Curve({(0, 2): 1, (1, 0): 1, (3, 0): -1}).homology
            >>> homology.genus, homology.marking
            (1, 'baker')
            >>> homology.intersection_form
            ((0, 1), (-1, 0))
        """
        return self._call(_operations.homology, _basis=True)

    def periods_kind_1(self):
        r"""Return first-kind half-periods and the normalized Riemann matrix.

        Generalizes elliptic period computation to genus :math:`g`. Returns a
        ``CurvePeriodsKind1`` record with half-period matrices ``omega``,
        ``omega_prime``, and the normalized Riemann matrix ``tau`` (the symmetric
        part of ``omega**-1 * omega_prime``). Full periods are ``2*omega`` and
        ``2*omega_prime``.

        For a holomorphic basis :math:`du_1,\ldots,du_g` and canonical
        cycles :math:`a_j,b_j`, the matrix entries are

        .. math::

            \begin{aligned}
            2\omega_{ij} &= \oint_{a_j} du_i, \\
            2\omega'_{ij} &= \oint_{b_j} du_i, \\
            \tau &= \omega^{-1}\omega'.
            \end{aligned}

        Rows index differentials and columns index cycles. For the elliptic
        curve :math:`y^2=x^3-x`, the automatic basis has just one form,
        :math:`du=dx/y`, so these matrices reduce to scalars:

        .. math::

            \begin{aligned}
            2\omega &= \oint_a \frac{dx}{y}, \\
            2\omega' &= \oint_b \frac{dx}{y}.
            \end{aligned}

        Uses the ordered first-kind basis bound at construction. If no basis
        was supplied, hyperelliptic curves get ``x**k dx/z`` in Baker marking;
        other supported curves get a basis from Newton polygon interior points.
        Supplying either basis at construction selects geometric-polygon marking,
        using an automatic geometric first-kind basis if necessary.
        See :ref:`custom-differential-bases` for definitions and a worked example.

        A non-positive-definite period matrix raises ``ValueError``.

        The lemniscatic curve has normalized period matrix ``tau = i``::

            >>> from generapy import Curve
            >>> from mpmath import mp
            >>> mp.dps = 15
            >>> curve = Curve({(0, 2): 1, (1, 0): 1, (3, 0): -1})
            >>> data = curve.periods_kind_1()
            >>> mp.re(data.tau[0, 0]), mp.im(data.tau[0, 0])
            (mpf('0.0'), mpf('1.0'))

        A supplied basis uses geometric-polygon marking::

            >>> custom = Curve(curve.polynomial, differentials_kind_1=(lambda x, y: 1 / y,))
            >>> data = custom.periods_kind_1()
            >>> data.marking
            'geometric-polygon'
            >>> curve.validate(data).passed
            True
        """
        self._check_precision()
        state = _operations._curve_cache_state(self.ctx)
        cached = self._first_kind_periods.get(state)
        if cached is not None:
            return self._copy_first_kind_periods(cached)
        result = self._call(_operations.periods, self._differentials_kind_1, _basis=True)
        self._first_kind_periods[state] = self._copy_first_kind_periods(result)
        return result

    def periods_kind_2(self):
        r"""Return second-kind half-periods and kappa.

        The ``CurvePeriodsKind2`` record contains ``eta``, ``eta_prime``
        and ``kappa``, with ``2*eta = -integral_a(dr)`` and
        ``2*eta_prime = -integral_b(dr)``. The returned ``kappa`` is the symmetric
        part of ``eta * omega**-1`` for the compatible first-kind half-periods.

        Second-kind differentials :math:`dr_i` are meromorphic forms with
        zero residues. Their cycle integrals define the half-periods with a
        minus sign:

        .. math::

            \begin{aligned}
            2\eta_{ij} &= -\oint_{a_j} dr_i, \\
            2\eta'_{ij} &= -\oint_{b_j} dr_i, \\
            \kappa &= \tfrac12\bigl(\eta\omega^{-1}
                       +(\eta\omega^{-1})^{\mathsf T}\bigr).
            \end{aligned}

        For :math:`y^2=x^3-x`, the automatic BEL form is
        :math:`dr=x\,dx/(4y)`. It has a double pole at infinity with zero
        residue, whereas :math:`du=dx/y` is holomorphic. In this basis,

        .. math::

            \begin{aligned}
            2\eta &= -\oint_a \frac{x\,dx}{4y}, \\
            2\eta' &= -\oint_b \frac{x\,dx}{4y}.
            \end{aligned}

        These periods describe the changes of a second-kind integral around
        cycles, analogous to the quasi-periods of Weierstrass zeta.

        Uses the second-kind basis bound through ``Curve(differentials_kind_2=...)``.
        With neither basis supplied, recognized hyperelliptic models construct
        the automatic BEL basis in Baker marking. Otherwise, a missing
        second-kind basis raises ``ValueError`` explaining how to supply it.
        Each supplied basis must contain one coefficient-of-``dx`` callable
        per genus. See :ref:`custom-differential-bases` for a complete example.

        Supplied second-kind forms must have zero residues and be regular on
        the integration paths; numerical convergence checks do not establish
        those properties. ``engine`` and ``marking`` describe the cycle basis.
        Compatible first-kind results are cached for later calls to
        :meth:`periods_kind_1`.
        """
        self._require_second_kind_basis("periods_kind_2()")
        first, second = self._call(
            _operations.periods, self._differentials_kind_1,
            second_kind=True, second_differentials=self._differentials_kind_2,
            _return_first=True, _basis=True)
        state = _operations._curve_cache_state(self.ctx)
        self._first_kind_periods[state] = self._copy_first_kind_periods(first)
        return second

    def riemann_matrix(self):
        r"""Return the normalized Riemann matrix of a plane algebraic curve.

        This is a convenience wrapper returning
        ``curve.periods_kind_1().tau``; see
        :meth:`Curve.periods_kind_1` for the input conventions.

            >>> from generapy import Curve
            >>> from mpmath import mp
            >>> mp.dps = 15
            >>> tau = Curve({(0, 2): 1, (1, 0): 1, (3, 0): -1}).riemann_matrix()
            >>> mp.im(tau[0, 0])
            mpf('1.0')
        """
        return self.periods_kind_1().tau

    def riemann_constant(self, *, base_place=None):
        r"""Return the vector of Riemann constants for a plane curve.

        The returned ``CurveRiemannConstant`` contains a direct representative of
        the normalized Jacobian vector ``value`` in mpmath's additive convention
        ``theta(A(D) + value, tau) = 0``, its literal ``(a, b)`` coordinates
        ``value = tau*a + b`` modulo the period lattice, the requested
        ``base_place`` (``None`` denotes the engine's natural base), and the
        maximum sheet residual of the direct contour integrations.

        Uses the same first-kind basis and marking as :meth:`periods_kind_1`
        and :meth:`abel_map_kind_1`, bound at construction. No second-kind
        forms are required. For the geometric engine, the value is computed
        from the canonical polygon and level-two contour integrals; theta
        functions and characteristic searches are not used.
        ``base_place`` may be a regular finite or chart-backed place. Changing
        the base uses ``K_Q = K_P + (g-1) A_P(Q)`` in normalized coordinates.
        The companion Abel map uses the same base for compatible coordinates.
        See :ref:`custom-differential-bases` for a worked example.

        In genus one the answer is the odd half-period ``(1+tau)/2``::

            >>> from generapy import Curve
            >>> from mpmath import mp
            >>> mp.dps = 15
            >>> curve = Curve({(0, 2): 1, (1, 0): 1, (3, 0): -1})
            >>> forms = (lambda x, y: 1 / y,)
            >>> custom = Curve(curve.polynomial, differentials_kind_1=forms)
            >>> constant = custom.riemann_constant()
            >>> periods = custom.periods_kind_1()
            >>> mp.almosteq(constant.value[0], (1 + periods.tau[0, 0]) / 2)
            True

        The direct level-two integrations are substantially more expensive than
        ordinary periods, although their cost does not include an exponential
        characteristic enumeration.
        """
        return self._call(
            _operations.riemann_constant, self._differentials_kind_1,
            base_place=base_place, _basis=True,
        )

    def validate(self, result):
        r"""Validate a result record returned by the curve functions.

        ``result`` is one of ``CurveBranchLocus``, ``CurveMonodromy``,
        ``CurveGenus``, ``CurveHomology``, ``CurvePeriodsKind1``,
        ``CurvePeriodsKind2`` or ``CurveRiemannConstant``. The returned
        ``CurveValidation`` record contains one named ``CurveCheck`` per
        invariant, the largest numerical residual among them, and whether
        every check passed.

        Structural checks are recomputed from the record; numerical
        integration residuals use recorded values. Tolerances use the
        current context precision; curve data is not recomputed.

            >>> from generapy import Curve
            >>> curve = Curve({(0, 2): 1, (1, 0): 1, (3, 0): -1})
            >>> report = curve.validate(curve.periods_kind_1())
            >>> report.passed
            True
            >>> report.checks[0]
            CurveCheck(name='tau_symmetry_residual', value=mpf('0.0'), passed=True)
        """
        self._check_precision()
        return _operations.validate(self.ctx, result)

    def fibre(self, x):
        r"""Return the labelled fibre of a plane algebraic curve over x.

        ``x`` must be a finite regular value of the ``x`` projection:
        it must not be a branch value or a value over which the projection
        drops degree. The returned tuple contains one ``CurvePlace`` per
        sheet, ordered deterministically by the real and imaginary parts
        of ``y``. Sheet labels refer to this particular fibre. At the
        monodromy base point they agree with the labels used by
        :attr:`Curve.monodromy`.

        >>> from generapy import Curve
            >>> from mpmath import mp
        >>> mp.dps = 15
        >>> curve = Curve({(0, 2): 1, (1, 0): 1, (3, 0): -1})
        >>> [mp.nstr(place.y, 6) for place in curve.fibre(2)]
        ['-2.44949', '2.44949']
        """
        return self._call(_operations.fibre, x)

    def path(self, start, end):
        r"""Return a lifted path between two regular finite places.

        ``start`` and ``end`` are regular finite places, each given as a
        ``(x, y)`` pair, a ``CurvePlace`` from :meth:`Curve.fibre`, or
        a chart-backed place from :meth:`Curve.chart_place`.  A guarded
        polyline in the x-plane avoids the branch values and is lifted by
        numerical continuation between the places' affine junction points;
        chart tails are joined at those junctions.  The returned ``CurvePath``
        record contains an opaque curve and numerical-context identity, the
        endpoint places, the sheet index reached, the continuation record
        carrying the numerical routing data used by
        :meth:`Curve.integral`, and any joined chart tails.

        Both junction points must have distinct ``x`` values, and ``end`` must
        lie on the sheet reached by continuation; otherwise ``ValueError`` is
        raised.

        >>> from generapy import Curve
            >>> from mpmath import mp
        >>> mp.dps = 15
        >>> curve = Curve({(0, 2): 1, (1, 0): -1})
        >>> path = curve.path((1, 1), (4, 2))
        >>> mp.nstr(path.start.y, 6), mp.nstr(path.end.y, 6)
        ('1.0', '2.0')
        """
        return self._call(_operations.path, start, end)

    def integral(self, differentials, path):
        r"""Integrate one differential or a differential basis along a path.

        ``differentials`` is either a single callable ``f(x, y)`` returning the
        coefficient of ``dx``, or a sequence of such callables; ``path`` is a
        ``CurvePath`` from :meth:`Curve.path`.  A single differential
        gives a scalar ``values`` entry, a sequence gives one entry per form.
        Chart tails joined to the path are integrated through their coordinate
        maps with the same differentials.  The returned ``CurveIntegral`` record
        also carries the maximum curve-equation residual encountered on the
        integration nodes and the number of path segments.  A path is bound to
        the curve and working precision at which it was constructed and cannot
        be reused with a different curve or precision.

        Unlike period and Abel-map bases, the sequence need not contain one
        holomorphic form per genus; its integrals only need to converge along
        this path. See :ref:`custom-differential-bases` for the callable convention.

        >>> from generapy import Curve
            >>> from mpmath import mp
        >>> mp.dps = 15
        >>> curve = Curve({(0, 2): 1, (1, 0): -1})
        >>> path = curve.path((1, 1), (4, 2))
        >>> integral = curve.integral(lambda x, y: 1 / y, path)
        >>> mp.nstr(integral.values, 12)
        '2.0'
        """
        return self._call(_operations.integral, differentials, path)

    def abel_map_kind_1(self, target, *, base_place=None,
                 reduce=False):
        r"""Evaluate the Abel map of a place or divisor on a plane curve.

        Generalizes elliptic integrals to genus :math:`g`, integrating :math:`g`
        differentials simultaneously from a base place to a target. Returns a
        ``CurveAbelMapKind1`` record whose ``value`` is an unnormalized vector
        in :math:`\mathbb{C}^g`, defined modulo its period lattice.

        For a base place :math:`Q` and an effective divisor
        :math:`D=P_1+\cdots+P_n`, the components are

        .. math::

            u_i(D;Q)=\sum_{\ell=1}^{n}\int_Q^{P_\ell}du_i,
            \qquad i=1,\ldots,g.

        Changing the paths by cycles changes this vector by
        :math:`2\omega m+2\omega'n`, with :math:`m,n\in\mathbb Z^g`.
        For a single point on :math:`y^2=x^3-x`, this is the elliptic
        integral

        .. math::

            u(P;Q)=\int_Q^P\frac{dx}{y}.

        The returned coordinates use this differential basis. Normalized
        Jacobian coordinates are :math:`(2\omega)^{-1}u`, with lattice
        :math:`\mathbb Z^g+\tau\mathbb Z^g`.

        ``target`` is a single place (as ``(x, y)`` pair or ``CurvePlace``),
        a chart-backed place, or a sequence of places (effective divisor). An
        empty sequence returns a record with a zero ``value`` vector.
        Uses the first-kind basis bound at construction, with exactly one
        holomorphic coefficient-of-``dx`` callable per genus. It shares this
        basis and marking with :meth:`periods_kind_1` and :meth:`riemann_constant`.
        For chart-backed endpoints on a recognized hyperelliptic curve, bind
        a custom basis to select geometric-polygon marking.
        See :ref:`custom-differential-bases` for a worked example.

        ``base_place`` defaults to the selected engine's natural base: the
        point at infinity for odd-degree hyperelliptic models, the first
        ordered finite branch point for even-degree hyperelliptic models, or
        sheet zero over the computational base point for the geometric-polygon
        engine.
        With ``reduce=True``, the result is reduced modulo the period lattice
        ``[2*omega, 2*omega_prime]`` of the selected first-kind basis.
        ``reduction_shift`` records the integer cycle shift, or is ``None``
        when reduction is disabled. ``engine`` and ``marking`` identify the
        calculation route and cycle convention.

        >>> from generapy import Curve
            >>> from mpmath import mp
        >>> mp.dps = 15
        >>> curve = Curve({(0, 2): 1, (1, 0): 1, (3, 0): -1})
        >>> point = (mp.mpf(2), mp.sqrt(6))
        >>> forms = (lambda x, y: 1 / y,)
        >>> custom = Curve(curve.polynomial, differentials_kind_1=forms)
        >>> value = custom.abel_map_kind_1(point, base_place=point).value
        >>> mp.nstr(mp.norm(value), 3)
        '0.0'
        """
        return self._call(
            _operations.abel_map, target, self._differentials_kind_1,
            _basis=True,
            base_place=base_place,
            reduce=reduce,
        )

    def abel_map_kind_2(
            self, target, *, base_place=None, reduce=False):
        r"""Evaluate second-kind Abelian integrals of a place or divisor.

        ``target`` and ``base_place`` have the same meanings as in :meth:`abel_map_kind_1`.
        Returns a ``CurveAbelMapKind2`` record with the second-kind ``value``.

        For a second-kind basis :math:`dr_1,\ldots,dr_g`, a base place
        :math:`Q`, and an effective divisor :math:`D=P_1+\cdots+P_n`,
        the components are

        .. math::

            r_i(D;Q)=\sum_{\ell=1}^{n}\int_Q^{P_\ell}dr_i,
            \qquad i=1,\ldots,g.

        For :math:`y^2=x^3-x`, the automatic BEL form gives, along a path
        between finite regular places avoiding its pole at infinity,

        .. math::

            r(P;Q)=\int_Q^P\frac{x\,dx}{4y}.

        This is a meromorphic Abelian integral rather than a holomorphic
        coordinate on the Jacobian. Adding cycles changes it by
        :math:`-2\eta m-2\eta'n`, using the signs in
        :meth:`periods_kind_2`. If first-kind reduction subtracts the cycle
        periods indexed by :math:`m,n`, the corresponding values satisfy

        .. math::

            \begin{aligned}
            u_{\mathrm{reduced}} &= u-2\omega m-2\omega'n, \\
            r_{\mathrm{reduced}} &= r+2\eta m+2\eta'n.
            \end{aligned}

        Uses the first- and second-kind bases bound at construction. With
        neither supplied, recognized hyperelliptic models use automatic BEL
        construction. If no compatible second-kind basis is available, raises
        ``ValueError`` explaining how to construct a curve with one.
        Supplied forms must have zero residues and be regular along the path.
        With ``reduce=True``, the bound first-kind Abel map determines the cycle
        shift applied to the second-kind value, recorded as ``reduction_shift``.
        See :ref:`custom-differential-bases` for definitions and an example.

        Chart-backed endpoints use the geometric polygon engine. Chart integration
        paths must converge; pole regularization is not performed.
        """
        self._require_second_kind_basis("abel_map_kind_2()")
        return self._call(
            _operations.abel_map, target, self._differentials_kind_1,
            _basis=True,
            second_kind=True,
            second_differentials=self._differentials_kind_2,
            base_place=base_place,
            reduce=reduce,
        )

    def lattice_reduce(self, value, periods):
        r"""Reduce a Jacobian vector modulo the period lattice.

        ``value`` is a genus-length column vector of Abelian coordinates. If
        ``periods`` is a ``CurvePeriodsKind1`` record, ``value`` is in the
        original differential basis and is reduced by the full lattice
        ``[2*omega, 2*omega_prime]``.  If ``periods`` is a normalized Riemann
        matrix ``tau``, ``value`` is in normalized coordinates and is reduced
        by ``[I, tau]``.  The returned ``CurveLatticeReduction`` record contains
        the equivalent vector and the integer lattice shift ``(m, n)``.

        >>> from generapy import Curve
            >>> from mpmath import mp
        >>> mp.dps = 15
        >>> curve = Curve({(0, 2): 1, (1, 0): 1, (3, 0): -1})
        >>> tau = curve.riemann_matrix()
        >>> reduced = curve.lattice_reduce(mp.matrix([2 + 1j]), tau)
        >>> reduced.shift
        (2, 1)
        >>> mp.nstr(reduced.value[0, 0], 3)
        '0.0'
        """
        self._check_precision()
        return _operations.lattice_reduce(self.ctx, value, periods)

    def chart(self, chart_curve, coordinate_map):
        r"""Return a user-supplied local chart of a plane algebraic curve.

        ``chart_curve`` gives the local curve as a sparse mapping from
        ``(t_power, w_power)`` pairs to coefficients.  ``coordinate_map(t, w)``
        must return the ambient triple ``(x, y, dx/dt)``, where ``x`` depends
        on ``t`` alone.  The returned ``CurveChart`` is bound to the ambient
        curve and working precision, and is accepted by the other chart
        methods and by :meth:`Curve.chart_place`.
        """
        return self._call(charts.chart, chart_curve, coordinate_map)

    def monomial_chart(self, x_power, y_power, *, source=None):
        r"""Return the monomial chart ``x = t**x_power, y = t**y_power*w``.

        ``source`` defaults to this curve; pass another ``CurveChart`` to
        compose coordinate maps. Negative powers describe places above
        infinity. Repeated factors are cleared to give a polynomial in
        ``t`` and ``w``; this does not normalize a singular chart.
        """
        self._check_precision()
        if source is None:
            source = self._polynomial
        return charts.monomial_chart(
            self.ctx, source, x_power, y_power)

    def chart_fibre(self, chart, t):
        r"""Return the ordered fibre of chart ``w`` values over ``t``.

        The values are ordered by real and imaginary part, like
        :meth:`Curve.fibre`.  A fibre whose values do not separate
        indicates that the chart does not resolve the requested place and is
        rejected.
        """
        self._check_precision()
        return charts.chart_fibre(self.ctx, chart, t)

    def chart_place(self, chart, seed, cutoff):
        r"""Return the chart-backed place reached by a local branch.

        ``seed`` is the branch value of ``w`` at ``t = 0``, for example from
        :meth:`Curve.chart_fibre`; the branch is continued along the
        straight chart path from ``t = 0`` to ``t = cutoff``.  The returned
        place is represented by its finite affine cutoff point together with a
        chart tail describing the local branch, and is bound to ``curve`` and
        the working precision.  The chart must parametrize ``curve``: the
        cutoff point is checked to lie on the curve.

        >>> from generapy import Curve
            >>> from mpmath import mp
        >>> mp.dps = 15
        >>> curve = Curve({(0, 2): 1, (1, 0): 1, (3, 0): -1})
        >>> chart = curve.monomial_chart(-2, -3)
        >>> [mp.nstr(value, 3) for value in curve.chart_fibre(chart, 0)]
        ['(-1.0 + 0.0j)', '(1.0 + 0.0j)']
        >>> place = curve.chart_place(chart, 1, mp.mpf("0.05"))
        >>> mp.nstr(place.x, 6)
        '400.0'
        """
        return self._call(
            charts.chart_place, chart, seed, cutoff)

    def chart_integral(self, chart, differentials, t_path, seed):
        r"""Integrate ambient differentials along a local chart branch.

        ``differentials`` are ambient ``f(x, y)`` callables returning the
        coefficient of ``dx``; they are pulled back through the chart's
        coordinate map, so a single callable gives a scalar and a sequence
        gives one entry per form.  ``t_path`` is a sequence of finite ``t``
        values along which the branch is continued from ``seed`` at
        ``t_path[0]``.  Closed chart loops therefore compute residues of
        pulled-back forms at the place.

        Continuation refines segments using local sheet separation and
        full-step/half-step consistency checks. Unresolved continuation or
        inconsistent quadrature branch labels raise ``ctx.NoConvergence``.
        These numerical checks do not certify that the supplied path avoids
        all critical values of the chart projection.
        """
        self._check_precision()
        return charts.chart_integral(
            self.ctx, chart, differentials, t_path, seed)

    def __repr__(self):
        return (
            f"Curve(x_degree={self.x_degree}, "
            f"y_degree={self.y_degree}, "
            f"ctx.prec={self._creation_state[0]})"
        )


CurveBranchLocus = _records.CurveBranchLocus
CurvePeriodsKind1 = _records.CurvePeriodsKind1
CurvePeriodsKind2 = _records.CurvePeriodsKind2
CurveAbelMapKind1 = _records.CurveAbelMapKind1
CurveAbelMapKind2 = _records.CurveAbelMapKind2
CurveMonodromy = _records.CurveMonodromy
CurveGenus = _records.CurveGenus
CurveHomology = _records.CurveHomology
CurveRiemannConstant = _records.CurveRiemannConstant
CurveValidation = _records.CurveValidation
CurveCheck = _records.CurveCheck
CurvePlace = _records.CurvePlace
CurveChart = _records.CurveChart
CurvePath = _records.CurvePath
CurveIntegral = _records.CurveIntegral
CurveLatticeReduction = _records.CurveLatticeReduction


__all__ = [
    "Curve",
    "CurveBranchLocus",
    "CurveChart",
    "CurveCheck",
    "CurvePeriodsKind1",
    "CurveGenus",
    "CurveHomology",
    "CurveIntegral",
    "CurveLatticeReduction",
    "CurveMonodromy",
    "CurvePath",
    "CurvePlace",
    "CurveRiemannConstant",
    "CurveAbelMapKind1",
    "CurveAbelMapKind2",
    "CurvePeriodsKind2",
    "CurveValidation",
]
