"""Shared arbitrary-precision Runge--Kutta integration for examples."""


def _add_scaled(left, right, scale):
    """Return ``left + scale*right`` componentwise."""
    return tuple(value + scale * increment
                 for value, increment in zip(left, right))


def rk4_step(rhs, state, step, *args):
    """Advance one classical fourth-order Runge--Kutta step.

    ``rhs(state, *args)`` may return any finite iterable. The returned state is
    a tuple and retains the numeric types supplied by the caller, including
    mpmath arbitrary-precision real and complex values.
    """
    state = tuple(state)
    k1 = tuple(rhs(state, *args))
    k2 = tuple(rhs(_add_scaled(state, k1, step / 2), *args))
    k3 = tuple(rhs(_add_scaled(state, k2, step / 2), *args))
    k4 = tuple(rhs(_add_scaled(state, k3, step), *args))
    return tuple(
        value + step * (d1 + 2 * d2 + 2 * d3 + d4) / 6
        for value, d1, d2, d3, d4 in zip(state, k1, k2, k3, k4)
    )


def rk4_steps(rhs, state, step, count, *args):
    """Advance ``count`` equal RK4 steps and return the final state."""
    if count < 0:
        raise ValueError("count must be nonnegative")
    state = tuple(state)
    for unused in range(count):
        state = rk4_step(rhs, state, step, *args)
    return state


def rk4_trajectory(rhs, initial, sample_points, *, substeps=1, args=()):
    """Integrate between ordered sample points with equal substeps per gap."""
    sample_points = tuple(sample_points)
    if not sample_points:
        return ()
    if substeps <= 0:
        raise ValueError("substeps must be positive")
    state = tuple(initial)
    result = [state]
    for left, right in zip(sample_points, sample_points[1:]):
        step = (right - left) / substeps
        state = rk4_steps(rhs, state, step, substeps, *args)
        result.append(state)
    return tuple(result)
