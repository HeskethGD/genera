"""Helpers for numerical functions implemented against an mpmath context."""

from functools import lru_cache, wraps
from threading import RLock
from weakref import WeakKeyDictionary, ref

from mpmath import mp


def resolve_context(ctx=None):
    """Return *ctx*, or the standard mpmath context when it is omitted."""
    return mp if ctx is None else ctx


def _context_state(ctx):
    """Return the numerical settings that can affect a cached calculation."""
    # Fixed-precision and interval contexts do not necessarily expose rounding
    # or trap_complex, so stable placeholders are used for missing settings.
    return (
        ctx.prec,
        getattr(ctx, "rounding", None),
        getattr(ctx, "trap_complex", None),
    )


def ctx_lru_cache(maxsize=128):
    """Cache a ``ctx``-first function separately for each numerical context.

    Cache entries include the context's numerical state, preventing a result
    calculated at one precision or rounding mode from being reused at another.

    Notes
    -----
    The mpmath development branch contains an equivalent context-aware cache,
    but it is not part of a released public API. Generapy temporarily owns this
    implementation. Once mpmath releases a supported cache utility, replace
    this implementation with the mpmath version.
    """
    if maxsize is not None and (not isinstance(maxsize, int) or maxsize < 0):
        raise ValueError("maxsize must be a nonnegative integer or None")

    def decorate(function):
        caches = WeakKeyDictionary()
        lock = RLock()

        def cache_for(ctx):
            with lock:
                cached = caches.get(ctx)
                if cached is not None:
                    return cached

                ctx_ref = ref(ctx)

                @lru_cache(maxsize=maxsize)
                def cached(_ctx_state, /, *args, **kwargs):
                    context = ctx_ref()
                    if context is None:  # pragma: no cover - guarded by call
                        raise RuntimeError("the numerical context no longer exists")
                    return function(context, *args, **kwargs)

                caches[ctx] = cached
                return cached

        @wraps(function)
        def wrapper(ctx, /, *args, **kwargs):
            return cache_for(ctx)(_context_state(ctx), *args, **kwargs)

        def cache_info(ctx):
            """Return cache statistics for *ctx*."""
            return cache_for(ctx).cache_info()

        def cache_clear(ctx=None):
            """Clear one context's cache, or every context cache if omitted."""
            with lock:
                if ctx is not None:
                    cached = caches.get(ctx)
                    if cached is not None:
                        cached.cache_clear()
                    return
                for cached in tuple(caches.values()):
                    cached.cache_clear()
                caches.clear()

        wrapper.cache_info = cache_info
        wrapper.cache_clear = cache_clear
        return wrapper

    return decorate
