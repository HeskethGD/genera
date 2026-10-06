from mpmath import mp

from generapy._context import ctx_lru_cache, resolve_context


class FakeContext:
    def __init__(self, prec=53, rounding=None, trap_complex=False):
        self.prec = prec
        self.rounding = rounding
        self.trap_complex = trap_complex


def test_resolve_context_uses_mpmath_by_default():
    assert resolve_context() is mp


def test_resolve_context_preserves_an_explicit_context():
    ctx = FakeContext()
    assert resolve_context(ctx) is ctx


def test_ctx_lru_cache_reuses_a_result_in_the_same_numerical_state():
    calls = []

    @ctx_lru_cache(maxsize=4)
    def evaluate(ctx, value, scale=1):
        calls.append((ctx.prec, value, scale))
        return ctx.prec + value * scale

    ctx = FakeContext()
    assert evaluate(ctx, 3, scale=2) == 59
    assert evaluate(ctx, 3, scale=2) == 59
    assert calls == [(53, 3, 2)]
    assert evaluate.cache_info(ctx).hits == 1
    assert evaluate.cache_info(ctx).misses == 1


def test_ctx_lru_cache_accepts_the_standard_mpmath_context():
    @ctx_lru_cache(maxsize=4)
    def evaluate(ctx, value):
        return ctx.sqrt(value)

    assert evaluate(mp, 4) == 2
    assert evaluate(mp, 4) == 2
    assert evaluate.cache_info(mp).hits == 1


def test_ctx_lru_cache_separates_numerical_states():
    calls = []

    @ctx_lru_cache()
    def evaluate(ctx):
        result = (ctx.prec, ctx.rounding, ctx.trap_complex)
        calls.append(result)
        return result

    ctx = FakeContext()
    assert evaluate(ctx) == (53, None, False)
    ctx.prec = 80
    assert evaluate(ctx) == (80, None, False)
    ctx.rounding = "up"
    assert evaluate(ctx) == (80, "up", False)
    ctx.trap_complex = True
    assert evaluate(ctx) == (80, "up", True)
    assert len(calls) == 4


def test_ctx_lru_cache_is_separate_for_each_context():
    calls = []

    @ctx_lru_cache()
    def evaluate(ctx, value):
        calls.append(ctx)
        return value

    first = FakeContext()
    second = FakeContext()
    assert evaluate(first, 1) == 1
    assert evaluate(second, 1) == 1
    assert evaluate(first, 1) == 1
    assert calls == [first, second]
    assert evaluate.cache_info(first).hits == 1
    assert evaluate.cache_info(second).hits == 0


def test_ctx_lru_cache_separates_internal_mpmath_rounding_modes():
    calls = []

    @ctx_lru_cache()
    def evaluate(ctx):
        calls.append(ctx._prec_rounding[1])
        return ctx._prec_rounding[1]

    ctx = mp.clone()
    assert evaluate(ctx) == 'n'
    ctx._prec_rounding[1] = 'f'
    assert evaluate(ctx) == 'f'
    ctx._prec_rounding[1] = 'n'
    assert evaluate(ctx) == 'n'
    assert calls == ['n', 'f']


def test_ctx_lru_cache_clear_can_target_one_context_or_all_contexts():
    calls = []

    @ctx_lru_cache()
    def evaluate(ctx):
        calls.append(ctx)
        return ctx.prec

    first = FakeContext()
    second = FakeContext()
    unused = FakeContext()
    evaluate(first)
    evaluate(second)

    evaluate.cache_clear(unused)
    evaluate.cache_clear(first)
    evaluate(first)
    evaluate(second)
    assert calls == [first, second, first]

    evaluate.cache_clear()
    evaluate(first)
    evaluate(second)
    assert calls == [first, second, first, first, second]


def test_ctx_lru_cache_obeys_its_size_bound():
    @ctx_lru_cache(maxsize=1)
    def evaluate(ctx, value):
        return value

    ctx = FakeContext()
    evaluate(ctx, 1)
    evaluate(ctx, 2)
    assert evaluate.cache_info(ctx).currsize == 1
    assert evaluate.cache_info(ctx).maxsize == 1


def test_ctx_lru_cache_rejects_an_invalid_size():
    try:
        ctx_lru_cache(-1)
    except ValueError as error:
        assert str(error) == "maxsize must be a nonnegative integer or None"
    else:
        raise AssertionError("negative cache sizes must be rejected")
