"""Exact, self-contained algebraic-curve catalog; no numerical imports."""
from fractions import Fraction as Q
from math import comb


def packed(terms):
    return [[i, j, str(Q(c))] for (i, j), c in sorted(terms.items()) if c]


def polynomial(roots):
    c = [Q(1)]
    for r in map(Q, roots):
        out = [Q(0)]*(len(c)+1)
        for i, a in enumerate(c):
            out[i] -= r*a
            out[i+1] += a
        c = out
    return c


def case(name, terms, genus, group, note='', numerators=None, expected_rejection=False):
    return dict(name=name, terms=packed(terms), genus=genus, group=group,
                note=note, numerators=numerators, expected_rejection=expected_rejection)


def hyper(name, roots, group='established'):
    c = polynomial(roots)
    return case(name, {(0, 2): 1, **{(i, 0): -v for i, v in enumerate(c)}},
                (len(roots)-1)//2, group, 'Specialized hyperelliptic dispatch')


def trig(mu):
    a,b,c,d,e,f,g,h,k = mu
    return {(0,3):1,(1,2):a,(0,2):d,(2,1):b,(1,1):e,(0,1):g,
            (4,0):-1,(3,0):-c,(2,0):-f,(1,0):-h,(0,0):-k}


def x_transform(terms, shift, scale):
    """Replace x by (X-shift)/scale using exact rational arithmetic."""
    out = {}
    for (i,j), c in terms.items():
        for k in range(i+1):
            out[k,j] = out.get((k,j), Q(0))+Q(c)*comb(i,k)*Q(-shift)**(i-k)/Q(scale)**i
    return out


def catalog():
    klein = {(3,1):1,(0,3):1,(1,0):1}
    general = trig((1,2,-3,4,-5,6,-7,8,-9))
    kova = {(2,4):1,(3,2):-4,(2,2):6,(1,2):-2,
            (2,0):Q(27,5),(1,0):-Q(26,5),(0,0):1}
    rows = [
        hyper('hyper-g1-lemniscatic', [-1,0,1]),
        hyper('hyper-g2-symmetric', [-2,-1,0,1,2]),
        hyper('hyper-g3-symmetric', [-3,-2,-1,0,1,2,3]),
        hyper('hyper-g2-even', [-2,-1,0,1,2,3]),
        hyper('hyper-g2-clustered', [-4,1,Q(9,8),3,6]),
        hyper('hyper-g3-generic', [2,3,5,7,11,13,17]),
        hyper('hyper-g4-symmetric', list(range(-4,5)), 'heavy'),
        case('general-g1-lemniscatic', {(0,2):1,(3,0):-1,(1,0):1},1,'diagnostic',
             'Explicit-basis override uses general geometry; automatic path stays hyperelliptic', [packed({(0,0):2})]),
        case('general-trig-g3',general,3,'established','Close branch pair; formerly hit panel depth limit'),
        case('general-trig-pure-g3',trig((0,0,2,0,0,3,0,5,7)),3,'established'),
        case('general-trig-dense-g3',trig((-3,-1,-9,6,-5,-4,-1,6,5)),3,'established',
             'Distinct from general-trig-g3'),
        case('general-klein-quartic-g3',klein,3,'established'),
        case('general-fermat-quartic-g3',{(4,0):1,(0,4):1,(0,0):-1},3,'established'),
        case('general-super-g4',{(0,3):1,(5,0):-1,(0,0):1},4,'heavy',
             'Actual equation y^3=x^5-1; fixes inconsistent prose in older catalog'),
        case('general-kovalevskaya-g3',kova,3,'established','Known basis required by generapy',
             [packed({(1,0):1}),packed({(1,1):1}),packed({(1,2):1,(0,0):-1})]),
        case('nonmonic-elliptic',{(1,2):1,(2,0):-1,(0,0):-1},1,'expanded',
             'w=x*y gives w^2=x^3+x'),
        case('general-trig-test-g3',{(0,3):1,(4,0):-1,(1,0):1,(0,0):-1},3,'expanded'),
        case('general-trig-infinity-g3',{(0,3):1,(4,0):-1,(0,0):1},3,'expanded'),
        case('klein-x-translated',x_transform(klein,10,1),3,'expanded','x -> X-10'),
        case('klein-x-scaled',x_transform(klein,0,8),3,'expanded','x -> X/8'),
        case('general-mu-x-translated',x_transform(general,5,1),3,'expanded'),
        hyper('hyper-g2-close-001',[-2,0,1,Q(101,100),3],'expanded'),
        hyper('hyper-g2-close-00001',[-2,0,1,Q(10001,10000),3],'stress'),
        case('hyper-g2-complex',{(0,2):1,(5,0):-1,(0,0):1},2,'expanded','Non-real branch points'),
        case('hyper-g1-linear-y',{(0,2):1,(1,1):2,(2,0):1,(3,0):-1,(1,0):1},1,
             'expanded','z=y+x gives z^2=x^3-x; specialized dispatch'),
    ]
    c = polynomial([1,2,3,4,4])
    rows.append(case('singular-g1-double-point',
                     {(0,2):1,**{(i,0):-v for i,v in enumerate(c)}},1,'negative',
                     'generapy should reject; Sage may normalize. Never a speed comparison.',
                     expected_rejection=True))
    return {r['name']: r for r in rows}


SMOKE = ('hyper-g1-lemniscatic','general-klein-quartic-g3')
QUICK = SMOKE + ('hyper-g2-symmetric','general-trig-g3','general-kovalevskaya-g3','nonmonic-elliptic')
