from sympy import Add, Integral, Poly, Pow, pprint, sqrt, Symbol, simplify, \
     diff, solve, sin, cos, tan, cot
from sympy.abc import x
from sympy.core.numbers import Half


def euler(*args, **kwargs) -> Add:
    """Integrate irrational function using Euler substitution"""
    integral = Integral(*args, **kwargs)
    _sqrt = []  # Искомый корень для замены

    def __sqrt_parse(functions_tuple):
        for f in functions_tuple:
            if type(f) is Pow and type(f.args[0]) is Add \
                    and type(f.args[1]) is Half:
                _sqrt.append(f)
                exit
            else:
                __sqrt_parse(f.args)

    __sqrt_parse(integral.args)
    if not _sqrt:  # Корень не найден
        raise Warning("Square root not found in expression")
    _sqrt = _sqrt[0]                                              # 0  1  2
    coeffs = Poly(_sqrt**2).all_coeffs()  # Коэфиценты многочлена: [a, b, c]
    assert(len(coeffs) == 3)
    x = integral.args[1][0]
    _t = Symbol('t')

    if coeffs[0] > 0:  # a > 0
        _x = (coeffs[2] - _t**2)/(2 * sqrt(coeffs[0])
                                    * _t - coeffs[1])  # x через t
        tx = _sqrt - sqrt(coeffs[0]) * x  # t через x
        dx = simplify(diff(_x, _t))  # dx через t
        sqrtt = sqrt(coeffs[0]) * _x + _t  # корень через t
    elif coeffs[2] > 0:  # c > 0
        _x = (coeffs[1]-2*sqrt(coeffs[2])*_t)/(_t**2 - coeffs[0])  # x через t
        tx = (_sqrt-sqrt(coeffs[2]))/x   # t через x
        dx = simplify(diff(_x, _t))  # dx через t
        sqrtt = _x*_t + sqrt(coeffs[2])  # корень через t
    else:  # third case   # sqrt(a*(x-x1)*(x-x2))
        sols = solve(_sqrt**2)
        if len(sols) == 2:
            _x = (_t**2 * sols[0] - coeffs[0] *
                  sols[1])/(_t**2 - coeffs[0])  # x через t
            tx = _sqrt/(x - sols[0])  # t через x
            dx = simplify(diff(_x, _t))  # dx через t
            sqrtt = _t*(_x - sols[0])  # корень через t
        else:
            raise Warning("Euler substitution is impossible. "
                          "Try different integration methods.")
    return Integral(simplify(integral.args[0].subs({_sqrt: sqrtt, x: _x})*dx),
                    (_t, integral.args[1][1:]), **kwargs).doit().subs(_t, tx)


def tan_ha(*args, **kwargs) -> Add:
    """Integrate function using Tangent half-angle substitution"""
    integral = Integral(*args, **kwargs)
    x = integral.args[1][0]
    _t = Symbol('t')
    _sin = (2 * _t)/(1 + _t**2)
    _cos = (1 - _t**2)/(1 + _t**2)
    _tg = _sin/_cos
    _ctg = 1/_tg
    _dx = 2/(1 + _t**2)
    _f = simplify(args[0].subs({sin(x): _sin, cos(x): _cos, tan(x): _tg,
                                cot(x): _ctg}) * _dx)
    return simplify(Integral(_f, (_t, integral.args[1][1:]),
                             **kwargs).doit().subs(_t, tan(x/2)))


def diff_binomial(*args, **kwargs) -> Add:
    """Integrate function as Differential binomial"""
    ...  # TODO: Реализация


Integral.euler = euler

# a = Integral(1/(2*sin(x) - cos(x) + 5), x)
# a = Integral((sin(x)**2/(1 + sin(x)**2)), x)
# a = Integral(1/(1 + sqrt(-1 + 50*x - x**2)), x)
a = Integral(1/(x + sqrt(x**2 + x + 1)), x)
# pprint(a, use_unicode=False)
print()
# pprint(tan_ha(*a.args))
pprint(a.euler(), use_unicode=False)
