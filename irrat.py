from sympy import Add, Integral, Poly, Pow, pprint, sqrt, Symbol, simplify, \
     diff, solve
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
                if f not in _sqrt:
                    _sqrt.append(f)
                exit
            else:
                __sqrt_parse(f.args)

    __sqrt_parse(integral.args)
    if len(_sqrt) != 1:  # Корень не найден или их несколько
        raise Warning("None or more than one square roots found in expression")
    _sqrt = _sqrt[0]                                              # 0  1  2
    coeffs = Poly(_sqrt**2).all_coeffs()  # Коэфиценты многочлена: [a, b, c]
    assert(len(coeffs) == 3)
    x = integral.args[1][0]
    t = Symbol('t')

    if coeffs[0] > 0:  # a > 0
        _x = (coeffs[2] - t**2)/(2 * sqrt(coeffs[0])
                                   * t - coeffs[1])  # x через t
        tx = _sqrt - sqrt(coeffs[0]) * x  # t через x
        dx = simplify(diff(_x, t))  # dx через t
        sqrtt = sqrt(coeffs[0]) * _x + t  # корень через t
    elif coeffs[2] > 0:  # c > 0
        _x = None
        tx = None
        dx = None
        sqrtt = None
        ...
    else:  # third case
        sols = solve(_sqrt**2)
        if sols:
            _x = None
            tx = None
            dx = None
            sqrtt = None
            ...
        else:
            raise Warning("Euler substitution is impossible. "
                          "Try different integration methods.")
    return Integral(simplify(integral.args[0].subs({_sqrt: sqrtt, x: _x}) *
                    dx), (t, integral.args[1][1:])).doit().subs(t, tx)


def tan_ha(*args, **kwargs) -> Add:
    """Integrate function using Tangent half-angle substitution"""
    ...  # TODO: Реализация


def diff_binomial(*args, **kwargs) -> Add:
    """Integrate function as Differential binomial"""
    ...  # TODO: Реализация


Integral.euler = euler

a = Integral(x*sqrt(x**2 - 2*x + 2), x)
# pprint(a, use_unicode=False)
print()
pprint(a.euler(), use_unicode=False)
