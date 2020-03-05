from sympy import Add, Integral, Poly, Pow, pprint, sqrt, Symbol
from sympy.abc import x
from sympy.core.numbers import Half


def euler(self: Integral) -> Add:
    """Implements Euler substitution for irrational integrals"""
    _sqrt = []  # Искомый корень для замены

    def __sqrt_parse(functions_tuple):
        for f in functions_tuple:
            if type(f) is Pow and type(f.args[0]) is Add \
                    and type(f.args[1]) is Half:
                _sqrt.append(f)
                exit
            else:
                __sqrt_parse(f.args)

    __sqrt_parse(self.args)
    if not _sqrt:  # Корень не найден
        raise Warning("Square root not found in expression")
    _sqrt = _sqrt[0]
    _coeffs = Poly(_sqrt**2).all_coeffs()  # Коэфиценты многочлена под корнем
    assert(len(_coeffs) == 3)
    _t = Symbol('t')

    if _coeffs[0] > 0:  # a > 0
        ...
    elif _coeffs[2] > 0:  # c > 0
        ...
    else:  # third case
        ...
    # return Integral(self.args[0].subs(x, _t),
    #                 (_t, self.args[1][1:]))  # .doit().subs(_t, x)
    return _sqrt  # TODO: Реализация подстановок


def diff_binomial(self: Integral) -> Add:
    ...  # TODO: Реализация


Integral.euler = euler

a = Integral(1/(x + sqrt(x**2 + x + 1)), x)
pprint(a)
print()
pprint(a.euler())
