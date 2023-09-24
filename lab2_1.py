import numpy as np


class Poly():
    def __init__(self, coefs: list):
        self.coefs = coefs

    def diff(self):
        b = np.array([len(self.coefs) - i for i in range(1, len(self.coefs) + 1)])
        return (b * self.coefs).tolist()

    def call(self, x0):
        res = 0
        for index, coef in enumerate(self.coefs[::-1]):
            res += coef * x0 ** index
        return res

    # Метод Ньютона, если C=1
    # Метод Ньютона-Бройдена, если C = const
    def newton(self, x0, e, C):
        polinom = Poly(self.coefs)
        xn = x0
        while True:
            f_xn = polinom.call(xn)
            df_xn = Poly(polinom.diff()).call(xn)
            if abs(f_xn) < e:
                return xn
            xn = xn - C * (f_xn / df_xn)

    def simple_newton(self, x0, e):
        polinom = Poly(self.coefs)
        f_x0 = polinom.call(x0)
        df_x0 = Poly(polinom.diff()).call(x0)
        x1 = x0 - f_x0 / df_x0
        xn = x1
        while True:
            f_xn = polinom.call(xn)
            if abs(f_xn) < e:
                return xn
            xn = xn - f_xn / df_x0

    def secant(self, x0, e, delta):
        polinom = Poly(self.coefs)
        f_x0 = polinom.call(x0)
        f_x00 = polinom.call(x0 - delta)
        df_x0 = (f_x0 - f_x00) / delta
        xn_ = x0
        xn = x0 - f_x0 / df_x0
        while True:
            f_xn = polinom.call(xn)
            f_xn_ = polinom.call(xn_)
            if abs(xn - xn_) < e:
                return xn
            xn_plus = xn - f_xn * (xn - xn_) / (f_xn - f_xn_)
            xn_ = xn
            xn = xn_plus


def main():
    coefs = [3, -4, -8, 10, -7]
    print('Метод Ньютона:', Poly(coefs).newton(x0=0.9, e=0.001, C=0.001))
    coefs = [1, 2, -3]
    print('Упрощенный метод Ньютона:', Poly(coefs).simple_newton(x0=0.8, e=0.001))
    coefs = [1, 0, -1, 1]
    print('Метод секущих:', Poly(coefs).secant(x0=-2, e=0.001, delta=0.1))
#   output:
#   Метод Ньютона: 2.0993841138798923
#   Упрощенный метод Ньютона: 1.0002295383516655
#   Метод секущих: -1.324717987082104


if __name__ == '__main__':
    main()
