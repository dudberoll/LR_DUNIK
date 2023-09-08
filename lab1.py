# Функция, для которой надо найти корень
def f(x):
    res = (x - 2) * (x - 7) * (x - 3) ** 2
    return res


# Задать левую и правую границы поиска, точность, максимальное число итераций
def dih(x0, x1, e, max_iter):
    x_left = x0
    x_right = x1
    for i in range(max_iter):
        # корень находится в интервале с заданной точностью
        if abs((x_right - x_left)) < e:
            break
        # перемещаем границу интервала
        # если находится в левом интервале, то меняем правую границу
        # елси находится справа, то меняем левую границу
        xn = (x_right + x_left) / 2
        if f(xn) * f(x_left) < 0:
            x_right = xn
        if f(xn) * f(x_right) < 0:
            x_left = xn
    return (x_right + x_left) / 2


def main():
    print(dih(0, 2.5, 0.001, 20))


if __name__ == '__main__':
    main()
