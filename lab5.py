import numpy as np

def lagrange(x_values, y_values , x):
    value = 0.0
    for i in range(len(x_values)):
        L = y_values[i]
        for j in range(len(x_values)):
            if j != i:
                L = L * (x - x_values[j]) / (x_values[i] - x_values[j])
        value += L

    return value

def linear_lagrange(x_values,y_values,x):
    i = 0
    while x_values[i] < x:
        i += 1
    i = i - 1

    P1 = (x-x_values[i+1])/(x_values[i]-x_values[i+1])
    P2 = (x - x_values[i])/(x_values[i+1] - x_values[i])

    value = P1*y_values[i] + P2*y_values[i+1]
    return value


def parabolic_lagrange(x_values,y_values,x):
    # окно - O[Xi, X+1]
    i = 0
    while x_values[i] < x:
        i += 1
    i = i - 1

    # O1 - первое окно - O[Xi-1, Xi+1]
    n, m = i-1, i+2
    x_values_1 = x_values[n-1:m+1]
    y_values_1 = y_values[n-1:m+1]
    L1 = lagrange(x_values_1,y_values_1,x)
    # O2 - второе окно - O[Xi, Xi+2]
    k, s = i, i+2
    x_values_2 = x_values[k-1:s+1]
    y_values_2 = y_values[k-1:s+1]
    L2 = lagrange(x_values_2,y_values_2,x)

    return (L1+L2)/2

def main():
    # Многочлены Лагранжа
    # Задайте сетку
    x_values = np.array([2,3,4,5])
    y_values = np.array([7,5,8,7])
    print('Значение в точке:',lagrange(x_values,y_values,2.5))

    # Линейная интерполяция
    x_values = np.array([-1, 0, 1, 3, 4])
    y_values = np.array([-1, 0, 1, 27, 64])
    print('Значение в точке:', linear_lagrange(x_values, y_values, 2))

    # Параболическая интерполяция
    print('Значение в точке:', parabolic_lagrange(x_values, y_values, 2))
if __name__ == '__main__':
    main()