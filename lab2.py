def next_value(xn, a):
    x_next = (xn + a / xn) * 1 / 2
    return x_next


def simple_iteration(x0, e, max_iter, a):
    x1 = next_value(x0, a)
    for i in range(max_iter):
        x1 = next_value(x1, a)
    return x1


def main():
    print(simple_iteration(2, 0.01, 10, 20))


#     output : 4.47213595499958

if __name__ == '__main__':
    main()
