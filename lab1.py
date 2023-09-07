import math

def f(x: int):
    res = ((x - 3) ** 2 * (x - 2)) * (x - 7)
    return res
def iter(x1,x0):
    return (x1-x0)/2
# def start_approx():
#
def dih(x0, x1, e, max_iter):

    xn = (x1 + x0)/2
    while max_iter != 0 and f(xn) != 0:
        max_iter -= 1
        if math.fabs(f(xn)) >= e:
            x1 = math.fabs((x1-xn))/2
            x0 = math.fabs((x0-xn))/2
            xn = iter(x1,x0)
        else:
            return xn
        return xn


def main():
    print(dih(0,2.5,0.01,30))

if __name__ == '__main__':
    main()
