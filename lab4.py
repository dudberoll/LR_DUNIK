import numpy as np


def iterations(A,X,e):
    lymda = 0
    k = 0
    while True:
        k += 1
        X_next = A.dot(X)
        lymda_next = X_next[0]/X[0]

        if np.linalg.norm(lymda_next - lymda) < e:
            return X_next/max(X_next), lymda_next

        X = X_next
        lymda = lymda_next

def jacobi_rotation(A, tol=1e-9, max_iter=100):
    n = len(A)
    eigenvectors = np.identity(n)
    for _ in range(max_iter):
        max_val = 0
        for i in range(n):
            for j in range(i + 1, n):
                if abs(A[i, j]) > max_val:
                    max_val = abs(A[i, j])
                    p, q = i, j

        if max_val < tol:
            break

        theta = 0.5 * np.arctan2(2 * A[p, q], A[q, q] - A[p, p])
        c, s = np.cos(theta), np.sin(theta)

        J = np.identity(n)
        J[p, p] = c
        J[p, q] = -s
        J[q, p] = s
        J[q, q] = c

        A = np.dot(np.dot(J.T, A), J)
        eigenvectors = np.dot(eigenvectors, J)

    eigenvalues = np.diag(A)
    return eigenvalues, eigenvectors
def main():
    # A = np.matrix([[5,1,2],[1,4,1],[2,1,3]])
    # X = np.matrix([[1],[1],[1]])
    #
    # eigenvector, eigenvalue = iterations(A,X,0.1)
    # print('Собственный вектор:',eigenvector, '\n'
    #       'Собственное значение:',eigenvalue)

    A = np.array([[4, 1, 2],
                  [1, 5, 3],
                  [2, 3, 6]])
    X = np.matrix([[1], [1], [1]])

    eigenvector, eigenvalue = iterations(A, X, 0.1)
    print('Собственный вектор:', eigenvector, '\n'
                                              'Собственное значение:', eigenvalue)

    A = np.array([[5, 1, 2],
                  [1, 4, 1],
                  [2, 1, 3]])
    eigenvalues, eigenvectors = jacobi_rotation(A)
    print("Собственные значения:", eigenvalues)
    print("Собственные векторы:")
    print(eigenvectors)
# Собственные значения: [6.8950141  2.86564459 2.23934132]
# Собственные векторы:
# [[ 0.75423962 -0.64320093 -0.13196649]
#  [ 0.43301255  0.63833949 -0.63641404]
#  [ 0.49358152  0.42286554  0.75997501]]

if __name__ == '__main__':
    main()