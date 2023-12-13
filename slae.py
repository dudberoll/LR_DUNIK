import numpy as np

def gauss(mat1,mat2):
    mat = np.concatenate((mat1,mat2), axis=1)
    num_rows = mat.shape[0]
    # Прямой ход
    for i in range(num_rows):
        # Деление на ведущие элементы
        value = mat[i, i]
        mat[i] = mat[i] / value
        # исключение элементов aij столбцов, лежащих ниже опорных строк матрицы А*
        for k in range(i + 1, num_rows):
            mat[k] -= mat[i] * mat[k, i]
    # Разбиваем на матрицу при иксах и матрицу ответов
    mat1 = mat[:, :-1]
    mat2 = mat[:, -1]
    # Решаем систему уравнений
    X = np.linalg.solve(mat1, mat2)
    return X

def rectangle(mat1, mat2):
    mat = np.concatenate((mat1,mat2),axis=1)
    num_rows = mat.shape[0]
    num_cols = mat.shape[1]
    # mat_next = np.zeros((num_rows,num_cols))

    for k in range(num_rows):
        # Делим о.с на ведущий элемент
        mat[k] = mat[k]/mat[k,k]
        mat_next = mat.copy()
        # Считаем новые значения матрица по правилу прямоугольника
        for i in range(k+1,num_rows):
            for j in range(k,num_cols):
                mat_next[i,j] = mat[i,j] - (mat[i,k] * mat[k,j])/mat[k,k]
        mat = mat_next

    mat1 = mat[:, :-1]
    mat2 = mat[:, -1]
    X = np.linalg.solve(mat1, mat2)
    return mat,X


def LeadingElement(mat1,mat2):
    matrix = np.concatenate((mat1,mat2), axis=1)
    n = len(matrix)

    for i in range(n):
        max_element = abs(matrix[i][i])
        max_row = i
        for k in range(i + 1, n):
            if abs(matrix[k][i]) > max_element:
                max_element = abs(matrix[k][i])
                max_row = k

        matrix[i], matrix[max_row] = matrix[max_row], matrix[i]

        for k in range(i + 1, n):
            factor = matrix[k][i] / matrix[i][i]
            for j in range(i, n + 1):
                if i == j:
                    matrix[k][j] = 0
                else:
                    matrix[k][j] -= factor * matrix[i][j]

    mat1 = matrix[:, :-1]
    mat2 = matrix[:, -1]
    solution = np.linalg.solve(mat1,mat2)

    # solution = [0 for i in range(n)]
    # for i in range(n - 1, -1, -1):
    #     solution[i] = matrix[i][n] / matrix[i][i]
    #     for k in range(i - 1, -1, -1):
    #         matrix[k][n] -= matrix[k][i] * solution[i]

    return solution
def simple_iterations(A, X , e):
    max_iter  = 10
    X_next = X
    x = X
    for _ in range(max_iter):
        X_next = A.dot(X_next) + x

        if np.linalg.norm(X_next - X) <= e:
            return X_next.T

        X = X_next

def gauss_seidel(A, b, eps=1e-4, max_iter=1000):
    n = len(A)
    x = np.zeros(n, dtype=np.float64)

    for _ in range(max_iter):
        x_next = np.copy(x)
        for i in range(n):
            s1 = np.dot(A[i, :i], x_next[:i])
            s2 = np.dot(A[i, i + 1:], x[i + 1:])
            x_next[i] = (b[i] - s1 - s2) / A[i, i]

        if np.linalg.norm(x_next - x) < eps:
            return x_next
        x = x_next

def gmres(A, b, x0=None, tol=1e-5, max_iter=100):
    n = len(b)
    if x0 is None:
        x0 = np.zeros(n)

    x = np.copy(x0)
    r = b - np.dot(A, x)
    beta = np.linalg.norm(r)
    q = [r / beta]

    H = np.zeros((max_iter + 1, max_iter))
    g = np.zeros(max_iter + 1)

    for k in range(max_iter):
        # Arnoldi ортогонализация
        v = np.dot(A, q[k])
        for j in range(k + 1):
            H[j, k] = np.dot(q[j], v)
            v = v - H[j, k] * q[j]

        H[k + 1, k] = np.linalg.norm(v)
        q.append(v / H[k + 1, k])

        # Решение подзадачи наименьших квадратов
        e1 = np.zeros(k + 1)
        e1[0] = 1.0
        y = np.linalg.lstsq(H[:k + 1, :k + 1], beta * e1, rcond=None)[0]
        g[:k + 1] = y

        # Обновление приближения
        x = x0 + np.dot(q[0], g[0])
        for j in range(1, k + 1):
            x = x + np.dot(q[j], g[j])

        # Проверка на сходимость
        if np.linalg.norm(b - np.dot(A, x)) < tol:
            return x
    return x


def lu_decomposition(A):
    n = len(A)
    L = [[0.0] * n for _ in range(n)]
    U = [[0.0] * n for _ in range(n)]

    for i in range(n):
        # Верхнетреугольная матрица U
        for k in range(i, n):
            total = sum(L[i][j] * U[j][k] for j in range(i))
            U[i][k] = A[i][k] - total

        # Нижнетреугольная матрица L
        for k in range(i, n):
            if i == k:
                L[k][i] = 1.0
            else:
                total = sum(L[k][j] * U[j][i] for j in range(i))
                L[k][i] = (A[k][i] - total) / U[i][i]

    return L, U

def solve_lu(L, U, b):
    n = len(L)
    y = [0.0] * n
    x = [0.0] * n

    # Решение Ly = b
    for i in range(n):
        y[i] = b[i] - sum(L[i][j] * y[j] for j in range(i))

    # Решение Ux = y
    for i in range(n - 1, -1, -1):
        x[i] = (y[i] - sum(U[i][j] * x[j] for j in range(i + 1, n))) / U[i][i]

    return x


def jacobi_method(A, b, max_iterations=100, tolerance=1e-6):
    n = len(b)
    x = np.zeros(n)  # начальное приближение

    for iteration in range(max_iterations):
        x_new = np.zeros_like(x)

        for i in range(n):
            sigma = sum(A[i, j] * x[j] for j in range(n) if j != i)
            x_new[i] = (b[i] - sigma) / A[i, i]

        if np.linalg.norm(x_new - x, ord=np.inf) < tolerance:
            return x_new

        x = x_new.copy()

    return x





