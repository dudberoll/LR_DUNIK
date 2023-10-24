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




def main():
    # Задание 1
    # Введите матрицу при иксах
    mat1 = np.matrix([[5, 0, 1], [2, 6, -2], [-3, 2, 10]]).astype(float)
    # Введите матрицу ответов
    mat2 = np.matrix([[11], [8], [6]]).astype(float)

    print(gauss(mat1,mat2))
#     output:
#           [[2.]
#            [1.]
#            [1.]]
    # Задание 2
    # Введите матрицу при иксах
    mat1 = np.matrix([[2, 1, 4], [3, 2, 1], [1, 3, 3]]).astype(float)
    # Введите матрицу ответов
    mat2 = np.matrix([[16], [10], [16]]).astype(float)

    final_matrix, solution = rectangle(mat1, mat2)
    print('Финальная матрица: ',final_matrix, '\n ', 'Ответ:' ,solution)
#     output:
#     Финальная матрица:  [[  1.    0.5   2.    8. ]
#  [  0.    1.  -10.  -28. ]
#  [  0.    0.    1.    3. ]]
#   Ответ: [[1.]
#  [2.]
#  [3.]]

    # Задание 3
    # Введите матрицу при иксах
    A = np.array([[10, 2, 1], [1, 5, 1], [2, 3, 10]], dtype=np.float64)
    # Введите матрицу ответов
    b = np.array([[7], [-8], [6]], dtype=np.float64)
    print('Ответ:',LeadingElement(A, b))
#     output:
#     Ответ: [ 1. -2.  1.]


#   Задание 4
#   Задайте матрицу
    A = np.matrix([[0,-0.1,-0.1],[-0.2,0,-0.1],[-0.2,-0.2,0]]).astype(float)
#   Задайте первое приближение
    X = np.matrix([[1.2],[1.3],[1.4]])

    print('Ответ:',simple_iterations(A,X,0.01))
# output:
# Ответ: [[0.999568 0.99946  0.999316]]


    # Задание 5
    #   Задайте матрицу
    A = np.array([[10, 2, 1], [1, 5, 1], [2, 3, 10]], dtype=np.float64)
    # Введите матрицу ответов
    b = np.array([7, -8, 6], dtype=np.float64)

    print('Ответ:',gauss_seidel(A, b))
# output:
# Ответ: [ 1.00000314 -2.00000037  0.99999948]

if __name__ == '__main__':
    main()