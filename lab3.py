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
def main():
    # Введите матрицу при иксах
    mat1 = np.matrix([[5, 0, 1], [2, 6, -2], [-3, 2, 10]]).astype(float)

    # Введите матрицу ответов
    mat2 = np.matrix([[11], [8], [6]]).astype(float)

    print(gauss(mat1,mat2))
#     output:
#           [[2.]
#            [1.]
#            [1.]]
if __name__ == '__main__':
    main()