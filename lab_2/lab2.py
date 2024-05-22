from copy import deepcopy
from math import log, ceil


def print_matrix(matrix: list[list]) -> print:
    for row in matrix:
        print(*row)
    print()


def add_vector_to_num(vec: list, num: int) -> list:
    new_vec = []
    for elem in vec:
        new_vec.append(elem + num)
    return new_vec


def add_vectors(vec1, vec2: list) -> list:
    res_vec = []
    for index in range(len(vec1)):
        res_vec.append(round(vec1[index] + vec2[index], 4))
    return res_vec


def mul_vector_to_num(vec: list, num: int) -> list:
    res_vec = []
    for index in range(len(vec)):
        if vec[index] != 0:
            res_vec.append(vec[index] * num)
        else:
            res_vec.append(vec[index])
    return res_vec


def div_vector_to_num(vec: list, num: int) -> list:
    res_vec = []
    for index in range(len(vec)):
        res_vec.append(round(vec[index] / num, 4))
    return res_vec


def reverse_matrix(matrix: list[list]) -> list[list]:
    res = []
    len_m = len(matrix)
    for j in range(len(matrix[0])):
        buff = []
        for i in range(len_m):
            buff.append(matrix[i][j])
        res.append(buff)
    return res


def gauss_method(matrix: list[list]) -> float:
    # Приведение к треугольному виду
    for i in range(len(matrix)):
        for j in range(i + 1, len(matrix)):
            mul = mul_vector_to_num(matrix[i], - matrix[j][i] / matrix[i][i])
            matrix[j] = add_vectors(mul, matrix[j])
    roots = []
    # Получение корней системы уранений
    for i in range(len(matrix) - 1, -1, -1):
        count = 0
        k = -1
        for j in range(i + 1, len(matrix)):
            count += matrix[i][j] * roots[k]
            k -= 1
        roots.append(round((matrix[i][-1] - count) / matrix[i][i], 4))
    roots.reverse()
    return roots


def get_matrix_system_reverse(matrix: list[list]) -> list[list]:
    res_matrix = []
    for i in range(len(matrix)):
        matrix_copy = deepcopy(matrix)
        for j in range(len(matrix)):
            if i == j:
                matrix_copy[j][-1] = 1
            else:
                matrix_copy[j][-1] = 0
        res_matrix.append(gauss_method(matrix_copy))
    # Создаем из исходной системы матрицу (те исключаем свободные члены)
    print("Исходная матрица:")
    m = [i[:len(i) - 1] for i in matrix]
    print_matrix(m)
    print("Обратная матрица:")
    print_matrix(reverse_matrix(res_matrix))
    # A * A^-1 = E 
    print("Результат их перемножения:")
    print_matrix(mul_matrix(m, reverse_matrix(res_matrix)))
    return reverse_matrix(res_matrix)


def cubic_norm(matrix: list[list], system=True, vec=False) -> float:
    if vec:
        return max(matrix)
    max_a = 0
    for row in range(len(matrix)):
        if system:
            sum_row = abs(sum(matrix[row][:-1]))
        else:
            sum_row = abs(sum(matrix[row]))
        if sum_row > max_a:
            max_a = sum_row
    if system:
        max_b = 0
        for row in range(len(matrix)):
            if matrix[row][-1] > max_b:
                max_b = matrix[row][-1]
        return max_a, max_b
    return max_a


def abs_and_rel_err(matrix: list[list], reversed_matrix: list[list]) -> tuple[float, float]:
    delta = 0.001
    norm_a, norm_b = cubic_norm(matrix)
    norm_rev_a = cubic_norm(reversed_matrix, system=False)
    abs_err_1 = delta / norm_b
    abs_err_2 = abs_err_1 * norm_rev_a
    rel_err = norm_a * norm_rev_a * abs_err_1
    return (round(abs_err_2, 6), round(rel_err, 6))


def taransform_for_iteration(matrix: list[list]) -> tuple[list, list[list]]:
    vec_c = []
    matrix_B = []
    for i in range(len(matrix)):
        vec = div_vector_to_num(matrix[i], matrix[i][i])
        matrix_B.append(vec[:-1])
        matrix_B[i] = mul_vector_to_num(matrix_B[i], -1)
        matrix_B[i][i] = 0
        vec_c.append(vec[-1])
    if cubic_norm(matrix_B, system=False) < 1:
        return vec_c, matrix_B
    return None


def mul_matrix(matrix1, matrix2: list[list]) -> list[list]:
    res_matrix = []
    for i in range(len(matrix1)):
        row = []
        for k in range(len(matrix2[0])):
            count = 0
            for j in range(len(matrix2)):
                count += matrix1[i][j] * matrix2[j][k]
            row.append(round(count, 4))
        res_matrix.append(row)
    return res_matrix


def add_matrix(matrix1, matrix2: list[list]) -> list[list]:
    res_matrix = []
    for i in range(len(matrix1)):
        row = []
        for j in range(len(matrix1[i])):
            row.append(round(matrix1[i][j] + matrix2[i][j], 4))
        res_matrix.append(row)
    return res_matrix


def iteration_method(matrix: list[list]) -> list[list]:
    eps = 0.01
    val = taransform_for_iteration(matrix)
    if val != None:
        c, b = val
    norm_c = cubic_norm(c, system=False, vec=True)
    norm_b = cubic_norm(b, system=False)
    count_interation = ceil(log(round(eps * (1 - norm_b) / norm_c, 4)) / log(norm_b))
    c = reverse_matrix([c])
    main_c = deepcopy(c)
    arr_steps = [c]
    for _ in range(count_interation):
        m = mul_matrix(b, c)
        buff = add_matrix(m, main_c)
        arr_steps.append(buff)
        c = buff
    return arr_steps[-1]


def transform_for_jkobi(A: list[list]):
    D, L, R, B, b = [], [], [], [], []
    for i in range(len(A)):
        row_D, row_L, row_R = [], [], []
        b.append(A[i][-1])
        for j in range(len(A)):
            if i == j:
                row_D.append(A[i][j])
                row_L.append(0)
                row_R.append(0)
            if j > i:
                row_D.append(0)
                row_L.append(0)
                row_R.append(A[i][j])
            if j < i:
                row_D.append(0)
                row_L.append(A[i][j])
                row_R.append(0)
        D.append(row_D)
        L.append(row_L)
        R.append(row_R)
    rev_D = []
    for i in range(len(D)):
        row_rev_D = []
        for j in range(len(D)):
            if i == j:
                row_rev_D.append(round(1 / D[i][j], 4))
            else:
                row_rev_D.append(0)
        rev_D.append(row_rev_D)
    for row in mul_matrix(rev_D, add_matrix(L, R)):
        B.append(mul_vector_to_num(row, -1))
    b = reverse_matrix([b])
    c = mul_matrix(rev_D, b)
    return D, L, R, B, b, c


def jkobi_method(A: list[list]):
    # Получаем преобразованные матрицы
    D, L, R, B, b, c = transform_for_jkobi(A)
    # Проверим сходмость Якоби
    for i in range(len(A)):
        sum_elem = 0
        for j in range(len(A)):
            if i != j:
                sum_elem += A[i][j]
        if A[i][i] < sum_elem:
            print("Нарушена сходимость")
            break
    # Вычисляем корни
    roots = []
    eps = 0.001
    start = reverse_matrix([[1 for _ in range(len(A))]])
    # a = abs(cubic_norm(start, system=False) - cubic_norm(c, system=False))
    # aa = round(((1 - cubic_norm(B, system=False)) / cubic_norm(B, system=False)) * eps, 4)
    while True:
        copy_start = deepcopy(start)
        start = add_matrix(mul_matrix(B, start), c)
        roots.append(start)
        if copy_start == start:
            break
    return roots[-1]


def check_roots(system, vector):
    for row in range(len(system)):
        row_res = 0
        for index1 in range(len(system)):
            row_res += system[row][index1] * vector[index1][0]

        row_res -= system[row][index1 + 1]
        print(f"x{row + 1}:", "{:0<20.15f}".format(row_res))


if __name__ == "__main__":
    # Тестовые примеры

    # Система уранений из 19 варианта
    matrix_1 = [[3.910, 0.129, 0.283, 0.107, 0.395],
                [0.217, 4.691, 0.279, 0.237, 0.432],
                [0.201, 0.371, 2.987, 0.421, 0.127],
                [0.531, 0.196, 0.236, 5.032, 0.458]]

    # Система уравнений из примера
    matrix_2 = [[5.526, 0.305, 0.887, 0.037, 0.774],
                [0.658, 2.453, 0.678, 0.192, 0.245],
                [0.398, 0.232, 4.957, 0.567, 0.343],
                [0.081, 0.521, 0.192, 4.988, 0.263]]

    # Система уранений из 16 варианта
    matrix_16 = [[2.923, 0.220, 0.159, 0.328, 0.605],
                 [0.363, 4.123, 0.268, 0.327, 0.496],
                 [0.169, 0.271, 3.906, 0.295, 0.590],
                 [0.241, 0.319, 0.257, 3.862, 0.896]]

    # Тестовые матрицы для тестирования операций
    matrix_3 = [[1, 2, 2],
                [3, 1, 2],
                [3, 3, 1]]

    matrix_4 = [[0, ],
                [2, ],
                [0, ]]

    matrix_5 = [[1, 2, 2]]

    matrix_6 = [[0, 2, 0]]

    # Создаем копии систем уравнений для тестирвоания
    copy_matrix_1 = deepcopy(matrix_1)
    copy_matrix_2 = deepcopy(matrix_1)
    copy_matrix_3 = deepcopy(matrix_1)
    copy_matrix_21 = deepcopy(matrix_2)
    copy_matrix_22 = deepcopy(matrix_2)
    copy_matrix_23 = deepcopy(matrix_2)
    copy_matrix_31 = deepcopy(matrix_16)
    copy_matrix_32 = deepcopy(matrix_16)
    copy_matrix_33 = deepcopy(matrix_16)
    copy_matrix_41 = deepcopy(matrix_1)

    # Решение системы линейных уравнений методом Гауса
    # -------------------------------------------------

    print("Метод Гауса")
    print("--------------------------------------------")
    print("19 вариант:", gauss_method(copy_matrix_1))
    print("16 вариант:", gauss_method(copy_matrix_31))
    print("Тестовый вариант", gauss_method(copy_matrix_21))
    print("---------------------------------------------")

    print(reverse_matrix([gauss_method(copy_matrix_31)]))

    # -------------------------------------------------

    # Решение системы линейных уранений методом Итераций
    # -------------------------------------------------

    print("Метод итераций")
    print("--------------------------------------------")
    print("19 вариант:", iteration_method(copy_matrix_2))
    print("16 вариант:", iteration_method(copy_matrix_32))
    print("Тестовый вариант", iteration_method(copy_matrix_22))
    print("---------------------------------------------")

    # -------------------------------------------------

    # Решение системы линейных уранений методом Якоби
    # -------------------------------------------------

    print("Метод Якоби")
    print("--------------------------------------------")
    print("19 вариант:", jkobi_method(copy_matrix_3))
    print("16 вариант:", jkobi_method(copy_matrix_33))
    print("Тестовый вариант", jkobi_method(copy_matrix_23))
    print("---------------------------------------------")

    # -------------------------------------------------

    # Проверка найденных корней для системы уравнений
    # -------------------------------------------------
    print("Подстановка корней методом Гауса")
    check_roots(matrix_1, reverse_matrix([gauss_method(matrix_1)]))
    print("--------------------------------------------")

    print("Подстановка корней методом итераций")
    check_roots(matrix_1, iteration_method(matrix_1))
    print("--------------------------------------------")

    print("Подстановка корней методом Якоби")
    check_roots(matrix_1, jkobi_method(matrix_1))
    print("--------------------------------------------")

    # -------------------------------------------------

    # Высчитывание абсолютной и относительной погрешности
    # -------------------------------------------------

    reversed_system = get_matrix_system_reverse(copy_matrix_41)

    abs_, rel_ = abs_and_rel_err(copy_matrix_41, reversed_system)
    print(f"Абсолютная погрешность: {abs_}")
    print(f"Относительная погрешность: {rel_}")
