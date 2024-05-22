import sympy as sp
from copy import deepcopy

# Получаем необходимые данные для текущей итерации методом Ньютона
def get_data_for_newton(J, params, system, start):
    # Создаем словарь приближенных корней вида {x1: valume1, x2: valume2 ...}
    lst_starts = {params[i]: start[i] for i in range(len(params))}

    # Находим значение якбиана в точках x1, x2...
    J_x = []
    for index1 in range(len(J)):
        buff = [] 
        for index2 in range(len(J)):
            val = J[index1][index2].subs(lst_starts).evalf()
            buff.append(round(val, 7))
        J_x.append(buff)

    # Обратная матрица
    J_x = sp.Matrix(J_x)
    J_x_inv = J_x.inv()

    # Находим значение системы функций в точках x1, x2 ....
    F_x = []
    for i in range(len(system)):
        F_x.append(system[i].subs(lst_starts).evalf())

    F_x = sp.Matrix(F_x)

    return J_x_inv, F_x

# Метод Ньютона нахождения корней 
def method_newton(params, system, start):
    eps = 0.000001

    # Строим якобиан от параметров x1, x2
    J = []
    for fun in system:
        buff = []
        for param in params:
            buff.append(sp.diff(fun, param))
        J.append(buff)

    # Получаем необходимые матрицы для расчета
    J_x_inv, F_x = get_data_for_newton(J, params, system, start)

    approximate_roots = []
    
    old_roots = sp.Matrix(deepcopy(start))
    new_roots = old_roots - J_x_inv * F_x

    approximate_roots.append((list(old_roots), None))

    while get_norm_vec(new_roots - old_roots) > eps:

        J_x_inv, F_x = get_data_for_newton(J, params, system, new_roots)

        old_roots = new_roots
        new_roots = old_roots - J_x_inv * F_x

        approximate_roots.append((list(old_roots), 
                                 "{:0<20.15f}".format(get_norm_vec(new_roots - old_roots))))
    

    return new_roots, approximate_roots

# Находим норму вектора
def get_norm_vec(vec):
    max_val = abs(vec[0])
    for val in vec[1:]:
        abs_val = abs(val)
        if abs_val > max_val:
            max_val = abs_val
    return max_val

# Получаем необходимые данные для текущей итерации методом спуска
def get_data_for_speed_down(params, start, diff_list):

    # Создаем словарь приближенных корней вида {x1: valume1, x2: valume2 ...}
    lst_starts = {params[i]: start[i] for i in range(len(params))}
    
    # Вычисление значения функции от точек x1, x2, ....
    val_func = []
    for index in range(len(diff_list)):
        val = diff_list[index].subs(lst_starts).evalf()
        val_func.append(round(val, 7))
    
    val_func = sp.Matrix(val_func)

    return val_func

# Метод спуска
def speed_damn(params, system, start):
    
    eps = 0.000001

    # Строим новую функцию путем возведения в квадрат исходных функций
    # и складываем их
    new_fun = 0
    for fun in system:
        new_fun += fun ** 2
    
    # Находим частные производные по x1, x2....
    diff_list = []
    for par in params:
        diff_list.append(sp.diff(new_fun, par))
    
    # Берем фиксированный шаг 
    alpha = 0.3

    # Будем собирать приближенные корни на каждом шаге
    approximate_roots = []

    # Получаем необходимые данные для итерации 
    val_func = get_data_for_speed_down(params, start, diff_list)

    # Стартовые значения
    old_root = sp.Matrix(start)
    new_root = old_root - alpha * val_func

    # Список приближенных корней в процессе итераций
    approximate_roots.append((list(old_root), None))

    while get_norm_vec(new_root - old_root) > eps:
        
        val_func = get_data_for_speed_down(params, new_root, diff_list)
        old_root = new_root
        new_root = sp.Matrix(old_root) - alpha * val_func

        approximate_roots.append((list(new_root), (get_norm_vec(new_root - old_root))))

    return new_root, approximate_roots

# Функция для проверки найденных корней
def check_roots(funcs_list, lst_roots, params):
    # Связываем параметры x1, x2 со значениями из params
    roots_dict = {params[i]: lst_roots[i] for i in range(len(params))}
    for i, func in enumerate(funcs_list):
        # Если значение = 0, то format не сможет преобразовать такое значение
        # поэтому обрабатываем отдельно
        if func.subs(roots_dict) != 0:
            print(f"Значение в функции {i + 1}:", "{:0<20.15f}".format(func.subs(roots_dict)))
        else:
            print(f"Значение в функции {i + 1}:", func.subs(roots_dict))
    

if __name__ == "__main__":
    # Тестовые примеры от двух систем с двумя переменными
    x1 = sp.Symbol("x1")
    x2 = sp.Symbol("x2")
    x3 = sp.Symbol("x3")

    # Система из примера
    system_1 = [
                sp.sin(x1 + 1.5) - x2 + 2.9, 
                sp.cos(x2 - 2) + x1
               ]
    
    # Система варианта 19
    system_2 = [
                sp.sin(x1 + x2) - x2 - 1.5,
                x1 + sp.cos(x2 - 0.5) - 0.5
               ]
    
    system_3 = [
                0.351*x1 + x2 + 0.572,
                0.867*x1 + x2 + 2.015
               ]
    
    system_4 = [
                0.867*x1 + x2 + 2.015,
                3.315*x1 + x2 + 3.342
               ]
    
    system_5 = [
                0.123*x1 + 0.351*x2 + x3 + 0.572,
                0.752*x1 + 0.867*x2 + x3 + 2.015,
                10.989*x1 + 3.315*x2 + x3 + 3.342
               ]
    
    system_6 = [
                sp.tan(x1 * x2 + 0.3) - x1 ** 2,
                0.9 * x1 ** 2 + 2 * x2 ** 2 - 1,
               ]
    
    #=========================================================
    print("Методо Ньютона")
    print("----------------------------------------------------")
    
    root1, data1 = method_newton([x1, x2], system_1, [5, 10])
    print("Из примера:", list(root1))
    root2, data2 = method_newton([x1, x2], system_2, [5, 10])
    print("Из 19 варианта:", list(root2))
    root3, data3 = method_newton([x1, x2], system_6, [0, 1])
    print("Вариант 11:", list(root3))


    print("----------------------------------------------------\n")
    #=========================================================

    #=========================================================
    print("Методо спуска")
    print("----------------------------------------------------")
    
    root3, data3 = speed_damn([x1, x2], system_1, [5, 10])
    print("Из примера:", list(root3))
    root4, data4 = speed_damn([x1, x2], system_2, [5, 10])
    print("Из 19 варианта:", list(root4))

    print("----------------------------------------------------\n")
    #==========================================================

    #==========================================================
    # Проверка полученных значений
    print("Проверка корней метода Ньютона")
    print("Для системы 1:")
    check_roots(system_1, root1, [x1, x2])
    print("Для системы 2:")
    check_roots(system_2, root2, [x1, x2])
    print("Проверка корней метода спуска")
    print("Для системы 1:")
    check_roots(system_1, root3, [x1, x2])
    print("Для системы 2:")
    check_roots(system_2, root4, [x1, x2])
    #==========================================================


    root, _ = method_newton([x1, x2], [0.351*x1 + x2 + 0.572, 0.867*x1 + x2 + 2.015], [10, 10])
    print(list(root))