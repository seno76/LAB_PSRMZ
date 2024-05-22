import sympy as sp
import matplotlib.pyplot as plt
import numpy as np
from copy import deepcopy
from lab4 import method_newton


def check_roots(funnc, arr_params):
    x = sp.Symbol("x")
    roots = []
    for val in arr_params:
        val_func = funnc.subs({x: val})
        roots.append(val_func)
    return roots

def lagrange_interpol_polinom(system):
    # system => [[x1, x2, x3...], [y1, y2, y3...]]

    # Вычисление базисных полиномов
    x = sp.Symbol("x")
    base_polinoms_lst = []
    for i in range(len(system[0])): 
        numerator, denominator = 1, 1
        for j, x_i in enumerate(system[0]):
            if j != i:
                numerator *= (x - x_i)
                denominator *= (system[0][i] - x_i)
        base_polinoms_lst.append(numerator / denominator)
    
    # Вычисление многочлена лагранжа
    L = 0
    for i, val in enumerate(system[1]):
        L += val * base_polinoms_lst[i]

    return L

def finite_differences(system):
    # Вычисление таблицы конечных разностей

    row_y = deepcopy(system[1])
    arr_rows_y = []
    count_variables = len(row_y)
    
    for _ in range(count_variables - 1):
        buff_row = []
        for i in range(len(row_y) - 1):
            buff_row.append(round(row_y[i + 1] - row_y[i], 3))
        arr_rows_y.append(buff_row)
        row_y = buff_row
    
    return arr_rows_y

def divided_differences(system):

    row_x = deepcopy(system[0])
    row_y = deepcopy(system[1])
    arr_rows_x = []
    for j in range(len(system[0]) - 1):
        buff_row = []
        for i in range(len(row_y) - 1):
            buff_row.append(round((row_y[i + 1] - row_y[i]) / (row_x[i + 1 + j] - row_x[i]), 3))
        arr_rows_x.append(buff_row)
        row_y = buff_row
    
    return arr_rows_x

def newtown_polinom(system):
    # Вычисление многочлена ньютона
    x = sp.Symbol("x")

    # Выбираем начальное значение функции => y0
    polinom_newton = system[1][0]
    
    # Получаем таблицу разделенных разностей
    arr_rows_x = divided_differences(system)

    for i in range(len(system[0]) - 1):
        koef = arr_rows_x[i][0]
        j = 0
        while i >= 0:
            koef *= (x - system[0][j])
            i -= 1
            j += 1
        polinom_newton += koef

    # Вычисление значение полинома в точке x1 + x2
    x1 = system[0][1]
    x2 = system[0][2]

    val_func = polinom_newton.subs({x: x1 + x2})

    return polinom_newton, val_func

def build_linear_spline(system):
    # Построим интерполяционные сплайн линейный
    x = sp.Symbol("x")
    
    # Определяем количество параметров
    params = []
    for i in range(len(system[0]) - 1):
        params_for_system = []
        params_for_system.append(sp.Symbol(f"a{i + 1}"))
        params_for_system.append(sp.Symbol(f"b{i + 1}"))
        params.append(params_for_system)

    # Составляем системы 
    systems = []
    j, i = 1, 0
    for _ in range(len(params)):
        buff_sys = []
        j -= 1
        for _ in range(len(params[0])):
            buff_sys.append(system[0][j] * params[i][0] + params[i][1] - system[1][j])
            j += 1
        i += 1
        systems.append(buff_sys)

    # Находим корни систем
    roots = []

    # Берем произвольные стартовые приближения
    intervals = [10 for _ in range(len(params[0]))]

    for i in range(len(systems)):
        root, _ = method_newton(params[i], systems[i], intervals)
        roots.append(list(root))
    
    # Находим итоговую систему линейного сплайна и диапозоны для каждой из строк системы

    linear_spline = []

    for pair_root in roots:
        linear_spline.append(pair_root[0] * x + pair_root[1])

    # Определяем диапозоны
    diaposon = []
    i = 1
    for _ in range(len(system[0]) - 1):
        buff = []
        i -= 1
        for _ in range(2):
            buff.append(system[0][i])
            i += 1
        diaposon.append(buff)

    # Объединяем наши уравнения и диапозоны в кортеж
    # (уравнение, диапозон принимаемых значений x)

    res_system = []
    for i in range(len(diaposon)):
        res_system.append((linear_spline[i], diaposon[i]))

    return res_system

def build_kvadro_spline(system):
    # Построим сплайн квадратичный
    x = sp.Symbol("x")

    # Определяем параметры
    params_kvadro = []

    # Определяем количество систем состоящих из 3 уравнений
    count_system = int(len(system[0]) / 3 + 0.99)

    for i in range(count_system):
        params_for_system = []
        params_for_system.append(sp.Symbol(f"a{i + 1}"))
        params_for_system.append(sp.Symbol(f"b{i + 1}"))
        params_for_system.append(sp.Symbol(f"c{i + 1}"))
        params_kvadro.append(params_for_system)

    # Составляем систему уравнений

    kvadro_system = []
    k = 0
    j = 0
    index_system = 0
    count_itr = len(params_kvadro[0])
    buff_row = None
    for i in range(count_system):
        if buff_row == None:
            syst = []
        else:
            syst = [buff_row]
        for _ in range(count_itr):
            last_row = pow(system[0][k], 2) * params_kvadro[index_system][0] \
                           + system[0][k] * params_kvadro[index_system][1] \
                           + params_kvadro[index_system][2] - system[1][j]
            syst.append(last_row)
            k += 1
            j += 1
        index_system += 1 
        count_itr -= 1
        kvadro_system.append(syst)
        if index_system < len(params_kvadro):
            buff_row = pow(system[0][k - 1], 2) * params_kvadro[index_system][0] \
                           + system[0][k - 1] * params_kvadro[index_system][1] \
                           + params_kvadro[index_system][2] - system[1][j - 1]

    # Вычисление корней систем

    roots = []
    intervals = [10 for i in range(len(params_kvadro[0]))]
    for i in range(len(params_kvadro)):
        root, _ = method_newton(params_kvadro[i], kvadro_system[i], intervals)
        roots.append(list(root))

    # Вычисление итоговой системы уравнений

    res_func = []
    for row in roots:
        buff = 0
        for i in range(len(row)):
            buff += row[i] * x ** (2 - i)
        res_func.append(buff)
    
    # Объединяем с диапозоном значений x

    res = []
    start = 0
    end = len(res_func)
    for i in range(len(res_func)):
        res.append((res_func[i], system[0][start:end + 1]))
        start = end
        end = end + len(res_func)
    
    return res

def draw_func_system(system, name):
    # Определение переменной x как символ из библиотеки SymPy
    x = sp.Symbol('x')

    # Построение графиков
    plt.figure(figsize=(8, 6))

    for equation, x_range in system:
        # Создание массива значений x в указанном диапазоне
        x_values = np.linspace(x_range[0], x_range[-1], 100)
        
        # Преобразование уравнения в функцию
        equation_func = sp.lambdify(x, equation, 'numpy')
        
        # Вычисление значений y для каждого значения x
        y_values = equation_func(x_values)
        
        # Построение графика уравнения
        plt.plot(x_values, y_values)

    # Настройка осей координат
    plt.axhline(0, color='black',linewidth=0.5)
    plt.axvline(0, color='black',linewidth=0.5)

    # Отображение четвертей координатной плоскости
    plt.axhline(0, color='black',linewidth=0.5)
    plt.axvline(0, color='black',linewidth=0.5)
    plt.xlim(-25, 25)
    plt.ylim(-25, 25)

    # Добавление подписей к осям
    plt.xlabel('x')
    plt.ylabel('y')
    # Добавление заголовка
    plt.title(name)

    # Добавление легенды
    plt.legend()

    # Отображение графика
    plt.grid(True)
    plt.show()

def draw_func(func, name):
    
    x = sp.Symbol('x')
    # Создание массива значений x
    x_values = np.linspace(-10, 10, 1000)

    # Вычисление значений функции для каждого значения x
    y_values = [func.subs({x: i})for i in x_values]

    # Построение графика функции
    plt.figure(figsize=(8, 6))
    plt.plot(x_values, y_values, label=f"{name}")

    # Настройка осей координат
    plt.axhline(0, color='black',linewidth=0.5)
    plt.axvline(0, color='black',linewidth=0.5)

    # Отображение четвертей координатной плоскости
    plt.axhline(0, color='black',linewidth=0.5)
    plt.axvline(0, color='black',linewidth=0.5)
    plt.xlim(-25, 25)
    plt.ylim(-25, 25)

    # Добавление подписей к осям
    plt.xlabel('x')
    plt.ylabel('y')

    # Добавление заголовка
    plt.title('График функции')

    # Добавление легенды
    plt.legend()

    # Отображение графика
    plt.grid(True)
    plt.show()

if __name__ == "__main__":
    # Система из варианта 19
    system_1 = [[0.234, 0.649, 1.382, 2.672, 2.849],
                [0.511, 0.982, 2.411, 3.115, 4.184]]
    
    # Система из примера
    system_2 = [[0.351, 0.867, 3.315, 5.013, 6.432],
                [-0.572, -2.015, -3.342, -5.752, -6.911]]
    
    test_system = deepcopy(system_1)
    
    # Нахождение интерполяционного многочлена Лагранжа
    print("Полученный многочлен Лагранжа:")
    print("----------------------------------------------")
    polinom1 = lagrange_interpol_polinom(test_system)
    print(polinom1)
    draw_func(polinom1, "Многочлен Лагранажа")
    print("----------------------------------------------")

    print(check_roots(polinom1, system_1[0]))

    # Построение таблицы конечных разностей
    print("Полученная таблица конечных разностей:")
    print("----------------------------------------------")
    print(finite_differences(test_system))
    print("----------------------------------------------")

    # Построение таблицы разделенных разностей
    print("Полученная таблица разделенных разностей:")
    print("----------------------------------------------")
    print(divided_differences(test_system))
    print("----------------------------------------------")

    # Построение полинома Ньютона
    print("----------------------------------------------")
    polinom2, val = newtown_polinom(test_system)
    print("Полученный многочлен ньютона:")
    print(polinom2)
    draw_func(polinom2, "Многочлен ньютона")
    print("Вычисление значения полинома в точке x1 + x2:")
    print(val)
    print("----------------------------------------------")

    print(check_roots(polinom2, system_1[0]))

    # Построение линейного сплайна
    print("Система линейного сплайна и значения параметров:")
    print("----------------------------------------------")
    system1 = build_linear_spline(test_system)
    print(system1)
    draw_func_system(system1, "Линейный сплайн")
    print("----------------------------------------------")

    for fun in system1:
        print("----------", check_roots(fun[0], system_1[0]))
    
    
    # Построение квадратичного сплайна
    print("Система квадратичного сплайна и значения параметров:")
    print("----------------------------------------------")
    system2 = build_kvadro_spline(test_system)
    print(system2)
    draw_func_system(system2, "Квадратичный сплайн")
    print("----------------------------------------------")

    for fun in system2:
        print(check_roots(fun[0], system_1[0]))