import sympy as sp
   
def method_iteration(func, x, range_roots):
    diff_1 = sp.diff(func, x)
    lst = {diff_1.subs({x: root}).evalf():root for root in range_roots}
    Q = abs(list(lst.keys())[0])
    flag = False
    for i in list(lst.keys())[1:]:
        if abs(i) > Q:
            Q = abs(i)
        if i < 0:
            flag = True 
    if flag: k = -int(Q / 2 + 0.99)
    else: k = int(Q / 2 + 0.99)
    phi_x = x - func / k
    diff_phi_x = sp.diff(phi_x, x)
    dict_roots = {abs(diff_phi_x.subs({x: root}).evalf()):root for root in range_roots}
    q = max(dict_roots.keys())
    eps = round((q / (1 - q) * 0.0000001).evalf(), 4)
    roots = []
    end = dict_roots[q]
    start = (phi_x.subs({x: end})).evalf()
    i = 0
    while abs(end - start) > eps:
        start, end = (phi_x.subs({x: start})).evalf(), start
        i += 1
    roots.append(round(start, 6))
    return roots

def method_newton(func, x, range_roots):
    eps = 0.000001
    diff_1 = sp.diff(func, x)
    roots = []
    for border in range_roots:
        start = border - 1
        end = border + 1
        i = 0
        while abs(end - start) > eps:
            start, end = start - (func.subs({x: start}) / diff_1.subs({x: start})).evalf(), start
            i += 1
        roots.append(round(start, 6))
    return roots

def check_roots(x, func, roots_lst):
    dict_roots = {}
    for root in roots_lst:
        dict_roots[root] = "{:0<20.15f}".format(func.subs({x: root}))
    return dict_roots

def method_lotsman(func, x, range_roots):
    eps = 0.000001
    diff_1 = sp.diff(func, x)
    diff_2 = sp.diff(diff_1, x)
    roots = []
    for border in range_roots:
        start = border - 1
        end = border + 1
        i = 0
        while abs(end - start) > eps:
            up = 2 * func.subs({x: start}) * diff_1.subs({x: start})
            down = 2 * diff_1.subs({x: start})**2 - func.subs({x: start})*diff_2.subs({x: start})
            start, end = (start - up / down).evalf(), start
            i += 1
        roots.append(round(start, 6))
    return roots

if __name__ == "__main__":

    # Создание объектов функций

    x = sp.Symbol("x")
    func = 2 - 0.5*x**2 - 0.5*x**-1*sp.sin(x) - x
    func_2 = 2*x**2 - x**3 - sp.exp(x)
    range_roots = [-100, 100]


    # В алгоритмах Ньютона и Лоцмана мы задаем некоторую 
    # область поиска корней в нашем примере от -100 до 100, 
    # если корень только один то он и будет выдан,
    # а если есть два корня > 0 и < 0 то они будут найдены.

    print()
    print("Нахождение корней методом Ньютона")
    print("----------------------------------------")

    roots = method_newton(func_2, x, range_roots)
    print("Вариант из примера:", check_roots(x, func_2, roots))
    roots = method_newton(func, x, range_roots)
    print("Вариант 19:", check_roots(x, func, roots))

    print("----------------------------------------\n")

    print("Нахождение корней методом Лоцмана")
    print("----------------------------------------")


    
    roots = method_lotsman(func_2, x, range_roots)
    print("Вариант из примера:", check_roots(x, func_2, roots))
    roots = method_lotsman(func, x, range_roots)
    print("Вариант 19:", check_roots(x, func, roots))

    print("----------------------------------------\n")

    # В алгоритме итераций необходимо первоначально знать 
    # область в которой необходимо осуществлять поиск
    # Поэтому чтобы получить несколько корней 
    # надо знать диапозон значений. Их мы получаем
    # путем граффического анализа функции.
    
    print("Нахождение методом простых итераций")
    print("----------------------------------------")

    range_roots = [-1, 2]
    roots = method_iteration(func_2, x, range_roots)
    print("Вариант из примера:", check_roots(x, func_2, roots))


    range_roots = [-4, -3]
    roots = method_iteration(func, x, range_roots)
    print("Вариант 19:", check_roots(x, func, roots))


    range_roots = [1, 2]
    roots = method_iteration(func, x, range_roots)
    print("Вариант 19:", check_roots(x, func, roots))

    print("----------------------------------------")