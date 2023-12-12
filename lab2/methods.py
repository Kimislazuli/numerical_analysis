import math

epsilon = 0.5 * 10 ** (-5)


def f(x):
    return math.cos(x) - 2 * math.e ** (-x)


def derivative(x):
    return 2 * math.e ** (-x) - math.sin(x)


def dichotomy(a, b):
    n = 0
    x_n = 0
    while (b - a) > epsilon:
        x_n = (b + a) / 2
        if f(x_n) * f(a) < 0:
            b = x_n
        else:
            a = x_n
        n += 1
    return x_n, n


def newton(x_0, m, M):
    n = 0
    x_n = x_0
    x_n_prev = 0
    while M / (2 * m) * (x_n_prev - x_n) ** 2 > epsilon:
        x_n_prev = x_n
        x_n = x_n - f(x_n) / derivative(x_n)
        n += 1
    return x_n, n


def modified_newton(x_0, m, M):
    n = 0
    x_n = x_0
    x_n_prev = 0
    while (1 - m / M) * abs(x_n_prev - x_n) > epsilon:
        x_n_prev = x_n
        x_n = x_n - f(x_n) / derivative(x_0)
        n += 1
    return x_n, n


def secant(x_0, x_1, m):
    n = 1
    x_n = x_1
    while abs(f(x_n)) / m > epsilon:
        x_n = x_n - f(x_n) * (x_n - x_0) / (f(x_n) - f(x_0))
        n += 1
    return x_n, n


def moving_secant(x_0, x_1, m):
    n = 1
    x_prev = x_0
    x_n = x_1
    while abs(f(x_n)) / m > epsilon:
        x_new = x_n - f(x_n) * (x_n - x_prev) / (f(x_n) - f(x_prev))
        x_prev, x_n = x_n, x_new
        n += 1
    return x_n, n


def phi(x):
    return -math.acos(2 * math.e ** (-x)) + 2 * math.pi


def relax_phi(x):
    return x - 0.98494 * f(x)


def simple_iteration(x_0, q):
    n = 1
    x_n = x_1 = phi(x_0)
    while abs(x_1 - x_0) * q ** n / (1 - q) > epsilon:
        x_n = phi(x_n)
        n += 1
    return x_n, n


def relaxation(x_0, q):
    n = 1
    x_n = x_1 = phi(x_0)
    while abs(x_1 - x_0) * q ** n / (1 - q) > epsilon:
        x_n = relax_phi(x_n)
        n += 1
    return x_n, n


print(dichotomy(3 * math.pi / 2, 4.8))
print(newton(3 * math.pi / 2, 1.01262, 0.10396))
print(modified_newton(3 * math.pi / 2, 1.01262, 1.01797))
print(secant(3 * math.pi / 2, 4.8, 1.01262))
print(moving_secant(3 * math.pi / 2, 4.8, 1.01262))
print(simple_iteration(3 * math.pi / 2, 0.01796))
print(relaxation(3 * math.pi / 2, 0.00264))
