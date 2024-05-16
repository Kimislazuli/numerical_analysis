import math
import matplotlib.pyplot as plt


def f(x, y):
    return 3 / (4 * math.sqrt(x))


def analytic(x):
    return 3 * math.sqrt(x) / 2 + 5 / 2


def forward_euler(x_0, y_0, n):
    h = 1 / n
    y_prev = y_0
    x_prev = x_0
    xs = [x_prev]
    ys = [y_prev]
    for i in range(n):
        y_prev -= h * f(x_prev, y_prev)
        ys.append(y_prev)
        x_prev -= h
        x_prev = float('{:.4f}'.format(x_prev))
        xs.append(x_prev)
    xs.reverse()
    ys.reverse()
    return xs, ys


def runge_kutta(x_0, y_0, n):
    h = 1 / n
    y_prev = y_0
    x_prev = x_0
    xs = [x_prev]
    ys = [y_prev]
    for i in range(n):
        if (abs(x_prev - h) == 0.0):
            x_prev += 0.25 * h
        k_1 = h * f(x_prev, y_prev)
        k_2 = h * f(x_prev - h / 2, y_prev - k_1 / 2)
        k_3 = h * f(x_prev - h / 2, y_prev - k_2 / 2)
        k_4 = h * f(x_prev - h, y_prev - k_3)
        y_prev -= 1 / 6 * (k_1 + 2 * k_2 + 2 * k_3 + k_4)
        ys.append(y_prev)
        x_prev -= h
        x_prev = float('{:.4f}'.format(x_prev))
        xs.append(x_prev)
    xs.reverse()
    ys.reverse()
    return xs, ys

def adams_bashforth(x_0, y_0, x_1, y_1, x_2, y_2, n):
    h = 1 / n
    xs = [x_0, x_1, x_2]
    ys = [y_0, y_1, y_2]
    for i in range(n - 2):
        new_y = ys[-1] - h / 12 * (23 * f(xs[-1], ys[-1]) - 16 * f(xs[-2], ys[-2]) + 5 * f(xs[-3], ys[-3]))
        new_x = float('{:.4f}'.format(xs[-1] - h))
        ys.append(new_y)
        xs.append(new_x)
    xs.reverse()
    ys.reverse()
    return xs, ys


x = [x / 1000 for x in range(0, 1001)]
y = [analytic(x) for x in x]
plt.plot(x, y, label="analytic", linewidth=0.5)

x_euler_10, y_euler_10 = forward_euler(1, 4, 10)
# plt.plot(x_euler_10, y_euler_10, label="euler 10", linewidth=0.5)

x_euler_20, y_euler_20 = forward_euler(1, 4, 20)
# plt.plot(x_euler_20, y_euler_20, label="euler 20", linewidth=0.5)

x_euler_30, y_euler_30 = forward_euler(1, 4, 30)
plt.plot(x_euler_30, y_euler_30, label="euler 30", linewidth=0.5)

x_rk_10, y_rk_10 = runge_kutta(1, 4, 10)
# plt.plot(x_rk_10, y_rk_10, label="runge-kutta 10", linewidth=0.5)

x_rk_20, y_rk_20 = runge_kutta(1, 4, 20)
# plt.plot(x_rk_20, y_rk_20, label="runge-kutta 20", linewidth=0.5)

x_rk_30, y_rk_30 = runge_kutta(1, 4, 30)
plt.plot(x_rk_30, y_rk_30, label="runge-kutta 30", linewidth=0.5)

x_ab_10, y_ab_10 = adams_bashforth(1, 4,
                                   x_rk_10[-2], y_rk_10[-2],
                                   x_rk_10[-3], y_rk_10[-3],
                                   10)
# plt.plot(x_ab_10, y_ab_10, label="adams-bashforth 10", linewidth=0.5)

x_ab_20, y_ab_20 = adams_bashforth(1, 4,
                                   x_rk_20[-2], y_rk_20[-2],
                                   x_rk_20[-3], y_rk_20[-3],
                                   20)
# plt.plot(x_ab_20, y_ab_20, label="adams-bashforth 20", linewidth=0.5)

x_ab_30, y_ab_30 = adams_bashforth(1, 4,
                                   x_rk_30[-2], y_rk_30[-2],
                                   x_rk_30[-3], y_rk_30[-3],
                                   30)
plt.plot(x_ab_30, y_ab_30, label="adams-bashforth 30", linewidth=0.5)


plt.legend()
plt.show()
