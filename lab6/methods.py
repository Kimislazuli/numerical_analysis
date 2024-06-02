import math
import matplotlib.pyplot as plt

def analytic(x):
    return 4.9 * x ** 2 - 4.9 * x + math.e ** (-x)+ math.e ** x - 2

def f(x, y):
    return y + 11.8+4.9 * x * (1 - x)

def y_prime(mu):
    return mu - 4.9

def runge_kutta_method(h, mu, x_0, y_0, iter):
    xs = [x_0]
    ys = [y_0]
    vs = [y_prime(mu)]
    for i in range(iter):
        K_1 = h * vs[i]
        L_1 = h * f(xs[i], ys[i])
        K_2 = h * (vs[i] + L_1 / 2)
        L_2 = h * f(xs[i] + h / 2, ys[i] + K_1 / 2)
        K_3 = h * (vs[i] + L_2 / 2)
        L_3 = h * f(xs[i] + h / 2, ys[i] + K_2 / 2)
        K_4 = h * (vs[i] + L_3 / 2)
        L_4 = h * f(xs[i] + h, ys[i] + K_3)
        xs.append(xs[i] + h)
        ys.append(ys[i] + 1 / 6 * (K_1 + 2 * (K_2 + K_3) + K_4))
        vs.append(vs[i] + 1 / 6 * (L_1 + 2 * (L_2 + L_3) + L_4))
    return xs, ys, vs

def adams_bashforth(h, mu, x_0, y_0):
    n = int(1 / h)
    xs, ys, vs = runge_kutta_method(h, mu, x_0, y_0, 2)
    for i in range(2, n):
        xs.append(xs[i] + h)
        ys.append(
            ys[i] + h / 12 * (
                    23 * vs[i]
                    - 16 * vs[i - 1]
                    + 5 * vs[i - 2]))
        vs.append(
            vs[i] + h / 12 * (
                    23 * f(xs[i], ys[i])
                    - 16 * f(xs[i - 1], ys[i - 1])
                    + 5 * f(xs[i - 2], ys[i - 2])
            )
        )

    return xs, ys, vs

def shooting(h, x_0, mu):
    new_mu = 0
    while True:
        xs, ys, vs = adams_bashforth(h, mu, x_0, mu)
        t = new_mu
        new_mu = mu - (ys[-1] - math.e - 1 / math.e + 2) / math.e
        mu = t
        if abs(new_mu - mu) < 1e-5:
            break
    xs, ys, vs = adams_bashforth(h, new_mu, x_0, new_mu)
    print(h, new_mu)
    return xs, ys

def tridiagonal(h, x_0):
    xs = [x_0]
    n = int(1 / h)
    lambdas = [1 / (1 + h)]
    mus = [4.9 * h / (1 + h)]
    for i in range(n - 1):
        xs.append(xs[i] + h)
        d_i = h * h * (11.8 + 4.9 * xs[i + 1] * (1 - xs[i + 1]))

        lambdas.append(- 1 / (lambdas[i] - 2 - h ** 2))
        mus.append((d_i - mus[i]) / (lambdas[i] - 2 - h ** 2))
    xs.append(1)
    ys = [0] * (n + 1)
    ys[-1] = math.e + 1 / math.e - 2
    for i in range(n, 0, -1):
        ys[i - 1] = lambdas[i - 1] * ys[i] + mus[i-1]
    return xs, ys


adams_xs_20, adams_ys_20 = shooting(0.05, 0, 0)
adams_xs_10, adams_ys_10 = shooting(0.1, 0, 0)

x = [x / 1000 for x in range(0, 1001)]
y = [analytic(x) for x in x]

tridiagonal(0.1, 0)
tridxs10, tridys10 = tridiagonal(0.1 ,0)
tridxs20, tridys20 = tridiagonal(0.05 ,0)
# plt.plot(tridxs20, tridys20, label="прогонка 0.05")
# plt.plot(adams_xs_20, adams_ys_20, label="стрельба 0.05", linewidth=0.55)
plt.plot(x, y, label="аналитическое", linewidth=0.55, color='magenta', linestyle='dashed')
plt.plot(adams_xs_10, adams_ys_10, label="стрельба 0.1", linewidth=0.7)
plt.plot(tridxs10, tridys10, label="прогонка 0.1")

plt.legend()
plt.show()
