import math


def f(x):
    return math.sin(2 + x ** 3)


def trapeze(a, b, h):
    m = math.ceil((b - a) / h)
    x = [a + h * i for i in range(0, m + 2)]
    sum_result = f(a) + f(b)
    for_sum = []
    for i in range(1, m):
        for_sum.append(f(x[i]))
        for_sum.append(f(x[i + 1]))
        sum_result += 2 * f(a + h * i)
    return sum_result * h / 2


def three_over_eight(a, b, h):
    m = math.ceil((b - a) / h)
    sum_result = 0
    curr_x = a
    next_x = a + h
    for i in range(m):
        sum_result += f(curr_x) + 3 * (f((2 * next_x + curr_x) / 3) + f((2 * curr_x + next_x) / 3)) + f(next_x)
        curr_x = next_x
        next_x += h

    return sum_result * h / 8


print(trapeze(1, 3, 0.1))
print(trapeze(1, 3, 0.05))
print(trapeze(1, 3, 0.025))
print()
print(three_over_eight(1, 3, 0.1))
print(three_over_eight(1, 3, 0.05))
print(three_over_eight(1, 3, 0.025))
