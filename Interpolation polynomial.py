import numpy as np
import matplotlib.pyplot as plt


def f(x):
    return np.sin(x)

def df(x):
    return np.cos(x)

def g(x):
    return 1 / (1 + 2 * x**2)

def dg(x):
    return -4 * x / (1 + 2 * x**2)**2


# L_k(x) recursive version
def L(k, j, xn, x):
    if j == len(xn):
        return 1
    if j == k:
        return L(k, j + 1, xn, x)
    return (x - xn[j]) / (xn[k] - xn[j]) * L(k, j + 1, xn, x)

# L_k'(x_k)
def dL(k, j, xn):
    if j == len(xn):
        return 0
    if j == k:
        return dL(k, j + 1, xn)
    return 1 / (xn[k] - xn[j]) + dL(k, j + 1, xn)


def lagrange(xn, yn, x, k=0):
    if k == len(xn):
        return 0
    return yn[k] * L(k, 0, xn, x) + lagrange(xn, yn, x, k + 1)


# divided difference 
def divided_diff(xn, yn):
    if len(xn) == 1:
        return yn[0]
    return (divided_diff(xn[1:], yn[1:]) - divided_diff(xn[:-1], yn[:-1])) / (xn[-1] - xn[0])

def newton(xn, yn, x, k=0, prod=1):
    if k == len(xn):
        return 0
    a = divided_diff(xn[:k + 1], yn[:k + 1])
    return a * prod + newton(xn, yn, x, k + 1, prod * (x - xn[k]))


def hermite(xn, yn, dyn, x, k=0):
    if k == len(xn):
        return 0
    Lk = L(k, 0, xn, x)
    h = (1 - 2 * (x - xn[k]) * dL(k, 0, xn)) * Lk**2
    hbar = (x - xn[k]) * Lk**2
    return yn[k] * h + dyn[k] * hbar + hermite(xn, yn, dyn, x, k + 1)


def run(func, dfunc, a, b, n, title):
    xn = np.linspace(a, b, n)
    yn = func(xn)
    dyn = dfunc(xn)

    x = np.linspace(a, b, 500)
    y = func(x)

    yl = lagrange(xn, yn, x)
    yn_ = newton(xn, yn, x)
    yh = hermite(xn, yn, dyn, x)

    print("Errors for", title)
    print("Lagrange:", np.max(np.abs(y - yl)))
    print("Newton:  ", np.max(np.abs(y - yn_)))
    print("Hermite: ", np.max(np.abs(y - yh)))
    print()

    plt.figure(figsize=(10, 5))
    plt.plot(x, y, 'k--', label='f(x)', linewidth=2)
    plt.plot(x, yl, label='Lagrange')
    plt.plot(x, yn_, ':', label='Newton')
    plt.plot(x, yh, label='Hermite')
    plt.plot(xn, yn, 'ro', label='Nodes')
    plt.title(title)
    plt.xlabel("x")
    plt.ylabel("y")
    plt.legend()
    plt.grid(True)
    plt.show()


run(f, df, -1, 1, 4, "f(x) = sin(x)")
run(g, dg, -1, 1, 4, "g(x) = 1/(1 + 2x^2)")