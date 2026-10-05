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


# L_k(x) computed recursively over j
def L(k, j, xn, x):
    if j == len(xn):
        return 1
    if j == k:
        return L(k, j + 1, xn, x)
    return (x - xn[j]) / (xn[k] - xn[j]) * L(k, j + 1, xn, x)

# L_k'(x_k) computed recursively over j
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


# divided difference f[x_0, ..., x_k]
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


# ---------- increasing n ----------
def g25(x):
    return 1 / (1 + 25 * x**2)

def dg25(x):
    return -50 * x / (1 + 25 * x**2)**2


def poly(method, func, dfunc, n, x):
    xn = np.linspace(-1, 1, n + 1)
    yn = func(xn)
    if method == "Lagrange":
        return lagrange(xn, yn, x)
    if method == "Newton":
        return newton(xn, yn, x)
    return hermite(xn, yn, dfunc(xn), x)


def error_table(func, dfunc, ns):
    x = np.linspace(-1, 1, 2000)
    print("n    Lagrange    Newton      Hermite")
    for n in ns:
        e = []
        for m in ["Lagrange", "Newton", "Hermite"]:
            e.append(np.max(np.abs(func(x) - poly(m, func, dfunc, n, x))))
        print(n, "  %.3e  %.3e  %.3e" % tuple(e))
    print()


def plot_n(func, dfunc, method, ns, title):
    x = np.linspace(-1, 1, 1000)
    plt.figure(figsize=(10, 7))
    for i in range(len(ns)):
        n = ns[i]
        xn = np.linspace(-1, 1, n + 1)
        plt.subplot(2, 2, i + 1)
        plt.plot(x, func(x), 'k--', label='original')
        plt.plot(x, poly(method, func, dfunc, n, x), label=method)
        plt.plot(xn, func(xn), 'ro', markersize=3)
        plt.title("n = " + str(n))
        plt.grid(True)
        plt.legend()
    plt.suptitle(title)
    plt.show()


error_table(g, dg, [3, 7, 11, 15, 20])
plot_n(g, dg, "Lagrange", [3, 7, 11, 15], "g(x) = 1/(1 + 2x^2), Lagrange")

error_table(g25, dg25, [5, 10, 15, 20])
plot_n(g25, dg25, "Lagrange", [5, 10, 15, 20], "g(x) = 1/(1 + 25x^2), Lagrange")
plot_n(g25, dg25, "Hermite", [5, 10, 15, 20], "g(x) = 1/(1 + 25x^2), Hermite")
