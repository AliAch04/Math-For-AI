import numpy as np
import matplotlib.pyplot as plt

# -------------------------------------------------------------------
# 1. Define Functions
# -------------------------------------------------------------------
def f(x):
    return np.sin(x)

def df(x): # Derivative needed for Hermite
    return np.cos(x)

def g(x):
    return 1 / (1 + 16 * x**2)

def dg(x):
    return -32 * x / (1 + 16 * x**2)**2


# -------------------------------------------------------------------
# 2. Interpolation Methods
# -------------------------------------------------------------------
def lagrange_interpolation(x_nodes, y_nodes, x_eval):
    """Calculates Lagrange Interpolation Polynomial at x_eval points."""
    n = len(x_nodes)
    p = np.zeros_like(x_eval)
    for i in range(n):
        L_i = np.ones_like(x_eval)
        for j in range(n):
            if i != j:
                L_i *= (x_eval - x_nodes[j]) / (x_nodes[i] - x_nodes[j])
        p += y_nodes[i] * L_i
    return p


def newton_interpolation(x_nodes, y_nodes, x_eval):
    """Calculates Newton Interpolation Polynomial using Divided Differences."""
    n = len(x_nodes)
    coef = np.zeros((n, n))
    coef[:, 0] = y_nodes
    
    for j in range(1, n):
        for i in range(n - j):
            coef[i, j] = (coef[i + 1, j - 1] - coef[i, j - 1]) / (x_nodes[i + j] - x_nodes[i])
            
    # Evaluate polynomial using Horner-like scheme
    p = coef[0, n - 1]
    for k in range(n - 2, -1, -1):
        p = coef[0, k] + (x_eval - x_nodes[k]) * p
    return p


def hermite_interpolation(x_nodes, y_nodes, dy_nodes, x_eval):
    """Calculates Hermite Interpolation Polynomial."""
    n = len(x_nodes)
    p = np.zeros_like(x_eval)
    
    for i in range(n):
        # Calculate L_i(x)
        L_i = np.ones_like(x_eval)
        dL_i_xi = 0.0  # Derivative L_i'(x_i)
        
        for j in range(n):
            if i != j:
                L_i *= (x_eval - x_nodes[j]) / (x_nodes[i] - x_nodes[j])
                dL_i_xi += 1.0 / (x_nodes[i] - x_nodes[j])
                
        # Basis polynomials h_i(x) and h_bar_i(x)
        h_i = (1 - 2 * (x_eval - x_nodes[i]) * dL_i_xi) * (L_i**2)
        h_bar_i = (x_eval - x_nodes[i]) * (L_i**2)
        
        p += y_nodes[i] * h_i + dy_nodes[i] * h_bar_i
        
    return p


# -------------------------------------------------------------------
# 3. Plotting & Error Calculation Function
# -------------------------------------------------------------------
def plot_and_analyze(func, dfunc, a, b, n_points, title):
    x_nodes = np.linspace(a, b, n_points)
    y_nodes = func(x_nodes)
    dy_nodes = dfunc(x_nodes)
    
    x_dense = np.linspace(a, b, 500)
    y_true = func(x_dense)
    
    y_lagrange = lagrange_interpolation(x_nodes, y_nodes, x_dense)
    y_newton = newton_interpolation(x_nodes, y_nodes, x_dense)
    y_hermite = hermite_interpolation(x_nodes, y_nodes, dy_nodes, x_dense)
    
    # Calculate Maximum Errors
    err_lagrange = np.max(np.abs(y_true - y_lagrange))
    err_newton = np.max(np.abs(y_true - y_newton))
    err_hermite = np.max(np.abs(y_true - y_hermite))
    
    print(f"--- Max Errors for {title} (n_points={n_points}) ---")
    print(f"Lagrange Error: {err_lagrange:.6e}")
    print(f"Newton Error:   {err_newton:.6e}")
    print(f"Hermite Error:  {err_hermite:.6e}\n")
    
    # Plotting
    plt.figure(figsize=(10, 5))
    plt.plot(x_dense, y_true, 'k--', label='Original f(x)', linewidth=2)
    plt.plot(x_dense, y_lagrange, label='Lagrange', alpha=0.7)
    plt.plot(x_dense, y_newton, ':', label='Newton', alpha=0.7)
    plt.plot(x_dense, y_hermite, label='Hermite', alpha=0.7)
    plt.plot(x_nodes, y_nodes, 'ro', label='Nodes')
    
    plt.title(f"{title} (n={n_points-1}, nodes={n_points})")
    plt.xlabel("x")
    plt.ylabel("y")
    plt.legend()
    plt.grid(True)
    plt.show()

# Run example for n = 3 (4 points) on [-1, 1]
plot_and_analyze(g, dg, -1, 1, 4, "g(x) = 1/(1+16x^2)")