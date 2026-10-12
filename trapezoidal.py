import math


def evaluate_function(x, expr):
  # fonctions mathematiques autorisees
  fct_dict = {
      "sin": math.sin,
      "cos": math.cos,
      "tan": math.tan,
      "exp": math.exp,
      "log": math.log,
      "sqrt": math.sqrt,
      "x": x,
      "pi": math.pi,
      "e": math.e,
  }
  # Evaluation
  return eval(expr, {"__builtins__": {}}, fct_dict)


def trapezoidal(f_expr, a, b, n):
  h = (b - a) / n
  S = evaluate_function(a, f_expr) + evaluate_function(b, f_expr)

  for i in range(1, n):
    x = a + (i * h)
    S = S + (2 * evaluate_function(x, f_expr))

  I = (h / 2) * S
  return I


def main():
  print("--- Numerical Integration: Trapezoidal Rule ---\n")

  print("Suggested functions:")
  print("1. sin(x)")
  print("2. cos(x)")
  print("3. x**2 ")
  print("4. Personalized")

  choix = input("Make your choice (1-4):")

  if choix == "1":
    f_expr = "sin(x)"
  elif choix == "2":
    f_expr = "cos(x)"
  elif choix == "3":
    f_expr = "x**2"
  elif choix == "4":
    f_expr = input(
        "Enter your expression (e.g., x**3 + 2*x or sqrt(x)) : "
    )
  else:
    print("Invalid choice! Defaulting to sin(x).")
    f_expr = "sin(x)"

  try:
    a = float(input("\nEnter the lower bound(a) : "))
    b = float(input("Enter the upper bound (b) : "))
    n = int(input("Enter the number of subdivisions (n) : "))

    resultat = trapezoidal(f_expr, a, b, n)
    print(
        f"\nLe résultat {f_expr} de {a} à {b} est : {resultat}"
    )

  except ValueError:
    print(
        "\nErreur : Veuillez entrer des valeurs numériques valides pour a, b et"
        " n."
    )


if __name__ == "__main__":
  main()
be
import time
t = time.time()
trapezoidal("sin(x)", 0, math.pi, 100000)
print(f"Time for n=100000: {time.time() - t:.3f} s")