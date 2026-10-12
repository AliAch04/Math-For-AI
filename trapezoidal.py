import math

def f(x) :
  return math.sin(x) 

def trapezoidal(a, b, n) :
  h = (b-a) / n
  S = f(b) - f(a)

  for i in range(1, n-1) :
    x = a + (i*h)
    S = S +(2 * f(x))

  I = (h/2) * S
  return I


print("result is : ", trapezoidal(0, 1, 20000))

