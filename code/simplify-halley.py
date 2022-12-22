import sympy as sp

fx = sp.symbols('fx')       # f(x)
dx = sp.symbols('dx')       # f'(x)
x = sp.symbols('x')         # x
dy = sp.symbols('dy')       # f'(y)

y = x - fx / dx             # y
ddx = (dy - dx) / (y - x)   # f''(x)
halley = x - (2 * fx * dx) / (2 * dx * dx - fx * ddx)

print(sp.simplify(halley - x))