import sympy as sp

e = sp.symbols('e')
e1 = e
e2 = e ** 2
e3 = e ** 3
e4 = e ** 4
c2, c3 = sp.symbols('c2 c3')

div = (e1 + c2*e2 + c3*e3 + e4) * (1 - 2*c2*e1 + (4*(c2**2) - 3*c3)*e2 + 8*(c2**3)*e3 + e4)

print(sp.latex(sp.simplify(sp.expand(div))))

