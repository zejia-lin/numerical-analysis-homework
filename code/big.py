# %%
import sympy as sp

# %%
c1, c2, c3, c4, c5, c6 = sp.symbols('c1 c2 c3 c4 c5 c6')
alpha = sp.symbols('alpha')
e = sp.symbols('e')
x = alpha + e
e1 = e
e2 = e ** 2
e3 = e ** 3
e4 = e ** 4
e5 = e ** 5
e6 = e ** 6

# %%
def taylor_1plus_theta(the, order):
    inv = 1 / (1 + the)
    return inv.series(e, 0, order)

def getO(the, order):
    return taylor_1plus_theta(the, order) - taylor_1plus_theta(the, order).removeO()

o3 = getO(e, 3)
o4 = getO(e, 4)
o6 = getO(e, 6)
o7 = getO(e, 7)
o6, o7

# %%
def df(err, expand=False):
    eee = lambda x: x
    if expand:
        eee = sp.expand
    return eee(1 + 2*c2*err + 3*c3*err**2 + 4*c4*err**3 + 5*c5*err**4 + 6*c6*err**5 + o6)


def f(err):
    return c1*err + c2*err**2 + c3*err**3 + c4*err**4 + c5*err**5 + c6*err**6 + o7


def taylor_inv(fenmu, order):
    inv = 1 / fenmu
    return 1 / inv.series(e, 0, order)

# %%
fx = f(e)
fx

# %%
dx = df(e)
dx

# %%
1 / taylor_inv(df(e), 6)

# %%
fx_div_dx = sp.collect(sp.expand(fx / taylor_inv(df(e), 6)), e)
fx_div_dx

# %%
y = x - fx_div_dx
y

# %%
dy = df(y - alpha, True)

# %%
sp.collect(dy, e)

# %%
z = x - 2 * fx / taylor_inv(dx + dy, 6)
sp.collect(sp.expand(z), e)

# %%
fz = f(z - alpha)

# %%
newx = z - 2 * fz / taylor_inv(dx + dy, 6)

# %%
newx


