# %%
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import sympy as sym
import decimal
from decimal import Decimal

plt.rcParams['text.usetex'] = False
plt.rcParams["font.family"] = "Arial Unicode MS"

# %%
mdl = np
roots = [
    1.6319808055660636,
    0.2575302854398608,
    -1.207647827130919,
    1.4044916482153411,
    1.1525907367571583,
    -1,
    1,
    1
]
funcs = [
    lambda x: x**3 + 4*x**2 - 15,
    lambda x: x**2 - mdl.exp(x) - 3*x + 2,
    lambda x: x*mdl.exp(x**2) - mdl.sin(x)**2 + 3*mdl.cos(x) + 5,
    lambda x: mdl.sin(x)**2 - x**2 + 1,
    lambda x: mdl.log(x**2 + 7*x + 14) - x - 2,
    lambda x: (x - 4) * (x + 1)**4 / (mdl.exp(x)),
    lambda x: mdl.exp(x**2 + 11*x - 12) - 1,
    lambda x: mdl.arctan(mdl.exp(x + 3) - 1) * (x - 1)**2
]

# %%
fig, axes = plt.subplots(2, 4, figsize=(10, 5), 
            sharex=True, sharey=True)
sns.set_context('talk')
xmin, xmax = -2, 2
mdl = np
for i in range(8):
    ax = axes.flatten()[i]
    x = np.linspace(xmin, xmax, 1000)
    y = funcs[i](x)
    sns.lineplot(x=x, y=y, ax=ax, linewidth=3, label=rf'$f_{i+1}$')
    ax.legend(loc='upper left')
    ax.hlines(0, xmin, xmax, color='k', linestyles='--')
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(-5, 5)
    # ax.xaxis.set_minor_locator(plt.MultipleLocator(0.5))
    ax.xaxis.set_major_locator(plt.MultipleLocator(1))

plt.subplots_adjust(left=0.05, right=0.99, bottom=0.07, top=0.98)
plt.savefig('../figure/eval-func.pdf')

# %%
mdl = sym
setattr(sym, 'arctan', sym.atan)
func_d1 = []
func_d2 = []
for f in funcs:
    x = sym.symbols('x')
    symf = f(x)
    func_d1.append(sym.diff(symf, x, 1))
    func_d2.append(sym.diff(symf, x, 2))

# for f in func_d1:
#     print(f"lambda x: {f},")
# print()
# for f in func_d2:
#     print(f"lambda x: {f},")

# %%
mdl = np
func_d1 = [
    lambda x: 3*x**2 + 8*x,
    lambda x: 2*x - np.exp(x) - 3,
    lambda x: 2*x**2*np.exp(x**2) + np.exp(x**2) - 2*np.sin(x)*np.cos(x) - 3*np.sin(x),
    lambda x: -2*x + 2*np.sin(x)*np.cos(x),
    lambda x: (2*x + 7)/(x**2 + 7*x + 14) - 1,
    lambda x: -(x - 4)*(x + 1)**4*np.exp(-x) + 4*(x - 4)*(x + 1)**3*np.exp(-x) + (x + 1)**4*np.exp(-x),
    lambda x: (2*x + 11)*np.exp(x**2 + 11*x - 12),
    lambda x: (x - 1)**2*np.exp(x + 3)/((np.exp(x + 3) - 1)**2 + 1) + (2*x - 2)*np.atan(np.exp(x + 3) - 1)
]

func_d2 = [
    lambda x: 2*(3*x + 4),
    lambda x: 2 - np.exp(x),
    lambda x: 4*x**3*np.exp(x**2) + 6*x*np.exp(x**2) + 2*np.sin(x)**2 - 2*np.cos(x)**2 - 3*np.cos(x),
    lambda x: 2*(-np.sin(x)**2 + np.cos(x)**2 - 1),
    lambda x: (-(2*x + 7)**2/(x**2 + 7*x + 14) + 2)/(x**2 + 7*x + 14),
    lambda x: (x + 1)**2*(20*x + (x - 4)*(x + 1)**2 - 8*(x - 4)*(x + 1) - 2*(x + 1)**2 - 40)*np.exp(-x),
    lambda x: ((2*x + 11)**2 + 2)*np.exp(x**2 + 11*x - 12),
    lambda x: (x - 1)**2*(np.exp(x + 3) - 2*(np.exp(x + 3) - 1)*np.exp(2*x + 6)/((np.exp(x + 3) - 1)**2 + 1))/((np.exp(x + 3) - 1)**2 + 1) + 4*(x - 1)*np.exp(x + 3)/((np.exp(x + 3) - 1)**2 + 1) + 2*np.atan(np.exp(x + 3) - 1)
]

# %%
np.seterr(all='raise')

def newton(f, fd, _, x):
    return x - f(x) / fd(x)


def halley(f, fd1, fd2, x):
    return x - (2 * f(x) * fd1(x)) / (2 * fd1(x)**2 - f(x)*fd2(x))


def neta(f, fd, _, x):
    w = x - f(x) / fd(x)
    z = w - f(x) / fd(x) * (f(x) - f(w)/2) / (f(x) - 5*f(w)/2)
    xn = z - f(z) / fd(x) * (f(x) - f(w)) / (f(x) - 3*f(w))
    return xn


def grau(f, fd, _, x):
    y = x - f(x) / fd(x)
    z = y - (y - x) / (2*f(y) - f(x)) * f(y)
    xn = z - (y - x) / (2*f(y) - f(x)) * f(z)
    return xn


def linz(f, fd, _, x):
    y = x - f(x) / fd(x)
    z = x - 2 * f(x) / (fd(x) + fd(y))
    xn = z - f(z) / fd(z)
    return xn


name_map = {
    newton.__name__: '牛顿',
    halley.__name__: 'Halley',
    neta.__name__: 'NM',
    grau.__name__: 'GM',
    linz.__name__: '\\textbf{本文}',
}


def coc(x0, x1, x2, x3):
    x0 = np.float128(x0)
    x1 = np.float128(x1)
    x2 = np.float128(x2)
    x3 = np.float128(x3)
    return np.log((x3 - x2) / (x2 - x1)) / np.log((x2 - x1) / (x1 - x0))


def compute(f, fd1, fd2, x, callback, root, tol, maxiter=100):
    x0, x1, x2, x3 = x, x, x, x
    for i in range(1, maxiter):
        try:
            tmp = x3
            x3 = callback(f, fd1, fd2, x3)
            x0 = x1
            x1 = x2
            x2 = tmp
            if abs(x3 - x2) + abs(f(x3)) < tol:
                return i, (x0, x1, x2, x3)
        except FloatingPointError as e:
            print(x3, x2, e, '\n')
            return -i, (x0, x1, x2, x3)
    return -1, (x0, x1, x2, x3)

# %%
mdl = np
np.seterr('ignore')
init_vals = [
    [1, 2],
    [0, 1],
    [-2, -1],
    [1, 2],
    [1, 2],
    [-1.5, -0.5],
    [0.5, 1.5],
    [0.5, 1.5]
]
# for i in [0]:# range(len(funcs)):
# print('\multirow{10}{*}{$f_' + str(i) + '$}', end='')
for initx_idx in range(2):
    for method in [newton, halley, neta, grau, linz]:
        for i in [3, 7]:
            inix = init_vals[i][initx_idx]
            iters, xs = compute(funcs[i], func_d1[i], func_d2[i], inix, method, roots[i], 1e-12, 1000)
            err = float(Decimal("{:.3g}".format(xs[3] - roots[i])))
            if abs(err) < 1e-16:
                err = '$< 10^{-16}$'
            else:
                err = f'${sym.latex(err)}$'
            if i > 3:
                print('&', end='\t')
            if method == newton:
                print(' ', '\multirow{5}{*}{' + str(inix) + '}', name_map[method.__name__], iters, err, sep='  & \t', end='')
            else:
                print(' ', ' ', name_map[method.__name__], iters, err, sep='  & \t', end='')
        print('\\\\\n')
        if method == linz:
            print('\\cline{2-5}\\cline{7-10}')
print('\\hline')

# %%
mdl = np
np.seterr(all='raise')
init_vals = [
    [1, 2],
    [0, 1],
    [-2, -1],
    [1, 2],
    [1, 2],
    [-1.5, -0.5],
    [0.5, 1.5],
    [0.5, 1.5]
]
for i in [7]:
    for inix in init_vals[i]:
        for method in [newton, halley, neta, grau, linz]:
            # print(inix, method.__name__)
            iters, xs = compute(funcs[i], func_d1[i], func_d2[i], inix, method, roots[i], 1e-12, 1000)
            err = "{:.3g}".format(xs[3] - roots[i])
            print(i, inix, name_map[method.__name__], iters, err, sep='\t')



