import numpy as np
import sympy
from scipy import integrate
import matplotlib.pyplot as plt

def one(x):
    return 1

def inner_product(f,g, sigma, a=-1, b=1):
    func = lambda x: f(x) * g(x) * sigma(x)

    return integrate.quad(func, a, b)[0]

def get_ortho_base(n, a, b, sigma=one):
    bases = [one]
    for i in range(1, n):
        phi_i = lambda x: x**i
        sub = lambda x: sum([((inner_product(phi_i, f, sigma, a, b) / inner_product(f, f, sigma, a, b)) * f(x)) for f in bases])
        bases.append(lambda x: phi_i(x) - sub(x))
    
    return bases

class PolyBase():
    def __init__(self, n):
        self.bases = [one] + [lambda x: x**i for i in range(1,n)]

    def orthogonalize(self, a, b, sigma):
        new_bases = []
        for ix, base in enumerate(self.bases):
            if ix == 0:
                new_bases.append(base)
            else:
                coefs = [(inner_product(base, f, sigma, a, b) / inner_product(f, f, sigma, a, b)) for f in new_bases]
                new_func = base
                for cx, cn in enumerate(coefs):
                    new_func = lambda x: new_func(x) - (cn * new_bases[cx](x))
                new_bases.append(new_func)
        self.ortho_bases = new_bases


def plot_functions_quick(func_list, a=-1, b=1, num_points=200):
    """
    Plots a list of functions over the interval [a, b].
    """
    # Generate the x values
    x_vals = np.linspace(a, b, num_points)
    
    plt.figure(figsize=(8, 6))
    
    for i, func in enumerate(func_list):
        # We use a list comprehension here instead of just func(x_vals).
        # Why? If one of your polynomials is just a constant (like P_0(x) = 1), 
        # func(x_vals) might return a single integer '1' instead of an array of 1s,
        # which will make matplotlib crash. This safely forces it into a list!
        y_vals = [func(x) for x in x_vals]
        
        plt.plot(x_vals, y_vals, label=f"Polynomial {i}")
        
    # Formatting to make it easy to read
    plt.axhline(0, color='black', linewidth=1) # Emphasize the x-axis
    plt.title(f"Polynomial Sanity Check on [{a}, {b}]")
    plt.xlabel("x")
    plt.ylabel("y")
    plt.grid(True, linestyle=':', alpha=0.7)
    
    # Place legend outside if you have a lot of functions, or just inside
    plt.legend(loc='best')
    
    # Render the plot
    plt.show()

poly_base = PolyBase(5)
poly_base.orthogonalize(-1,1,one)

plot_functions_quick(poly_base.bases)
plot_functions_quick(poly_base.ortho_bases)