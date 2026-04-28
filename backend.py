from sympy import integrate, Symbol, factor, exp, oo, expand, simplify, collect

x = Symbol('x', real=True)

# def get_orthonormal_expressions(n,a,b,sigma):
#     expressions = [x**i for i in range(n + 1)] 
#     ortho_expressions = [x**0]

#     for expression in expressions[1:]:
#         new_expression = expression
#         for orth_exp in ortho_expressions:
#             new_expression = new_expression - get_coef(expression, orth_exp, a, b, sigma) * orth_exp
#         ortho_expressions.append(new_expression)

#     normed_expressions = [ortho_expressions[0]]
#     for orth_exp in ortho_expressions[1:]:
#         val = orth_exp.subs({x: 1})
#         normed_expressions.append(factor(orth_exp * (1/val)))
    
#     return normed_expressions

def get_ortho_expressions(n,a,b,sigma, start=None):
    expressions = [x**i for i in range(n)]
    if start == None:
        ortho_expressions = [x**0]
    else:
        ortho_expressions = start

    basis_norms = [simplify(integrate(express * express * sigma, (x, a, b))) for express in ortho_expressions]
    index = len(ortho_expressions)
    for expression in expressions[len(ortho_expressions):]:
        print(expression)
        new_expression = expression
        shortcut_possible = a == (-b)
        for ox, orth_exp in enumerate(ortho_expressions):
            # Shortcut if interval is symmetric over 0 and function is odd:
            # This works because in these cases, integrate(expression * orth_exp * sigma, (x, a, b)) = 0
            if shortcut_possible and (index + ox) % 2 == 1:
                continue
            numer = integrate(expression * orth_exp * sigma, (x, a, b)) 
            new_expression = new_expression - ((numer / basis_norms[ox]) * orth_exp)
            print(new_expression)
            
        new_expression = collect(expand(new_expression), x)
        ortho_expressions.append(new_expression)
        basis_norms.append(simplify(integrate(new_expression * new_expression * sigma, (x, a, b))))
        index += 1

    return ortho_expressions


if __name__ == "__main__":
    # normed_expressions = get_ortho_expressions(5, 0, oo, exp(-x))
    # for ne in normed_expressions:
    #     print(ne)
    import pickle

    with open("orthogonal_polys_cache.pkl", "rb") as f:
        cache = pickle.load(f)
    
    for key in cache:
        print(f"{key}:")
        for val in cache[key]:
            print(val)
        