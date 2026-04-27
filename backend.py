from sympy import integrate, Symbol, factor, exp, oo, expand, simplify

x = Symbol('x', real=True)

def get_coef(expression1, expression2, a, b, sigma=1):
    coef = integrate(expression1 * expression2 * sigma, (x, a, b)) / integrate(expression2 * expression2 * sigma, (x, a, b))
    return simplify(coef)

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

def get_ortho_expressions(n,a,b,sigma):
    expressions = [x**i for i in range(n)]
    ortho_expressions = [x**0]

    for expression in expressions[1:]:
        new_expression = expression
        for orth_exp in ortho_expressions:
            new_expression = new_expression - get_coef(expression, orth_exp, a, b, sigma) * orth_exp
            
        new_expression = expand(new_expression)
        ortho_expressions.append(new_expression)
    

    return ortho_expressions

def approximate(ortho_base, func, a, b, sigma):
    coefs = [get_coef(func, ob, a, b, sigma) for ob in ortho_base]
    expression = ortho_base[0] * coefs[0]
    for cx in range(1, len(coefs)):
        expression = expression + ortho_base[cx] * coefs[cx]

    return expression

if __name__ == "__main__":
    normed_expressions = get_ortho_expressions(5, 0, oo, exp(-x))
    for ne in normed_expressions:
        print(ne)