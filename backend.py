from sympy import *

x = Symbol('x', real=True)

def get_coef(expression1, expression2, a, b, sigma=1):
    return integrate(expression1 * expression2 * sigma, (x, a, b)) / integrate(expression2 * expression2 * sigma, (x, a, b))

expressions = [x**i for i in range(5)]

ortho_expressions = [x**0]

for expression in expressions[1:]:
    new_expression = expression
    for orth_exp in ortho_expressions:
        new_expression = new_expression - get_coef(expression, orth_exp, -1, 1) * orth_exp
    ortho_expressions.append(new_expression)

normed_expressions = [ortho_expressions[0]]
for orth_exp in ortho_expressions[1:]:
    val = orth_exp.subs({x: 1})
    normed_expressions.append(orth_exp * (1/val))

for norm_exp in normed_expressions:
    print(factor(norm_exp))
    print(norm_exp.subs({x: 1}))

def get_orthonormal_expressions(n,a,b,sigma):
    expressions = [x**i for i in range(n)]
    ortho_expressions = [x**0]

    for expression in expressions[1:]:
        new_expression = expression
        for orth_exp in ortho_expressions:
            new_expression = new_expression - get_coef(expression, orth_exp, a, b, sigma) * orth_exp
        ortho_expressions.append(new_expression)

    normed_expressions = [ortho_expressions[0]]
    for orth_exp in ortho_expressions[1:]:
        val = orth_exp.subs({x: 1})
        normed_expressions.append(orth_exp * (1/val))
    
    return normed_expressions


