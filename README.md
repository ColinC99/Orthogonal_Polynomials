# Orthogonal Polynomial Visualizer

This is a simple GUI application which allows for visualizing different orthogonal polynomial bases as well as the approximation of any function that can be represented by a sympy expression in the basis of any of these orthogonal bases.
It uses matplotlib to render a 2D graph displaying the selected option. It also includes the ability to customize the x and y limits of the graph. 
The process for approximating functions using orthogonal polynomials is very similar to Fourier series approximation, except using a polynomial basis ($\{1, x, x^2, x^3, \dots\}$) rather than a trigonometric one (\dots, e^{-i2x}, e^{-ix}, 1, e^{ix}, e^{i2x}, \dots).

## Mathematical Background

To understand how this program generates and uses polynomials to approximate functions (much like a Fourier series does with sines and cosines), we need to define three core concepts: the inner product, orthogonality, and the Gram-Schmidt algorithm.

### 1. The Inner Product
When dealing with standard vectors, we use the dot product to determine how much two vectors "overlap." For continuous functions, we use an **inner product**, which is essentially a dot product for functions with an infinite number of points. 

Given two functions $f(x)$ and $g(x)$, defined on an interval $[a,b]$ with a weight function $\sigma(x)$, the inner product is defined as:

$$\langle f,g \rangle = \int_{a}^{b} f(x) g(x) \sigma(x) dx$$

The simplest option for the weight function is $\sigma(x)=1$. These inner products are essentially fully determined by the 3 parameters $a$, $b$, and $\sigma$. 
Keep in mind that the inner product parameters used for the approximation of a function using an orthogonal polynomial basis will always be the same as the inner product parameters used to generate the orthogonal polynomials in the first place.

### 2. Orthogonality
In standard geometry, two vectors are orthogonal (perpendicular) if their dot product is zero. In function spaces, two functions are orthogonal if their inner product is zero (the dot product is also an inner product).

A set of polynomials $P_0(x), P_1(x), P_2(x) \dots$ is considered an orthogonal basis if every pair of distinct polynomials satisfies this condition:

$$\langle P_n,P_m \rangle = 0 \quad \text{for} \quad n \neq m$$

Having an orthogonal basis is important because it allows us to calculate the coefficients of a functions representation in that basis independently for each basis function.

### 3. The Gram-Schmidt Algorithm
We know we want a set of orthogonal polynomials, but how do we build them? We start with the standard, non-orthogonal polynomial basis of simple powers: $\{1, x, x^2, x^3, \dots\}$. 

The Gram-Schmidt process takes this simple set and "straightens them out" one by one by subtracting any overlapping parts (this is the same Gram-Schmidt process used for orthogonalizing a vector basis). 

1. We define the first polynomial as our starting point:

$$P_0(x) = 1$$

2. For the next power ($x$), we subtract its "projection" onto the previous polynomial to make it orthogonal:

$$P_1(x) = x - \frac{\langle x,P_0 \rangle}{\langle P_0,P_0 \rangle} P_0(x)$$

3. We repeat this process for every subsequent power $x^n$, subtracting its projection onto *all* previously calculated polynomials:

$$P_n(x) = x^n - \sum_{k=0}^{n-1} \frac{\langle x^n,P_k \rangle}{\langle P_k,P_k \rangle} P_k(x)$$

Depending on the choice of $a, b, \sigma$ we generate a different set of orthogonal polynomials.
This program uses this exact recursive formula to generate orthogonal polynomial bases using any choice of $a$, $b$, and $\sigma$. 
When selecting a historical orthogonal polynomial basis, keep in mind that this program only displays the un-normalized version of these orthogonal polynomials. 
This is because the historical bases often use their own unique normaliztion method and I wanted to keep things very general. 
Furthermore, the only thing affected by not normalizing the basis is the way the polynomials look when they are themselves looked at. 
Approximation will be exactly identical to the approximation using the normalized basis vectors.

## Function Approximation
Once we have our orthogonal basis, we can approximate any target function $f(x)$ by calculating a coefficient $c_k$ for each polynomial:

$$f(x) \approx \sum_{k=0}^{n} c_k P_k(x) \quad \text{where} \quad c_k = \frac{\langle f,P_k \rangle}{\langle P_k,P_k \rangle}$$

This program calculates each coefficient (although using a numnerical approximation for greater efficiency) and then approximates the function using the above sum to the $n$ selected using the slider in the GUI. 
