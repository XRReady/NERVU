# The "NERVU Algebra"

This algebra augments the real numbers with two additional components that track information normally lost under multiplication and division by zero. It forms a three-dimensional commutative, non-associative baric algebra in which the effects of zero are recorded rather than erased. This requires two additional spaces I have called the Vanished space $\mathbb{V}$ and the Undefined space $\mathbb{U}$.

## Algebra Basis

Basis: e0 (r), e1 (rv), e2 (ru)

```
a = (a0,a1,a2)
b = (b0,b1,b2)
```

r is the real numbers
rv is the vanished space with a real coefficient
ru is the undefined space with a real coefficient

## "Pure" Zero (sting literal) Mapping

a * 0 = (a2, a0+a1, 0)
a / 0 = (a1, 0, a0+a2)

mutliplying by `0` transforms $\mathbb{R} \to \mathbb{V}$ AND $\mathbb{U} \to \mathbb{R}$
mutliplying by $0^{-1}$ (or "divide by 0") transforms $\mathbb{R} \to \mathbb{U}$ AND $\mathbb{V} \to \mathbb{R}$

Note:

when `a0 * 0` happens it leaves behind `0` in the real space just like you would expect for the reals. A division by `0` is a bit trickier. It can leave anything behind but supplants the "undefined" with the transform from $\mathbb{R} \to \mathbb{U}$

```
a' = (a0,0,0)
a' * 0 = (0,a0,0) # notice the reals still have a0 * 0 = 0
```


## RV identity

av * 0 = av

## RU identity

au / 0 = au

## Mutliplication

Multiplication:
```
                         r    rv      ru
        Row/Col ->      e0    e1      e2    
        r   e0       | +e0   +e1     +e2    
        rv  e1       | +e1   +e1     +e0    
        ru  e2       | +e2   +e0     +e2 
```

bilinear product:

```
c0 = a0b0 + a2b1 + b1a2
c1 = a0b1 + a1b0 + a1b1
c2 = a0b2 + a2b0 + a2b2
``` 

This squares like spaces and converts r to the greedy e1 or e2 space.

## Division

Solve the linear system:

$$
M(b) x = a
$$

where the matrix (M(b)) determined by the right-multiplication action is:

$$
M(b)=
\begin{bmatrix}
b_0 & b_2 & b_1\
b_1 & b_0+b_1 & 0\
b_2 & 0 & b_0+b_2
\end{bmatrix}.
$$

Division is possible **if and only if** (M(b)) is invertible. The determinant factorizes as:

$$
\det M(b) = (b_0^2 - b_1 b_2) (b_0 + b_1 + b_2).
$$

So division $a / b$ will fail (no unique solution) exactly when:

$$
\boxed{b_0^2 = b_1 b_2 \quad \text{or} \quad b_0+b_1+b_2 = 0.}
$$


## Bilinear form

I haven't given much thought to the inner product. This is an area for development

### Gram Matrix (inner product)

The algebra is $V = \mathbb{R} e_0 \oplus \mathbb{R} e_1 \oplus \mathbb{R} e_2$

I believe this should produce the following Gram Matrix:

```
G = [[1, 0, 0],
     [0, 1, 0],
     [0, 0, 1]]

```

They should all be orthogonal.

### ChatGPT5 Feedback:

#### Recommendation

* If your primary goal is an **orthonormal** coordinate system and standard geometric intuition, keep (G=I). This aligns with your “three orthogonal directions” intent.
* If you need **Frobenius-style invariance** for proofs or derivations, add the degenerate associative form alongside (G=I) and use it only where that invariance is required.

## Addition behaves normally

```
a + b:

c0 = a0 + b0
c1 = a1 + b1
c2 = a2 + b2
```

## Subtraction behaves normally

```
a - b:

c0 = a0 - b0
c1 = a1 - b1
c2 = a2 - b2
```



## Name and Interpretation

The NERVU Algebra consists of elements called nervu, each represented as a baric triple $(r,v,u)$. It is implemented by the `N3RVU` class.

When I created this algebra I did it purely over an obssession of the information lost when mutliplying by `0`. I realized that if you consider a pure `0` that was never the result of a previous operation, then we could create a meaningful space of numbers that preserve what was multiplied by that pure zero `0`. I further realized two dimensions were needed for both the vanished numbers $\mathbb{V}$ as a result of this operation and the undefined numbers $\mathbb{U}$ as a result of dividing by `0`.

I was obsessed with this idea for many years and took some time to develop a coherent algebra. I was disappointed to realize that it was non-associative because I wanted a "complete number" that didn't sacrifice information for self consistency. In hindsight the non-associativty was obvious. 

What does this have to do with the name of the algebra? Well, I could no longer call it a number. But it is interestingly enough possibly the most direct mathematical encoding of time-dependent neural processing in a minimal algebra. The algebra is N3 (3 dimensions) and encodes RVU (real, vanished, undefined). This is the N3RVU class defining the baric triples and NERVU algebra. The `3` was replaced with an `E` to be reminiscent of a nerve or neuron.

This reflects the algebra’s structural analogy to biphasic synaptic interaction in neural systems. The v-component acts as an inhibitory channel. The u-component as an excitatory channel. The r-component carries the baseline signal state. The baric weight (w(r,v,u)=r+v+u) expresses a conserved global activation measure, paralleling excitation–inhibition balance in neural population models. Because multiplication in the algebra is commutative but non-associative, the order of interaction affects outcomes, mirroring sequential neural state updating in biological signaling pathways. Thus the term Nervu denotes an algebraic element encoding a dynamic, balanced, and directionally responsive signal state.

## NERVU Algebra Property Summary

```bash
=== RVU Algebraic Property Summary ===

Core Structural Properties:
    • Unital (e₀ is identity):        True
    • Associative:                    False
    • Commutative:                    True
    • Flexible:                       True
    • Alternative:                    False
    • Power-associative:              True

Jordan / Lie Structure:
    • Jordan identity holds:          False
    • Jacobi identity (commutator):   True

Baric Structure:
    • Weight map multiplicative:      True

Zero Divisors:
    • Zero divisors exist (sampled):  True


=== Left Multiplication Operators L_{eᵢ} ===

L(ee0 ) =
    [ 1 0 0 ]
    [ 0 1 0 ]
    [ 0 0 1 ]

L(ee1 ) =
    [ 0 0 1 ]
    [ 1 1 0 ]
    [ 0 0 0 ]

L(ee2 ) =
    [ 0 1 0 ]
    [ 0 0 0 ]
    [ 1 0 1 ]


=== Center of the Algebra ===
    Basis: e₀, e₁, e₂   (the center is the whole algebra)


=== Idempotents (solutions to x² = x) ===
    (0, 0, 0)
    (0, 0, 1)
    (0, 1, 0)
    (1/3, 1/3, 1/3)
    (2/3, -1/3, -1/3)
    (1, -1, 0)
    (1, 0, -1)
    (1, 0, 0)


=== Sample Units (non-zero elements with det(Lₓ) ≠ 0) ===
    (-2, -2, -1), (-2, -2, 1), (-2, -2, 2), (-2, -1, -2)
    (-2, -1, -1), (-2, -1, 1), (-2, -1, 2), (-2, 1, -2)
    (-2, 1, -1)


=== Determinant of Left Multiplication ===
    det(Lₓ) = x0**3 + x0**2*x1 + x0**2*x2 - x0*x1*x2 - x1**2*x2 - x1*x2**2
    Factored: (x0**2 - x1*x2)*(x0 + x1 + x2)


=== Jordan Identity Counterexample ===
    x = (1, 1, -1)
    y = (0, 1, 0)
    ( (x²·y)·x ) − ( x²·(y·x) ) = (-2, 2, 0)


=== Literal Zero Operators (Spec Check) ===
    Z(e₁) = (0, 1, 0)
    U(e₂) = (0, 0, 1)
```

## Attribution Notice

```python
"""
`n3rvu_algebra.md` — Attribution Notice

This work and the conceptual algebra presented here were developed by
Dale Spencer over several years of independent study and refinement.

If you use, copy, reference, extend, or fork this material—whether in
academic writing, software, teaching, or online publication—please include
appropriate credit to the author.

If you redistribute modified versions of this code or publish derivative
works (including forks), please retain this attribution notice so that the
origin of the algebraic formulation remains clear.

Informal citation:
    Dale Spencer, creator of the NERVU Algebra (N3RVU baric triple formulation).

Formal citation example:
    Spencer, D. (2025). The NERVU Algebra (N3RVU): A baric triple algebra
    modeling vanished and undefined quantities.

— Dale Spencer | TexasDataSafe | texasdatasafe@gmail.com
"""
```