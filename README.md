
# The NERVU Algebra (N3RVU)

[![License: Attribution Only](https://img.shields.io/badge/License-Attribution--Only-blue.svg)](#attribution-notice)

## Repository Structure

```
N3RVU/
├── .gitignore
├── CITATION.cff                 # Formal citation metadata
├── README.md                    # Project overview (this file)
├── n3rvu.py                     # Core N3RVU algebra implementation
├── N3RVU_Intro.ipynb            # Interactive notebook introduction
├── n3rvu_algebra.md             # Formal algebra description and derivations
├── n3rvu_sympy.py               # Symbolic algebra property checker (SymPy)
├── n3rvu_sympy_out.txt          # Saved console output of property report

```

## What is the NERVU Algebra?

It augments the real numbers with two additional components that track information normally lost under multiplication and division by zero. It forms a three-dimensional commutative, non-associative baric algebra in which the effects of zero are recorded rather than erased. This requires two additional spaces I have called the Vanished space $\mathbb{V}$ and the Undefined space $\mathbb{U}$.

The same process could be applied to the complex numbers by replacing the base type in the `N3RVU` class. I have kept it simple so that it can be easily understoood.

## Why and How Did I Create It?

I was obssessed with the idea of losing information when multiplying by `0` in favor of consistency within the ring axioms. I feel like there is a tension between the $0 \in \mathbb{R}$ and the idealization of `0` as a pure concept beyond a number system like the reals. I had no real purpose other than my obsession that has endured for many years.

I eventually realized that the key step was to formalize the concept of a pure `0`. I attempted to write a paper back in 2024 when I made this realization, which constructed the set of vanished reals and bridged that set to the reals from a pure zero $0_z$ set. But it had several flaws. The most troubling flaw was that I tried to do everything through this new set of vanished reals which proved to ultimately cause problems with the algebraic properties and axioms. 

I realized this year (2025) that I actually needed two extra dimensions to smoothly create a consistent and meaningful algebra. So the set of vanished numbers $\mathbb{V}$ and the new set of undefined numbers were born $\mathbb{U}$. Within the set themselves associativity holds, but it is a misnomer to call them numbers since they do not meet the requirements of a field. The word *number* refers to the fact that they are constructed from the reals as coefficients $\text{Vanished} = \mathbb{R} \mathbb{V}$ and $\text{Undefined} = \mathbb{R} \mathbb{U}$.

I decided not to to fix my paper (I will one day) because the computational path has a construct that represents my idea of a pure zero both philosophically and symbolically, namely the literal `0`. It is guaranteed to never have been the result of any previous arithmetic and quite literally pulls from the idea of my singleton set ${0_z}$ in the 2024 paper.

## ELI5

When we multiply by a `0` it absorbs any other number like this $a * 0 = 0$. Now imagine that instead of absorbing it we store it in a safe by sending it to a new space called $\mathbb{V}$ (for vanished) and it becomes $a * 0 = av$. If you are familiar with imaginary numbers we use a similar type of notation $a * \sqrt{-1} = ai$.

A similar thing happens when we divide by zero (mutliply by $0^{-1}$) like this $a * 1/0 = au$. It goes to a new space called $\mathbb{U}$. It is really just that simple. I am embarassed to say that it took me a long time to realize this and come up with a clean notation. The break through was the concept of a pure zero $0_z$.

We get the values back out of the space with their inverse operation $a * 0 * 0^{-1} = a$. This first multiplication puts $a$ into $au$ and then $au * 0^{-1}$ pulls it back out. This actually creates a push-pull tension between $\mathbb{V}$ and $\mathbb{U}$. The operation  that would put a real number into $\mathbb{V}$ would by definition pull anything out of $\mathbb{U}$ and vice versa. This fact will be important later.

This explanation is actually an oversimplification and special case when we think of $a$ starting with no $v$ or $u$ components, which is a good place to start since we are augmenting the reals.

## Important Properties Quick Reference

We can now create the core properties of the algebra just from the simple statements above. You can go to a [Detailed Algebra Description](n3rvu_algebra.md) for more information. The first element represents the real numbers, then the vanished, and then the undefined.

```bash
Basis: e0 (r), e1 (rv), e2 (ru)

a = (a0,a1,a2)              

a * 0 = (a2, a0+a1, 0) # notice a0 is going into the the V-safe and a2 is coming out of the U-safe

a / 0 = (a1, 0, a0+a2) # notice a0 is going into the the U-safe and a1 is coming out of the V-safe
```

As a result of how we construct the vanished and undefined sets, we have these interesting identities. The above and the follwing identies will result in the algebra being non-associative for the baric triple $(r,rv,ru)$.

### Vanished Set Identity

$av * 0 = av$

### Undefined Set Identity

$au * 1/0 = au$

### Inherited Scalar Identity

An interesting feature of both sets is that they "inherit" the identity of the reals since they are both $\text{Vanished} = \mathbb{R} \mathbb{V}$ and $\text{Undefined} = \mathbb{R} \mathbb{U}$ when you apply a scalar from the reals.

$av * 1 = av$

$au * 1 = au$


## Why is N3RVU Non-Associative?

This is by design, and there is no way around it. It is because once you pack something into $\mathbb{V}$ or $\mathbb{U}$ the order of operations matters. You can't ever get them to "talk again" to anything in the real space. Any operation that would pull it out of $\mathbb{V}$ or $\mathbb{U}$ puts what is currently in the real space to $\mathbb{U}$ or $\mathbb{V}$. This represents the push-pull I mentioned earlier. More precisely the NERVU algebra says to know the final state, you must know the order in which information vanished or became undefined.

Here is an example sequence that illustrates my point. This can be shown compactly with only multiplication but let's begin by first understanding the basic tension where we only apply scalars to the `nervu` defined as `a`.

```bash
#First sequence (a*0+2/0)
a = (1,0,0)
a * 0 = (0,1,0)
a + 2 = (2,1,0)
a / 0 = (1,0,2)

#Second sequence (a/0+2*0)
a = (1,0,0)
a / 0 = (0,0,1)
a + 2 = (2,0,1)
a * 0 = (1,2,0)
```

Now consider a relevant example for associativty `(xy)z = (xy)z`.

```bash
x     = (0,1,0)
y     = (0,1,0)
z     = (0,0,1)

xy    = x * y   = (0,1,0)
yz    = y * z   = (1,0,0)

(xy)z = (x*y)*z = (1,0,0)
x(yz) = x*(y*z) = (0,1,0)
```

To fully understand why this happens you would need to go through my bilinear product formulation from the matrix multiplication table. But I have already provided the intuition.

See [n3rvu_sympy.py](n3rvu_sympy.py) for verification of associativity and other properties.

## Is the NERVU Algebra Useful?

To be honest, I created it because I was obsessed with the idea. I simply didn't like the loss of information that was happening and I wanted to do something about it. After I had fully completed the formalism of the algebra an unexpected possible use emerged. I originally wanted it to be associative, but in hindsight that was impossible due to my construction. There may be an application in neural processing. It might be the most direct mathematical encoding of time-dependent neural processing in a minimal algebra. This arises from the push-pull tension I mentioned earlier, and creates a structural analogy to biphasic synaptic interaction in neural systems.

The `v`-component acts as an inhibitory channel. The `u`-component as an excitatory channel. The `r`-component carries the baseline signal state. The baric weight $w(r,v,u)=r+v+u$ expresses a conserved global activation measure, paralleling excitation–inhibition balance in neural population models. Because multiplication in the algebra is commutative but non-associative, the order of interaction affects outcomes, mirroring sequential neural state updating in biological signaling pathways.

I would not have made this connection without the impressive search capabilities of modern LLMs, which tie semantic statistical patterns together using cosine similarity within the contextualized embedding space. In other words, it identified a closely related topic because the underlying NERVU algebraic structure produced similar patterns in the embedding representation.

## What Does the Name Mean?

`N3RVU` is the class and represents a 3-dimensional (N3) algebra. It consists of the baric triple $(r,v,u)$.

Due to the potential connection to neural processing we can turn the `3` into and `E` yielding `NERVU`. This is supposed to be reminiscent of the word for nerve or neuron. It might also be nerve-wracking to some or nerve you.

## Is Division by Zero in the Reals Undefined?

Yes. The NERVU algebra doesn't change that. It is a result and requirement of the axioms to be consistent for that number system.

## Can Division in NERVU Result in an Undefined Value?

Yes! This may seem strange because it is specifically designed to handle the scalar `0`, but that is not the same as the `nervu` `(0,0,0)`. The important point here is to notice that operating with scalars applied to a `nervu` is different than two `nervu`'s applied to each other through the algebra.

When doing division between two `nervu` we are actually solving for $x$ in $M(b) \cdot x = a$. This has no solution under the following conditions and will throw a `ZeroDivisionError` exception:

$$
\det M(b) = (b_0^2 - b_1 b_2) (b_0 + b_1 + b_2)
$$

$$
\boxed{b_0^2 = b_1 b_2 \quad \text{or} \quad b_0+b_1+b_2 = 0}
$$

The most trivial example of this is when $b = (0,0,0)$. How do you like them poetic appples?

```python
a = N3RVU(1,2,3)
b = N3RVU(0,0,0)
c = N3RVU(1,1,1)
d = N3RVU(2,8,0.5)

# all of these operations will throw ZeroDivisionError
a / b
a / c
a / d
```


## Interactive Notebook (Google Colab)

You can explore the N3RVU algebra live in your browser — no installation required.  
The notebook walks through the same examples as this README, including special zero actions, non-associativity, and interactive widgets.

**Open in Colab:**  
[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/drive/1gLETvw3aaPowiwKslyVYPDvJKEYICaBR?authuser=1#scrollTo=XgkJqhVXRxvF)

Once open, just run the setup cell to load [n3rvu.py](n3rvu.py) (it will auto-download from GitHub if needed) and start experimenting.

## Getting Started

You can try the algebra directly using the demo section of [n3rvu.py](n3rvu.py).
Run it as a standalone script:

```bash
python n3rvu.py
```

Example interactive usage:

```python
from n3rvu import N3RVU, ZERO

a = N3RVU(1, 2, 3)
b = N3RVU(2, 1, 1)

print(a * b)    # Bilinear product
print(a / b)    # Solves for M(b) * x = a, can result in ZeroDivisionError
print(a / 2)    # Scalar division
print(a * 2)    # Scalar multiplication
print(a * 0)    # Multiplication by zero (transforms reals to the V-space and U's to real-space)
print(a / 0)    # Division by zero (transforms reals to the U-space and V's to real-space)
print(a.Z())    # Multiply by pure zero
print(a.U())    # Divide by pure zero
```

You can also explore the full symbolic property verification using:

```bash
python n3rvu_sympy.py
```

The results will print a full algebraic property report like the one stored in [n3rvu_sympy_out.txt](n3rvu_sympy_out.txt).

---

## Attribution Notice

```
This project and the conceptual algebra it presents were developed by
Dale Spencer (Texas Data Safe). You are free to use, modify, reproduce, publish,
or distribute this material for any purpose, including commercial use, provided
that proper attribution is included and this notice is preserved in derivative
works.

To cite this work, please use the citation metadata provided in the
CITATION.cff file included in this repository. Most platforms and tools can
automatically format the citation from that file.

Informal reference:
    Dale Spencer, creator of the NERVU Algebra (N3RVU baric triple formulation).

— Dale Spencer | Texas Data Safe | texasdatasafe@gmail.com
```
