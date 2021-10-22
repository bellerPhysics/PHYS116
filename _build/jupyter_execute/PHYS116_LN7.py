#!/usr/bin/env python
# coding: utf-8

# # LN7: Fourier transforms
# 
# ## From Fourier series to transforms
Suppose we have a function that is *not* periodic, and the domain of interest 
is all real numbers $ (-\infty,\infty) $
Recall that on any *finite* interval of length $L$, we could represent any function by its Fourier series, and the series converges to the original function everywhere except the endpoints of the interval. 

So by taking $L \rightarrow \infty$, we can accurately represent our arbitrary function at all finite $x$.

Let's use the complex Fourier series:
# \begin{align*}
#     f(x) &= \sum_{n=-\infty}^\infty c_n e^{2n\pi i x /L} \\
#     c_n &= \frac{1}{L} \int_{-L/2}^{L/2} du f(u) e^{-2n\pi i u/L} \\
#     \Rightarrow f(x) &= \frac{1}{L} \sum_{n=-\infty}^\infty \int_{-L/2}^{L/2}  du f(u) e^{2n\pi i (x-u)/L}
# \end{align*}

# Let $\Delta q = 2\pi /L$. Then

# $$   
# f(x) = \frac{\Delta q}{2\pi} \sum_{n=-\infty}^\infty 
#         \int_{-L/2}^{L/2} du f(u) e^{n \Delta q i(x-u)}         
# $$

# Now, instead of summing over all integers $n$, we are going to divide the *real* number line into intervals of length $\Delta q$. The $n$th interval is centered at
# $ q_n \equiv n \Delta q$:

# $$   
# f(x) = \frac{1}{2\pi}  \sum_{n=-\infty}^\infty \Delta q 
#         \int_{-L/2}^{L/2} du f(u) e^{q_n i (x-u)}
# $$

# As $L\rightarrow \infty$, $\Delta q = 2\pi /L \rightarrow 0$ and so the sum 
# \begin{equation*} \sum_{n=-\infty}^\infty \Delta q (\dots\text{function of $q_n$}\dots) \end{equation*}
# becomes the integral

# $$
# \int_{-\infty}^\infty dq_n (\dots \text{function of $q_n$}\dots ) 
# $$

# We have no more need for $n$, so let's use $q$ in place of $q_n$. We'll also replace $\pm L/2$ with $\pm \infty$.

# \begin{align*}
#     f(x) &= \frac{1}{2\pi} \int_{-\infty}^\infty dq \int_{-\infty}^\infty du f(u) e^{iq(x-u)} \\
#     &= \frac{1}{2\pi} \int_{-\infty}^\infty dq e^{iqx} \int_{-\infty}^\infty du f(u) e^{-iqu}
# \end{align*}

# Define
# 
# $$ 
# \tilde f(q) \equiv \frac{1}{2\pi} \int_{-\infty}^\infty du f(u) e^{-iqu}. 
# $$ 
# 
# (We could use $x$ instead of $u$ here, since the variable being integrated over can have whatever name we like as long as it doesn't already appear in the equation.)
# 
# Then,
# 
# $$
# f(x) = \int_{-\infty}^\infty dq \tilde f(q) e^{iqx} 
# $$

# We have found

# 
# ```{important}
# 
# \begin{align*}
#         f(x) &= \int_{-\infty}^\infty dq \tilde f(q) e^{iqx} \\
#         \tilde f(q) &\equiv \frac{1}{2\pi} \int_{-\infty}^\infty dx f(x) e^{-iqx} 
#     \end{align*}
# 
# For any function $f$, $\tilde f$ is its ***Fourier transform***.
# 
# ```
# 

# ```{note}
# $f$ and $\tilde f$ are functions of different variables, $x$ and $q$ respectively. Since $iqx$ must be dimensionless, the units of $q$ must be the reciprocal of the units of $x$. For example,
# 
# * For a function of time $f(t)$, the Fourier transform is a function of frequency $\tilde f(\omega)$, where $\omega$ has units of (e.g.) $\mathrm{s}^{-1}=\mathrm{Hz}$.
# 
# * For a function of position $f(x)$, the Fourier transform is a function of ***wavenumber*** $q$, which has units of $1/\mathrm{length}$.
# ```

# ```{note}
#         
# "Fourier transform" can refer to the function $\tilde f(q)$ or to the act of taking $\frac{1}{2\pi}\int_{-\infty}^\infty dx f(x) e^{-iqx}$. Thus we might also "***take the inverse Fourier transform***" of any function:
# 
# $$ 
# \int_{-\infty}^\infty dq \tilde g(q) e^{-iqx}. 
# $$    
# ```

# ### Example 1

# 
# *Suppose our step function example from the Fourier Series section were not periodic, but instead*
# 
# $$ 
# f(x) = \begin{cases} 1 & -2 \leq x < 2 \\
#                         0 & \text{otherwise} 
#           \end{cases}
# $$
# 
# *Find its Fourier transform.*
# 
# 
# 
# 
# \begin{align*}
#     \tilde f(q) &= \frac{1}{2\pi} \int_{-\infty}^\infty dx f(x) e^{-iqx} 
#         = \frac{1}{2\pi} \int_{-2}^{2} dx (1) e^{-iqx} \\
#         &= \frac{1}{2\pi} 
#             \left[ \frac{e^{-iqx}}{-iq}\right]^2_{-2}
#         = -\frac{1}{q \pi}\frac{e^{-2qi}-e^{2qi}}{2i} = \frac{\sin(2q)}{q\pi}
# \end{align*}
# 
# 
# \begin{align*}
#     f(x) &= \int_{-\infty}^\infty dq \tilde f(q) e^{iqx}  \\
#     &= \int_{-\infty}^\infty dq \frac{\sin(2q)}{\pi q} (\cos (qx) + i \sin(qx))
# \end{align*}
# 
# 
# Now note that 
# 
# * $\sin(2q)$ is an odd function of $q$.
# 
# * $(1/(\pi q)) = (1/\pi) q^{-1}$ is also an odd function of $q$.
# 
# * the interval of integration is symmetric about $0$. 
# 
# Therefore, $\sin(2q)/(\pi q)$ is an even function, and its product with an odd function (such as $\sin(qx))$ will integrate to $0$, leaving:
# 
# $$
# f(x) = \int_{-\infty}^\infty dq \frac{\sin(2q)}{\pi q} \cos(qx) 
# $$
# 
# 

# ### Example 2

# \begin{align*}
# f(x) &= \cos(3x) \\
# \tilde f(q) &= \frac{1}{2\pi}\int_{-\infty}^\infty dx \cos(3x) e^{-iqx} \\
# &= \frac{1}{2\pi} \int_{-\infty}^\infty dx \cos(3x) (\cos(qx) - i \sin(qx)) 
# \end{align*}
# 
# The $i \sin(qx)$ term is odd and so contributes $0$ to the integral (since $\cos(3x)$ is even). 

# That leaves    
# 
# \begin{align*}
# \tilde f(q)    &= \frac{1}{2\pi} \int_{-\infty}^\infty dx \frac{e^{3ix}+e^{-3ix}}{2} 
#         \frac{e^{iqx} + e^{-iqx}}{2} \\
#         &= \frac{1}{8\pi} \int_{-\infty}^\infty dx
#             \left[ e^{(3+q)ix}+ e^{-(3+q)ix} + e^{(3-q)ix} + e^{-(3-q)ix} \right]
# \end{align*}

# 
# To evaluate this integral, we will need to introduce the ***Dirac delta function***, defined by
# 
# \begin{align*} 
# \delta(u)& \equiv \begin{cases} \infty & u=0 \\ 0 & u\neq 0 \end{cases} \\
# \int_a^b du \delta(u)& = \begin{cases} 1 & a < 0 < b \\ 0 & \text{otherwise} \end{cases} 
# \end{align*}
# 
# 
# It can be shown that a $\delta$ function in an integral "picks out" the integrand's value at the point where $\delta$'s argument is zero:
# 
# $$
# \int_a^c du \delta(u-b) g(u) = \begin{cases} g(b) & a < b < c \\ 0 & \text{otherwise} \end{cases}
# $$
# 
# Also, $\delta(u)$ is an even function of $u$. 

# Now, it's easy to show that the inverse Fourier transform of $\delta(q)$ is simply 1:
# 
# $$
# \int_{-\infty}^\infty dq \delta(q) e^{-iqx} = e^{-0ix} = 1.
# $$
# 
# $\delta(q)$ is therefore an expression of the Fourier transform of the function 1:
# 
# $$
#     \frac{1}{2\pi} \int_{-\infty}^\infty dq e^{iqx} = \delta(q) \\
#     \Rightarrow \int_{-\infty}^\infty dq e^{iqx} = 2\pi \delta(q).
# $$
# 

# Returning to our example, an integral like 
# 
# $$
#  \int_{-\infty}^\infty dx e^{(3+q)ix}
# $$
# 
# is equal to $2\pi$ times the Dirac delta function evaluated at $3+q$ instead of $q$. 
# 
# $$
# \tilde f(q) = \frac{2\pi}{8\pi} \left[ \delta(3+q) +\delta(-3-q) +\delta(3-q) +\delta(-3+q) \right] 
# $$
# 
# Using the evenness of $\delta(u)$:
# 
# $$
# \tilde f(q) = \frac{1}{2} \left[ \delta(q+3) + \delta(q-3) \right] 
# $$
# 

# We can recover the original function as
# \begin{align*}
#     f(x) &= \int_{-\infty}^\infty dq \tilde f(q) e^{-iqx} 
#     = \frac{1}{2} \left[ \int_{-\infty}^\infty dq \delta(q+3) e^{-iqx} + \int_{-\infty}^\infty dq \delta(q-3)e^{-iqx}\right] \\
#     &= \frac{1}{2} \left[ e^{-3ix} + e^{3ix} \right] = \cos(3x) \quad \checkmark
# \end{align*}    
# 

# What does this example show us? The Fourier transform of a simple harmonic oscillation, with periodicity $L=2\pi/3$, is an infinite spike at $q=3=2\pi/L$, and zero everywhere else. 
# 
# What about for more general $f(x)$? Note that the Fourier transform is *linear*:

# 
# \begin{align*}
#     f(x) &= g(x) + c h(x) \\
#     \tilde f(q) &= \frac{1}{2\pi} 
#         \int_{-\infty}^\infty dq f(x) e^{iqx} 
#     = \frac{1}{2\pi} \int_{-\infty}^\infty dq g(x) e^{iqx} 
#         + c \int_{-\infty}^\infty dq h(x) e^{iqx}   \\
#     &= \tilde g(q) + c \tilde h(q) \quad \checkmark
# \end{align*}

# For example, if $f(x) = \cos(3x) + 2 \cos(4 x)$ then
# 
# $$ 
# \tilde f(q) = \frac{1}{2} \left[ \delta(q+3) + \delta(q-3) + 2 \delta(q+4) + 2\delta(q-4) \right]
# $$

# So, we can think of any function $\tilde f(x)$ as a linear combination of SHO functions $\cos(qx)$ and $\sin(qx)$—or equivalently, $e^{iqx}$ and $e^{-iqx}$. 
# 
# Now $q$ takes on the value of all real numbers, not just integers, so the linear combination *sum* becomes a linear combination *integral*: the inverse Fourier transform!
# 
# Meanwhile, the Fourier transform itself just takes the inner product of $f(x)$ with basis element $e^{iqx}$, that is, assigns to $\tilde f(q)$ the coefficient of $e^{iqx}$ (what we used to call $c_n$ for Fourier series) when $f(x)$ is expanded in the $\{ e^{iqx} \}$ basis. 

# We therefore speak of the space of $q$-values as "frequency-space" (for functions of time) or "wavenumber space" (for functions of position). More generally, we can refer to $q$-space as ***reciprocal space*** or ***Fourier space***.
# 

# ## Parseval's Theorem for the Fourier transform

# 
# Let's find an analog of Parseval's Theorem,
# 
# $$ 
# \frac{1}{2\pi} \int_{-\pi}^\pi |f(x)|^2 = \sum_{n=-\infty}^\infty |c_n|^2
# $$
# 
# for Fourier transforms. 

# 
# First note that $|f(x)|^2 = f(x) (f(x))^*$. 
# We can take the complex conjugate of
# \begin{align*}
#     f(x) &= \int_{-\infty}^\infty dq \tilde f(q) e^{-iqx} 
# \end{align*}
# to obtain 
# \begin{align*}
#     (f(x))^* &= \int_{-\infty}^\infty dq (\tilde f(q))^* e^{+iqx} 
# \end{align*}
# So,
# \begin{align*}
#     \frac{1}{2\pi} \int_{-\infty}^\infty dx |f(x)|^2 &= 
#         \frac{1}{2\pi}  \int_{-\infty}^\infty dx  \int_{-\infty}^\infty dq \tilde f(q) e^{-iqx} 
#         \int_{-\infty}^\infty dq' (\tilde f(q'))^* e^{+iq'x} \\
#         &= \frac{1}{2\pi}  \int_{-\infty}^\infty dq \tilde f(q)
#                 \int_{-\infty}^\infty dq' (\tilde f(q'))^* 
#                 \int_{-\infty}^\infty dx e^{i(q'-q) x} \\
# \end{align*}                                            

# Now, from our discussion above about the Dirac delta function, we know 
# $ \dfrac{1}{2\pi}   \int_{-\infty}^\infty dx e^{i(q'-q)x} $ is just $\delta(q'-q)$.

# 
# Thus,
# \begin{align*}
#     \frac{1}{2\pi} \int_{-\infty}^\infty dx |f(x)|^2 &= 
#         \int_{-\infty}^\infty dq \tilde f(q) \int_{-\infty}^\infty dq' (\tilde f(q'))^* \delta(q'-q)\\      
# \end{align*}
# 

# Now integrating over $q'$, the $\delta$ function "picks out" the value of $(\tilde f(q'))^*$ where $q'-q=0$, that is, $q'=q$:                                                                                                     

# 
# 
# \begin{align*}
#     \frac{1}{2\pi} \int_{-\infty}^\infty dx |f(x)|^2
#         &= \int_{-\infty}^\infty dq \tilde f(q) (\tilde f(q))^* = \int_{-\infty}^\infty dq |\tilde f(q)|^2
# \end{align*}
# This is Parseval's Theorem for Fourier transforms. 

# ## Fourier transform of a derivative

# 
# Suppose we know the Fourier transform $\tilde f(q)$ of $f(x)$. What is the Fourier transform of $h(x) = \dfrac{df}{dx}$?
# 

# \begin{align*}
#     \tilde h(q) = \frac{1}{2\pi} \int_{-\infty}^\infty dx \dfrac{df}{dx} e^{iqx} = \frac{1}{2\pi} \left[ \int_{-\infty}^\infty dx \frac{d}{dx} \left(f(x) e^{iqx}\right) - \int_{-\infty}^\infty dx f(x) \frac{d}{dx} e^{iqx} \right]  \\
# \end{align*}
# Let's assume $f(x) \rightarrow 0$ as $x\rightarrow \pm \infty$:
# \begin{align*}
#        \tilde h(q) &= \frac{1}{2\pi} \left[ {\left. f(x) e^{iqx} \right|^\infty_{-\infty}}  - iq \int_{-\infty}^\infty dx f(x) e^{iqx} \right] \\
#        &= \frac{1}{2\pi} \left[ 0-0 - iq \int_{-\infty}^\infty dx f(x) e^{iqx} \right] \\
#        &= - i q \tilde f(q)
# \end{align*}

# ```{important}
# 
# $$
# h(x) = \dfrac{df}{dx} \Leftrightarrow \tilde h(q) = -iq \tilde f(q)
# $$ 
# 
# Differentiation in position-space becomes multiplication by $-iq$ in Fourier space!
# 
# ```

# We can iterate this result to get the Fourier transform of the $n$th derivative: Let $k(x) = dh/dx = d^2 f/dx^2$. Then
# 

# 
# $$
# \tilde k(q) = -iq \tilde h(q) = (-iq)^2 \tilde f(q) = -q^2 \tilde f(q)
# $$
# 
# 

# and more generally,

# ```{important}
# $$
# \text{ For $p(x) = d^n f/dx^n$, } \quad  \tilde p(q) = (-iq)^n \tilde f(q).
# $$
# ```

# This property can be very useful in solving differential equations "in Fourier space"... as long we can compute the inverse Fourier transform! 

# ### Example 1

# *Solve the differential equation $\dfrac{d^2 f}{dx^2} = - 9 f$.*
# 
# Taking the Fourier transforms of both sides, and exploiting linearity:

# 
# 
# \begin{align*}
#     (-iq)^2 \tilde f(q) &= - 9 \tilde f(q) \\
#     \Rightarrow (9-q^2)\tilde f(q) &= 0 
# \end{align*}

# Since this equation is true for *all* $q$, $\tilde f(q)$ must vanish everywhere except $q=\pm 3$. So (for general complex coefficients $a$ and $b$) the general solution satisfies

# \begin{align*}
# \tilde f(q) &= a \delta(q+3) + b \delta(q-3) \\
#     \Rightarrow f(x) &= \int_{-\infty}^\infty dq \left(a \delta(q+3) + b\delta(q-3) \right) e^{-iqx} = a e^{3ix} + b e^{-3ix} \\
#     &= A \cos(3x) + i B \sin(3x)
# \end{align*}
# where $A=a+b$, $B=a-b$. 

# ### Example 2

# 
# \begin{align*}
#     \frac{df}{dx} &= \cos(3x) \\
#     \Rightarrow -iq \tilde f(q) &= \frac{1}{2} \left[ \delta(q+3) + \delta(q-3) \right] \\
#     \Rightarrow \tilde f(q) &= - \frac{1}{2iq} \left[ \delta(q+3) + \delta(q-3) \right] \\
#     \Rightarrow f(x) &= \int_{-\infty}^\infty dq \left(-\frac{1}{2iq} \left[\delta(q+3) + \delta(q-3) \right] \right) e^{-iqx} \\
#     &= \left(-\frac{1}{2i} \right) \left[ \int_{-\infty}^\infty dq \frac{e^{-iqx}}{q} \delta(q+3) + \int_{-\infty}^\infty dq \frac{e^{-iqx}}{q} \delta(q-3) \right]  \\
#     &= \left(-\frac{1}{2i} \right) \left[ \frac{e^{3ix}}{-3} + \frac{e^{-3ix}}{3} \right] = \frac{1}{3} \sin(3x) \quad \checkmark 
# \end{align*}
# 

# ## Summary
# 
# * Fourier transforms generalize Fourier series (along with Parseval's Theorem) to not-necessarily-periodic functions on $(-\infty,\infty)$
# 
# * The Fourier transform of a function of time picks out the frequencies (maybe infinitely many!) of the SHOs that sum to give the function.
# 
# * Differentiation for $f(x)$ become multiplication by $(iq)$ in Fourier space.
# 
# * We have also introduced the very important Dirac delta function $\delta(u)$. 

# $\vdots$  
# $\vdots$   
# $\vdots$  
# $\vdots$  
