# Órbitas Periódicas

En un sistema dinámico, una órbita periódica es aquella donde la solución $x(t)$ a la ecuación diferencial $\tfrac{dx}{dt} = F(x(t), p)$ es una función periódica, es decir, $x(t) = x(t + T)$ para todo $t \in \mathbb{R}$ donde $T$ es el periodo de $x(t)$.

En otras palabras, se cumple que la integral 

$$ \int_{0}^{T} F(x(t), p)\ dt = 0 $$

# Implementación numérica

### Continuación de las ramas periódicas

Como en la continuación de las ramas de equilibrio, nuevamente implementaremos el método de Newton sobre un conjunto de ecuaciones que incluyen el método de Pseudo Longitud de Arco.

Las ecuaciones a resolver son las siguientes:

$$

G(x_{ini}, p, T) = \left\{ \begin{matrix} 

\int_{0}^{T} F(x(t), p)\ dt = x(T) - x(0) = x(x_{ini}, p, T) - x_{ini} \\ \\

(x_{ini} - x_{ini}^{0}) \cdot \dfrac{dx_{ini}^{0}}{dt} = (x_{ini} - x_{ini}^{0}) \cdot F(x_{ini}^{0}, p^{0}) \\ \\

(x_{ini} - x_{ini}^{0}) \cdot \dfrac{dx_{ini}^{0}}{ds} + (p - p^{0}) \dfrac{dp^{0}}{ds} + (T - T^{0})\dfrac{dT^{0}}{ds} - \Delta s

\end{matrix}\right. 

$$

$$ G : (D_{x} \times D_{p} \times D_{T}) \subset \mathbb{R}^{n + 2} \rightarrow \mathbb{R}^{n+2}$$

Si consideramos la siguiente notación:

$$s_{k} = k\Delta s$$

$$x_{ini}^{k} = x_{ini}(s_{k})$$

$$p^{k} = p(s_{k})$$

$$T^{k} = T(s_{k}) $$

$$k \in 0, 1, ..., N$$

donde $N$ es el número de órbitas periódicas, podemos encontrar cada una resolviendo las ecuaciones:

$$

\left\{ \begin{matrix} 

x(x_{ini}^{k}, p^{k}, T^{k}) - x^{k}_{ini} = 0 \\ \\

(x^{k}_{ini} - x_{ini}^{k-1}) \cdot \dfrac{dx_{ini}^{k-1}}{dt} = 0 \\ \\

(x_{ini}^{k} - x_{ini}^{k-1}) \cdot \dfrac{dx_{ini}^{k-1}}{ds} + (p^{k} - p^{k-1}) \dfrac{dp^{k-1}}{ds} + (T^{k} - T^{k-1})\dfrac{dT^{k-1}}{ds} - \Delta s = 0

\end{matrix}\right. 

$$

mediante el Método de Newton:

$$

\begin{pmatrix}

\Delta x_{ini}^{k^{\nu}} \\ \\ \\
\Delta p^{k^{\nu}} \\ \\ \\
\Delta T^{k^{\nu}}

\end{pmatrix} = - \begin{pmatrix}

\dfrac{\partial x}{\partial x_{ini}}(x_{ini}^{k^{\nu}}, p^{k^{\nu}}, T^{k^{\nu}}) - \dfrac{\partial x_{ini}}{\partial x_{ini}} && \dfrac{\partial x}{\partial p}(x^{k^{\nu}}_{ini}, p^{k^{\nu}}, T^{k^{\nu}}) && \dfrac{\partial x}{\partial T}(x^{k^{\nu}}_{ini}, p^{k^{\nu}},T^{k^{\nu}}) \\ \\

\dfrac{dx_{ini}^{k-1}}{dt}^{*} && 0 && 0  \\ \\

\dfrac{dx^{k-1}_{ini}}{ds}^{*} && \dfrac{dp^{k-1}}{ds} && \dfrac{dT^{k-1}}{ds}

\end{pmatrix} \begin{pmatrix} 

x(x_{ini}^{k^{\nu}}, p^{k^{\nu}}, T^{k^{\nu}}) - x^{k^{\nu}}_{ini} \\ \\

(x^{k^{\nu}}_{ini} - x_{ini}^{k-1}) \cdot \dfrac{dx_{ini}^{k-1}}{dt} \\ \\

(x_{ini}^{k^{\nu}} - x_{ini}^{k-1}) \cdot \dfrac{dx_{ini}^{k-1}}{ds} + (p^{k^{\nu}} - p^{k-1}) \dfrac{dp^{k-1}}{ds} + (T^{k^{\nu}} - T^{k-1})\dfrac{dT^{k-1}}{ds} - \Delta s

\end{pmatrix}

$$

$$

= - \begin{pmatrix}

\dfrac{\partial x}{\partial x_{ini}}(x_{ini}^{k^{\nu}}, p^{k^{\nu}}, T^{k^{\nu}}) - I_{n} && \dfrac{\partial x}{\partial p}(x_{ini}^{k^{\nu}}, p^{k^{\nu}}, T^{k^{\nu}}) && F(x(x_{ini}^{k^{\nu}}, p^{k^{\nu}}, T^{k^{\nu}}), p^{k^{\nu}}) \\ \\

F(x_{ini}^{k-1}, p^{k-1})^{*} && 0 && 0  \\ \\

\dfrac{dx^{k-1}_{ini}}{ds}^{*} && \dfrac{dp^{k-1}}{ds} && \dfrac{dT^{k-1}}{ds}

\end{pmatrix} \begin{pmatrix} 

x(x_{ini}^{k^{\nu}}, p^{k^{\nu}}, T^{k^{\nu}}) - x^{k^{\nu}}_{ini} \\ \\

(x^{k^{\nu}}_{ini} - x_{ini}^{k-1}) \cdot \dfrac{dx_{ini}^{k-1}}{dt} \\ \\

(x_{ini}^{k^{\nu}} - x_{ini}^{k-1}) \cdot \dfrac{dx_{ini}^{k-1}}{ds} + (p^{k^{\nu}} - p^{k-1}) \dfrac{dp^{k-1}}{ds} + (T^{k^{\nu}} - T^{k-1})\dfrac{dT^{k-1}}{ds} - \Delta s

\end{pmatrix}

$$

$$ x^{k^{\nu + 1}}_{ini} = x^{k^{\nu}}_{ini} + \Delta x^{k^{\nu}}_{ini} $$

$$ p^{k^{\nu + 1}} = p^{k^{\nu}} + \Delta p^{k^{\nu}} $$

$$ T^{k^{\nu + 1}} = T^{k^{\nu}} + \Delta T^{k^{\nu}} $$

$$ x^{k^{0}} = x^{k-1} + \dfrac{dx_{ini}^{k-1}}{ds}\Delta s $$

$$ p^{k^{0}} = p^{k-1} + \dfrac{dp^{k-1}}{ds}\Delta s $$

$$ T^{k^{0}} = T^{k-1} + \dfrac{dT^{k-1}}{ds}\Delta s $$

donde $\nu$ es el número de iteraciones.

Al iniciar en una bifurcación de Hopf $(x_{ini}^{0}, p^{0})$, tenemos que

$$ T^{0} = \dfrac{2\pi}{\text{Eigenval}(F_{x}(x_{ini}^{0}, p^{0}))} $$

$$ \dfrac{dx_{ini}^{0}}{dt} = \dfrac{d \phi}{dt}(0) $$

$$ \dfrac{dx_{ini}^{0}}{ds} = \dfrac{\phi(0)}{||\phi(0)||} $$

$$ \dfrac{dp^{0}}{ds} = 0 $$

$$ \dfrac{dT^{0}}{ds} = 0 $$

donde

$$ \dfrac{d\phi}{dt} = \dfrac{\partial F}{dx}(x_{ini}^{0}, p^{0})\ \phi (t) $$

Entonces las ecuaciones a resolver para encontrar la nueva dirección es:

$$

\begin{pmatrix}

\dfrac{\partial x}{\partial x_{ini}}(x_{ini}^{k}, p^{k}, T^{k}) - I_{n} && \dfrac{\partial x}{\partial p}(x_{ini}^{k}, p^{k}, T^{k}) && F(x(x_{ini}^{k}, p^{k}, T^{k}), p^{k}) \\ \\

F(x_{ini}^{k-1}, p^{k-1})^{*} && 0 && 0  \\ \\

\dfrac{dx^{k-1}_{ini}}{ds}^{*} && \dfrac{dp^{k-1}}{ds} && \dfrac{dT^{k-1}}{ds}

\end{pmatrix} \begin{pmatrix}

\dfrac{dx^{k}_{ini}}{ds} \\ \\ \dfrac{dp^{k}}{ds} \\ \\ \dfrac{dT^{k}}{ds}

\end{pmatrix} = \begin{pmatrix}

\\ 0_{n+1\times 1} \\ \\ \\ \\ 1 \\ \\

\end{pmatrix} $$

$$

\left|\left|\dfrac{dx^{k}_{ini}}{ds}\right|\right|^{2} + \left|\dfrac{dp_{k}}{ds}\right|^{2} + \left|\dfrac{dT_{k}}{ds}\right|^{2} = 1

$$