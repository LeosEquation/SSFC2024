# Ramas de equilibrio

En un sistema dinámico, una solución de equilibrio es aquella donde su derivada respecto al tiempo es nula, es decir, los valores de $x$ y $p$ donde la ecuación diferencial $F(x, p) = 0$.

El teorema de la función inversa asegura que es posible establecer una función diferenciable $x_{eq}: D_{p}\subset\mathbb{R}^{m} \rightarrow \mathbb{R}^{n}$ ($p \mapsto x_{eq}$) tal que $F(x_{eq}(p), p) = 0$ para todo $p \in D_{p}$, entonces la derivada de $x_{eq}$ respecto a $p$ es

$$ \dfrac{dx_{eq}}{dp} = -\left(\dfrac{\partial F}{\partial x}\right)^{-1} \dfrac{\partial F}{\partial p} $$

A esta función $x_{eq}$ se le conoce como rama de equilibrio.

Es posible que existan 2 o más soluciones de equilibrio para un mismo valor del parámetro $p$ en el sistema descrito por $F$. Cada una de estas soluciones pertenece a una rama de equilibrio distinta.

Existen puntos donde la solución de equilibrio pertenece a 2 ramas de equilibrio distintas. A esto se le conoce como Nodo de Saddle y forma parte de la Teoría de Bifurcaciones.

# Implementación numérica

### Continuación de las ramas de equilibrio

Para encontrar numéricamente las ramas de equilibrio, se implementa el método de Newton para encontrar los puntos de equilibrio de $F(x,p)$, sin embargo el Jacobiano no es cuadrado al considerar el parámetro $p$ por lo que no se puede implementar de manera eficiente el método. La solución a este problema es implementar una condición extra para cuadratizarlo. La condición extra que se implementará en este proyecto es el método de Pseudo Longitud de Arco, que controla la distancia entre cada punto de la continuación y permite implementarla a pesar de la existenca de bifurcaciones de Nodo de Saddle.

Entonces, al considerar el sistema dinámico y la condición de continuación, obtenemos la siguiente función:

$$ G(x, p) = \left\{ \begin{matrix} 

F(x,p) \\

(x - x_{0}) \cdot \tfrac{dx_{0}}{ds} + (p - p_{0}) \tfrac{dp_{0}}{ds} - \Delta s

\end{matrix}\right. $$

$$ G : (D_{x} \times D_{p}) \subset \mathbb{R}^{n + 1} \rightarrow \mathbb{R}^{n+1}$$

Las ramas de equilibrio se pueden parametrizar respecto a la longitud de arco $s$, $F(x_{eq}(s), p(s)) = 0$ para todo $s \in \mathbb{R}$, donde $x_{0} = x_{eq}(s = 0)$ y $p_{0} = p(s = 0)$ es un punto de equilibrio de referencia calculado analíticamente o encontrado numéricamente.

Como $F$ es constante para todo $s$, tenemos que

$$ \dfrac{dF}{ds}(x_{0}, p_{0}) = 0 $$

Es decir

$$ \dfrac{\partial F}{\partial x}\dfrac{dx_{0}}{ds} + \dfrac{\partial F}{\partial p}\dfrac{dp_{0}}{ds} = 0_{n\times 1} $$

$$ \begin{pmatrix}

\dfrac{\partial F}{\partial x} && \dfrac{\partial F}{\partial p}

\end{pmatrix}_{n \times n + 1} \begin{pmatrix}

\dfrac{dx_{0}}{ds} \\ \\

\dfrac{dp_{0}}{ds}

\end{pmatrix}_{n + 1 \times 1} = 0_{n \times 1}

$$

Entonces, el cálculo de la derivada respecto a la longitud de arco en $x_{0}$ y $p_{0}$ está dado por el espacio nulo del jacobiano completo de $F$ respecto a $x$ y $p$. Las derivadas deben estar normalizadas para que el método funcione, es decir:

$$

\left|\left|\dfrac{dx_{0}}{ds}\right|\right|^{2} + \left|\dfrac{dp_{0}}{ds}\right|^{2} = 1

$$

No es necesario encontrar el espacio nulo en cada paso de la continuación, basta con encontrar el primero y predecir los siguientes resolviendo un sistema de ecuaciones. Esto resulta más eficiente que calcular el espacio nulo en cada paso. Para ello obtenemos la derivada de $G$ respecto a $s$, recordando que $\Delta s$ puede verse como $s - s_{0}$:

$$ \dfrac{dG}{ds}(x_{eq},p) = \begin{pmatrix}

\dfrac{\partial F}{\partial x}(x_{eq}, p) && \dfrac{\partial F}{\partial p}(x_{eq}, p)  \\ \\

\dfrac{dx_{0}}{ds}^{*} && \dfrac{dp_{0}}{ds}

\end{pmatrix}_{n+1\times n+1} \begin{pmatrix}

\dfrac{dx_{eq}}{ds} \\ \\ \dfrac{dp}{ds}

\end{pmatrix}_{n+1\times 1} - \begin{pmatrix}

\\ 0_{n\times 1} \\ \\ 1 \\ \\

\end{pmatrix} $$

Como $G$ también debe ser $0$ para todo $s$, su derivada es también cero:

$$ \begin{pmatrix}

\dfrac{\partial F}{\partial x} (x_{eq},p) && \dfrac{\partial F}{\partial p}(x_{eq},p) \\ \\

\dfrac{dx_{0}}{ds}^{*} && \dfrac{dp_{0}}{ds}

\end{pmatrix}_{n+1\times n+1} \begin{pmatrix}

\dfrac{dx_{eq}}{ds} \\ \\ \dfrac{dp}{ds}

\end{pmatrix}_{n+1\times 1} - \begin{pmatrix}

\\ 0_{n\times 1} \\ \\ 1 \\ \\

\end{pmatrix} = 0_{n+1 \times 1}$$

Entonces las ecuaciones a resolver para encontrar la nueva dirección del método de continuación es:

$$\begin{pmatrix}

\dfrac{\partial F}{\partial x}(x_{eq},p) && \dfrac{\partial F}{\partial p}(x_{eq},p)  \\ \\

\dfrac{dx_{0}}{ds}^{*} && \dfrac{dp_{0}}{ds}

\end{pmatrix}_{n+1\times n+1} \begin{pmatrix}

\dfrac{dx_{eq}}{ds} \\ \\ \dfrac{dp}{ds}

\end{pmatrix}_{n+1\times 1} = \begin{pmatrix}

\\ 0_{n\times 1} \\ \\ 1 \\ \\

\end{pmatrix} $$

$$

\left|\left|\dfrac{dx_{eq}}{ds}\right|\right|^{2} + \left|\dfrac{dp}{ds}\right|^{2} = 1

$$

Si consideramos la siguiente notación:

$$s_{k} = k\Delta s$$

$$x_{k} = x_{eq}(s_{k})$$

$$p_{k} = p(s_{k})$$

$$k \in 0, 1, ..., N$$

donde $N$ es el número de puntos de la continuación, tenemos que a partir de $x_{0}$ y $p_{0}$ podemos encontrar numéricamente $N$ puntos resolviendo para todo $k \geq 1$ el siguiente sistema de ecuaciones:

$$

\left\{ \begin{matrix} 

F(x_{k},p_{k}) = 0 \\

(x_{k} - x_{k-1}) \cdot \tfrac{dx_{k-1}}{ds} + (p_{k} - p_{k-1}) \tfrac{dp_{k-1}}{ds} - \Delta s

\end{matrix}\right.

$$

mediante el Método de Newton:

$$

\begin{pmatrix}

\Delta x_{k}^{\nu} \\ \\ \\
\Delta p_{k}^{\nu}

\end{pmatrix} = - \begin{pmatrix}

\dfrac{\partial F}{\partial x}(x_{k}^{\nu}, p_{k}^{\nu}) && \dfrac{\partial F}{\partial p}(x_{k}^{\nu}, p_{k}^{\nu})  \\ \\

\dfrac{dx_{k-1}}{ds}^{*} && \dfrac{dp_{k-1}}{ds}

\end{pmatrix}_{n+1\times n+1}^{-1} \begin{pmatrix} 

F(x_{k}^{\nu},p_{k}^{\nu}) \\ \\ \\

(x_{k}^{\nu} - x_{k-1}) \cdot \tfrac{dx_{k-1}}{ds} + (p_{k}^{\nu} - p_{k-1}) \tfrac{dp_{k-1}}{ds} - \Delta s

\end{pmatrix}_{n+1\times 1}

$$

$$ x_{k}^{\nu + 1} = x_{k}^{\nu} + \Delta x_{k}^{\nu} $$

$$ p_{k}^{\nu + 1} = p_{k}^{\nu} + \Delta p_{k}^{\nu} $$

donde $\nu$ es el número de iteraciones y

$$ 

\begin{pmatrix}

\dfrac{\partial F}{\partial x}(x_{0}, p_{0}) && \dfrac{\partial F}{\partial p}(x_{0}, p_{0})

\end{pmatrix}_{n \times n + 1} \begin{pmatrix}

\dfrac{dx_{0}}{ds} \\ \\

\dfrac{dp_{0}}{ds}

\end{pmatrix}_{n + 1 \times 1} = 0_{n \times 1}

$$

$$

\left|\left|\dfrac{dx_{0}}{ds}\right|\right|^{2} + \left|\dfrac{dp_{0}}{ds}\right|^{2} = 1

$$

$$

\begin{pmatrix}

\dfrac{\partial F}{\partial x}(x_{k},p_{k}) && \dfrac{\partial F}{\partial p}(x_{k},p_{k})  \\ \\

\dfrac{dx_{k-1}}{ds}^{*} && \dfrac{dp_{k-1}}{ds}

\end{pmatrix}_{n+1\times n+1} \begin{pmatrix}

\dfrac{dx_{k}}{ds} \\ \\ \dfrac{dp_{k}}{ds}

\end{pmatrix}_{n+1\times 1} = \begin{pmatrix}

\\ 0_{n\times 1} \\ \\ 1 \\ \\

\end{pmatrix} $$

$$

\left|\left|\dfrac{dx_{k}}{ds}\right|\right|^{2} + \left|\dfrac{dp_{k}}{ds}\right|^{2} = 1

$$























