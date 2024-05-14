# Ramas de equilibrio en modelos de ecuaciones diferenciales

Sea $F: D_{F} \subseteq \mathbb{R}^{n}\times\mathbb{R} \rightarrow \mathbb{R}^{n}$ una función que define un modelo de ecuaciones diferenciales tal que $F(x,p) = \dfrac{d\textbf{x}}{dt}$ donde $\textbf{x}: D_{\textbf{t}}\subseteq \mathbb{R} \rightarrow \mathbb{R}^{n}$, $t\in D_{\textbf{t}}$ y $ p \in \mathbb{R}$ con $n \in \mathbb{N}$. 

Entonces, $(\textbf{x},p)$ es un punto de equilibrio si $F(x,p) = 0$.

Si existe una función $\textbf{x}:D_{\textbf{p}} \subset \mathbb{R} \rightarrow \mathbb{R}^{n}$ continuamente diferenciable tal que $F(\textbf{x}(p),p) = 0$ para toda $p \in D_{\textbf{p}}$, entonces $\textbf{x}$ será una rama de equilibrio de $F$.

El Teorema de la Función Implícita asegura que existen ramas de equilibrio y estas son únicas si la función $F$ satisface que [Doedel]

- $F(\textbf{x}_{0},p_{0}) = 0$ con $\textbf{x}_{0} \in \mathbb{R}^{n}$ y $p_{0} \in \mathbb{R}$

- La matriz jacobiana $F_{\textbf{x}}(\textbf{x}_{0},p_{0})$ debe tener una matriz inversa acotada, es decir, para algún $M > 0$

$$||F_{\textbf{x}}(\textbf{x}_{0},p_{0})^{-1}|| \leq M$$

- $F(\textbf{x}_{0},p_{0})$ y $F_{\textbf{x}}(\textbf{x}_{0},p_{0})$ son continuamente diferenciables en $D_{F}$

Además, este teorema asegura que

$$\dfrac{d\textbf{x}}{dp} = F_{\textbf{x}}(x(p),p)^{-1}F_{p}(x(p),p)$$

# Estabilidad en las ramas de equilibrio

Un punto estable sucede cuando en el punto $\textbf{x}_{0}$ y $p_{0}$ la parte real de los eigenvalores $\lambda$ del jacobiano $F_{\textbf{x}}(\textbf{x}_{0},p_{0})$ son negativas y un punto inestable sucede cuando al menos uno es positivo.

# Bifurcacion de punto límite en un modelo de ecuaciones diferenciales

Una bifurcación de punto límite (LP por sus siglas en inglés) sucede cuando en el punto $\textbf{x}_{0}$ y $p_{0}$ los eigenvalores $\lambda$ del jacobiano $F_{\textbf{x}}(\textbf{x}_{0},p_{0})$ son cero. En este puntp puede haber 2 ramas de equilibrio. 

Cuando $n = 1$ estas bifurcaciones se pueden clasificar [Caitlin] en:

- Nodo de Sadle: Sucede cuando $\tfrac{\partial F}{\partial p} \neq 0$ y $\tfrac{\partial^{2} F}{\partial x^{2}} \neq 0$ 

- Transcrítica: Sucede cuando $\tfrac{\partial F}{\partial p} = 0$ y $\tfrac{\partial^{2} F}{\partial x^{2}} \neq 0$ 

- Pitchfork: Sucede cuando $\tfrac{\partial F}{\partial p} = 0$, $\tfrac{\partial^{2} F}{\partial p^{2}} = 0$, $\tfrac{\partial^2 F}{\partial p \partial x} \neq 0$, $\tfrac{\partial^{3} F}{\partial p^2 \partial x} =\neq 0$ y $\tfrac{\partial^{3} F}{\partial x^3} =\neq 0$

# Ramas periódicas en un modelo de ecuaciones diferenciales

Una solución periódica sucede cuando existe un tiempo $T \in D_{\textbf{t}}$ que satisface:

$$\dfrac{d\textbf{x}}{dt} = F(\textbf{x},p)$$

$$\textbf{x}(t) = \textbf{x}(t + T)$$

Entonces una rama periódica, de manera similar a la de equilibrio, es una función continua $\textbf{x}: D_{p}\subseteq \mathbb{R}\rightarrow\mathbb{R}^{n}$ tal que la solución de $F(\textbf{x}(p),p)$ es periódica para toda $p\in D_{p}$

# Bifurcaciones de Hopf

Una bifurcación de Hopf sucede cuando en el punto $\textbf{x}_{0}$ y $p_{0}$ los eigenvalores $\lambda$ del jacobiano $F_{\textbf{x}}(\textbf{x}_{0},p_{0})$ son puramente imaginarios, es decir $\lambda = \pm i \beta$ con $\beta \in \mathbb{R^{n}}$.

Una bifurcación de Hopf se caracteriza por ser un punto de equilibrio y ser un punto con solución periódica, es decir, es la intersección entre una rama de equilibrio y una rama periódica.




