# Método de Newton

El método de Newton se utiliza para entontrar raices de una función de forma iterativa. Cuando la función va de $F:\mathbb{R}^{n+1}\rightarrow \mathbb{R}^{n}$ la iteración está dada por:

$$ \left(\begin{array}{c} x^{i+1}_{1} \\ \vdots \\ x^{i+1}_{n} \\ p_{i+1} \end{array}\right) = \left(\begin{array}{c} x^{i}_{1} \\ \vdots \\ x^{i}_{n} \\ p_{i} \end{array}\right) - \dfrac{||F(\textbf{x}_{i},p_{i})||}{\nabla||F(\textbf{x}_{i},p_{i})||\cdot \nabla||F(\textbf{x}_{i},p_{i})||}\nabla||F(\textbf{x}_{i},p_{i})|| $$

# Continuación por Pseudo Longitud de Arco.

Este método es un método de continuación que predice el avance en una rama de equilibrioa partir de un punto inicial. Este consta de resolver el siguiente sistema:

$$F(\textbf{x}_{i+1},p_{i+1}) = 0$$
$$(\textbf{x}_{i+1}-\textbf{x}_{i})\cdot \left.\dfrac{d\textbf{x}}{ds}\right|_{s_{i}} + (p_{i+1} - p_{i}) \left.\dfrac{dp}{ds}\right|_{s_{i}} - \Delta s = 0$$

donde $s$ representa la longitud de arco de la rama de equilibrio, $\textbf{x}_{i} = \textbf{x}(s_{i})$, $p_{i} = p(s_{i})$ y $\Delta s = \sqrt{(\textbf{x}_{i+1} - \textbf{x}_{i})\cdot (\textbf{x}_{i+1} - \textbf{x}_{i}) + (p_{i+1}-p_{i})^{2}}$

La solución a este sistema [Doedel] nos proporciona este sistema de ecuaciones que podemos resolver de forma automática en Julia:

$$
\left(\begin{array}{cc} 
F_{\textbf{x}}(\textbf{x}_{i+1},p_{i+1}) & F_{p}(\textbf{x}_{i+1},p_{i+1}) \\ \left.\dfrac{d\textbf{x}}{ds}\right|_{s_{i}}^{T} & \left.\dfrac{dp}{ds}\right|_{s_{i}}
\end{array}\right)
\left(\begin{array}{c} 
\left.\dfrac{d\textbf{x}}{ds}\right|_{s_{+1}} \\ \left.\dfrac{dp}{ds}\right|_{s_{i+1}}
\end{array}\right)
=
\left(\begin{array}{c} 
\textbf{0} \\ 1
\end{array}\right)
 $$

Con esto, podemos asociar $\Delta \textbf{x}_{i+1} = \left.\dfrac{d\textbf{x}}{ds}\right|_{s_{+1}}$ y $\Delta p_{i+1} = \left.\dfrac{dp}{ds}\right|_{s_{+1}}$ y de esta manera tener $\textbf{x}_{i+2} = \textbf{x}_{i+1} + \Delta \textbf{x}_{i+1}$ y $p_{i+2} = p_{i+1} + \Delta p_{i+1}$.

Para obtener el paso inicial $i = 0$, se utiliza la derivada de la función implícita explicada en el documento de antecedentes.

# Puntos Límite

Para encontrar puntos límite necesitamos resolver el sistema 

$$
S(\textbf{x},p,\nu) = \begin{cases}
F(\textbf{x},p) = 0\\
F_{\textbf{x}}(\textbf{x},p)\nu = 0\\
\nu_{0}\cdot \nu -1 = 0
\end{cases}
$$

donde $\nu $ es un eigenvalor de $F$ en $x$ y $p$.

Se puede hacer una aproximación al punto LP mediante el método de Newton, pero como ahora tenemos un sistema invertible:

$$
\begin{pmatrix}
\textbf{x}_{i+1}\\
p_{i+1}\\
\nu_{i+1}
\end{pmatrix} = 
\begin{pmatrix}
\textbf{x}_{i}\\
p_{i}\\
\nu_{i}
\end{pmatrix} - 
\left.\begin{pmatrix}
\dfrac{\partial S_{1}}{\partial \textbf{x}} & \dfrac{\partial S_{1}}{\partial p} & \dfrac{\partial S_{1}}{\partial \nu} \\
\dfrac{\partial S_{2}}{\partial \textbf{x}} & \dfrac{\partial S_{2}}{\partial p} & \dfrac{\partial S_{2}}{\partial \nu} \\
\dfrac{\partial S_{3}}{\partial \textbf{x}} & \dfrac{\partial S_{3}}{\partial p} & \dfrac{\partial S_{3}}{\partial \nu} \\
\end{pmatrix}\right|_{\textbf{x}_{i},p_{i},\nu_{i}}^{-1} \cdot S(\textbf{x}_{i},p_{i},\nu_{i})
$$

donde:

- $\frac{\partial S_{1}}{\partial \textbf{x}} (\textbf{x}_{i},p_{i},\nu_{i}) = F_{\textbf{x}}(\textbf{x}_{i},p_{i})$
- $\frac{\partial S_{1}}{\partial p} (\textbf{x}_{i},p_{i},\nu_{i}) = F_{p}(\textbf{x}_{i},p_{i})$
- $\frac{\partial S_{1}}{\partial \nu} (\textbf{x}_{i},p_{i},\nu_{i}) = 0_{n\times n}$

- $\frac{\partial S_{2}}{\partial \textbf{x}} (\textbf{x}_{i},p_{i},\nu_{i})_{i,j} = \frac{\partial}{\partial x_{j}} \sum_{k=1}^{n} \frac{\partial F_{i}}{\partial x_{k}}\nu_{k} = \sum_{k=1}^{n} \frac{\partial^{2} F_{i}}{\partial x_{j}\partial x_{k}}\nu_{k}  $
- $\frac{\partial S_{2}}{\partial p} (\textbf{x}_{i},p_{i},\nu_{i})_{i} = \sum_{k=1}^{n} \frac{\partial^{2} F_{i}}{\partial p\partial x_{k}}\nu_{k}$
- $\frac{\partial S_{2}}{\partial \nu} (\textbf{x}_{i},p_{i},\nu_{i})_{i,j} = \frac{\partial}{\partial \nu_{j}} \sum_{k=1}^{n} \frac{\partial F_{i}}{\partial x_{k}}\nu_{k} = \sum_{k=1}^{n} \frac{\partial F_{i}}{\partial x_{k}}\delta_{j,k} = \frac{\partial F_{i}}{\partial x_{j}}$

- $\frac{\partial S_{3}}{\partial \textbf{x}} (\textbf{x}_{i},p_{i},\nu_{i}) = 0_{1\times n}$
- $\frac{\partial S_{3}}{\partial p} (\textbf{x}_{i},p_{i},\nu_{i}) = 0$
- $\frac{\partial S_{3}}{\partial \nu} (\textbf{x}_{i},p_{i},\nu_{i}) = \nu_{0}^{T}$

# Bifurcaciones de Hopf

Para encontrar puntos límite necesitamos resolver el sistema 

$$
S(\textbf{x},p,\nu) = \begin{cases}
F(\textbf{x},p) = 0\\
(F_{\textbf{x}}(\textbf{x},p)^{2} + \kappa I_{n})v = 0\\
v \cdot v -1 = 0 \\
w \cdot v -1 = 0
\end{cases}
$$

donde $v $ es un eigenvalor de $F$ en $x$ y $p$.

Como en el caso anterior, usamos método de Newton

$$
\begin{pmatrix}
\textbf{x}_{i+1}\\
p_{i+1}\\
v_{i+1}\\
\kappa_{i+1}
\end{pmatrix} = 
\begin{pmatrix}
\textbf{x}_{i}\\
p_{i}\\
v_{i}\\
\kappa_{i}
\end{pmatrix} - 
\left.\begin{pmatrix}
\frac{\partial S_{1}}{\partial \textbf{x}} & \frac{\partial S_{1}}{\partial p} & \frac{\partial S_{1}}{\partial v} & \frac{\partial S_{1}}{\partial \kappa}\\
\frac{\partial S_{2}}{\partial \textbf{x}} & \frac{\partial S_{2}}{\partial p} & \frac{\partial S_{2}}{\partial v} & \frac{\partial S_{2}}{\partial \kappa} \\
\frac{\partial S_{3}}{\partial \textbf{x}} & \frac{\partial S_{3}}{\partial p} & \frac{\partial S_{3}}{\partial v} & \frac{\partial S_{3}}{\partial \kappa}\\
\frac{\partial S_{4}}{\partial \textbf{x}} & \frac{\partial S_{4}}{\partial p} & \frac{\partial S_{4}}{\partial v} & \frac{\partial S_{4}}{\partial \kappa} 
\end{pmatrix} \right|_{\textbf{x}_{i},p_{i},v_{i},\kappa_{i}}^{-1} * S(\textbf{x}_{i},p_{i},v_{i},\kappa_{i})
$$

donde:

- $\frac{\partial S_{1}}{\partial \textbf{x}} = F_{\textbf{x}}$
- $\frac{\partial S_{1}}{\partial p} = F_{p}$
- $\frac{\partial S_{1}}{\partial v} = 0_{n\times n}$
- $\frac{\partial S_{1}}{\partial \kappa} = 0$

- $\frac{\partial S_{2}}{\partial \textbf{x}} = \frac{\partial}{\partial x_{j}}\sum_{k}\left(\sum_{l}\frac{\partial F_{i}}{\partial x_{l}}\frac{\partial F_{l}}{\partial x_{k}} + \kappa\delta_{i,k})\right)v_{k} = \sum_{k}\sum_{l}\frac{\partial}{\partial x_{j}}\left(\frac{\partial F_{i}}{\partial x_{l}}\frac{\partial F_{l}}{\partial x_{k}}\right)v_{k}$
- $\frac{\partial S_{3}}{\partial \textbf{x}} = \sum_{k}\sum_{l}\frac{\partial}{\partial p}\left(\frac{\partial F_{i}}{\partial x_{l}}\frac{\partial F_{l}}{\partial x_{k}}\right)v_{k}$
- $\frac{\partial S_{2}}{\partial v} = \frac{\partial}{\partial v_{j}}\sum_{k}\left(\sum_{l}\frac{\partial F_{i}}{\partial x_{l}}\frac{\partial F_{l}}{\partial x_{k}} + \kappa\delta_{i,k})\right)v_{k} = \sum_{k}\left(\sum_{l}\frac{\partial F_{i}}{\partial x_{l}}\frac{\partial F_{l}}{\partial x_{k}} + \kappa\delta_{i,k}\right)\delta_{k,j} = \sum_{l}\frac{\partial F_{i}}{\partial x_{l}}\frac{\partial F_{l}}{\partial x_{j}} + \kappa \delta_{i,j}$
- $\frac{\partial S_{2}}{\partial \kappa} = \frac{\partial}{\partial \kappa}\sum_{k}\left(\sum_{l}\frac{\partial F_{i}}{\partial x_{l}}\frac{\partial F_{l}}{\partial x_{k}} + \kappa\delta_{i,k}\right)v_{k} = \frac{\partial}{\partial \kappa}\sum_{k}\kappa\delta_{i,k}v_{k} = v_{i}$

- $\frac{\partial S_{3}}{\partial \textbf{x}} = 0_{1\times n}$
- $\frac{\partial S_{3}}{\partial p} = 0$
- $\frac{\partial S_{3}}{\partial v} = \frac{\partial}{\partial v_{j}}\sum_{k}v_{k}^{2} = \sum_{k}2v_{k}\delta_{k,j} = 2v_{j} = 2v^{T}$
- $\frac{\partial S_{3}}{\partial \kappa} = 0$

- $\frac{\partial S_{1}}{\partial \textbf{x}} = 0_{1\times n}$
- $\frac{\partial S_{1}}{\partial p} = 0$
- $\frac{\partial S_{1}}{\partial v} = \frac{\partial}{\partial v_{j}}\sum_{k}w_{k}v_{k} = \sum_{k}w_{k}\delta_{k,j} = w_{j} = w^{T}$
- $\frac{\partial S_{1}}{\partial \kappa} = 0$



