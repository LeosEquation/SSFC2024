# Método de Newton

El método de Newton se utiliza para entontrar raices de una función de forma iterativa. Cuando la función va de $F:\mathbb{R}^{n+1}\rightarrow \mathbb{R}^{n}$ la iteración está dada por:

$$ \left(\begin{array}{c} x^{i+1}\_{1} \\\ \vdots \\\ x^{i+1}\_{n} \\\ p\_{i+1} \end{array}\right) = \left(\begin{array}{c} x^{i}\_{1} \\\ \vdots \\\ x^{i}\_{n} \\\ p\_{i} \end{array}\right) - \dfrac{||F(\textbf{x}\_{i},p\_{i})||}{\nabla||F(\textbf{x}\_{i},p\_{i})||\cdot \nabla||F(\textbf{x}\_{i},p\_{i})||}\nabla||F(\textbf{x}\_{i},p\_{i})|| $$

# Continuación por Pseudo Longitud de Arco.

Este método es un método de continuación que predice el avance en una rama de equilibrioa partir de un punto inicial. Este consta de resolver el siguiente sistema:

$$F(\textbf{x}\_{i+1},p\_{i+1}) = 0$$
$$(\textbf{x}\_{i+1}-\textbf{x}\_{i})\cdot \left.\dfrac{d\textbf{x}}{ds}\right|\_{s\_{i}} + (p\_{i+1} - p\_{i}) \left.\dfrac{dp}{ds}\right|\_{s\_{i}} - \Delta s = 0$$

donde $s$ representa la longitud de arco de la rama de equilibrio, $\textbf{x}\_{i} = \textbf{x}(s\_{i})$, $p\_{i} = p(s\_{i})$ y $\Delta s = \sqrt{(\textbf{x}\_{i+1} - \textbf{x}\_{i})\cdot (\textbf{x}\_{i+1} - \textbf{x}\_{i}) + (p\_{i+1}-p\_{i})^{2}}$

La solución a este sistema [Doedel] nos proporciona este sistema de ecuaciones que podemos resolver de forma automática en Julia:

$$
\begin{pmatrix}
F\_{\textbf{x}}(\textbf{x}\_{i+1},p\_{i+1}) & F\_{p}(\textbf{x}\_{i+1},p\_{i+1}) \\\ \left.\dfrac{d\textbf{x}}{ds}\right|\_{s\_{i}}^{T} & \left.\dfrac{dp}{ds}\right|\_{s\_{i}}
\end{pmatrix}
\begin{pmatrix} 
\left.\dfrac{d\textbf{x}}{ds}\right|\_{s\_{+1}} \\\ \left.\dfrac{dp}{ds}\right|\_{s\_{i+1}}
\end{pmatrix} = 
\begin{pmatrix} 
\textbf{0} \\\ 1
\end{pmatrix}
$$

Con esto, podemos asociar $\Delta \textbf{x}\_{i+1} = \left.\dfrac{d\textbf{x}}{ds}\right|\_{s\_{+1}}$ y $\Delta p\_{i+1} = \left.\dfrac{dp}{ds}\right|\_{s\_{+1}}$ y de esta manera tener $\textbf{x}\_{i+2} = \textbf{x}\_{i+1} + \Delta \textbf{x}\_{i+1}$ y $p\_{i+2} = p\_{i+1} + \Delta p\_{i+1}$.

Para obtener el paso inicial $i = 0$, se utiliza la derivada de la función implícita explicada en el documento de antecedentes.

# Puntos Límite

Para encontrar puntos límite necesitamos resolver el sistema 

$$
S(\textbf{x},p,\nu) = \begin{cases}
F(\textbf{x},p) = 0\\
F\_{\textbf{x}}(\textbf{x},p)\nu = 0\\
\nu\_{0}\cdot \nu -1 = 0
\end{cases}
$$

donde $\nu $ es un eigenvalor de $F$ en $x$ y $p$.

Se puede hacer una aproximación al punto LP mediante el método de Newton, pero como ahora tenemos un sistema invertible:

$$
\begin{pmatrix}
\textbf{x}\_{i+1}\\
p\_{i+1}\\
\nu\_{i+1}
\end{pmatrix} = 
\begin{pmatrix}
\textbf{x}\_{i}\\
p\_{i}\\
\nu\_{i}
\end{pmatrix} - 
\left.\begin{pmatrix}
\dfrac{\partial S\_{1}}{\partial \textbf{x}} & \dfrac{\partial S\_{1}}{\partial p} & \dfrac{\partial S\_{1}}{\partial \nu} \\
\dfrac{\partial S\_{2}}{\partial \textbf{x}} & \dfrac{\partial S\_{2}}{\partial p} & \dfrac{\partial S\_{2}}{\partial \nu} \\
\dfrac{\partial S\_{3}}{\partial \textbf{x}} & \dfrac{\partial S\_{3}}{\partial p} & \dfrac{\partial S\_{3}}{\partial \nu} \\
\end{pmatrix}\right|_{\textbf{x}\_{i},p\_{i},\nu\_{i}}^{-1} \cdot S(\textbf{x}\_{i},p\_{i},\nu\_{i})
$$

donde:

- $\frac{\partial S\_{1}}{\partial \textbf{x}} (\textbf{x}\_{i},p\_{i},\nu\_{i}) = F\_{\textbf{x}}(\textbf{x}\_{i},p\_{i})$
- $\frac{\partial S\_{1}}{\partial p} (\textbf{x}\_{i},p\_{i},\nu\_{i}) = F\_{p}(\textbf{x}\_{i},p\_{i})$
- $\frac{\partial S\_{1}}{\partial \nu} (\textbf{x}\_{i},p\_{i},\nu\_{i}) = 0\_{n\times n}$

- $\frac{\partial S\_{2}}{\partial \textbf{x}} (\textbf{x}\_{i},p\_{i},\nu\_{i})\_{i,j} = \frac{\partial}{\partial x\_{j}} \sum_{k=1}^{n} \frac{\partial F\_{i}}{\partial x\_{k}}\nu\_{k} = \sum_{k=1}^{n} \frac{\partial^{2} F\_{i}}{\partial x\_{j}\partial x\_{k}}\nu\_{k}  $
- $\frac{\partial S\_{2}}{\partial p} (\textbf{x}\_{i},p\_{i},\nu\_{i})\_{i} = \sum_{k=1}^{n} \frac{\partial^{2} F\_{i}}{\partial p\partial x\_{k}}\nu\_{k}$
- $\frac{\partial S\_{2}}{\partial \nu} (\textbf{x}\_{i},p\_{i},\nu\_{i})\_{i,j} = \frac{\partial}{\partial \nu\_{j}} \sum_{k=1}^{n} \frac{\partial F\_{i}}{\partial x\_{k}}\nu\_{k} = \sum_{k=1}^{n} \frac{\partial F\_{i}}{\partial x\_{k}}\delta\_{j,k} = \frac{\partial F\_{i}}{\partial x\_{j}}$

- $\frac{\partial S\_{3}}{\partial \textbf{x}} (\textbf{x}\_{i},p\_{i},\nu\_{i}) = 0\_{1\times n}$
- $\frac{\partial S\_{3}}{\partial p} (\textbf{x}\_{i},p\_{i},\nu\_{i}) = 0$
- $\frac{\partial S\_{3}}{\partial \nu} (\textbf{x}\_{i},p\_{i},\nu\_{i}) = \nu\_{0}^{T}$

# Bifurcaciones de Hopf

Para encontrar puntos límite necesitamos resolver el sistema 

$$
S(\textbf{x},p,\nu) = \begin{cases}
F(\textbf{x},p) = 0\\
(F\_{\textbf{x}}(\textbf{x},p)^{2} + \kappa I\_{n})v = 0\\
v \cdot v -1 = 0 \\
w \cdot v -1 = 0
\end{cases}
$$

donde $v $ es un eigenvalor de $F$ en $x$ y $p$.

Como en el caso anterior, usamos método de Newton

$$
\begin{pmatrix}
\textbf{x}\_{i+1}\\
p\_{i+1}\\
v\_{i+1}\\
\kappa\_{i+1}
\end{pmatrix} = 
\begin{pmatrix}
\textbf{x}\_{i}\\
p\_{i}\\
v\_{i}\\
\kappa\_{i}
\end{pmatrix} - 
\left.\begin{pmatrix}
\frac{\partial S\_{1}}{\partial \textbf{x}} & \frac{\partial S\_{1}}{\partial p} & \frac{\partial S\_{1}}{\partial v} & \frac{\partial S\_{1}}{\partial \kappa}\\
\frac{\partial S\_{2}}{\partial \textbf{x}} & \frac{\partial S\_{2}}{\partial p} & \frac{\partial S\_{2}}{\partial v} & \frac{\partial S\_{2}}{\partial \kappa} \\
\frac{\partial S\_{3}}{\partial \textbf{x}} & \frac{\partial S\_{3}}{\partial p} & \frac{\partial S\_{3}}{\partial v} & \frac{\partial S\_{3}}{\partial \kappa}\\
\frac{\partial S\_{4}}{\partial \textbf{x}} & \frac{\partial S\_{4}}{\partial p} & \frac{\partial S\_{4}}{\partial v} & \frac{\partial S\_{4}}{\partial \kappa} 
\end{pmatrix} \right|_{\textbf{x}\_{i},p\_{i},v\_{i},\kappa\_{i}}^{-1} * S(\textbf{x}\_{i},p\_{i},v\_{i},\kappa\_{i})
$$

donde:

- $\frac{\partial S\_{1}}{\partial \textbf{x}} = F\_{\textbf{x}}$
- $\frac{\partial S\_{1}}{\partial p} = F\_{p}$
- $\frac{\partial S\_{1}}{\partial v} = 0\_{n\times n}$
- $\frac{\partial S\_{1}}{\partial \kappa} = 0$

- $\frac{\partial S\_{2}}{\partial \textbf{x}} = \frac{\partial}{\partial x\_{j}}\sum_{k}\left(\sum_{l}\frac{\partial F\_{i}}{\partial x\_{l}}\frac{\partial F\_{l}}{\partial x\_{k}} + \kappa\delta\_{i,k})\right)v\_{k} = \sum_{k}\sum_{l}\frac{\partial}{\partial x\_{j}}\left(\frac{\partial F\_{i}}{\partial x\_{l}}\frac{\partial F\_{l}}{\partial x\_{k}}\right)v\_{k}$
- $\frac{\partial S\_{3}}{\partial \textbf{x}} = \sum_{k}\sum_{l}\frac{\partial}{\partial p}\left(\frac{\partial F\_{i}}{\partial x\_{l}}\frac{\partial F\_{l}}{\partial x\_{k}}\right)v\_{k}$
- $\frac{\partial S\_{2}}{\partial v} = \frac{\partial}{\partial v\_{j}}\sum_{k}\left(\sum_{l}\frac{\partial F\_{i}}{\partial x\_{l}}\frac{\partial F\_{l}}\{\partial x\_{k}} + \kappa\delta\_{i,k}\right)v\_{k} = \sum_{k}\left(\sum_{l}\frac{\partial F\_{i}}{\partial x\_{l}}\frac{\partial F\_{l}}{\partial x\_{k}} + \kappa\delta\_{i,k}\right)\delta\_{k,j} = \sum_{l}\frac{\partial F\_{i}}{\partial x\_{l}}\frac{\partial F\_{l}}{\partial x\_{j}} + \kappa \delta\_{i,j}$
- $\frac{\partial S\_{2}}{\partial \kappa} = \frac{\partial}{\partial \kappa}\sum_{k}\left(\sum_{l}\frac{\partial F\_{i}}{\partial x\_{l}}\frac{\partial F\_{l}}{\partial x\_{k}} + \kappa\delta\_{i,k}\right)v\_{k} = \frac{\partial}{\partial \kappa}\sum_{k}\kappa\delta\_{i,k}v\_{k} = v\_{i}$

- $\frac{\partial S\_{3}}{\partial \textbf{x}} = 0\_{1\times n}$
- $\frac{\partial S\_{3}}{\partial p} = 0$
- $\frac{\partial S\_{3}}{\partial v} = \frac{\partial}{\partial v\_{j}}\sum_{k}v\_{k}^{2} = \sum_{k}2v\_{k}\delta\_{k,j} = 2v\_{j} = 2v^{T}$
- $\frac{\partial S\_{3}}{\partial \kappa} = 0$

- $\frac{\partial S\_{1}}{\partial \textbf{x}} = 0\_{1\times n}$
- $\frac{\partial S\_{1}}{\partial p} = 0$
- $\frac{\partial S\_{1}}{\partial v} = \frac{\partial}{\partial v\_{j}}\sum_{k}w\_{k}v\_{k} = \sum_{k}w\_{k}\delta\_{k,j} = w\_{j} = w^{T}$
- $\frac{\partial S\_{1}}{\partial \kappa} = 0$
