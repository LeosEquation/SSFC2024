# Simulaciones numéricas en módelos de ecuaciones diferenciales.

Este proyecto forma parte de un programa de Servicio Social para obtener el grado de Licenciado en Física por la Universidad Nacional Autónoma de México. Estas simulaciones están programadas en [Julia](https://julialang.org/) y se basan en el uso de la paquetería [TaylorSeries.jl](https://github.com/JuliaDiff/TaylorSeries.jl).

#### Autor

- Leonel Mayorga López. Estudiante de Licenciatura en Física, Facultad de Ciencias, Universidad Nacional Autónoma de México (UNAM).

#### Asesor de Servicio Social

- Luis Benet Fernández. Doctor en Física Teórica y Profesor Asociado, Instituto de Ciencias Físicas, Universidad Nacional Autónoma de México (UNAM)

#### Lista de métodos

- Continuación por pseudo-longitud de arco

## Continuación por pseudo-longitud de arco

Este método (también usualmente conocido como *Pseudo Arc Length Continuation* (PACL) en inglés) es utilizado para encontrar las raíces del sistema de un ecuaciones diferenciales dado, partiendo de un parámetro inicial y sus respectivas raíces del sistema. Calcula numéricamente los pasos $\Delta \vec{x}$ y $\lambda$ a dar para encontrar la siguiente raíz usando el método de Newton a partir del paso anterior. La ventaja de este método es que cada paso cambia manteniendo que $||\Delta \vec{x}||^{2} + \Delta \lambda ^{2} \approxx \Delta s$ para todos los pasos dados y siempre el avance resulta en un punto cercano a una raíz, así se asegura que el método de Newton converga en cada paso.

Tenemos un sistema de ecuaciones diferenciales expresado de esta manera:

$$ \dfrac{d\vec{x}}{dt} = F(\vec{x},\lambda) $$

$$F:\mathbb{R}^{n+1} \rightarrow \mathbb{R}^{n}$$

$$ (\vec{x},\lambda) \xmapsto[]{}   F(\vec{x},\lambda) $$

Según el teorema de la función implícita (*IFT* por sus siglas en inglés) si $F$ es continuamente diferenciable y existe un conjunto $\Omega \subset \mathbb{R}^{n+1}$ tal que $F(\vec{x},\lambda) = 0$ para todo $(\vec{x},\lambda) \in \Omega$, entonces existe una función $G:D\subset\mathbb{R} \rightarrow \mathbb{R}^{n}$ tal que $(G(\lambda),\lambda)\in\Omega$ para todo $\lambda \in D$ y además se asegura que

$$\dfrac{dG}{d\lambda} = -J_{x}^{-1}J_{\lambda} = -\begin{bmatrix} \dfrac{\partial F_{1}}{\partial x_{1}} & \cdots & \dfrac{\partial F_{1}}{\partial x_{n}} \\\ \vdots & \ddots & \vdots \\\ \dfrac{\partial F_{n}}{\partial x_{1}} & \cdots & \dfrac{\partial F_{n}}{\partial x_{n}} \end{bmatrix}^{-1} \begin{bmatrix} \dfrac{\partial F_{1}}{\partial \lambda} \\\ \vdots \\\ \dfrac{\partial F_{n}}{\partial \lambda} \end{bmatrix}$$


