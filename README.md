# Simulaciones numéricas en módelos de ecuaciones diferenciales.

Este proyecto forma parte de un programa de Servicio Social para obtener el grado de Licenciado en Física por la Universidad Nacional Autónoma de México. Estas simulaciones están programadas en [Julia](https://julialang.org/) y se basan en el uso de la paquetería [TaylorSeries.jl](https://github.com/JuliaDiff/TaylorSeries.jl).

#### Autor

- Leonel Mayorga López. Estudiante de Licenciatura en Física, Facultad de Ciencias, Universidad Nacional Autónoma de México (UNAM).

#### Asesor de Servicio Social

- Luis Benet Fernández. Doctor en Física Teórica y Profesor Asociado, Instituto de Ciencias Físicas, Universidad Nacional Autónoma de México (UNAM)

#### Lista de métodos

- Continuación por pseudo-longitud de arco

## Continuación por pseudo-longitud de arco

Este método (también usualmente conocido como *Pseudo Arc Length Continuation* (PACL) en inglés) es utilizado para encontrar las raíces del sistema de un ecuaciones diferenciales dado, partiendo de un parámetro inicial y sus respectivas raíces del sistema. Para encontrar estas raíces se usa el método de Newton sobre la variable y el parámetro del sistema.

Tenemos un sistema de ecuaciones diferenciales expresado de esta manera:

$$ \dfrac{d**x**}{dt} = F(\vec{x},\lambda) $$

$$F:\mathbb{R}^{n+1} \rightarrow \mathbb{R}^{n}$$

$$ (\vec{x},\lambda) \xmapsto[]{}   F(\vec{x},\lambda) $$

Si $F$ es continuamente diferenciable y existe un conjunto $\Omega \subset \mathbb{R}^{n+1}$ tal que $F(\vec{x},\lambda) = 0$ para todo $(\vec{x},\lambda) \in \Omega$, entonces existe una función $G:D\subset\mathbb{R} \rightarrow \mathbb{R}^{n}$ tal que $(G(\lambda),\lambda)\in\Omega$ para todo $\lambda \in D$
