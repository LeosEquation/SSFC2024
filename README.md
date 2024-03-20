# Simulaciones numéricas en módelos de ecuaciones diferenciales.

Este proyecto forma parte de un programa de Servicio Social para obtener el grado de Licenciado en Física por la Universidad Nacional Autónoma de México. Estas simulaciones están programadas en [Julia](https://julialang.org/) y se basan en el uso de la paquetería [TaylorSeries.jl](https://github.com/JuliaDiff/TaylorSeries.jl).

#### Autor

- Leonel Mayorga López. Estudiante de Licenciatura en Física, Facultad de Ciencias, Universidad Nacional Autónoma de México (UNAM).

#### Asesor de Servicio Social

- Luis Benet Fernández. Doctor en Física Teórica y Profesor Asociado, Instituto de Ciencias Físicas, Universidad Nacional Autónoma de México (UNAM)

#### Lista de métodos

- Continuación por pseudo-longitud de arco

## Continuación por pseudo-longitud de arco

Este método (también usualmente conocido como *Pseudo Arc Length Continuation* (PACL) en inglés) es utilizado para encontrar las raíces del sistema de un ecuaciones diferenciales dado, partiendo de un parámetro inicial y sus respectivas raíces del sistema. Para encontrar estas raíces se usa el método de Newton sobre la variable y el parámetro del sistema. Si tenemos un sistema de ecuaciones diferenciales expresado de esta manera

$$ \dfrac{d**x**}{dt} = F(**x**,\lambda) $$

donde $F:\mathbb{R}^{n+1} \rightarrow \mathbb{R}^{n}$, $ \begin{pmatrix} **x** \\\ \lambda \end{pmatrix}   F(**x**,\lambda)$
