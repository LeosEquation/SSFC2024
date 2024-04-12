# Simulaciones numéricas en módelos de ecuaciones diferenciales.

Este proyecto forma parte de un programa de Servicio Social para obtener el grado de Licenciado en Física por la Universidad Nacional Autónoma de México. Estas simulaciones están programadas en [Julia](https://julialang.org/) y se basan en el uso de la paquetería [TaylorSeries.jl](https://github.com/JuliaDiff/TaylorSeries.jl).

#### Autor

- Leonel Mayorga López. Estudiante de Licenciatura en Física, Facultad de Ciencias, Universidad Nacional Autónoma de México (UNAM).

#### Asesor de Servicio Social

- Luis Benet Fernández. Doctor en Física Teórica y Profesor Asociado, Instituto de Ciencias Físicas, Universidad Nacional Autónoma de México (UNAM)

#### Para empezar...

En este repositorio se pueden consultar el análisis de puntos de equilibrio y bifurcaciones de 3 ejemplos en la carpeta `examples`: uno unidimensional (`example1D`), un análisis del modelo Depredador-Presa (`Predator-Prey model`) y un análisis a la reacción $A\rightarrow B\rightarrow C$ (`ABC-Reaction`). Se busca comparar lo obtenido con lo escrito por [Eusebius J. Doedel](https://users.encs.concordia.ca/~doedel/) en su presentación [An Introduction to Numerical Continuation Methods with Applications](https://users.encs.concordia.ca/~doedel/courses/comp-6361/slides.pdf).

El código para lograr esto consta de 6 archivos en la carpeta `src`.

- `diff_tools.jl`: Usando la paquetería [TaylorSeries.jl](https://github.com/JuliaDiff/TaylorSeries.jl), se desarrollan herramientas de diferenciación automática que se adapten a nuestras funciones, variables y parámetros de entrada. Se pueden hacer cálculos de gradientes y jacobianos.

- `newton.jl`: Se desarrollan variós métodos del Método de Newton que se adaptan a los valores de entrada que tenemos para hallar puntos de equilibrio de forma iterativa.

- `implicit_function.jl`: Mediante el uso del Teorema de la Función Implícita, se desarrollan varios métodos para hallar el valor de la derivada de la función implícita haciendo uso de nuestras herramientas de diferenciación automática.

- `psarc.jl`: Aquí se concentran las funciones y métodos para llevar a cabo el método de Continuación por Pseudo Longitud de Arco (PALC por sus siglas en inglés).

- `equilibrium.jl`: Haciendo uso de los métodos PALC y Newton, se desarrollan las funciónes que devolveran las ramas de equilibrio.

- `bifurcation.jl`: Aquí se encuentran las funciones para hallar los puntos de bifurcación en las ramas de equilibrio.

- `stability.jl`: Aquí se encuentran las funciones que determinaran los intervalos de estabilidad e inestabilidad de las ramas de equilibrio con un enfoque a la graficación simultánea.

