<h1 align="center">Modelos de EDO para el crecimiento tumoral y su tratamiento</h1>

<p align="center">
  <img src="https://img.shields.io/badge/Python-3776AB?style=for-the-badge&logo=python&logoColor=white" alt="Python">
  <img src="https://img.shields.io/badge/NumPy-013243?style=for-the-badge&logo=numpy&logoColor=white" alt="NumPy">
  <img src="https://img.shields.io/badge/SciPy-8CAAE6?style=for-the-badge&logo=scipy&logoColor=white" alt="SciPy">
  <img src="https://img.shields.io/badge/Matplotlib-11557C?style=for-the-badge&logo=plotly&logoColor=white" alt="Matplotlib">
</p>

<p align="center">
  <img src="figs/fig1_angiogenesis.jpg" alt="Ilustración de la angiogénesis tumoral: vasos nuevos que invaden un tumor" width="360">
</p>

<p align="center"><em>Figura 1. Angiogénesis tumoral: el tumor recluta vasos nuevos (rojo) que le llevan oxígeno y nutrientes, aumentan su capacidad de carga y abren la vía a la diseminación.</em></p>

## Descripción del proyecto

Este proyecto estudia, con ecuaciones diferenciales ordinarias, cómo crece un tumor sólido y cómo responde a dos tratamientos: los fármacos antiangiogénicos y la radioterapia. Se basa en el artículo de Sachs, Hlatky y Hahnfeldt (2001), *Simple ODE Models of Tumor Growth and Anti-Angiogenic or Radiation Treatment*.

La Figura 1 resume la idea central. Un tumor pequeño puede quedarse latente, pero para crecer necesita vasos sanguíneos nuevos que le lleven oxígeno y nutrientes. Esos vasos aumentan su capacidad de carga, es decir, el número máximo de células que el entorno puede sostener. Si se frena la formación de vasos, se frena también el tumor.

El trabajo avanza en tres capas de modelos, de la más simple a la más detallada:

1. **Una EDO para el tumor.** Los modelos logístico generalizado y de Gompertz describen por qué los tumores desaceleran al crecer.
2. **Dos EDO para tumor y vasos.** El sistema (N, K) explica esa desaceleración con un mecanismo: el tumor construye su propia capacidad de carga y a la vez la limita.
3. **Dos EDO para la radiación.** El modelo lineal-cuadrático (LQ) calcula cuántas células sobreviven a una dosis dada.

Todo el desarrollo es determinista y vale para tumores ya grandes, de 1 mm³ o más, donde las fluctuaciones por número pequeño de células se pueden despreciar. El mensaje del artículo es que pocos parámetros bastan para capturar lo esencial; lo difícil es estimarlos.

## Objetivos y alcance

El objetivo es entregar un reporte de máximo 10 páginas y una exposición presencial de 15 minutos que cubran los cinco puntos que pide la orientación del curso:

1. **Derivación del modelo:** hipótesis biológicas y deducción de los modelos de una EDO, de angiogénesis (ecs. 6–7) y LQ (ecs. 8–9).
2. **Solución analítica:** separación de variables, equilibrio del sistema (N, K) e integración del modelo LQ hasta ln(N/N₀) = −αD − βGD².
3. **Análisis cualitativo:** equilibrios, estabilidad, atractor global y bifurcaciones, incluido el caso no autónomo con tratamiento.
4. **Simulación numérica:** método de Euler contra soluciones exactas, reproducción de la Fig. 1 del artículo y comparación del modelo LQ con su fórmula.
5. **Interpretación biológica:** la angiogénesis como explicación mecanicista de la desaceleración y el modelo LQ como herramienta de planificación en radioterapia.

Se omite la Sección 5 del artículo (generalizaciones con McKendrick–von Foerster, PDEs y modelos estocásticos). Como extensiones propias, el equipo añade el tratamiento constante con sus umbrales de erradicación, la adimensionalización del sistema (N, K) y un simulador interactivo.

## Modelos de una EDO

Casi todos los tumores sólidos reducen su tasa de crecimiento específica a medida que crecen, y una EDO de una variable lo captura con la ecuación dN/dt = f(N). El artículo usa el logístico generalizado, con capacidad de carga K y exponente ν, que tiene un único atractor global en N = K:

$$
\frac{dN}{dt} = \frac{\mu N}{\nu}\left[1-\left(\frac{N}{K}\right)^{\nu}\right], \qquad \mu>0
$$

Con ν = 1 queda el logístico. En el límite ν → 0 queda la ecuación de Gompertz, dN/dt = −μ N ln(N/K). Ambos se resuelven por separación de variables, y el cambio u = (N/K)^(−ν) convierte el caso general en una EDO lineal:

$$
N(t) = K\left[1 + \left(\left(\tfrac{K}{N_0}\right)^{\nu} - 1\right)e^{-\mu t}\right]^{-1/\nu}, \qquad N(t) = K\exp\!\left(\ln\tfrac{N_0}{K}\,e^{-\mu t}\right)\ \text{(Gompertz)}
$$

La forma de la curva la fija ν. El punto de inflexión, donde el tumor crece más rápido, está en N = K(1+ν)^(−1/ν):

| Modelo | ν | Inflexión en N |
|---|---|---|
| Logístico | 1 | K/2 |
| Generalizado | 0.25 | ≈ 0.41 K |
| Gompertz | → 0 | K/e ≈ 0.37 K |

El artículo menciona que el mejor ajuste clínico a tumores de mama sugiere ν cercano a 0.25. La Figura 2 compara una familia amplia de modelos: los de arriba crecen sin saturar y los de abajo se estabilizan. En la imagen, Gompertz es el que más rápido se acerca a la capacidad de carga y las variantes de von Bertalanffy son las más lentas. La Figura 3 muestra las soluciones exactas de los tres modelos del cuadro.

<p align="center">
  <img src="figs/fig2_comparacion_modelos.png" alt="Comparación de modelos de crecimiento tumoral sin y con saturación" width="400">
</p>

<p align="center"><em>Figura 2. Modelos de crecimiento tumoral sin saturación (arriba) y con saturación (abajo). Imagen aportada por el equipo.</em></p>

<p align="center">
  <img src="figs/fig3_modelos_exactos.png" alt="Soluciones exactas de los modelos logístico, generalizado y de Gompertz" width="480">
</p>

<p align="center"><em>Figura 3. Soluciones exactas con K = 5000, N₀ = 300 y μ = 0.1 para ν = 1, ν = 0.25 y Gompertz. Los puntos marcan el punto de inflexión de cada curva.</em></p>

## Angiogénesis: el sistema tumor–vasos

El artículo trata la capacidad de carga como una variable K(t), proporcional a la neovascularización, que el tumor construye y a la vez frena. Esto da un sistema de dos EDO para las células tumorales N y la capacidad de carga K:

$$
\frac{dN}{dt} = N\,F\!\left(\frac{N}{K}\right), \qquad \frac{dK}{dt} = -\alpha\, c(t)\,K + \omega N - \gamma N^{2/3} K
$$

F es una función decreciente con F(1) = 0, por ejemplo la de Gompertz, F(x) = −μ ln x. En la ecuación de K, el término ωN es la estimulación de corto alcance que el tumor ejerce sobre las células endoteliales, y γN^(2/3)K es la inhibición de largo alcance. El término −αc(t)K representa un fármaco antiangiogénico con concentración c(t). El exponente 2/3 puede leerse como una relación entre superficie y volumen; esa lectura es del equipo, el artículo no la justifica.

Sin tratamiento (c = 0) el sistema tiene un único atractor global en el cuadrante positivo:

$$
N^{*} = K^{*} = \left(\frac{\omega}{\gamma}\right)^{3/2}, \qquad J = \begin{pmatrix} -\mu & \mu \\ \omega/3 & -\omega \end{pmatrix}, \quad \operatorname{tr}J = -(\mu+\omega) < 0, \quad \det J = \tfrac{2}{3}\mu\omega > 0
$$

La traza negativa y el determinante positivo confirman que el equilibrio es localmente estable. Esto da la explicación mecanicista de la desaceleración: al crecer, el tumor recluta vasos y puede seguir creciendo, pero la inhibición va ganando peso, K se acerca a un límite y el tumor no puede superarlo.

Al escalar N y K con (ω/γ)^(3/2) y el tiempo con 1/μ, el sistema sin tratamiento queda con un solo parámetro, ω/μ, y el equilibrio pasa a ser el punto (1, 1). La Figura 4 muestra el plano de fase en esas variables.

<p align="center">
  <img src="figs/fig4_plano_fase.png" alt="Plano de fase del sistema tumor–vasos con nulclinas y trayectorias" width="380">
</p>

<p align="center"><em>Figura 4. Plano de fase (n, k) sin tratamiento, con n = N/(ω/γ)^(3/2), k = K/(ω/γ)^(3/2), F(x) = −μ ln x y ω/μ = 1.5. Las líneas punteadas son las nulclinas (n = k y k = n^(1/3)) y las curvas continuas, trayectorias que convergen a (1, 1).</em></p>

## Tratamiento antiangiogénico

Un inhibidor como la endostatina frena o revierte la formación de vasos, y esto entra al modelo por el término −αc(t)K. El artículo supone una dosis diaria en la que la concentración de cada bolus administrado en el instante t₀ decae exponencialmente:

$$
c(t) = \sum_{i=10}^{19} c_0\, H(t-i)\, e^{-r\,(t-i)}
$$

Con tratamiento el sistema deja de ser autónomo, porque c depende del tiempo, y ya no hay un equilibrio fijo que estudiar: se analiza la trayectoria. El artículo describe lo que debe verse: el tumor crece y estimula una neovasculatura que crece rápido; al tratar, K baja y arrastra al tumor tras él; al parar el tratamiento el crecimiento se reanuda, más lento al principio porque K quedó por debajo de su valor previo.

El artículo no da los valores numéricos de su Fig. 1, así que los parámetros de la Figura 5 son supuestos del equipo: μ = 0.35 y ω = 0.45 (días⁻¹), γ elegido para que el equilibrio sin tratamiento sea 14 cc, α = 1, c₀ = 0.5, r = 1 día⁻¹, V₀ = 0.4 cc y K₀ = 1 cc. Con ellos, el tumor sin tratar llega a unos 11.3 cc el día 25, y con tratamiento alcanza un máximo cercano a 3.5 cc, baja a unos 1.9 cc y rebrota hasta 2.8 cc.

<p align="center">
  <img src="figs/fig5_tratamiento.png" alt="Simulación del tratamiento antiangiogénico diario" width="520">
</p>

<p align="center"><em>Figura 5. Tratamiento antiangiogénico diario en los días 10–19, con F(x) = −μ ln x. Se muestran el tumor V y la capacidad de carga K con tratamiento, y el tumor sin tratamiento como referencia. Parámetros supuestos, calculados con un integrador de paso pequeño.</em></p>

## Radioterapia: modelo lineal-cuadrático

La radiación ionizante daña sobre todo el ADN, en forma de roturas de doble cadena (DSB). La mayoría se repara en media hora, pero algunas se reparan mal, y una mala reparación puede bastar para matar a la célula en la siguiente mitosis. El artículo modela esto con dos EDO para el número medio de DSB por célula U(t) y el número de células N(t), con tasa de dosis R(t):

$$
\frac{dU}{dt} = \delta R - \frac{U}{\tau} - \gamma U^{2}, \qquad \frac{dN}{dt} = -\left(\alpha R + \kappa U^{2}\right)N
$$

δ es el número de DSB por unidad de dosis, τ el tiempo de reparación, γ la constante de reacción binaria entre DSB, α el número de lesiones letales no reparables por unidad de dosis y 2κ la tasa de mala reparación letal por par de DSB. Para dosis menores de unos 5 Gy se aproxima γ = 0, la primera ecuación se integra de forma explícita y, para tiempos mucho mayores que τ, se obtiene el modelo LQ:

$$
\ln\frac{N}{N_0} = -\alpha D - \beta\, G\, D^{2}, \qquad \beta = \frac{\delta^{2}\kappa\,\tau}{2}, \qquad G = 2\int_0^T f(t)\int_0^t e^{-(t-t')/\tau} f(t')\,dt'\,dt, \quad f(t) = \frac{R(t)}{D}
$$

α mide el daño que no se puede reparar y crece linealmente con la dosis; β mide la mala reparación y crece con el cuadrado. G es el factor de Lea–Catcheside: vale 1 para una dosis aguda y es menor que 1 si la irradiación se prolonga, porque durante ella la reparación actúa. Para una tasa de dosis constante durante un tiempo T, con x = T/τ, resulta G = 2(x − 1 + e⁻ˣ)/x², que tiende a 1 cuando x → 0 y a 2/x cuando x → ∞.

La Figura 6 muestra cómo la razón α/β cambia la forma de la curva: con α/β pequeño la curva se dobla antes. La Figura 7 muestra que cuanto más se prolonga la irradiación, menos pesa el término cuadrático.

<p align="center">
  <img src="figs/fig6_supervivencia_LQ.png" alt="Curvas de supervivencia del modelo lineal-cuadrático" width="440">
</p>

<p align="center"><em>Figura 6. Curvas de supervivencia LQ con α = 0.2 Gy⁻¹ y G = 1. Los puntos marcan D = α/β, donde los términos lineal y cuadrático pesan lo mismo.</em></p>

<p align="center">
  <img src="figs/fig7_funcion_G.png" alt="Factor G en función de la duración relativa de la irradiación" width="420">
</p>

<p align="center"><em>Figura 7. Factor G para una tasa de dosis constante, en función de x = T/τ. Para x pequeño G ≈ 1 (dosis aguda); para x grande decae como 2/x.</em></p>

## Tratamiento constante y bifurcaciones

Esta sección es una extensión del equipo: el artículo plantea el término de tratamiento con c(t) general, y aquí se fija c constante para estudiar cómo cambian los equilibrios con la dosis. Es donde aparecen las bifurcaciones que pide la orientación.

En la EDO de una variable con tratamiento, dN/dt = −αcN + f(N), el equilibrio positivo del logístico generalizado existe solo por debajo de un umbral. Al cruzarlo choca con N = 0 y le cede la estabilidad, una bifurcación transcrítica:

$$
N^{*} = K\left(1-\frac{\nu\alpha c}{\mu}\right)^{1/\nu}, \qquad c < c^{*} = \frac{\mu}{\nu\alpha}
$$

En el límite de Gompertz (ν → 0) el equilibrio es N\* = K·e^(−αc/μ) y el umbral tiende a infinito: con una dosis constante el tumor se reduce, pero este modelo nunca lo erradica.

En el sistema de angiogénesis (6)–(7) con c constante y F(x) = −μ ln x aparece también un umbral, ahora finito. El determinante del Jacobiano cambia de signo exactamente en él:

$$
N^{*} = K^{*} = \left(\frac{\omega-\alpha c}{\gamma}\right)^{3/2}, \qquad c < c^{*} = \frac{\omega}{\alpha}, \qquad \det J = \tfrac{2}{3}\mu(\omega-\alpha c)
$$

Dentro del modelo, una dosis constante por encima de ω/α colapsa el tumor, mientras que por debajo solo reduce su tamaño de equilibrio. Estos son resultados del modelo, no una recomendación clínica. La Figura 8 muestra ambos diagramas de bifurcación.

<p align="center">
  <img src="figs/fig8_bifurcacion.png" alt="Diagramas de bifurcación con tratamiento constante" width="640">
</p>

<p align="center"><em>Figura 8. Diagramas de bifurcación con dosis constante. Izquierda: equilibrio N*/K de la EDO de una variable para ν = 1, ν = 0.25 y Gompertz; los círculos son los umbrales αc/μ = 1/ν. Derecha: equilibrio del sistema (6)–(7) con umbral en αc/ω = 1. Línea continua: estable; línea discontinua: inestable.</em></p>

## Glosario de parámetros y variables

La letra α tiene dos significados según el modelo: en (7) es la fuerza del fármaco y en (9) es el daño letal directo de la radiación. En el reporte hay que aclararlo siempre.

| Símbolo | Significado | Modelo |
|---|---|---|
| N | Células tumorales (o volumen V) | Todos |
| K | Capacidad de carga, proporcional a la neovascularización | Una EDO y angiogénesis |
| μ, ν | Tasa de crecimiento y exponente que gobierna cómo se acerca N a K | Una EDO y angiogénesis |
| F(x) | Función decreciente con F(1) = 0, por ejemplo −μ ln x | Angiogénesis |
| ω | Estimulación de corto alcance del tumor sobre los vasos | Angiogénesis |
| γ | Inhibición de largo alcance (angiogénesis); constante de reacción binaria entre DSB (LQ) | Angiogénesis y LQ |
| α | Fuerza del fármaco (angiogénesis); lesiones letales no reparables por Gy (LQ) | Angiogénesis y LQ |
| c(t) | Concentración del agente en el tumor | Angiogénesis |
| U | Número medio de DSB por célula | LQ |
| δ | DSB producidas por unidad de dosis | LQ |
| τ | Tiempo de reparación | LQ |
| κ | 2κ es la tasa de mala reparación letal por par de DSB | LQ |
| β = δ²κτ/2 | Coeficiente del término cuadrático | LQ |
| R, D | Tasa de dosis y dosis total (Gy) | LQ |
| G | Factor de Lea–Catcheside, vale 1 para dosis aguda | LQ |

## Estructura del código y cómo reproducir las figuras

Las Figuras 3 a 8 de este documento se generaron con un script de Python (NumPy, SciPy y Matplotlib) que acompaña a esta documentación (`figuras_proyecto.py`). Las Figuras 1 y 2 son imágenes aportadas por el equipo y no las genera el script. Para reproducir las demás:

```bash
pip install numpy scipy matplotlib
python figuras_proyecto.py   # escribe las figuras 3 a 8 en figs/
```

## Referencias

- R. K. Sachs, L. R. Hlatky y P. Hahnfeldt. Simple ODE Models of Tumor Growth and Anti-Angiogenic or Radiation Treatment. *Mathematical and Computer Modelling* 33 (2001), 1297–1305.
- P. Hahnfeldt, D. Panigrahy, J. Folkman y L. Hlatky. Tumor development under angiogenic signaling: a dynamical theory of tumor growth, post-vascular dormancy and treatment response. *Cancer Research* (modelo de capacidad de carga dinámica; referencia 17 del artículo).
- S. B. Curtis. Lethal and potentially lethal lesions induced by radiation: a unified repair model. *Radiation Research* 106 (1986), 252–270 (origen de las ecs. 8–9; referencia 23 del artículo).
- H. D. Thames y J. H. Hendry. *Fractionation in Radiotherapy*. Taylor & Francis, Londres, 1987 (el modelo LQ; referencia 19 del artículo).
- Figuras 1 y 2: imágenes aportadas por el equipo. La Figura 1 viene de un banco de imágenes y trae marca de agua: confirmar su licencia o sustituirla antes de la entrega. Figuras 3 a 8: generadas por el equipo con Python.
