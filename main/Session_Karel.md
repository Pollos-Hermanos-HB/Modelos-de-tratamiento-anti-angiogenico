# Modelo LQ (Linear-Quadratic) en Radiobiología

Significado de Variables y Parámetros:

*   **\(U(t)\)**: Número promedio de rupturas de doble cadena (DSBs) en el ADN por célula en el tiempo \(t\).
*   **\(N(t)\)**: Número promedio de células viables (clonogénicas) en la población en el tiempo \(t\).
*   **\(\delta\)**: Número promedio de DSBs inducidos directamente por unidad de dosis (Gy).
*   **\(\tau\)**: Constante de tiempo de reparación del ADN. Indica la rapidez con la que la célula repara las DSBs.
*   **\(\gamma\)**: Constante de tasa de interacción binaria entre dos DSBs (cinética de acción de masas).
*   **\(\alpha\)**: Número promedio de lesiones letales irreparables producidas directamente por unidad de dosis (daño de "un solo impacto").
*   **\(\kappa\)**: Tasa de reparación errónea letal de DSBs por par de rupturas. Determina la probabilidad de que una interacción entre DSBs cause la muerte celular.

---

La ecuación es:
\[
 \frac{dU}{dt} = \delta R - \frac{U}{\tau} - \gamma U^2 
\]

Asumiendo **\(\gamma = 0\)**, la ecuación se reduce a una EDO lineal de primer orden:
\[
 \frac{dU}{dt} = \delta R - \frac{U}{\tau} 
\]

Para una **dosis aguda** (administrada en un instante muy corto, \(t=0\)), la tasa de dosis \(R(t)\) se modela como una función Delta de Dirac: \(R(t) = D \cdot \delta(t)\). Integrando en un intervalo infinitesimal alrededor de \(t=0\), obtenemos la condición inicial justo después de la irradiación:
\[
 U(0) = \delta D 
\]

Para \(t > 0\), la radiación ha cesado (\(R=0\)), por lo que la ecuación queda:
\[
 \frac{dU}{dt} = -\frac{U}{\tau} 
\]

Integrando con la condición inicial \(U(0) = \delta D\):
\[
 U(t) = \delta D \cdot e^{-t/\tau} 
\]

Calculamos la integral de \(U^2\) desde \(t=0\) hasta \(\infty\):
\[
 \int_{0}^{\infty} U^2 dt = \int_{0}^{\infty} (\delta D \cdot e^{-t/\tau})^2 dt = \delta^2 D^2 \int_{0}^{\infty} e^{-2t/\tau} dt 
\]
\[
 \int_{0}^{\infty} U^2 dt = \delta^2 D^2 \left[ -\frac{\tau}{2} e^{-2t/\tau} \right]_{0}^{\infty} = \frac{\delta^2 D^2 \tau}{2} 
\]

---

Obtención del Modelo LQ:

A partir de:
\[
 \frac{dN}{dt} = -(\alpha R + \kappa U^2) N 
\]

Separamos variables e integramos desde \(t=0\) hasta \(t=\infty\):
\[
 \int_{N_0}^{N} \frac{dN}{N} = -\alpha \int_{0}^{\infty} R dt - \kappa \int_{0}^{\infty} U^2 dt 
\]

Sabiendo que \(\int R dt = D\) (dosis total) y sustituyendo la integral de \(U^2\):
\[
 \ln\left(\frac{N}{N_0}\right) = -\alpha D - \kappa \left( \frac{\delta^2 D^2 \tau}{2} \right) 
\]

Agrupando constantes, definimos **\(\beta \equiv \frac{\delta^2 \kappa \tau}{2}\)**. Obtenemos la fórmula LQ:
\[
 \ln\left(\frac{N}{N_0}\right) = -\alpha D - \beta D^2 
\]

QED

---

Interpretaciones biológicas:

*   **El término \(\alpha D\) (Lineal):** Representa la muerte celular por **daño irreparable directo**. Un solo fotón o partícula impacta y rompe el ADN de forma tan compleja que la célula no puede arreglarlo. Es proporcional a la dosis.
*   **El término \(\beta D^2\) (Cuadrático):** Representa la muerte por **mala reparación**. La radiación crea muchas rupturas (\(D\)). La probabilidad de que *dos* rupturas interactúen entre sí para formar una aberración letal es proporcional al cuadrado de la dosis (\(D^2\)). Depende del tiempo (\(\tau\)): si la célula repara rápido (menor \(\tau\)), hay menos tiempo para que dos rupturas interactúen, reduciendo \(\beta\).
