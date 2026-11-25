# Bloque LQ: función G de Lea-Catcheside, simulación y razón α/β

**Responsable:** Brayan · **Ecuaciones del paper:** (10)–(11) · **Presupuesto:** ≈1.5 páginas (texto + 2 figuras principales; el resto de figuras va al repositorio/anexo).

## 1. Objetivo

Calcular la función $G$ para distintos perfiles de irradiación, verificar sus límites, comprobar numéricamente que la fórmula LQ $\ln(N/N_0)=-\alpha D-\beta G D^2$ reproduce la solución de (8)–(9), y explicar por qué prolongar o fraccionar la dosis aumenta la supervivencia celular.

## 2. Modelo y convención

Sea $U(t)$ el número medio de roturas de doble cadena (DSB) reparables por célula y $N(t)$ el número de células. Usamos la convención (**supuesta**: el extracto del plan no reproduce (8)–(9) y debe contrastarse con el paper):

$$\dot U=\delta\,\dot D(t)-\frac{U}{\tau}-\gamma U^{2}\qquad(8),\qquad\qquad \dot N=-N\big[\alpha\,\dot D(t)+\kappa\,U^{2}\big]\qquad(9)$$

con $U(0)=0$, $N(0)=N_0$ y $\dot D\ge 0$ la tasa de dosis (Gy/h), $D=\int\dot D\,dt$ la dosis total. Significado: $\delta$ = DSB por célula y por Gy; $\tau$ = tiempo característico de reparación; $\gamma$ = tasa de interacción de lesiones (término cuadrático que acelera su eliminación cuando $U$ es grande); $\alpha$ = tasa de muerte lineal por Gy (daño irreparable); $\kappa$ = tasa de muerte por pares de lesiones, $\propto U^2$.

**Aclaración de notación.** En (8)–(9), $\alpha$ [Gy$^{-1}$] y $\gamma$ [DSB$^{-1}$h$^{-1}$] **no** son los $\alpha,\gamma$ de (6)–(7) (potencia del fármaco antiangiogénico y tasa de inhibición vascular). En el código se usa `GAMMA_LQ` para evitar colisiones con el `params.py` del bloque de angiogénesis.

Con $\gamma=0$, (8) es lineal: $U(t)=\delta\int_{-\infty}^{t}\dot D(s)e^{-(t-s)/\tau}ds$ (para dosis aguda, $U=\delta D e^{-t/\tau}$ y $\int U^2dt=\delta^2D^2\tau/2$, como obtuvo Karel). Integrando (9):

$$\ln\frac{N}{N_0}=-\alpha D-\kappa\!\int_0^\infty\!U^2dt=-\alpha D-\beta\,G\,D^2,\qquad \beta=\frac{\delta^2\kappa\tau}{2}.\qquad(10)$$

## 3. La función de Lea-Catcheside $G$

**Definición.** Sustituyendo la expresión de $U$ y usando $\int_{\max(t_1,t_2)}^\infty e^{-2t/\tau}dt=\tfrac{\tau}{2}e^{-2\max(t_1,t_2)/\tau}$, se obtiene $\kappa\int U^2dt=\beta D^2 G$ con

$$G=\frac{1}{D^2}\iint \dot D(t_1)\dot D(t_2)\,e^{-|t_1-t_2|/\tau}dt_1dt_2=\frac{2}{D^2}\int dt_1\,\dot D(t_1)\!\int_{-\infty}^{t_1}\!dt_2\,\dot D(t_2)\,e^{-(t_1-t_2)/\tau}\qquad(11)$$

(la segunda forma integra sobre la región triangular $t_2<t_1$). Como $e^{-|t_1-t_2|/\tau}\le1$ y $\iint\dot D\dot D=D^2$, siempre $0<G\le1$.

**Dosis aguda.** $\dot D=D\,\delta_{\rm Dirac}(t)$ ⇒ $e^{-|t_1-t_2|/\tau}=1$ ⇒ $G=1$.

**Dosis constante en $[0,T]$**, $\dot D=D/T$, $x=T/\tau$:

$$G=\frac{2}{T^2}\int_0^T\!\!\int_0^{t_1}\!e^{-(t_1-t_2)/\tau}dt_2\,dt_1=\frac{2\tau}{T^2}\int_0^T\!\big(1-e^{-t_1/\tau}\big)dt_1=\boxed{\frac{2\,(x-1+e^{-x})}{x^{2}}}$$

*Límites.* Para $x\to0$: $x-1+e^{-x}=\tfrac{x^2}{2}-\tfrac{x^3}{6}+\dots$ ⇒ $G=1-\tfrac{x}{3}+O(x^2)\to1$ (se recupera la dosis aguda). Para $x\to\infty$: $G\sim 2(x-1)/x^2\to 2/x\to0$.

**Perfil exponencial**, $\dot D=(D/\theta)e^{-t/\theta}$, $t\ge0$ (dosis total $D$ entregada con constante de tiempo $\theta$, p. ej. una fuente que decae), $y=\theta/\tau$. La integral interna es $\int_0^{t_1}e^{-t_2/\theta}e^{-(t_1-t_2)/\tau}dt_2=\dfrac{e^{-t_1/\tau}-e^{-t_1/\theta}}{1/\theta-1/\tau}$; multiplicando por $e^{-t_1/\theta}$ e integrando en $t_1$ queda $\dfrac{\tau\theta^2}{2(\tau+\theta)}$, de modo que

$$G=\frac{2}{\theta^{2}}\cdot\frac{\tau\theta^{2}}{2(\tau+\theta)}=\boxed{\frac{1}{1+\theta/\tau}}$$

con $G\to1$ si $\theta\ll\tau$ (equivale a dosis aguda) y $G\to\tau/\theta\to0$ si $\theta\gg\tau$; en $\theta=\tau$, $G=1/2$.

**Extra (fraccionamiento).** Para $n$ fracciones agudas iguales separadas por $\Delta$: $G=\dfrac1n+\dfrac{2}{n^2}\sum_{j=1}^{n-1}(n-j)\,e^{-j\Delta/\tau}$; $G\to1$ si $\Delta\to0$ y $G\to1/n$ si $\Delta\to\infty$.

## 4. Por qué $\alpha$ es lineal, $\beta$ cuadrático y $G$ modula solo el segundo

El término $\alpha D$ cuenta lesiones letales producidas por una sola "pasada" de radiación y por tanto irreparables: es proporcional a la dosis y no depende de cuándo se entrega. El término $\beta G D^2$ proviene de $\kappa U^2$: la muerte requiere la **coincidencia** de dos lesiones reparables, y la probabilidad de coincidencia crece como el producto de ambas, de ahí $D^2$. Cuanto más se reparan las lesiones mientras se sigue irradiando, menos coinciden: $G$ mide la fracción de coincidencias que sobrevive a la reparación. Por eso $G$ multiplica a $\beta$ y no a $\alpha$.

## 5. Metodología numérica

Se integra $(U,L)$ con $L=\ln(N/N_0)$ (evita *underflow* cuando $S\ll1$): $\dot L=-(\alpha\dot D+\kappa U^2)$, con $U(0)=L(0)=0$. La dosis aguda se trata como un **salto exacto** $U(0^+)=\delta D$, $L(0^+)=-\alpha D$; la dosis constante como una tasa $D/T$ sobre $[0,T]$. Tras la irradiación se integra una cola de $15\tau$ (el residuo es $\sim e^{-30}$). Métodos: RK4 (principal) y Euler explícito (contraste), con paso $h=0.02\,\tau$. Se compara $|\Delta\ln S|=|\ln S_{\rm num}-\ln S_{\rm LQ}|$ y se estima $G$ numérica con $\gamma=0$ como $G_{\rm num}=-(L_f+\alpha D)/(\beta D^2)$. Como controles independientes: (i) $G$ por la doble integral (11) discretizada, (ii) `scipy.integrate.solve_ivp` (DOP853, `rtol=1e-12`) frente a nuestro RK4, (iii) orden de convergencia. Para $\gamma>0$ y dosis aguda hay solución cerrada (Bernoulli): $U(t)=U_0e^{-t/\tau}/[1+\gamma\tau U_0(1-e^{-t/\tau})]$, con $U_0=\delta D$; el parámetro adimensional relevante es $\gamma\tau\delta D$ (cociente entre el término cuadrático y el de reparación lineal).

**Tabla 1. Parámetros — TODOS SUPUESTOS** (el paper no da valores numéricos para ellos; órdenes de magnitud típicos en radiobiología celular, a justificar en el reporte).

| Parámetro | Valor | Unidades | Justificación |
|---|---|---|---|
| $\tau$ | 1 | h | reparación de DSB con vida media ≲1 h |
| $\delta$ | 35 | DSB/(célula·Gy) | del orden de decenas de DSB por Gy |
| $\alpha$ | 0.15 | Gy$^{-1}$ | valor típico en tejidos |
| $\beta$ | 0.05 | Gy$^{-2}$ | da $\alpha/\beta=3$ Gy (tejido de respuesta tardía) |
| $\kappa$ | $8.16\times10^{-5}$ | DSB$^{-2}$h$^{-1}$ | derivado: $\kappa=2\beta/(\delta^2\tau)$ |
| $\gamma$ | $2\times10^{-3}$ | DSB$^{-1}$h$^{-1}$ | elegido para que $\gamma\tau\delta D\sim0.1$–$0.4$ a 2–5 Gy |
| $\alpha/\beta$ (barrido) | 3, 10, 30 | Gy | $\alpha$ fijo, $\beta=\alpha/(\alpha/\beta)$ |
| $T/\tau$ (barrido) | 0, 1, 4, 16 | – | 0 = aguda |

## 6. Resultados

**Fig. 1** (`fig_G_vs_Ttau.pdf`). *Título:* "Función de Lea-Catcheside según el tiempo de irradiación". Ejes: $T/\tau$ (o $\theta/\tau$) en escala log vs. $G$ adimensional. Curvas: dosis constante, perfil exponencial, asíntota $2/x$, línea $G=1$ y círculos con la estimación numérica desde (8)–(9). *Pie:* "$G$ decrece monótonamente de 1 (dosis aguda) a 0 al prolongarse la irradiación; los marcadores numéricos ($\gamma=0$) coinciden con las fórmulas cerradas (diferencias $\lesssim10^{-6}$ salvo discretización)."

**Fig. 2** (`fig_LQ_vs_numerico.pdf`). *Título:* "Solución numérica de (8)–(9) frente a la fórmula LQ". Fila superior: $S=N/N_0$ vs. $D$ (Gy), semilog, para dosis aguda y constante con $T=4\tau$. Fila inferior: $|\Delta\ln S|$ de RK4, Euler y del caso $\gamma>0$. *Pie:* "Con $\gamma=0$ RK4 reproduce LQ a $10^{-8}$–$10^{-9}$; Euler ($h=0.02\tau$) tiene error $\sim10^{-2}$. Con $\gamma>0$ la supervivencia numérica queda por encima de LQ y la brecha crece con $D$."

**Fig. 3** (`fig_supervivencia_alpha_beta.pdf`). *Título:* "Supervivencia para distintas razones $\alpha/\beta$ y tiempos de irradiación". Panel izquierdo: dosis aguda, $\alpha/\beta=3,10,30$ Gy ($\alpha$ fijo). Paneles central y derecho: $\alpha/\beta=3$ y $10$ Gy con $T/\tau=0,1,4,16$. Círculos: numérico. *Pie:* "Cuanto menor es $\alpha/\beta$, más curvada es la curva y más sensible a $G$; con $\alpha/\beta$ alto domina el término lineal y el tiempo de irradiación importa poco."

Figuras adicionales para el repositorio/anexo: `fig_trayectorias_U_S.pdf` ($U(t)$ y $S(t)$ a 6 Gy: el pico de $U$ baja al prolongar la dosis), `fig_fraccionamiento.pdf` ($G$ vs. $\Delta/\tau$ para $n=2,3,5,10$), `fig_mapa_gamma.pdf` (desviación respecto a LQ en el plano $(D,\ \gamma\tau\delta)$) y `fig_convergencia_euler_rk4.pdf`.

**Tabla 2. Cifras obtenidas con los parámetros supuestos de la Tabla 1** (salida de `figures.py` y `check_scipy.py`).

| Caso | $\max|\Delta\ln S|$ RK4, $\gamma=0$ | idem Euler, $\gamma=0$ | $\gamma>0$: desviación de $\ln S$ a 2 / 5 / 10 Gy |
|---|---|---|---|
| Dosis aguda | $1.7\times10^{-8}$ | $5.1\times10^{-2}$ | +3.4 % / +11.6 % / +23.7 % |
| Constante, $T=4\tau$ | $2.1\times10^{-9}$ | $5.8\times10^{-3}$ | +0.9 % / +4.1 % / +10.6 % |

(Desviación $=(\ln S_{\rm num}-\ln S_{\rm LQ})/|\ln S_{\rm LQ}|$; positiva = más supervivencia que LQ. Orden observado: RK4 ≈ 4.06, Euler ≈ 1.07. Máx. $|\Delta|$ RK4–SciPy $\le1.2\times10^{-8}$.)

**Tabla 3. Efecto de $G$ con $\alpha/\beta=3$ Gy** (fórmula LQ).

| Irradiación | $G$ | $S(2\text{ Gy})$ | $S(10\text{ Gy})$ |
|---|---|---|---|
| Aguda | 1 | 0.607 | $1.50\times10^{-3}$ |
| Constante, $T=4\tau$ | 0.377 | 0.687 | $3.38\times10^{-2}$ |
| Constante, $T=16\tau$ | 0.117 | 0.724 | 0.124 |

A 10 Gy, prolongar la irradiación de aguda a $16\tau$ multiplica la supervivencia por ≈80; a 2 Gy el efecto es modesto (≈19 %), porque el término lineal domina.

## 7. Interpretación biológica

Con dosis aguda, las lesiones reparables aparecen todas a la vez, $U$ alcanza su máximo $\delta D$ y la tasa de muerte por coincidencias $\kappa U^2$ es máxima. Si la misma dosis se entrega en un tiempo $T\gtrsim\tau$, la reparación opera durante la irradiación: $U$ nunca llega a acumularse (en el régimen constante, $U\to\delta(D/T)\tau$, mucho menor que $\delta D$ cuando $T\gg\tau$) y el término cuadrático se reduce en la fracción $1-G$. El término $\alpha D$ no cambia, de modo que el efecto es mayor cuanto menor es $\alpha/\beta$. En radioterapia esto sustenta el "efecto de tasa de dosis" y la ventaja de fraccionar: los tejidos de respuesta tardía (bajo $\alpha/\beta$) se benefician más de fracciones pequeñas que los tumores de $\alpha/\beta$ alto.

## 8. ¿Por qué fraccionar o prolongar la dosis cambia la supervivencia celular?

Porque la muerte celular tiene dos componentes de naturaleza distinta. La parte lineal ($\alpha D$) depende solo de la dosis total y no se ve afectada por el calendario. La parte cuadrática ($\beta G D^2$) exige que dos lesiones reparables coexistan, y su coexistencia depende del tiempo: cada lesión desaparece con constante $\tau$. Si la dosis se reparte (en el tiempo o en fracciones separadas por $\Delta\gtrsim\tau$), las lesiones de una parte se reparan antes de que lleguen las de la siguiente, y los pares letales se destruyen antes de formarse. Matemáticamente, $G=\frac{1}{D^2}\iint\dot D\dot D\,e^{-|t_1-t_2|/\tau}$ pondera cada par de instantes de irradiación por la probabilidad de que sus lesiones aún coexistan; separar los instantes aumenta $|t_1-t_2|$ y reduce $G$. Casos extremos: $n$ fracciones muy espaciadas dan $G\to1/n$ (el término cuadrático es el de una sola fracción de dosis $D/n$, sumado $n$ veces), y la irradiación infinitamente prolongada da $G\to0$ (queda solo la parte lineal). Con la misma $D$ y el mismo $\alpha$, más fraccionamiento o más tiempo implica más supervivencia (Tabla 3); en la práctica, para tumores, esto se compensa con más dosis total.

## 9. Limitaciones

(i) La convención de (8)–(9) con $\gamma>0$ es una suposición a verificar contra el paper; (ii) todos los parámetros son supuestos y no se ajustaron a datos; (iii) el modelo es de campo medio: $U$ es un promedio por célula y $N$ evoluciona de forma determinista, sin variabilidad entre células ni fluctuaciones; (iv) no hay repoblación, redistribución en el ciclo celular, reoxigenación ni efectos de hipoxia, que importan en tratamientos de varias semanas, y $\tau$ y $\delta$ se toman constantes; (v) con $\gamma>0$ la fórmula LQ deja de ser exacta y sobrestima la muerte a dosis altas (aquí, +12 % en $\ln S$ a 5 Gy en dosis aguda), por lo que el umbral "$\gamma=0$ es razonable bajo ≈5 Gy" depende del valor real de $\gamma$ (criterio: $\gamma\tau\delta D\ll1$); (vi) en el código, la dosis aguda es un salto instantáneo y las "dosis constante/exponencial" se truncan (30$\theta$) o discretizan, con error numérico ≤$10^{-8}$ en RK4.

## 10. Conclusiones

(1) $G$ vale 1 para dosis aguda, $2(x-1+e^{-x})/x^2$ para dosis constante y $1/(1+\theta/\tau)$ para perfil exponencial, con los límites esperados. (2) La integración numérica de (8)–(9) con $\gamma=0$ reproduce LQ con el $G$ correspondiente a precisión de RK4. (3) Prolongar o fraccionar la dosis reduce el término cuadrático y aumenta la supervivencia, sobre todo con $\alpha/\beta$ bajo y dosis altas. (4) Con $\gamma>0$ LQ pierde exactitud a dosis altas.

## 11. Referencias

- Sachs, R. K., Hlatky, L. R. y Hahnfeldt, P. (2001). Simple ODE models of tumor growth and anti-angiogenic or radiation treatment. *Mathematical and Computer Modelling*, 33(12–13), 1297–1305.
- Lea, D. E. y Catcheside, D. G. (1942). The mechanism of the induction by radiation of chromosome aberrations in *Tradescantia*. *Journal of Genetics*, 44, 216–245.
- Fowler, J. F. (1989). The linear-quadratic formula and progress in fractionated radiotherapy. *British Journal of Radiology*, 62, 679–694.
- Dale, R. G. (1985). The application of the linear-quadratic dose-effect equation to fractionated and protracted radiotherapy. *British Journal of Radiology*, 58, 515–528.
- Hall, E. J. y Giaccia, A. J. *Radiobiology for the Radiologist* (Wolters Kluwer; usar la edición disponible).


---

# Anexo B. Lista de verificación para Brayan

- [ ] Confirmar en el paper la forma exacta de (8)–(9) y de (11); si difiere de la convención usada, ajustar `f` en `lq.py` y la Sección 2.
- [ ] Declarar en el reporte que $\tau,\delta,\alpha,\beta,\kappa,\gamma$ son **supuestos** (Tabla 1) y justificarlos en una línea.
- [ ] Notación: aclarar que $\alpha,\gamma$ de (8)–(9) no son los de (6)–(7); usar `GAMMA_LQ` en el código.
- [ ] Revisar que $G$: aguda = 1; constante = $2(x-1+e^{-x})/x^2$ con límites 1 y $2/x$; exponencial = $1/(1+\theta/\tau)$.
- [ ] Verificar numéricamente $\ln(N/N_0)=-\alpha D-\beta G D^2$ (Fig. 2) y $G$ (Fig. 1) con $\gamma=0$.
- [ ] Incluir $\gamma>0$ y dejar clara la condición $\gamma\tau\delta D\ll1$; ajustar la afirmación de "≈5 Gy" al $\gamma$ elegido.
- [ ] Figuras con título, ejes con unidades (Gy, h, $T/\tau$, $S$), leyenda legible y pie; en el reporte, solo las Figs. 1 y 3 (resto al repositorio).
- [ ] Texto ≤1.5 páginas y estilo/colores comunes (Adán); revisión cruzada con Roger.
- [ ] Tener lista la respuesta oral a "¿por qué fraccionar o prolongar la dosis cambia la supervivencia?" (Sección 8): lineal fijo, cuadrático depende del tiempo, $G\to1/n$ y $G\to0$.
- [ ] (Opcional) Animación Manim de la región triangular $t_2<t_1$ de (11) y de $e^{-(t_1-t_2)/\tau}$.
- [ ] Verificar referencias y que el código corra desde cero con el README.