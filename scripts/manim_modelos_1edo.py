"""
Animación de Manim (Community Edition) para el bloque "Modelos de una EDO".

Cómo correrla (desde la carpeta del repo):
    pip install manim
    manim -pql scripts/manim_modelos_1edo.py ModelosUnaEDO      # calidad baja, rápida
    manim -pqh scripts/manim_modelos_1edo.py ModelosUnaEDO      # calidad alta, para la entrega
Requiere LaTeX instalado (MathTex lo usa).

Parámetros: los mismos supuestos del informe (no vienen del paper).
"""
from manim import *
import numpy as np

MU, K, N0 = 0.5, 1.0, 0.01     # supuestos del equipo
T_MAX = 20.0


# ---------- Matemática (misma que el informe) ----------
def sol(t, nu):
    """Solución exacta del logístico generalizado (ec. 4); Gompertz si nu ~ 0."""
    if nu < 1e-6:
        return K * np.exp(np.log(N0 / K) * np.exp(-MU * t))
    D = (K / N0) ** nu - 1.0
    return K / (1.0 + D * np.exp(-MU * t)) ** (1.0 / nu)


def n_inf(nu):
    """Valor de N en la inflexión: K (1+nu)^(-1/nu)."""
    return K * (1.0 + nu) ** (-1.0 / nu)


def t_inf(nu):
    """Instante en que N(t) = N_inf (de u = 1/(1+nu))."""
    D = (K / N0) ** nu - 1.0
    return np.log(D / nu) / MU


def fit(group, max_w=12.5, max_h=5.4):
    """Escala un grupo para que quepa en pantalla."""
    if group.width > max_w:
        group.scale_to_fit_width(max_w)
    if group.height > max_h:
        group.scale_to_fit_height(max_h)
    return group


class ModelosUnaEDO(Scene):
    def clear_all(self):
        if self.mobjects:
            self.play(*[FadeOut(m) for m in self.mobjects])

    def construct(self):
        self.seccion_gancho()
        self.seccion_modelos()
        self.seccion_gompertz()
        self.seccion_generalizado()
        self.seccion_inflexion()
        self.seccion_cierre()

    # ------------------------------------------------------------------
    # 1. Gancho: exponencial vs. tumor real
    # ------------------------------------------------------------------
    def seccion_gancho(self):
        self.next_section("gancho")
        titulo = Text("¿Por qué los tumores se frenan solos?", font_size=40)
        titulo.to_edge(UP)
        self.play(Write(titulo))

        ejes = Axes(
            x_range=[0, T_MAX, 5], y_range=[0, 1.2, 0.5],
            x_length=9, y_length=4.6, tips=False,
            axis_config={"include_numbers": False},
        ).shift(DOWN * 0.7)
        et = Text("tiempo", font_size=24).next_to(ejes.x_axis, DOWN, buff=0.2)
        en = MathTex("N").next_to(ejes.y_axis, LEFT, buff=0.2)
        linea_K = DashedLine(ejes.c2p(0, K), ejes.c2p(T_MAX, K), color=GREY_B)
        tK = MathTex("K").next_to(linea_K, RIGHT, buff=0.15)
        self.play(Create(ejes), FadeIn(et), FadeIn(en))
        self.play(Create(linea_K), FadeIn(tK))

        t_exp = np.log(1.15 / N0) / MU           # hasta que pasa de K
        curva_exp = ejes.plot(lambda t: N0 * np.exp(MU * t),
                              x_range=[0, t_exp], color=ORANGE)
        l_exp = Text("Exponencial: sin límite", font_size=26, color=ORANGE)
        l_exp.next_to(ejes.c2p(t_exp, 1.12), RIGHT, buff=0.2)
        self.play(Create(curva_exp), run_time=2.5)
        self.play(FadeIn(l_exp))
        self.wait(0.5)

        curva_real = ejes.plot(lambda t: sol(t, 1.0),
                               x_range=[0, T_MAX], color=TEAL)
        l_real = Text("Tumor real: se frena cerca de K", font_size=26, color=TEAL)
        l_real.next_to(ejes.c2p(10, 0.5), RIGHT, buff=0.6).shift(DOWN * 0.3)
        self.play(Create(curva_real), run_time=3)
        self.play(FadeIn(l_real))
        self.wait(1.5)
        self.clear_all()

    # ------------------------------------------------------------------
    # 2. Una ecuación, varios modelos (y el límite a Gompertz)
    # ------------------------------------------------------------------
    def seccion_modelos(self):
        self.next_section("modelos")
        titulo = Text("Una ecuación, varios modelos", font_size=40).to_edge(UP)
        gen = MathTex(
            r"\frac{dN}{dt}=\frac{\mu N}{\nu}"
            r"\left[1-\left(\frac{N}{K}\right)^{\nu}\right]"
        ).scale(1.2)
        self.play(Write(titulo))
        self.play(Write(gen))
        self.wait(1)

        caso1 = MathTex(r"\nu=1:", r"\ \text{logístico}", font_size=40)
        caso2 = MathTex(r"\nu\to0:", r"\ \text{Gompertz}", font_size=40)
        caso3 = MathTex(r"\nu\approx0.25:", r"\ \text{mama}", font_size=40)
        casos = VGroup(caso1, caso2, caso3).arrange(RIGHT, buff=1.0)
        casos.next_to(gen, DOWN, buff=0.9)
        for c in casos:
            self.play(FadeIn(c, shift=UP * 0.3))
        self.wait(1)

        # Límite nu -> 0
        self.play(FadeOut(casos), gen.animate.to_edge(UP).shift(DOWN * 1.2).scale(0.8))
        paso1 = MathTex(r"x^{\nu}=e^{\nu\ln x}\approx 1+\nu\ln x", font_size=44)
        paso2 = MathTex(r"\frac{1}{\nu}\left(1-x^{\nu}\right)\;\longrightarrow\;-\ln x",
                        font_size=44)
        gomp = MathTex(r"\frac{dN}{dt}=-\mu N\ln\frac{N}{K}", font_size=52)
        pasos = VGroup(paso1, paso2).arrange(DOWN, buff=0.5).next_to(gen, DOWN, buff=0.7)
        self.play(Write(paso1))
        self.play(Write(paso2))
        self.wait(0.5)
        gomp.next_to(pasos, DOWN, buff=0.7)
        caja = SurroundingRectangle(gomp, color=YELLOW, buff=0.25)
        self.play(Write(gomp), Create(caja))
        self.wait(2)
        self.clear_all()

    # ------------------------------------------------------------------
    # 3. Gompertz: separación de variables
    # ------------------------------------------------------------------
    def seccion_gompertz(self):
        self.next_section("gompertz")
        titulo = Text("Gompertz: separación de variables", font_size=36).to_edge(UP)
        self.play(Write(titulo))

        pasos = VGroup(
            MathTex(r"\frac{dN}{N\ln(N/K)}=-\mu\,dt"),
            MathTex(r"y=\ln\frac{N}{K},\ dy=\frac{dN}{N}\ \Rightarrow\ "
                    r"\int\frac{dy}{y}=-\mu\int dt"),
            MathTex(r"\ln|y|=-\mu t+c\ \Rightarrow\ y=A\,e^{-\mu t}"),
            MathTex(r"N(0)=N_0\ \Rightarrow\ A=\ln\frac{N_0}{K}"),
        ).arrange(DOWN, buff=0.45, aligned_edge=LEFT)
        fit(pasos, max_w=11.5, max_h=4.2)
        pasos.next_to(titulo, DOWN, buff=0.5)

        for p in pasos:
            self.play(Write(p))
            self.wait(0.6)

        resultado = MathTex(
            r"N(t)=K\exp\!\left[\ln\!\left(\frac{N_0}{K}\right)e^{-\mu t}\right]"
        ).scale(1.1)
        resultado.to_edge(DOWN, buff=0.6)
        caja = SurroundingRectangle(resultado, color=YELLOW, buff=0.25)
        self.play(Write(resultado), Create(caja))
        self.wait(2.5)
        self.clear_all()

    # ------------------------------------------------------------------
    # 4. Logístico generalizado: cambio de variable
    # ------------------------------------------------------------------
    def seccion_generalizado(self):
        self.next_section("generalizado")
        titulo = Text("Logístico generalizado: un cambio de variable",
                      font_size=36).to_edge(UP)
        self.play(Write(titulo))

        pasos = VGroup(
            MathTex(r"u=\left(\frac{N}{K}\right)^{\nu}"),
            MathTex(r"\frac{du}{dt}=\nu u\cdot\frac{1}{N}\frac{dN}{dt}"
                    r"=\mu\,u(1-u)"),
            MathTex(r"\frac{1}{u(1-u)}=\frac{1}{u}+\frac{1}{1-u}"),
            MathTex(r"u=\frac{1}{1+D\,e^{-\mu t}},\quad D=\left(\frac{K}{N_0}\right)^{\nu}-1"),
        ).arrange(DOWN, buff=0.45, aligned_edge=LEFT)
        fit(pasos, max_w=11.5, max_h=4.2)
        pasos.next_to(titulo, DOWN, buff=0.5)

        for p in pasos:
            self.play(Write(p))
            self.wait(0.6)

        resultado = MathTex(
            r"N(t)=K\left[1+\left(\left(\tfrac{K}{N_0}\right)^{\nu}-1\right)"
            r"e^{-\mu t}\right]^{-1/\nu}"
        ).scale(1.0)
        fit(resultado, max_w=12.5, max_h=1.2)
        resultado.to_edge(DOWN, buff=0.6)
        caja = SurroundingRectangle(resultado, color=YELLOW, buff=0.25)
        self.play(Write(resultado), Create(caja))
        self.wait(2.5)
        self.clear_all()

    # ------------------------------------------------------------------
    # 5. Equilibrios e inflexión (con nu variable)
    # ------------------------------------------------------------------
    def seccion_inflexion(self):
        self.next_section("inflexion")
        equil = MathTex(
            r"f'(0)=\frac{\mu}{\nu}>0\ \text{(inestable)},\qquad "
            r"f'(K)=-\mu<0\ \text{(estable)}",
            font_size=36,
        ).to_edge(UP)
        infl = MathTex(
            r"\frac{d^2N}{dt^2}=f'(N)\frac{dN}{dt}=0\ \Rightarrow\ "
            r"N_{\mathrm{inf}}=K\,(1+\nu)^{-1/\nu}",
            font_size=36,
        ).next_to(equil, DOWN, buff=0.35)
        self.play(Write(equil))
        self.play(Write(infl))

        ejes = Axes(
            x_range=[0, T_MAX, 5], y_range=[0, 1.1, 0.5],
            x_length=9, y_length=3.8, tips=False,
            axis_config={"include_numbers": False},
        ).to_edge(DOWN, buff=0.6).shift(LEFT * 1.2)
        linea_K = DashedLine(ejes.c2p(0, K), ejes.c2p(T_MAX, K), color=GREY_B)
        tK = MathTex("K").next_to(linea_K, RIGHT, buff=0.15)
        self.play(Create(ejes), Create(linea_K), FadeIn(tK))

        nu = ValueTracker(1.0)
        curva = always_redraw(lambda: ejes.plot(
            lambda t: sol(t, nu.get_value()), x_range=[0, T_MAX], color=TEAL))
        punto = always_redraw(lambda: Dot(
            ejes.c2p(t_inf(nu.get_value()), n_inf(nu.get_value())),
            color=RED, radius=0.09))
        etiqueta = always_redraw(lambda: VGroup(
            MathTex(rf"\nu={nu.get_value():.2f}", font_size=40),
            MathTex(rf"N_{{\mathrm{{inf}}}}={n_inf(nu.get_value()):.2f}\,K",
                    font_size=40, color=RED),
        ).arrange(DOWN, aligned_edge=LEFT, buff=0.3)
         .next_to(ejes, RIGHT, buff=0.5))

        self.play(Create(curva), FadeIn(punto), FadeIn(etiqueta))
        self.wait(1)
        # nu = 1 -> K/2 ; nu = 0.25 -> 0.41 K ; nu -> 0 -> K/e
        self.play(nu.animate.set_value(0.25), run_time=3)
        self.wait(1)
        self.play(nu.animate.set_value(0.01), run_time=3)
        nota = Text("ν → 0: Gompertz, inflexión en K/e", font_size=24, color=YELLOW)
        nota.next_to(etiqueta, DOWN, buff=0.4)
        self.play(FadeIn(nota))
        self.wait(2.5)
        self.clear_all()

    # ------------------------------------------------------------------
    # 6. Cierre: ¿de dónde sale K?
    # ------------------------------------------------------------------
    def seccion_cierre(self):
        self.next_section("cierre")
        t1 = Text("Los tres modelos suponen un techo K fijo", font_size=38)
        t1.to_edge(UP, buff=1.0)
        self.play(Write(t1))

        ejes = Axes(
            x_range=[0, T_MAX, 5], y_range=[0, 1.6, 0.5],
            x_length=9, y_length=3.8, tips=False,
            axis_config={"include_numbers": False},
        ).shift(DOWN * 0.4)
        curva = ejes.plot(lambda t: sol(t, 1.0), x_range=[0, T_MAX], color=TEAL)
        linea_K = DashedLine(ejes.c2p(0, K), ejes.c2p(T_MAX, K), color=GREY_B)
        tK = MathTex("K").next_to(linea_K, RIGHT, buff=0.15)
        self.play(Create(ejes), Create(curva), Create(linea_K), FadeIn(tK))
        self.wait(0.5)

        flecha = Arrow(ejes.c2p(14, K), ejes.c2p(14, 1.45), color=YELLOW, buff=0)
        self.play(Create(flecha), linea_K.animate.shift(UP * ejes.y_axis.unit_size * 0.4),
                  tK.animate.shift(UP * ejes.y_axis.unit_size * 0.4))
        pregunta = Text("¿Y si K lo construye el propio tumor?",
                        font_size=40, color=YELLOW)
        pregunta.to_edge(DOWN, buff=0.5)
        self.play(Write(pregunta))
        self.wait(3)
        self.clear_all()
