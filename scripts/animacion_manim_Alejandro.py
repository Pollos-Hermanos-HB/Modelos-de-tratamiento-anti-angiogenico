"""
Animacion Manim Community: particulas (N, K) convergiendo al equilibrio.

Ejecutar (calidad baja, rapido):
    manim -pql animacion_manim.py AngiogenesisConvergence
Calidad alta:
    manim -pqh animacion_manim.py AngiogenesisConvergence

No usa LaTeX (solo Text), asi que funciona aunque no tengas TeX instalado.
"""
import numpy as np
from scipy.integrate import solve_ivp
from manim import *

MU, OM, GA = 0.5, 1.0, 1.0
NS = (OM / GA) ** 1.5          # equilibrio N* = K*
LIM = 3.0                      # rango de los ejes
T_FIN = 20.0                   # tiempo simulado


def campo(N, K):
    N = max(N, 1e-9)
    K = max(K, 1e-9)
    return np.array([-MU * N * np.log(N / K),
                     OM * N - GA * N ** (2 / 3) * K])


def trayectoria(N0, K0):
    sol = solve_ivp(lambda t, y: campo(*y), [0, T_FIN], [N0, K0],
                    t_eval=np.linspace(0, T_FIN, 800), rtol=1e-8, atol=1e-10)
    return sol.t, sol.y


class AngiogenesisConvergence(Scene):
    def construct(self):
        titulo = Text("Angiogenesis: convergencia al equilibrio", font_size=30)
        titulo.to_edge(UP)

        axes = Axes(
            x_range=[0, LIM, 1], y_range=[0, LIM, 1],
            x_length=6.5, y_length=5.0,
            axis_config={"include_numbers": False, "include_tip": True},
        ).shift(0.4 * DOWN + 1.5 * LEFT)
        lx = Text("N (celulas)", font_size=22).next_to(axes.x_axis, DOWN, buff=0.15)
        ly = Text("K (capacidad)", font_size=22).next_to(axes.y_axis, LEFT, buff=0.15)
        ly.rotate(PI / 2)

        self.play(Write(titulo), Create(axes), FadeIn(lx), FadeIn(ly))

        # ---- campo vectorial (flechas normalizadas)
        flechas = VGroup()
        for N in np.linspace(0.15, LIM - 0.15, 14):
            for K in np.linspace(0.15, LIM - 0.15, 12):
                v = campo(N, K)
                n = np.linalg.norm(v)
                if n < 1e-9:
                    continue
                v = v / n * 0.18
                p0 = np.array([N, K])
                flechas.add(Arrow(axes.c2p(*(p0 - v / 2)),
                                  axes.c2p(*(p0 + v / 2)),
                                  buff=0, stroke_width=2,
                                  max_tip_length_to_length_ratio=0.5,
                                  color=BLUE_E))
        self.play(FadeIn(flechas, lag_ratio=0.02), run_time=2)

        # ---- nulclinas
        nul_N = axes.plot(lambda x: x, x_range=[0.02, LIM], color=RED)
        nul_K = axes.plot(lambda x: (OM / GA) * x ** (1 / 3),
                          x_range=[0.02, LIM], color=GREEN)
        leyenda = VGroup(
            VGroup(Line(ORIGIN, 0.5 * RIGHT, color=RED),
                   Text("dN/dt = 0  (K = N)", font_size=20)).arrange(RIGHT),
            VGroup(Line(ORIGIN, 0.5 * RIGHT, color=GREEN),
                   Text("dK/dt = 0", font_size=20)).arrange(RIGHT),
            VGroup(Dot(color=YELLOW),
                   Text("equilibrio (N*, K*)", font_size=20)).arrange(RIGHT),
        ).arrange(DOWN, aligned_edge=LEFT).to_corner(UR).shift(0.2 * DOWN)
        self.play(Create(nul_N), Create(nul_K), FadeIn(leyenda))

        # ---- equilibrio
        eq = Dot(axes.c2p(NS, NS), color=YELLOW, radius=0.11)
        halo = Circle(radius=0.3, color=YELLOW).move_to(eq)
        self.play(FadeIn(eq), Create(halo))
        self.play(FadeOut(halo))

        # ---- particulas
        inicios = [(0.2, 0.2), (2.8, 0.3), (0.3, 2.8), (2.8, 2.8),
                   (1.6, 0.2), (0.2, 1.8), (2.2, 1.2), (0.8, 2.4)]
        colores = [ORANGE, PINK, TEAL, PURPLE, MAROON, GOLD, LIGHT_BROWN, WHITE]

        datos, puntos, estelas = [], [], []
        for (N0, K0), col in zip(inicios, colores):
            t, y = trayectoria(N0, K0)
            datos.append((t, y))
            d = Dot(axes.c2p(N0, K0), color=col, radius=0.07)
            tr = TracedPath(d.get_center, stroke_color=col, stroke_width=3)
            puntos.append(d)
            estelas.append(tr)
        self.add(*estelas)
        self.play(*[FadeIn(d) for d in puntos])

        def mover(i):
            t, y = datos[i]

            def f(m, alpha):
                tt = alpha * T_FIN
                N = np.interp(tt, t, y[0])
                K = np.interp(tt, t, y[1])
                m.move_to(axes.c2p(N, K))
            return UpdateFromAlphaFunc(puntos[i], f)

        self.play(*[mover(i) for i in range(len(inicios))],
                  run_time=9, rate_func=linear)

        msg = Paragraph(
            "Todas las trayectorias",
            "llegan al equilibrio",
            "(atractor global)",
            font_size=24, color=YELLOW, alignment="left", line_spacing=1.1
        )
        msg.next_to(axes, RIGHT, buff=0.6).shift(UP * 0.5)
        self.play(Write(msg))
        self.wait(2)


class CambioGamma(Scene):
    """Extra: muestra como el equilibrio se mueve al cambiar gamma (N* = (om/ga)^1.5)."""

    def construct(self):
        axes = Axes(x_range=[0, 3, 1], y_range=[0, 3, 1],
                    x_length=6, y_length=5).shift(1.5 * LEFT)
        self.add(axes, Text("N", font_size=22).next_to(axes.x_axis, DOWN),
                 Text("K", font_size=22).next_to(axes.y_axis, LEFT))
        diag = axes.plot(lambda x: x, x_range=[0.02, 3], color=RED)
        self.add(diag)

        g = ValueTracker(0.6)
        nul_K = always_redraw(lambda: axes.plot(
            lambda x: (OM / g.get_value()) * x ** (1 / 3),
            x_range=[0.02, 3], color=GREEN))
        eq = always_redraw(lambda: Dot(
            axes.c2p(*(2 * [(OM / g.get_value()) ** 1.5])), color=YELLOW))
        txt = always_redraw(lambda: Text(
            f"gamma = {g.get_value():.2f}   N* = {(OM / g.get_value()) ** 1.5:.2f}",
            font_size=26).to_corner(UR))
        self.add(nul_K, eq, txt)
        self.play(g.animate.set_value(3.0), run_time=6, rate_func=smooth)
        self.play(g.animate.set_value(0.6), run_time=6, rate_func=smooth)
        self.wait()
