"""Modelos de una EDO: logístico generalizado y Gompertz.

Genera la Figura 1 del informe y la tabla de errores de Euler.
Uso (desde la carpeta scripts/):
    python modelos_1edo.py              # guarda la figura del informe y abre 3 ventanas interactivas
    python modelos_1edo.py --no-show    # solo guarda (sin ventana, p. ej. en un servidor)
Salida: ../figuras/fig_modelos_1edo.pdf (vectorial) y .png

Notación (igual que en el informe y que el paper):
    N  número de células        K  capacidad de carga
    mu, nu  parámetros del modelo (nu = 0 significa Gompertz)
    f(N) = N * F(N/K)           lado derecho de dN/dt
"""
import sys
from pathlib import Path
import numpy as np
import matplotlib
if "--no-show" in sys.argv:      # sin ventana: backend que solo escribe archivos
    matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ---------- Parámetros ----------
# Ninguno viene del paper: son SUPUESTOS del equipo (se justifican en el informe).
MU = 0.5          # velocidad de crecimiento (tiempo en unidades arbitrarias)
K = 1.0           # capacidad de carga; N se normaliza, así N/K está en [0, 1]
N0 = 0.01         # tamaño inicial del tumor (N0 < K)
T = 20.0          # tiempo total de simulación
NU_MAMA = 0.25    # nu típico para tumores de mama (Sachs et al., 2001)
H_FIG = 1.0                  # paso grande de Euler, para la Fig. 1(a)
H_LIST = [0.5, 0.1, 0.01]    # pasos para estudiar el error (Fig. 1(b) y tabla)
T_LARGO = 200.0              # horizonte de los datos en las ventanas (largo plazo)
T_VISTA = 60.0               # tramo que se ve al abrir la ventana (se puede ampliar)

CARPETA_FIG = Path(__file__).resolve().parent.parent / "figuras"


# ---------- Modelo ----------
def F(x, mu=MU, nu=1.0):
    """F(x) con x = N/K. nu = 0 es Gompertz."""
    x = np.asarray(x, dtype=float)
    if nu == 0:
        return -mu * np.log(x)
    return (mu / nu) * (1.0 - x**nu)


def f(N, mu=MU, nu=1.0, K=K):
    """Lado derecho de dN/dt = f(N) = N F(N/K)."""
    return N * F(N / K, mu, nu)


def exacta(t, N0=N0, mu=MU, nu=1.0, K=K):
    """Solución exacta: ec. (3) (Gompertz) o ec. (4) (generalizado)."""
    t = np.asarray(t, dtype=float)
    if nu == 0:
        return K * np.exp(np.log(N0 / K) * np.exp(-mu * t))
    Q = (K / N0) ** nu - 1.0
    return K / (1.0 + Q * np.exp(-mu * t)) ** (1.0 / nu)


def n_inflexion(nu, K=K):
    """N en el punto de inflexión: K (1+nu)^(-1/nu); K/e si nu = 0."""
    return K / np.e if nu == 0 else K * (1.0 + nu) ** (-1.0 / nu)


def tasa_especifica(N, mu=MU, nu=1.0, K=K):
    """(1/N) dN/dt = F(N/K)."""
    return F(np.asarray(N) / K, mu, nu)


# ---------- Euler ----------
def euler(nu, h, N0=N0, mu=MU, T=T, K=K):
    """Euler explícito: N_{k+1} = N_k + h f(N_k). Devuelve (t, N)."""
    n = int(round(T / h))
    t = np.linspace(0.0, n * h, n + 1)
    N = np.empty(n + 1)
    N[0] = N0
    for k in range(n):
        N[k + 1] = N[k] + h * f(N[k], mu, nu, K)
    return t, N


def error_maximo(nu, h):
    t, N = euler(nu, h)
    return np.max(np.abs(N - exacta(t, nu=nu)))


# ---------- Figura 1 ----------
def figura():
    """Figura 1 del informe (compacta, t en [0, T]); se guarda en ../figuras/."""
    plt.rcParams.update({"font.size": 8, "axes.grid": True, "grid.alpha": 0.3})
    fig, ax = plt.subplots(1, 3, figsize=(7.2, 2.35))
    tt = np.linspace(0, T, 600)
    casos = [(1.0, "C0", "ν=1 (logístico)"), (NU_MAMA, "C2", f"ν={NU_MAMA}"),
             (0.0, "C3", "Gompertz")]

    # (a) exactas, Euler con paso grande y puntos de inflexión
    for nu, color, etiqueta in casos:
        ax[0].plot(tt, exacta(tt, nu=nu), color, label=etiqueta)
        t, N = euler(nu, H_FIG)
        ax[0].plot(t, N, ".", color=color, ms=3)
        Ni = n_inflexion(nu)
        ax[0].plot(tt[np.argmin(np.abs(exacta(tt, nu=nu) - Ni))], Ni, "k^", ms=4)
    ax[0].set_xlabel("t (u.a.)"); ax[0].set_ylabel("N/K")
    ax[0].set_title(f"(a) exacta (línea), Euler h={H_FIG:g} (puntos)", fontsize=7.5)
    ax[0].legend(fontsize=6.5, loc="lower right")

    # (b) error de Euler (logístico y Gompertz, tres pasos)
    for nu, color, nombre in [(1.0, "C0", "log."), (0.0, "C3", "Gomp.")]:
        for h, ls in zip(H_LIST, ["-", "--", ":"]):
            t, N = euler(nu, h)
            err = np.abs(N - exacta(t, nu=nu)) + 1e-16
            ax[1].semilogy(t, err, color, ls=ls, lw=1.1, label=f"{nombre} h={h:g}")
    ax[1].set_ylim(1e-8, 1)
    ax[1].set_xlabel("t (u.a.)")
    ax[1].set_title("(b) error |Euler − exacta|", fontsize=7.5)
    ax[1].legend(fontsize=5.5, ncol=2, loc="lower left")

    # (c) tasa específica
    x = np.linspace(0.005, 1.5, 400)
    for nu, color, _ in casos:
        ax[2].plot(x, tasa_especifica(x * K, nu=nu), color)
    ax[2].set_ylim(-0.5, 4)
    ax[2].set_xlabel("N/K"); ax[2].set_ylabel("(1/N)dN/dt")
    ax[2].set_title("(c) tasa específica", fontsize=7.5)

    plt.tight_layout()
    CARPETA_FIG.mkdir(exist_ok=True)
    fig.savefig(CARPETA_FIG / "fig_modelos_1edo.pdf")
    fig.savefig(CARPETA_FIG / "fig_modelos_1edo.png", dpi=300)
    print(f"Figura guardada en {CARPETA_FIG}")
    plt.close(fig)


# ---------- Ventanas interactivas (una por gráfica, comportamiento a largo plazo) ----------
def _zoom_con_rueda(fig, ax, paso=1.25):
    """Rueda del ratón = zoom alrededor del cursor. Las herramientas de la barra
    (cruz para mover, lupa para recuadro) y 'inicio' (restablecer) siguen disponibles."""
    def en_rueda(evento):
        if evento.inaxes is not ax or evento.xdata is None:
            return
        f = 1.0 / paso if evento.button == "up" else paso
        x0, x1 = ax.get_xlim()
        ax.set_xlim(evento.xdata - (evento.xdata - x0) * f, evento.xdata + (x1 - evento.xdata) * f)
        y0, y1 = ax.get_ylim()
        if ax.get_yscale() == "log":
            ly, l0, l1 = np.log10(evento.ydata), np.log10(y0), np.log10(y1)
            ax.set_ylim(10 ** (ly - (ly - l0) * f), 10 ** (ly + (l1 - ly) * f))
        else:
            ax.set_ylim(evento.ydata - (evento.ydata - y0) * f, evento.ydata + (y1 - evento.ydata) * f)
        fig.canvas.draw_idle()
    fig.canvas.mpl_connect("scroll_event", en_rueda)


def _ventana(titulo):
    fig, ax = plt.subplots(figsize=(10, 6.2), constrained_layout=True)
    try:
        fig.canvas.manager.set_window_title(titulo)
    except Exception:
        pass
    ax.grid(True, alpha=0.3)
    fig.get_layout_engine().set(rect=(0, 0.03, 1, 0.97))   # deja espacio abajo para la ayuda
    fig.text(0.5, 0.008, "Rueda: zoom · cruz de la barra: mover · lupa: recuadro · casa: restablecer",
             ha="center", va="bottom", fontsize=8, color="gray")
    _zoom_con_rueda(fig, ax)
    return fig, ax


def figuras_interactivas():
    """Abre tres ventanas (a, b, c) con datos hasta T_LARGO para explorar libremente."""
    plt.rcParams.update({"font.size": 11})
    casos = [(1.0, "C0", "ν = 1 (logístico)"), (NU_MAMA, "C2", f"ν = {NU_MAMA}"),
             (0.0, "C3", "Gompertz")]
    tt = np.linspace(0.0, T_LARGO, 6000)

    # (a) exactas, Euler con paso grande, inflexiones y la capacidad de carga K
    fig, ax = _ventana("(a) Crecimiento N(t)")
    ax.axhline(1.0, color="gray", ls="--", lw=1, label="K")
    for nu, color, etiqueta in casos:
        ax.plot(tt, exacta(tt, nu=nu) / K, color, label=etiqueta + " (exacta)")
        t, N = euler(nu, H_FIG, T=T_LARGO)
        ax.plot(t, N / K, ".", color=color, ms=3, label=f"{etiqueta} (Euler h={H_FIG:g})")
        Ni = n_inflexion(nu)
        ti = tt[np.argmin(np.abs(exacta(tt, nu=nu) - Ni))]
        ax.plot(ti, Ni / K, "k^", ms=8)
    ax.set_xlim(0, T_VISTA); ax.set_ylim(-0.03, 1.12)
    ax.set_xlabel("t (u.a.)"); ax.set_ylabel("N/K")
    ax.set_title("(a) Exacta y Euler; ▲ = punto de inflexión")
    ax.legend(loc="lower right", fontsize=8, ncol=2)

    # (b) error de Euler, escala logarítmica
    fig, ax = _ventana("(b) Error de Euler")
    for nu, color, nombre in [(1.0, "C0", "logístico"), (0.0, "C3", "Gompertz")]:
        for h, ls in zip(H_LIST, ["-", "--", ":"]):
            t, N = euler(nu, h, T=T_LARGO)
            ax.semilogy(t, np.abs(N - exacta(t, nu=nu)) + 1e-16, color, ls=ls, lw=1.4,
                        label=f"{nombre}, h = {h:g}")
    ax.set_xlim(0, T_VISTA); ax.set_ylim(1e-16, 1)
    ax.set_xlabel("t (u.a.)"); ax.set_ylabel("|Euler − exacta|")
    ax.set_title("(b) Error absoluto de Euler (1e-16 = cero numérico)")
    ax.legend(loc="lower left", fontsize=8, ncol=2)

    # (c) tasa específica; dominio amplio (N/K de 0 a 3)
    fig, ax = _ventana("(c) Tasa específica")
    x = np.linspace(0.001, 3.0, 4000)
    for nu, color, etiqueta in casos:
        ax.plot(x, tasa_especifica(x * K, nu=nu), color, label=etiqueta)
    ax.axhline(0.0, color="gray", lw=1); ax.axvline(1.0, color="gray", ls="--", lw=1)
    ax.set_xlim(0, 1.6); ax.set_ylim(-1.0, 4.0)
    ax.set_xlabel("N/K"); ax.set_ylabel("(1/N) dN/dt")
    ax.set_title("(c) Tasa específica de crecimiento (se anula en N = K)")
    ax.legend(loc="upper right")

    plt.show()                   # abre las tres ventanas; el programa sigue al cerrarlas


if __name__ == "__main__":
    print("Error máximo de Euler (nu, h, error):")
    for nu in (1.0, 0.0):
        for h in H_LIST:
            print(f"  nu={nu:g}  h={h:<5g}  {error_maximo(nu, h):.2e}")
    print("Puntos de inflexión / K:",
          {nu: round(n_inflexion(nu) / K, 3) for nu in (1.0, NU_MAMA, 0.0)})
    figura()                     # guarda la figura del informe
    if "--no-show" not in sys.argv:
        figuras_interactivas()   # ventanas para explorar