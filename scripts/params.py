import numpy as np

# ---------- Parámetros ----------

# Seccion 1
MU = 0.5          # velocidad de crecimiento (tiempo en unidades arbitrarias)
K = 1.0           # capacidad de carga; N se normaliza, así N/K está en [0, 1]
N0 = 0.01         # tamaño inicial del tumor (N0 < K)
T = 20.0          # tiempo total de simulación
NU_MAMA = 0.25    # nu típico para tumores de mama (Sachs et al., 2001)
H_FIG = 1.0                  # paso grande de Euler, para la Fig. 1(a)
H_LIST = [0.5, 0.1, 0.01]    # pasos para estudiar el error (Fig. 1(b) y tabla)
T_LARGO = 200.0              # horizonte de los datos en las ventanas (largo plazo)
T_VISTA = 60.0               # tramo que se ve al abrir la ventana (se puede ampliar)


# Seccion 5
# --- Semilla fija (el bloque es determinista; se fija por convención del equipo)
SEED = 20260101
np.random.seed(SEED)

# --- Unidades: tiempo en horas (h), dosis en Gy, U en DSB/célula
TAU      = 1.0      # h            tiempo característico de reparación (SUPUESTO)
DELTA    = 35.0     # DSB/(cél·Gy) DSBs inducidos por Gy               (SUPUESTO)
ALPHA    = 0.15     # 1/Gy         daño letal lineal                   (SUPUESTO)
BETA     = 0.05     # 1/Gy^2       => alpha/beta = 3 Gy                (SUPUESTO)
KAPPA    = 2.0*BETA/(DELTA**2*TAU)   # 1/(DSB^2 h), de beta = delta^2*kappa*tau/2
GAMMA_LQ = 2.0e-3   # 1/(DSB h)    interacción de lesiones, ec. (8)    (SUPUESTO)

# --- Numérico
H_REL  = 0.02       # paso h = H_REL*TAU
TAIL   = 15.0       # tiempo de cola tras la irradiación, en unidades de TAU

# --- Barridos para las figuras
DOSES        = np.linspace(0.0, 10.0, 41)         # Gy
AB_RATIOS    = [3.0, 10.0, 30.0]                  # Gy (alpha fijo, beta = alpha/ratio)
T_OVER_TAU   = [0.0, 1.0, 4.0, 16.0]              # 0 = dosis aguda
FIG_DIR      = "figs"

# --- Estilo común de figuras (vector PDF)
RC = {
    "figure.figsize": (6.4, 4.2), "font.size": 10, "axes.grid": True,
    "grid.alpha": 0.3, "axes.titlesize": 11, "legend.fontsize": 8.5,
    "savefig.format": "pdf", "pdf.fonttype": 42, "lines.linewidth": 1.8,
}