""""
Modelo 1-EDO con tratamiento constante (logístico generalizado)
Método de Euler puro
"""

import numpy as np
import matplotlib.pyplot as plt

# ============================================================
# Parámetros del modelo
# ============================================================
mu = 0.1
nu = 1.0  # logístico clásico (nu=1)
K = 1.0
alpha = 1.0
c_star = mu / (nu * alpha)  # = 0.1

dt = 0.005
t_max = 25000
rel_tol = 1e-4  # tolerancia relativa


# ============================================================
# Función g(N) y tiempo de convergencia con Euler
# ============================================================
def g(N, c):
    return -alpha * c * N + (mu * N / nu) * (1.0 - (N / K) ** nu)


def tiempo_convergencia(c):
    """Parte desde una perturbación relativa fija (5%) por debajo de N*."""
    if c >= c_star - 1e-14:
        return np.nan
    N_star = K * (1.0 - c / c_star) ** (1.0 / nu)
    if N_star < 1e-12:
        return np.nan

    N = N_star * 0.95
    t = 0.0
    while t < t_max:
        if abs(N - N_star) < rel_tol * N_star:
            return t
        N = N + dt * g(N, c)
        if N < 0.0:
            N = 0.0
        t += dt
    return t


# ============================================================
# Barrido de valores de c
# ============================================================
cs_far = np.linspace(0.25 * c_star, 0.90 * c_star, 18)
cs_near = np.linspace(0.90 * c_star, 0.997 * c_star, 22)
cs = np.unique(np.concatenate([cs_far, cs_near]))

tiempos_num = []
taus_teo = []
Nstars = []

for i, c in enumerate(cs, start=1):
    print(f"{i}/{len(cs)} c={c:.5f}")
    t_conv = tiempo_convergencia(c)
    tiempos_num.append(t_conv)
    g_prime_abs = mu * (1.0 - c / c_star)
    taus_teo.append(1.0 / g_prime_abs)
    Nstars.append(K * (1.0 - c / c_star) ** (1.0 / nu))

tiempos_num = np.array(tiempos_num)
taus_teo = np.array(taus_teo)
Nstars = np.array(Nstars)
delta = c_star - cs

# ============================================================
# Gráfica
# ============================================================
plt.rcParams.update(
    {
        "font.size": 11,
        "axes.labelsize": 12,
        "axes.titlesize": 13,
        "legend.fontsize": 10,
        "figure.facecolor": "white",
        "axes.facecolor": "white",
        "axes.grid": True,
        "grid.alpha": 0.3,
        "lines.linewidth": 2,
    }
)

fig, ax = plt.subplots(figsize=(11, 6))

# --- Tiempo de convergencia vs c ---
mask = np.isfinite(tiempos_num)
ax.plot(
    cs[mask],
    tiempos_num[mask],
    "o-",
    color="#1f77b4",
    markersize=5,
    label="Tiempo numérico (Euler)",
)
ax.plot(cs, taus_teo, "--", color="#ff7f0e", lw=2, label=r"Teórico $1/|g'(N^*)|$")
ax.axvline(c_star, color="gray", ls=":", lw=1.5, label=r"$c^*$")
ax.set_xlabel(r"Dosis $c$")
ax.set_ylabel("Tiempo de convergencia")
ax.legend(loc="upper left", framealpha=0.95)
ax.set_xlim(cs[0] * 0.95, c_star * 1.02)
ax.grid(True, alpha=0.3)

fig.text(
    0.5,
    0.01,
    r"Parámetros: $\mu=0.1$, $\nu=1$, $K=1$, $\alpha=1$  $\Rightarrow$  $c^*=0.1$.  "
    r"Perturbación inicial: $5\%$ por debajo de $N^*$.  Tolerancia relativa: $10^{-4}$.",
    ha="center",
    fontsize=9,
    style="italic",
    color="#444444",
)

plt.tight_layout(rect=[0, 0.04, 1, 1])
plt.savefig("euler-func_combergence.png", dpi=150)
plt.show()
