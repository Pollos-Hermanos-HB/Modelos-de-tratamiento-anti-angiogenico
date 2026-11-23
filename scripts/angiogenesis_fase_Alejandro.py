"""
Modelo de angiogenesis (Hahnfeldt et al. / Sachs et al. 2001), c(t) = 0:

    dN/dt = -mu * N * ln(N/K)            (Gompertz con K dinamico)
    dK/dt = omega*N - gamma*N^(2/3)*K

Equilibrio: N* = K* = (omega/gamma)^(3/2)
Jacobiano en el equilibrio: J = [[-mu, mu], [omega/3, -omega]]

Genera:
  fase.png        plano de fase (campo, nulclinas, trayectorias)
  series.png      N(t), K(t) para varias condiciones iniciales
  comparacion.png Euler vs RK45
  y imprime: convergencia, Jacobiano analitico vs numerico, tabla de autovalores.
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp


# ---------------------------------------------------------------- modelo
def campo(N, K, mu, om, ga):
    return (-mu * N * np.log(N / K),
            om * N - ga * N ** (2 / 3) * K)


def equilibrio(om, ga):
    return (om / ga) ** 1.5


def jacobiano(mu, om):
    return np.array([[-mu, mu], [om / 3, -om]])


def jacobiano_numerico(mu, om, ga, h=1e-6):
    Ns = equilibrio(om, ga)
    p = np.array([Ns, Ns])
    f = lambda v: np.array(campo(v[0], v[1], mu, om, ga))
    J = np.zeros((2, 2))
    for j in range(2):
        e = np.zeros(2)
        e[j] = h
        J[:, j] = (f(p + e) - f(p - e)) / (2 * h)
    return J


def resolver(N0, K0, mu, om, ga, T=80, **kw):
    def rhs(t, y):
        N, K = np.maximum(y, 1e-9)
        return campo(N, K, mu, om, ga)
    return solve_ivp(rhs, [0, T], [N0, K0], rtol=1e-8, atol=1e-10, **kw)


def euler(N0, K0, mu, om, ga, T=80, h=0.01):
    n = int(T / h)
    t = np.linspace(0, T, n + 1)
    Y = np.zeros((n + 1, 2))
    Y[0] = (N0, K0)
    for i in range(n):
        N, K = np.maximum(Y[i], 1e-9)
        dN, dK = campo(N, K, mu, om, ga)
        Y[i + 1] = Y[i] + h * np.array([dN, dK])
    return t, Y


# ---------------------------------------------------------------- figuras
def plano_fase(mu=0.5, om=1.0, ga=1.0, archivo="fase.png"):
    Ns = equilibrio(om, ga)
    lim = 3 * Ns
    fig, ax = plt.subplots(figsize=(7, 6))

    X, Y = np.meshgrid(np.linspace(0.05, lim, 22), np.linspace(0.05, lim, 22))
    U, V = campo(X, Y, mu, om, ga)
    n = np.hypot(U, V) + 1e-12
    ax.quiver(X, Y, U / n, V / n, color="lightsteelblue", angles="xy")

    x = np.linspace(0.05, lim, 300)
    ax.plot(x, x, "r--", label="dN/dt = 0  (K = N)")
    ax.plot(x, (om / ga) * x ** (1 / 3), "g--", label="dK/dt = 0")

    inicios = [(0.1, 0.1), (lim, 0.2), (0.2, lim), (lim, lim),
               (0.5 * lim, 0.1), (0.1, 0.6 * lim), (0.8 * Ns, 2 * Ns)]
    print("Convergencia (condicion inicial -> estado final):")
    for N0, K0 in inicios:
        s = resolver(N0, K0, mu, om, ga, max_step=0.1)
        ax.plot(s.y[0], s.y[1], lw=1.5)
        ax.plot(N0, K0, "k.", ms=4)
        print(f"  ({N0:.2f},{K0:.2f}) -> ({s.y[0, -1]:.4f},{s.y[1, -1]:.4f})")

    ax.plot(Ns, Ns, "ko", ms=8, label=f"equilibrio ({Ns:.2f}, {Ns:.2f})")
    ax.set(xlim=(0, lim), ylim=(0, lim), xlabel="N (celulas tumorales)",
           ylabel="K (capacidad de carga)",
           title=f"Plano de fase: mu={mu}, omega={om}, gamma={ga}")
    ax.legend(loc="upper left")
    fig.savefig(archivo, dpi=150, bbox_inches="tight")
    plt.close(fig)


def series_tiempo(mu=0.5, om=1.0, ga=1.0, archivo="series.png"):
    Ns = equilibrio(om, ga)
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(10, 4), sharey=True)
    for N0, K0 in [(0.1, 0.1), (3 * Ns, 0.3), (0.3, 3 * Ns), (3 * Ns, 3 * Ns)]:
        s = resolver(N0, K0, mu, om, ga, T=40, max_step=0.1)
        a1.plot(s.t, s.y[0])
        a2.plot(s.t, s.y[1])
    for a, nm in [(a1, "N(t)"), (a2, "K(t)")]:
        a.axhline(Ns, color="k", ls=":", label="equilibrio")
        a.set(xlabel="t", title=nm)
    a1.legend()
    fig.savefig(archivo, dpi=150, bbox_inches="tight")
    plt.close(fig)


def comparar_euler(mu=0.5, om=1.0, ga=1.0, archivo="comparacion.png"):
    N0, K0 = 0.2, 2.0
    s = resolver(N0, K0, mu, om, ga, T=40, max_step=0.05)
    fig, ax = plt.subplots(figsize=(6, 5))
    ax.plot(s.y[0], s.y[1], "k", lw=2, label="RK45 (referencia)")
    for h, c in [(0.5, "tab:red"), (0.1, "tab:orange"), (0.01, "tab:green")]:
        t, Y = euler(N0, K0, mu, om, ga, T=40, h=h)
        ax.plot(Y[:, 0], Y[:, 1], "--", color=c, label=f"Euler h={h}")
    ax.plot(equilibrio(om, ga), equilibrio(om, ga), "ko")
    ax.set(xlabel="N", ylabel="K", title="Euler vs RK45")
    ax.legend()
    fig.savefig(archivo, dpi=150, bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------- main
if __name__ == "__main__":
    mu, om, ga = 0.5, 1.0, 1.0

    plano_fase(mu, om, ga)
    series_tiempo(mu, om, ga)
    comparar_euler(mu, om, ga)

    print("\nJ analitico:\n", jacobiano(mu, om))
    print("J numerico:\n", jacobiano_numerico(mu, om, ga).round(5))
    J = jacobiano(mu, om)
    print(f"traza = {np.trace(J):.4f}  (esperado {-(mu + om):.4f})")
    print(f"det   = {np.linalg.det(J):.4f}  (esperado {2 * mu * om / 3:.4f})")
    print(f"discriminante = {mu**2 - 2*mu*om/3 + om**2:.4f} (>0: autovalores reales)")

    print("\n  mu    om    ga |   N*     lambda1   lambda2")
    for m, o, g in [(0.5, 1, 1), (0.5, 2, 1), (0.5, 4, 1),
                    (0.5, 1, 2), (0.5, 1, 4), (1.0, 1, 1)]:
        lam = np.sort(np.linalg.eigvals(jacobiano(m, o)).real)
        print(f"{m:5.2f} {o:5.2f} {g:5.2f} | {equilibrio(o, g):6.3f} "
              f"{lam[0]:9.4f} {lam[1]:9.4f}")
