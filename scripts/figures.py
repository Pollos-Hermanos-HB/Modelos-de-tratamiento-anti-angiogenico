"""
figures.py -- Genera TODAS las figuras (PDF vectorial) y la tabla de resultados.
Uso:  python figures.py
"""
import os, csv
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import params as P
import lq

os.makedirs(P.FIG_DIR, exist_ok=True)
plt.rcParams.update(P.RC)
def save(fig, name):
    fig.tight_layout(); fig.savefig(os.path.join(P.FIG_DIR, name)); plt.close(fig)

def kappa_of(beta):               # beta = delta^2 kappa tau / 2
    return 2.0*beta/(P.DELTA**2*P.TAU)

# ======================================================= Fig. 1: G vs T/tau
x = np.logspace(-2, 2, 400)
xn = np.array([0.1, 0.5, 1, 2, 5, 10, 30])
Gn = [lq.G_numeric(lq.schedule_const(1.0, xi*P.TAU), 1.0) for xi in xn]
yn = np.array([0.1, 0.5, 1, 3, 10])
Gen = [lq.G_numeric(lq.schedule_exp(1.0, yi*P.TAU), 1.0) for yi in yn]
fig, ax = plt.subplots()
ax.loglog(x, lq.G_constant(x), label=r"dosis constante: $2(x-1+e^{-x})/x^2$")
ax.loglog(x, lq.G_exponential(x), "--", label=r"perfil exponencial: $1/(1+\theta/\tau)$")
ax.loglog(x[x > 3], 2/x[x > 3], ":", c="k", label=r"asíntota $2/x$")
ax.axhline(1.0, c="gray", ls="-.", lw=1, label="dosis aguda: $G=1$")
ax.plot(xn, Gn, "o", mfc="none", c="C0", label=r"numérico (8)-(9), $\gamma=0$")
ax.plot(yn, Gen, "s", mfc="none", c="C1")
ax.set_xlabel(r"$T/\tau$  (o $\theta/\tau$ para el perfil exponencial)"); ax.set_ylabel(r"$G$ (adimensional)")
ax.set_title("Función de Lea-Catcheside según el tiempo de irradiación"); ax.legend()
save(fig, "fig_G_vs_Ttau.pdf")

# ============================================ Fig. 2: numérico vs LQ (aguda / constante)
Dg = np.linspace(0.5, 10, 20)
cases = {"Dosis aguda": (lambda D: lq.schedule_acute(D), lambda: lq.G_acute()),
         r"Dosis constante, $T=4\tau$": (lambda D: lq.schedule_const(D, 4*P.TAU), lambda: lq.G_constant(4.0)[0])}
fig, axs = plt.subplots(2, 2, figsize=(9, 6.2), sharex=True)
rows = []
for j, (name, (sch, Gf)) in enumerate(cases.items()):
    G = Gf()
    Dfine = np.linspace(0, 10, 200)
    ref = lq.lnS_LQ(Dfine, G)
    n0 = np.array([lq.simulate(sch(D), gamma=0.0)[0] for D in Dg])
    n1 = np.array([lq.simulate(sch(D), gamma=P.GAMMA_LQ)[0] for D in Dg])
    ne = np.array([lq.simulate(sch(D), gamma=0.0, method="euler")[0] for D in Dg])
    ax = axs[0, j]
    ax.semilogy(Dfine, np.exp(ref), c="k", label=r"fórmula LQ ($\gamma=0$)")
    ax.semilogy(Dg, np.exp(n0), "o", mfc="none", c="C0", label=r"RK4, $\gamma=0$")
    ax.semilogy(Dg, np.exp(n1), "--", c="C3", label=fr"RK4, $\gamma={P.GAMMA_LQ:g}$ DSB$^{{-1}}$h$^{{-1}}$")
    ax.set_title(name); ax.set_ylabel(r"$S=N/N_0$")
    if j == 0: ax.legend()
    ax = axs[1, j]
    lqD = lq.lnS_LQ(Dg, G)
    ax.semilogy(Dg, np.abs(n0 - lqD) + 1e-18, "o-", c="C0", ms=3, label=r"|error| RK4, $\gamma=0$")
    ax.semilogy(Dg, np.abs(ne - lqD) + 1e-18, "s-", c="C2", ms=3, label=r"|error| Euler, $\gamma=0$")
    ax.semilogy(Dg, np.abs(n1 - lqD), "^--", c="C3", ms=3, label=r"|desviación| $\gamma>0$")
    ax.set_xlabel("Dosis total $D$ (Gy)"); ax.set_ylabel(r"$|\Delta \ln S|$")
    if j == 0: ax.legend()
    for D, a, b, c, d in zip(Dg, lqD, n0, ne, n1):
        rows.append([name.replace("$", "").replace("\\", ""), D, a, b, c, d])
save(fig, "fig_LQ_vs_numerico.pdf")
with open("resultados_lnS.csv", "w", newline="") as fh:
    w = csv.writer(fh); w.writerow(["caso", "D_Gy", "lnS_LQ", "lnS_RK4_g0", "lnS_Euler_g0", "lnS_RK4_gpos"]); w.writerows(rows)

# ================================== Fig. 3: supervivencia, alpha/beta y tiempo de irradiación
Dc = np.linspace(0, 12, 200)
fig, axs = plt.subplots(1, 3, figsize=(11, 3.8), sharey=True)
for r, c in zip(P.AB_RATIOS, ["C0", "C1", "C2"]):
    beta = P.ALPHA/r
    axs[0].semilogy(Dc, np.exp(lq.lnS_LQ(Dc, 1.0, beta=beta)), c=c, label=fr"$\alpha/\beta={r:g}$ Gy")
    for D in (2, 6, 10):                                   # verificación numérica
        L = lq.simulate(lq.schedule_acute(D), kappa=kappa_of(beta))[0]
        axs[0].semilogy(D, np.exp(L), "o", mfc="none", c=c)
axs[0].set_title(r"Dosis aguda ($G=1$), $\alpha$ fijo"); axs[0].set_ylabel(r"$S=N/N_0$")
for k, r in enumerate([3.0, 10.0]):
    ax = axs[k+1]; beta = P.ALPHA/r
    for xt, c in zip(P.T_OVER_TAU, ["k", "C0", "C1", "C3"]):
        G = 1.0 if xt == 0 else lq.G_constant(xt)[0]
        lab = "aguda" if xt == 0 else fr"$T/\tau={xt:g}$"
        ax.semilogy(Dc, np.exp(lq.lnS_LQ(Dc, G, beta=beta)), c=c, label=lab)
        if xt > 0:
            for D in (6, 10):
                L = lq.simulate(lq.schedule_const(D, xt*P.TAU), kappa=kappa_of(beta))[0]
                ax.semilogy(D, np.exp(L), "o", mfc="none", c=c)
    ax.set_title(fr"$\alpha/\beta={r:g}$ Gy, dosis constante en $[0,T]$")
for ax in axs:
    ax.set_xlabel("Dosis total $D$ (Gy)"); ax.set_ylim(1e-8, 1.5); ax.legend()
save(fig, "fig_supervivencia_alpha_beta.pdf")

# ===================================== Fig. 4: trayectorias U(t) y S(t)
D = 6.0
fig, axs = plt.subplots(1, 2, figsize=(9, 3.6))
for xt, c in [(0.0, "k"), (1.0, "C0"), (4.0, "C1"), (16.0, "C3")]:
    sch = lq.schedule_acute(D) if xt == 0 else lq.schedule_const(D, xt*P.TAU)
    _, tr = lq.simulate(sch, store=True)
    lab = "aguda" if xt == 0 else fr"$T/\tau={xt:g}$"
    axs[0].plot(tr[:, 0], tr[:, 1], c=c, label=lab) # pyright: ignore[reportOptionalSubscript]
    axs[1].semilogy(tr[:, 0], np.exp(tr[:, 2]), c=c, label=lab) # pyright: ignore[reportOptionalSubscript]
axs[0].set_xlim(0, 20); axs[1].set_xlim(0, 20)
axs[0].set_xlabel("Tiempo $t$ (h)"); axs[0].set_ylabel(r"$U$ (DSB/célula)"); axs[0].set_title(r"Daño reparable $U(t)$, $D=6$ Gy")
axs[1].set_xlabel("Tiempo $t$ (h)"); axs[1].set_ylabel(r"$S=N/N_0$"); axs[1].set_title("Supervivencia en el tiempo"); axs[1].legend()
save(fig, "fig_trayectorias_U_S.pdf")

# ===================================== Fig. 5: fraccionamiento
r = np.logspace(-2, 1.5, 300)
fig, ax = plt.subplots()
for n, c in zip([2, 3, 5, 10], ["C0", "C1", "C2", "C3"]):
    ax.semilogx(r, lq.G_fractions(n, r), c=c, label=f"$n={n}$ fracciones")
    for ri in (0.3, 2.0):
        Gn_ = lq.G_numeric(lq.schedule_frac(1.0, n, ri*P.TAU), 1.0)
        ax.plot(ri, Gn_, "o", mfc="none", c=c)
ax.set_xlabel(r"Intervalo entre fracciones $\Delta/\tau$"); ax.set_ylabel(r"$G$ (adimensional)")
ax.set_title(r"$G$ para $n$ fracciones agudas iguales (círculos: numérico)"); ax.legend()
save(fig, "fig_fraccionamiento.pdf")

# ===================================== Fig. 6: mapa de desviación por gamma > 0
Dm = np.linspace(0.5, 12, 24)
gam = np.logspace(-4, -1.5, 22)                       # DSB^-1 h^-1
fig, axs = plt.subplots(1, 2, figsize=(10, 3.9), sharey=True)
for ax, (name, sch, G) in zip(axs, [("Dosis aguda", lq.schedule_acute, 1.0),
                                    (r"Dosis constante $T=4\tau$", lambda D: lq.schedule_const(D, 4*P.TAU), lq.G_constant(4.0)[0])]):
    Z = np.zeros((len(gam), len(Dm)))
    for i, g in enumerate(gam):
        for k, Dk in enumerate(Dm):
            L = lq.simulate(sch(Dk), gamma=g, h=0.05*P.TAU)[0]
            ref = lq.lnS_LQ(Dk, G)
            Z[i, k] = 100*(L - ref)/abs(ref)
    pc = ax.pcolormesh(Dm, gam*P.DELTA*P.TAU, Z, shading="nearest", cmap="viridis")
    ax.set_yscale("log"); ax.set_xlabel("Dosis total $D$ (Gy)"); ax.set_title(name)
    fig.colorbar(pc, ax=ax, label=r"$(\ln S-\ln S_{LQ})/|\ln S_{LQ}|$ (%)")
axs[0].set_ylabel(r"$\gamma\,\tau\,\delta$ (adimensional)")
save(fig, "fig_mapa_gamma.pdf")

# ===================================== Fig. 7: convergencia Euler vs RK4
hs = np.array([0.4, 0.2, 0.1, 0.05, 0.025, 0.0125])*P.TAU
D = 6.0; Gc = lq.G_constant(4.0)[0]; ref = lq.lnS_LQ(D, Gc)
fig, ax = plt.subplots()
for m, c in [("euler", "C2"), ("rk4", "C0")]:
    err = [abs(lq.simulate(lq.schedule_const(D, 4*P.TAU), method=m, h=h)[0] - ref) for h in hs]
    ax.loglog(hs/P.TAU, err, "o-", c=c, label=m.upper())
    print(f"orden observado {m}: {np.polyfit(np.log(hs), np.log(err), 1)[0]:.2f}")
ax.set_xlabel(r"Paso $h/\tau$"); ax.set_ylabel(r"$|\ln S_{num}-\ln S_{LQ}|$"); ax.legend()
ax.set_title(r"Convergencia ($\gamma=0$, $D=6$ Gy, $T=4\tau$)")
save(fig, "fig_convergencia_euler_rk4.pdf")

# ===================================== Resumen en consola (para la tabla del reporte)
print("\nResumen: máximo |Δ ln S| respecto a la fórmula LQ (D en 0.5-10 Gy)")
for name in cases:
    sel = [r_ for r_ in rows if r_[0] == name.replace("$", "").replace("\\", "")]
    a = np.array([[r_[2], r_[3], r_[4], r_[5]] for r_ in sel])
    print(f"  {name:32s} RK4 g=0: {np.max(abs(a[:,1]-a[:,0])):.2e} | Euler g=0: {np.max(abs(a[:,2]-a[:,0])):.2e} | "
          f"g>0: {np.max(abs(a[:,3]-a[:,0])):.2e}  (a 10 Gy: {100*(a[-1,3]-a[-1,0])/abs(a[-1,0]):+.1f} %)")
print("\nS(2 Gy) y S(10 Gy):")
for name, G in [("aguda", 1.0), ("T=4tau", lq.G_constant(4.0)[0]), ("T=16tau", lq.G_constant(16.0)[0])]:
    print(f"  {name:8s} G={G:.4f}  S(2)={np.exp(lq.lnS_LQ(2.0,G)):.4f}  S(10)={np.exp(lq.lnS_LQ(10.0,G)):.3e}")
print("Figuras guardadas en", P.FIG_DIR)
