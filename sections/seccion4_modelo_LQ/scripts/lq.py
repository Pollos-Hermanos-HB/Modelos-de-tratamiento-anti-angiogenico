"""
lq.py -- Función G de Lea-Catcheside e integración de (8)-(9).

Convención (SUPUESTA; verificar contra las ecs. (8)-(9) del paper):
  (8)  dU/dt = delta*Ddot(t) - U/tau - gamma*U^2
  (9)  dN/dt = -N*( alpha*Ddot(t) + kappa*U^2 )
Se integra L = ln(N/N0) en lugar de N (evita underflow): dL/dt = -(alpha*Ddot + kappa*U^2).
"""
import numpy as np
import params as P

# ------------------------------------------------------------------ G analítica
def G_acute():
    """Dosis aguda (delta de Dirac): G = 1."""
    return 1.0

def G_constant(x):
    """Dosis constante en [0,T], x = T/tau:  G = 2(x - 1 + e^{-x})/x^2."""
    x = np.atleast_1d(np.asarray(x, dtype=float))
    xs = np.where(x < 1e-4, 1.0, x)                       # evita 0/0
    exact = 2.0*(xs + np.expm1(-xs))/xs**2
    series = 1.0 - x/3.0 + x**2/12.0                      # desarrollo en x -> 0
    return np.where(x < 1e-4, series, exact)

def G_exponential(y):
    """Tasa Ddot=(D/theta)e^{-t/theta}, t>=0, y = theta/tau:  G = 1/(1+y)."""
    return 1.0/(1.0 + np.asarray(y, dtype=float))

def G_fractions(n, dt_over_tau):
    """n fracciones agudas iguales separadas Delta: G = 1/n + (2/n^2) sum_{j=1}^{n-1}(n-j)e^{-j Delta/tau}."""
    r = np.atleast_1d(np.asarray(dt_over_tau, dtype=float))
    out = np.full_like(r, 1.0/n)
    for j in range(1, n):
        out += (2.0/n**2)*(n - j)*np.exp(-j*r)
    return out

def G_from_profile(rate, t_end, tau=P.TAU, m=20000):
    """Verificación independiente: G = (1/D^2) ∬ Ddot(t1)Ddot(t2) e^{-|t1-t2|/tau} dt1 dt2
    (regla del punto medio, recursión O(m)). rate(t) en Gy/h."""
    h = t_end/m
    t = (np.arange(m) + 0.5)*h
    r = rate(t)*h                      # dosis por celda
    D = r.sum()
    e = np.exp(-h/tau)
    S, acc = 0.0, 0.0
    for ri in r:                       # S_i = sum_{j<i} r_j e^{-(i-j)h/tau}
        acc += ri*S
        S = e*(S + ri)
    return (2.0*acc + np.sum(r**2))/D**2

# ------------------------------------------------------------------ cronogramas
def schedule_acute(D):                 return [("imp", D)]
def schedule_const(D, T):              return [("rate", lambda s: D/T, T)]
def schedule_frac(D, n, Delta):
    s = []
    for i in range(n):
        s.append(("imp", D/n))
        if i < n-1: s.append(("wait", Delta))
    return s
def schedule_exp(D, theta, n_theta=30.0):
    return [("rate", lambda s: (D/theta)*np.exp(-s/theta), n_theta*theta)]

# ------------------------------------------------------------------ integración
def simulate(schedule, gamma=0.0, method="rk4", h=None, tail=P.TAIL*P.TAU,
             alpha=P.ALPHA, delta=P.DELTA, tau=P.TAU, kappa=P.KAPPA, store=False):
    """
    Integra (U, L=ln N/N0) con U(0)=0, L(0)=0 sobre un cronograma de eventos:
      ("imp", D)        pulso agudo: U += delta*D, L -= alpha*D (salto exacto)
      ("rate", f, T)    tasa Ddot=f(s) durante T horas (s = tiempo local)
      ("wait", T)       sin radiación
    Se añade una cola sin radiación de duración `tail`. method: "rk4" | "euler".
    Devuelve (lnS_final, trayectoria|None); trayectoria = columnas (t, U, ln N/N0).
    """
    h = P.H_REL*tau if h is None else h
    U, L, t = 0.0, 0.0, 0.0
    traj = [(0.0, 0.0, 0.0)] if store else None

    def f(s, U, L, rate):
        r = rate(s)
        return delta*r - U/tau - gamma*U*U, -(alpha*r + kappa*U*U)

    zero = lambda s: 0.0
    events = list(schedule) + [("wait", tail)]
    for ev in events:
        if ev[0] == "imp":
            U += delta*ev[1]; L -= alpha*ev[1]
            if store: traj.append((t, U, L))
            continue
        if ev[0] == "rate":
            rate, T = ev[1], ev[2]
        else:
            rate, T = zero, ev[1]
        n = max(1, int(np.ceil(T/h))); hh = T/n; s = 0.0
        for _ in range(n):
            if method == "euler":
                dU, dL = f(s, U, L, rate)
                U, L = U + hh*dU, L + hh*dL
            else:
                a1, b1 = f(s, U, L, rate)
                a2, b2 = f(s+hh/2, U+hh/2*a1, L+hh/2*b1, rate)
                a3, b3 = f(s+hh/2, U+hh/2*a2, L+hh/2*b2, rate)
                a4, b4 = f(s+hh,   U+hh*a3,   L+hh*b3,   rate)
                U += hh/6*(a1+2*a2+2*a3+a4); L += hh/6*(b1+2*b2+2*b3+b4)
            s += hh; t += hh
            if store: traj.append((t, U, L))
    return L, (np.array(traj) if store else None)

def lnS_LQ(D, G, alpha=P.ALPHA, beta=P.BETA):
    """ln(N/N0) = -alpha D - beta G D^2."""
    return -alpha*np.asarray(D) - beta*G*np.asarray(D)**2

def G_numeric(schedule, D, **kw):
    """G estimada desde (8)-(9) con gamma=0:  G = -(L_final + alpha D)/(beta D^2)."""
    L, _ = simulate(schedule, gamma=0.0, **kw)
    return -(L + P.ALPHA*D)/(P.BETA*D**2)