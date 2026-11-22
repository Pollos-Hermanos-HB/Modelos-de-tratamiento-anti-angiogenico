"""
check_scipy.py -- Referencia independiente con scipy.integrate.solve_ivp (DOP853, tolerancias estrictas)
para validar el RK4 propio de lq.py, con gamma = 0 y gamma > 0.
"""
import numpy as np
from scipy.integrate import solve_ivp
import params as P
import lq

def lnS_scipy(D, T, gamma):
    """Dosis constante D/T en [0,T] (T>0) o aguda (T=0); cola de P.TAIL*tau sin radiación."""
    def rhs(rate):
        return lambda t, y: [P.DELTA*rate - y[0]/P.TAU - gamma*y[0]**2,
                             -(P.ALPHA*rate + P.KAPPA*y[0]**2)]
    opts = dict(method="DOP853", rtol=1e-12, atol=1e-14)
    if T == 0:
        y0 = [P.DELTA*D, -P.ALPHA*D]                      # salto exacto del pulso
    else:
        y0 = solve_ivp(rhs(D/T), (0, T), [0.0, 0.0], **opts).y[:, -1]
    return solve_ivp(rhs(0.0), (0, P.TAIL*P.TAU), y0, **opts).y[1, -1]

if __name__ == "__main__":
    for gamma in (0.0, P.GAMMA_LQ):
        for D, T in [(6.0, 0.0), (6.0, 4.0*P.TAU), (10.0, 16.0*P.TAU)]:
            sch = lq.schedule_acute(D) if T == 0 else lq.schedule_const(D, T)
            a = lq.simulate(sch, gamma=gamma)[0]; b = lnS_scipy(D, T, gamma)
            print(f"gamma={gamma:7.1e} D={D:4.1f} T/tau={T/P.TAU:5.1f}  RK4={a:+.8f}  SciPy={b:+.8f}  |dif|={abs(a-b):.1e}")