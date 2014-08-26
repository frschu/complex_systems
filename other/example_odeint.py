from scipy.integrate import odeint

def f(x, t):
    return x

def df(x, t):
    return 1


x0 = 0.1
dt = 0.1
t_max = 1
n_runs = int(t_max / dt)
ts = [dt * i for i in range(n_runs)]
xs = odeint(f, x0, ts)

