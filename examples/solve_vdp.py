
import numpy as np

from scipy.integrate import solve_ivp

import matplotlib.pyplot as plt

MU = 10.0

def f(t,y):

    q = np.empty_like(y)

    q[0] = y[1]
    q[1] = MU*(1.0 - y[0]**2)*y[1] - y[0]

    return q


def main():

    y0 = [1.0,1.0]

    sol = solve_ivp(f,[0.0,100.0], y0, method='BDF')

    print(sol.t.shape)
    print(sol.y.shape)

    res = np.vstack((sol.t, sol.y)).T
    print(res.shape)

    np.savetxt("vdp_python.txt",res)

    plt.plot(sol.t,sol.y.T)
    plt.show()

if __name__ == '__main__':
    main()