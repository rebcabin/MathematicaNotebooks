import numpy as np
from typing import List, Tuple


def arow(m: int, x: float) -> List[float]:
    return [x ** i for i in range(m)]


def ground_truth(x: float) -> float:
    return np.sin(2 * np.pi * x)


Z_SIGMA = 1.0
X_MIN = 0.0
X_MAX = 1.0


def zfake(x: float) -> float:
    return np.random.normal(ground_truth(x), Z_SIGMA)


import pprint
pp = pprint.PrettyPrinter(indent=2)


def my_reduce(bin_fun, seq, a_priori):
    result = a_priori
    for s in seq:
        result = bin_fun(result, s)
    return result


def mk_data(n: int, m: int):
    xs = np.linspace(0, 1.0, n)
    yf = [zfake(x) for x in xs]
    a = [arow(m, x) for x in xs]
    zs = yf

    return zs, a

def mc_normal(n: int, m: int):
    zs, a = mk_data(n, m)

    at = np.transpose(a)
    ata = np.dot(at, a)
    atzs = np.dot(at, zs)
    atainv = np.linalg.inv(ata)

    result = np.dot(atainv, atzs)
    return result


def mc_rls(n: int, m: int):
    def rls(xiLambda, zArow):
        xi, Lambda = xiLambda
        z, aVec = zArow

        a = np.array([aVec])
        at = np.transpose(a)
        Pi = Lambda + np.dot(at, a)

        predxi = np.dot(at, z)
        attenxi = np.dot(Lambda, xi)

        newXi = np.linalg.solve(Pi, predxi + attenxi)
        # 60.seconds
        # newXi = np.dot(np.linalg.inv(Pi), predxi + attenxi)
        # 75.seconds
        # linalg.inv is eliminated for this benchmark
        newLambda = Pi

        return (newXi, newLambda)

    zs, a = mk_data(n, m)

    seq = [(zs[i], a[i]) for i in range(len(zs))]
    a_priori = [[0.0] * m, np.diag([1e-6] * m)]

    result = my_reduce(rls, seq, a_priori)
    return result
