import numpy as np
import numpy.random as rd
DIM = 32

def spin(v):
    return lat[v[0]][v[1]]

def flip(v):
    lat[v[0]][v[1]] *= -1

def neigh(v):
    i, j = v[0], v[1]
    n = (i + 1) % DIM;
    e = (j + 1) % DIM;
    w = (j - 1 + DIM) % DIM;
    s = (i - 1 + DIM) % DIM;
    return [(n, j), (i, e), (i, w), (s, j)];

def elat():
    e = 0
    for i in range(DIM):
        for j in range(DIM):
            e -= lat[i][j] * (lat[i][(j + 1) % DIM] + lat[(i + 1) % DIM][j]);
    return e

def wolff(steps, P):
    mag = 0
    for i in range(steps):
        A = [(rd.randint(0, DIM - 1), rd.randint(0, DIM - 1))]
        C = []
        m = np.sum(lat)
        for v in A:
            C.append(v)
            for w in neigh(v):
                if spin(v) == spin(w):
                    if rd.random() < P and w not in A:
                        A.append(w)
            if A == C:
                break

        for v in A:
            flip(v)
            m += 2 * spin(v)
        mag += float(m) / (DIM * DIM)
        print i, elat()
    return mag

steps = int(1e3)
for T in np.arange(2, 3, 0.05):
    lat = rd.choice([-1, 1], (DIM, DIM))
    #lat = np.ones((DIM, DIM))
    P = 1 - np.exp(-1.0 / 2.6)
    wolff(steps, P)
    break
