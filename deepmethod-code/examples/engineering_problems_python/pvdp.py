import sys
from math import pi

def objective(X: list):
    x1, x2, x3, x4 = X
    return 0.6224 * x1 * x3 * x4 + 1.7781 * x2 * x3**2 + 3.1661 * x1**2 * x4 + 19.84 * x1**2 * x3

def g1_func(X: list):
    return 0.0193 * X[2] - X[0]

def g2_func(X: list):
    return 0.00954 * X[2] - X[1]

def g3_func(X: list):
    return 1296000 - pi * X[2]**2 * X[3] - (4/3) * pi * X[2]**3

def g4_func(X: list):
    return X[3] - 240

def penalty_func(X: list):
    alpha = 1e6
    obj_value = objective(X)
    g1 = g1_func(X)
    g2 = g2_func(X)
    g3 = g3_func(X)
    g4 = g4_func(X)

    penalty = (
        max(0, g1)**2 +
        max(0, g2)**2 +
        max(0, g3)**2 +
        max(0, g4)**2
    )

    return obj_value + alpha * penalty


if __name__ == "__main__":
    X = list()
    for i in sys.argv[1:]:
        X.append(float(i))
    print(penalty_func(X))