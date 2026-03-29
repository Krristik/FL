#!/usr/bin/env python3
import sys
from math import sin, sqrt


def schwefel(xx: list):
    sum = 0
    for x in xx:
        sum += x * sin(sqrt(abs(x)))

    return 418.9829 * len(xx) - sum


xx = list()
for i in sys.argv[1:]:
    xx.append(float(i))

print(schwefel(xx))
