"""
Calculate a trajectory by brute force.

still broken.
"""

interval = 0.00001

t1 = 0.0
x2 = 0.0
while x2 < 2.67 / 2:
    t1 += interval
    x1 = 0 * t1 + 0 * t1 ** 2 / 2 + 1.25 * t1 ** 3 / 6
    v1 = 0 + 0 * t1 + 1.25 * t1 ** 2 / 2
    a1 = 0 + 1.25 * t1
    x2 = x1 + v1 * t1 + a1 * t1 ** 2 / 2 + -1.25 * t1 ** 3 / 6
print t1, x1, x2, v1, a1

x2 = 0.0
t2 = 0.0
while x2 < 2.67 / 2:
    t2 += interval
    x2 = x1 + v1 * t2 + a1 * t2 ** 2 / 2 + -1.25 * t2 ** 3 / 6
    v2 = v1 + a1 * t2 + -1.25 * t2 ** 2 / 2
    a2 = a1 + -1.25 * t2
print t1 + t2, x2, v2, a2

t3 = 0.0
x4 = 0.0
while x4 < 2.67:
    t3 += interval
    x3 = x2 + v2 * t3 + a2 * t3 ** 2 / 2 + -1.25 * t3 ** 3 / 6
    v3 = v2 + a2 * t3 + -1.25 * t3 ** 2 / 2
    a3 = a2 + -1.25 * t3
    x4 = x3 + v3 * t3 + a3 * t3 ** 2 / 2 + 1.25 * t3 ** 3 / 6
print t1 + t2 + t3, x3, x4, v3, a3

t4 = 0.0
x4 = 0.0
while x4 < 2.67:
    t4 += interval
    x4 = x3 + v3 * t4 + a3 * t4 ** 2 / 2 + 1.25 * t4 ** 3 / 6
    v4 = v3 + a3 * t4 + 1.25 * t4 ** 2 / 2
    a4 = a3 + 1.25 * t4
print t1 + t2 + t3 + t4, x4, v4, a4

print 0.0, 0.0, 0.0, 1.25
print t1, x1, v1, a1, -1.25
print t1 + t2, x2, v2, a2, -1.25
print t1 + t2 + t3, x3, v3, a3, 1.25
print t1 + t2 + t3 + t4, x4, v4, a4, 0.0
