# Use numbers:
vx, vy, ax, ay, arx, ary, bx, by, brx, bry = map(sympify,
        [1, 2, 3, 4, 5, 6, 7, 8, 9, 10])

theta = atan((ax - vx) / (ay - vy))
A1 = arx ** 2 * sin(theta) ** 2 + ary ** 2 * cos(theta) ** 2
A2 = brx ** 2 * sin(theta) ** 2 + bry ** 2 * cos(theta) ** 2
B1 = 2 * (ary ** 2 - arx ** 2) * sin(theta) * cos(theta)
B2 = 2 * (bry ** 2 - brx ** 2) * sin(theta) * cos(theta)
C1 = arx ** 2 * cos(theta) ** 2 + ary ** 2 * sin(theta) ** 2
C2 = brx ** 2 * cos(theta) ** 2 + bry ** 2 * sin(theta) ** 2
D1 = -2 * A1 * ax - B1 * ay
D2 = -2 * A2 * bx - B2 * by
E1 = -B1 * ax - 2 * C1 * ay
E2 = -B2 * bx - 2 * C2 * by
F1 = A1 * ax ** 2 + B1 * ax * ay + C1 * ay ** 2 - arx ** 2 * ary ** 2
F2 = A2 * bx ** 2 + B2 * bx * by + C2 * by ** 2 - brx ** 2 * bry ** 2
k, b = symbols("k b")
eq1 = Eq((B1 * b + 2 * C1 * k * b + D1 + E1 * k) ** 2 - 4 * (A1 + B1 * k + C1 * k ** 2) * (C1 * b **2 + E1 * b + F1), 0)
eq2 = Eq((B2 * b + 2 * C2 * k * b + D2 + E2 * k) ** 2 - 4 * (A2 + B2 * k + C2 * k ** 2) * (C2 * b **2 + E2 * b + F2), 0)
result = solve((eq1, eq2), (k, b))
