import numpy as np
import matplotlib.pyplot as plt


def f1(A, B, C, q, r, w0, a32, a33):  # dp/dt
    return (C - B) * (3 * w0 ** 2 * a32 * a33 - q * r) / A


def f2(A, B, C, r, p, w0, a31, a33):  # dq/dt
    return (A - C) * (3 * w0 ** 2 * a31 * a33 - p * r) / B


def f3(A, B, C, p, q, w0, a31, a32):  # dr/dt
    return (B - A) * (3 * w0 ** 2 * a32 * a31 - q * p) / C


def f4(q, r, w0, a22, a23, gama):  # d(alpha)/dt
    return (q * np.cos(gama) - r * np.sin(gama)) / (a22 * np.cos(gama) - a23 * np.sin(gama)) - w0


def f5(p, q, r, a21, a22, a23, gama):  # d(gama)/dt
    return p - a21 * (q * np.cos(gama) - r * np.sin(gama)) / (a22 * np.cos(gama) - a23 * np.sin(gama))


def f6(q, r, a22, a23, gama):  # d(beta)/dt
    return (q * a23 - r * a22) / (a23 * np.sin(gama) - a22 * np.cos(gama))


[alpha, beta, gama, A, B, C, p, q, r] = [float(i) for i in input().split()]


def A21(beta):
    return np.sin(beta)


def A22(beta, gama):
    return np.cos(beta) * np.cos(gama)


def A23(beta, gama):
    return -np.cos(beta) * np.sin(gama)


def A31(alpha, beta):
    return -np.sin(alpha) * np.cos(beta)


def A32(alpha, beta, gama):
    return np.cos(alpha) * np.sin(gama) + np.sin(alpha) * np.sin(beta) * np.cos(gama)


def A33(alpha, beta, gama):
    return np.cos(alpha) * np.cos(gama) - np.sin(alpha) * np.sin(beta) * np.sin(gama)


w0 = np.sqrt(9.81 / 6400000)

t = 0
h = 1
q = 1.1*w0



Alpha = np.array([])
Alpha_v = np.array([])

Beta = np.array([])
Beta_v = np.array([])

Gama = np.array([])
Gama_v = np.array([])


while (t < 100000):
    k11 = f1(A, B, C, q, r, w0, A32(alpha, beta, gama), A33(alpha, beta, gama))
    k12 = f2(A, B, C, r, p, w0, A31(alpha, beta), A33(alpha, beta, gama))
    k13 = f3(A, B, C, p, q, w0, A31(alpha, beta), A32(alpha, beta, gama))
    k14 = f4(q, r, w0, A22(beta, gama), A23(beta, gama), gama)
    k15 = f5(p, q, r, A21(beta), A22(beta, gama), A23(beta, gama), gama)
    k16 = f6(q, r, A22(beta, gama), A23(beta, gama), gama)

    Alpha = np.append(Alpha, alpha)
    Alpha_v = np.append(Alpha_v, k14)

    Beta = np.append(Beta, beta)
    Beta_v = np.append(Beta_v, k16)

    Gama = np.append(Gama, gama)
    Gama_v = np.append(Gama_v, k15)

    k21 = f1(A, B, C, q + k12 * h / 2, r + k13 * h / 2, w0,
             A32(alpha + k14 * h / 2, beta + k16 * h / 2, gama + k15 * h / 2),
             A33(alpha + k14 * h / 2, beta + k16 * h / 2, gama + k15 * h / 2))
    k22 = f2(A, B, C, r, p, w0, A31(alpha + k14 * h / 2, beta + k16 * h / 2),
             A33(alpha + k14 * h / 2, beta + k16 * h / 2, gama + k15 * h / 2))
    k23 = f3(A, B, C, p, q, w0, A31(alpha + k14 * h / 2, beta + k16 * h / 2),
             A32(alpha + k14 * h / 2, beta + k16 * h / 2, gama + k15 * h / 2))
    k24 = f4(q, r, w0, A22(beta + k16 * h / 2, gama + k15 * h / 2), A23(beta + k16 * h / 2, gama + k15 * h / 2), gama)
    k25 = f5(p, q, r, A21(beta + k16 * h / 2), A22(beta + k16 * h / 2, gama + k15 * h / 2),
             A23(beta + k16 * h / 2, gama + k15 * h / 2), gama)
    k26 = f6(q, r, A22(beta + k16 * h / 2, gama + k15 * h / 2), A23(beta + k16 * h / 2, gama + k15 * h / 2), gama)

    k31 = f1(A, B, C, q + k22 * h / 2, r + k23 * h / 2, w0,
             A32(alpha + k24 * h / 2, beta + k26 * h / 2, gama + k25 * h / 2),
             A33(alpha + k24 * h / 2, beta + k26 * h / 2, gama + k25 * h / 2))
    k32 = f2(A, B, C, r, p, w0, A31(alpha + k24 * h / 2, beta + k26 * h / 2),
             A33(alpha + k24 * h / 2, beta + k26 * h / 2, gama + k25 * h / 2))
    k33 = f3(A, B, C, p, q, w0, A31(alpha + k24 * h / 2, beta + k26 * h / 2),
             A32(alpha + k24 * h / 2, beta + k26 * h / 2, gama + k25 * h / 2))
    k34 = f4(q, r, w0, A22(beta + k26 * h / 2, gama + k25 * h / 2), A23(beta + k26 * h / 2, gama + k25 * h / 2), gama)
    k35 = f5(p, q, r, A21(beta + k26 * h / 2), A22(beta + k26 * h / 2, gama + k25 * h / 2),
             A23(beta + k26 * h / 2, gama + k25 * h / 2), gama)
    k36 = f6(q, r, A22(beta + k26 * h / 2, gama + k25 * h / 2), A23(beta + k26 * h / 2, gama + k25 * h / 2), gama)

    k41 = f1(A, B, C, q + k32 * h, r + k33 * h, w0,
             A32(alpha + k34 * h, beta + k36 * h, gama + k35 * h),
             A33(alpha + k34 * h, beta + k36 * h, gama + k35 * h))
    k42 = f2(A, B, C, r, p, w0, A31(alpha + k34 * h, beta + k36 * h),
             A33(alpha + k34 * h, beta + k36 * h, gama + k35 * h))
    k43 = f3(A, B, C, p, q, w0, A31(alpha + k34 * h, beta + k36 * h),
             A32(alpha + k34 * h, beta + k36 * h, gama + k35 * h))
    k44 = f4(q, r, w0, A22(beta + k36 * h, gama + k35 * h), A23(beta + k36 * h, gama + k35 * h), gama)
    k45 = f5(p, q, r, A21(beta + k36 * h), A22(beta + k36 * h, gama + k35 * h),
             A23(beta + k36 * h, gama + k35 * h), gama)
    k46 = f6(q, r, A22(beta + k36 * h, gama + k35 * h), A23(beta + k36 * h, gama + k35 * h), gama)

    p += h * (k11 + 2 * k21 + 2 * k31 + k41) / 6
    q += h * (k12 + 2 * k22 + 2 * k32 + k42) / 6
    r += h * (k13 + 2 * k23 + 2 * k33 + k43) / 6
    alpha += h * (k14 + 2 * k24 + 2 * k34 + k44) / 6
    gama += h * (k15 + 2 * k25 + 2 * k35 + k45) / 6
    beta += h * (k16 + 2 * k26 + 2 * k36 + k46) / 6

    t += h
    #print(alpha)


fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(16, 10), dpi=150)




ax1.scatter(Alpha, Alpha_v, s=2)
ax2.scatter(Beta, Beta_v, s=2)
ax3.scatter(Gama, Gama_v, s=2)

ax1.grid(True)
ax2.grid(True)
ax3.grid(True)

ax1.set_title('Alpha')
ax2.set_title("Beta")
ax3.set_title("Gama")


plt.savefig("Eq.png")
plt.show()
#0.1 0 0 20 40 10 0 0 0