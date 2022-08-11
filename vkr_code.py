import math
from prettytable import PrettyTable


def calculations(x1, a_L1, mode='old') -> float:
    # Пустые массивы
    t_T = []  # массив времени
    x = []  # перемещение
    x_zvezd = []  # первая производная по перемещению
    x_zvezd_2 = []  # вторая производная по перемещению
    beta = []  # угол поворота
    beta_zvezd = []  # первая производная по углу поворота
    beta_zvezd_2 = []  # вторая производная по углу поворота
    X = []  # гидродинамическая составляющая силы
    Xi = []  # инерционная составляющая силы
    P = []  # упор
    dA = []  # работа

    # Исходные данные
    J = x1[0]  # поступь
    k = 2 / 0.93  # коэффициент
    v0 = 1  # отнормированная скорость набегающего потока
    e_L = 0.1  # относительная толщина профиля
    a_L = a_L1  # относительный разбег
    koefficent = x1[1]
    f_L = 0  # кривизна
    dt_T = 0.01  # относительный шаг по времени
    Re = 1e6  # Число Рейнольдса
    Sh = math.pi / J  # число Струхаля
    mu = (0.77 * e_L) / (1 - 0.6 * e_L)  # функция относительной толщины профиля
    A_L = 1 / (2 * (1 + mu ** 2))
    alpha = math.atan(2 * f_L * (1 + mu ** 2))
    R = (1 + mu) / math.cos(alpha)
    r = R / (1 + 2 * mu)

    # Расчет присоединенных масс
    mu_11 = (math.pi / 4) * (A_L ** 2) * (r ** 2 + R ** 2 - 2 * math.cos(2 * alpha))
    mu_22 = (math.pi / 4) * (A_L ** 2) * (r ** 2 + R ** 2 + 2 * math.cos(2 * alpha))
    mu_12 = (math.pi / 4) * (A_L ** 2) * math.sin(alpha)
    mu_16 = (math.pi / 8) * (A_L ** 3) * (r ** 2 + R ** 2 + 4 * (r + R) * math.cos(alpha)) * math.sin(alpha)
    mu_26 = (math.pi / 8) * (A_L ** 3) * (
            r ** 3 + R ** 3 + (r ** 2 + R ** 2) * math.cos(alpha) + 2 * (r + R) * math.cos(2 * alpha))

    # Заполнение массива по времени
    for i in range(0, 105):
        t_T.append(float("{0:.2f}".format(i * dt_T)))

    # Заполнение массива по перемещению
    n = 1.5
    for t in t_T:
        if mode == 'old':
            c = 0
            b = 1
        elif mode == 'new':
            c = x1[2]
            b = x1[3]
        x.append((a_L / 2) * (1 - ((c + math.cos(2 * math.pi * t + ((2 * c) / n) + n)) * math.cos(n) + (
                b * math.sin(2 * math.pi * t + ((2 * c) / n) + n)) * math.sin(n)) / (math.sqrt(
            (c + math.cos(2 * math.pi * t + ((2 * c) / n) + n)) ** 2 + (
                    b * math.sin(2 * math.pi * t + ((2 * c) / n) + n)) ** 2))))

    # Заполнение массива первых производных по перемещению
    for j in range(len(x) - 1):
        if j == 0:
            x_zvezd.append(0.0)
        elif j != 0:
            x_zvezd.append((x[j + 1] - x[j - 1]) / (2 * dt_T))

    # Заполнение массива вторых производных по перемещению
    for j in range(len(x_zvezd) - 1):
        if j == 0:
            x_zvezd_2.append(0.0)
        elif j != 0:
            x_zvezd_2.append((x_zvezd[j + 1] - x_zvezd[j - 1]) / (2 * dt_T))

    # Заполнение массива углов поворота
    for i in range(len(t_T) - 1):
        beta.append(0.5 * ((math.pi / 2) + math.atan(v0 - koefficent * x_zvezd[i])))

    # Заполнение массива первых производных по углу поворота
    for j in range(len(beta) - 1):
        if j == 0:
            beta_zvezd.append(0.0)
        elif j != 0:
            beta_zvezd.append((beta[j + 1] - beta[j - 1]) / (2 * dt_T))

    # Заполнение массива вторых производных по углу поворота
    for j in range(len(beta_zvezd) - 1):
        if j == 0:
            beta_zvezd_2.append(0.0)
        elif j != 0:
            beta_zvezd_2.append((beta_zvezd[j + 1] - beta_zvezd[j - 1]) / (2 * dt_T))

    # Расчет гидродинамической составляющей силы
    for i in range(len(x_zvezd)):
        # Коэффициент нормальной силы на пластине (Формула Релея)
        C_n = 2 * math.pi * ((math.sin(abs(beta[i]))) / (4 + math.pi * (math.sin(abs(beta[i]))))) * (
                1 - 0.5 + (1.6 / (abs(beta[i]) * (6 / math.pi)) ** 0.7))
        X.append(k * (-(1 + Sh * x_zvezd[i]) * abs(1 + Sh * x_zvezd[i]) * (C_n * math.sin(beta[i]) + (
                ((2 + 2.4 * e_L + 17 * e_L ** 3) * 0.455) / ((math.log10(Re * abs(1 - x_zvezd[i]))) ** 2.58)))))

    # Расчет инерционной составляющей силы
    for i in range(len(beta_zvezd_2)):
        Xi.append(-1 * (-2 * Sh ** 2 * (
                x_zvezd_2[i] * (mu_12 - mu_11 * (math.cos(beta[i])) ** 2 - mu_22 * (math.sin(beta[i])) ** 2) + x_zvezd[
            i] * beta_zvezd[i] * (math.sin(2 * beta[i]) * (mu_11 - mu_22) + 2 * math.cos(2 * beta[i]) * mu_12) -
                beta_zvezd_2[i] * (math.cos(beta[i]) * mu_16 - math.sin(beta[i]) * mu_26) + beta_zvezd[i] ** 2 * (
                        math.sin(beta[i]) * mu_16 + math.cos(beta[i]) * mu_26))))
    Xi[0] = Xi[-2]
    Xi[1] = Xi[-1]

    # Расчет упора
    for i in range(len(Xi)):
        P.append(Xi[i] + X[i])

    # Расчет КПД
    for i in range(len(P)):
        dA.append(x_zvezd[i] * dt_T * P[i])  # работа

    A = -sum(dA)  # суммарная работа
    P_average = sum(P) * dt_T  # средний по периоду упор
    kpd = P_average / (Sh * A)  # КПД

    return kpd


machineAcc = 0.000000001


# Исследующий поиск
def utilSearch(b, h, f, a_L, mode):
    bres = b[:]
    fb = f(bres, a_L, mode)
    for i in range(0, len(bres)):
        bn = bres
        bn[i] = bn[i] + h[i]
        fc = f(bn, a_L, mode)
        if (fc + machineAcc >= fb):
            bres = bn
            fb = fc
        else:
            bn[i] = bn[i] - 2 * h[i]
            fc = f(bn, a_L, mode)
            if (fc + machineAcc >= fb):
                bres = bn
                fb = fc
    print(f"bres = {bres}, h = {h} KPD = {calculations(a_L1=a_L, x1=bres, mode='new')}")
    return bres


Path1 = []
Path2 = []
Path3 = []
Path4 = []


# Метод конфигураций Хука-Дживса
def HJ(b1, h, e, f, a_L, mode):
    z = 0.01
    runOuterLoop = True
    while (runOuterLoop):
        runOuterLoop = False
        runInnerLoop = True
        xk = b1  # step1
        b2 = utilSearch(b1, h, f, a_L, mode)  # step2
        Path1.append(b1)
        Path2.append(b2)
        Path3.append(xk)
        while (runInnerLoop):
            Path1.append(b1)
            Path2.append(b2)
            runInnerLoop = False
            for i in range(len(b1)):  # step3
                xk[i] = b1[i] + 2 * (b2[i] - b1[i])
            Path3.append(xk)
            x = utilSearch(xk, h, f, a_L, mode)  # step4
            Path4.append(x)
            b1 = b2  # step5

            fx = f(x, a_L, mode)
            fb1 = f(b1, a_L, mode)
            if (fx + machineAcc < fb1):  # step6
                b2 = x
                runInnerLoop = True  # to step3
            elif (fx - machineAcc > fb1):  # step7
                runOuterLoop = True  # to step1
                break
            else:
                s = 0
                for i in range(len(h)):
                    s += h[i] * h[i]
                if (e * e + machineAcc > s):  # step8
                    break  # to step10
                else:
                    for i in range(len(h)):  # step9
                        h[i] = h[i] * z
                    runOuterLoop = True  # to step1
    return b1  # step10


def main():
    table = PrettyTable()
    table.field_names = ['№', 'a/L', 'J', 'Ja', 'koef', 'c', 'b', 'KPD']
    a_L_list = [2, 3, 4, 5, 10, 20, 50, 100]  # [10, 20, 50, 100]
    h_list = [[0.5, 0.2, 0.01, 0.05], [0.1, 0.2, 0.01, 0.1], [0.5, 0.25, 0.01, 0.1], [0.1, 0.4, 0.01, 0.1],
              [0.5, 1, 0.01, 0],
              [1, 1, 0.01, 0.1], [0.5, 1, 0.01, 0.1], [0.5, 1, 0.01, 0.1]]

    b_list = [[3, 1, 0.02, 0.85], [8, 1, 0.01, 0.8], [14, 1, 0.02, 0.8], [20, 1, 0.015, 0.85],
              [41.3, 5.252, 0.02, 0.81],
              [65, 20, 0, 0.81], [100, 100, 0, 0.81], [400, 100, 0, 0.81]]
    for i in range(len(a_L_list)):
        res = HJ(f=calculations, h=h_list[i], b1=b_list[i], e=1e-3, a_L=a_L_list[i], mode='new')
        table.add_row([i + 1, a_L_list[i], res[0], res[0] / (math.pi * a_L_list[i]), res[1], res[2], res[3],
                       calculations(res, a_L_list[i], mode='new')])
        print(f"[INFO] Calculated: {i + 1}/{len(a_L_list)}")
        print(
            f"a/L = {a_L_list[i]}: J = {res[0]}, Ja = {res[0] / (math.pi * a_L_list[i])}, koef = {res[1]}, c = {res[2]}, b = {res[3]} KPD = {calculations(res, a_L_list[i], mode='new')}")
    print(table)

if __name__ == '__main__':
    main()
