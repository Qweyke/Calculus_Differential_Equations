""" ЧИСЛЕННЫЕ МЕТОДЫ РЕШЕНИЯ ОБЫКНОВЕННЫХ ДИФФЕРЕНЦИАЛЬНЫХ УРАВНЕНИЙ
Решение диф. уров итерационным методом (Эйлер) (самый простой)

Метод Эйлера — это простой численный метод для приближенного решения обыкновенных дифференциальных уравнений (ОДУ)
с начальными условиями, известной как задача Коши. Работает для диф. уров первого порядка,
 если диф. ур выше необходимо понизить его порядок.

Изначально имеются заданные условия, пример x(t) = const, где t - начальный момент времени, к примеру.
Выбираем значение шага h, (с каждым шагом отклонение повышается) на каждой итерации прибавляем
значение шага к (изначальному моменту времени + предыдущие шаги)

x' = f(t, x)
Формула для решения Коши: x_i = x_i-1 + h * f(t_i-1, x_i-1)


Бояршинов МГ, 2006 -- Численные методы.ч4
"""

from math import *
import matplotlib.pyplot as plt

c1, c2 = 1, 1  # начальные условия, рассчитанные на основе x(0) = 1 и x'(0) = 1
t, t_e = 0, 10*pi  # левый и правый концы интервала значений


def equat1(u, v, t_h):  # разбиение диф. ура 2 порядка на систему уравнений первого
    # x" + x = 0
    # u = x
    # v = x'
    du_dt = v
    dv_dt = -u
    return du_dt, dv_dt


def equat1_ac(t_h):
    return sin(t_h) + cos(t_h)  # аналитическое решение


def equat2(u, v, t_h):
    # x" + x = e^t
    # u = x
    # v = x'
    du_dt = v
    dv_dt = e**t_h - u
    return du_dt, dv_dt


def equat2_ac(t_h):
    return (sin(t_h) / 2) + (cos(t_h) / 2) + ((exp(1)**t_h) / 2)  # аналитическое решение


def equat3(u, v, t_h):
    # x" + x = sin(t)
    # u = x
    # v = x'
    du_dt = v
    dv_dt = sin(t_h) - u
    return du_dt, dv_dt


def equat3_ac(t_h):
    return (t_h * sin(t_h) / 4) + cos(t_h) + sin(t_h) - ((t**2 * cos(t_h)) / 4)  # аналитическое решение


def euler(eq, eq_ac, t_end, h, ax, title):

    arr_v = []
    arr_u = []
    arr_ac = []
    arr_t = []

    index = 1

    arr_u.append(1)  # показатели, данные в условии (значение x(0) = 1) u = x
    arr_v.append(1)  # показатели, данные в условии (значение x'(0) = 1) v = x'
    arr_ac.append(eq_ac(0))  # аналит. решение для начального t
    arr_t.append(0)  # начальное t

    while arr_t[index - 1] < t_end:

        arr_t.append(arr_t[index - 1] + h)  # увеличиваем значение t с каждым шагом
        du_dt, dv_dt = eq(arr_u[index - 1], arr_v[index - 1], arr_t[index])  # производные от U "x" и V, для подстановки

        arr_u.append(arr_u[index - 1] + h * du_dt)  # итерационное решение для U
        arr_v.append(arr_v[index - 1] + h * dv_dt)  # итерационное решение для V
        arr_ac.append(eq_ac(arr_t[index]))  # аналитическое решение

        index += 1

    ax.plot(arr_t, arr_v, color='red', ls='--', marker='*', label='V(t)')
    ax.plot(arr_t, arr_u, color='blue', ls='-', marker='^', label='U(t) "x"')
    ax.plot(arr_t, arr_ac, color='green', ls='dotted', marker='8', label='Аналитическое решение')
    ax.set_title(title)
    ax.set_xlabel('Время')
    ax.set_ylabel('Решения Коши')
    ax.legend()
    ax.grid(True)


fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(8, 6))
euler(equat1, equat1_ac, t_e, 0.1, ax1, 'x" + x = 0')
euler(equat2, equat2_ac, t_e, 0.1, ax2, 'x" + x = e^t')
euler(equat3, equat3_ac, t_e, 0.1, ax3, 'x" + x = sin(t)')

plt.tight_layout()
plt.show()
