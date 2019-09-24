# main code for labo 3 equipe 1 Lucka Barbeau Matthew Coffey

import numpy as np
from sympy import *
import matplotlib.pyplot as plt

# define variables

rho_0 = 1
u_0 = 800
v_0 = 800
p_0 = 1 * 10 ** 5

rho_x = 0.15
u_x = 50
v_x = -75
p_x = 0.2 * 10 ** 5

rho_y = -0.1
u_y = -30
v_y = 40
p_y = 0.5 * 10 ** 5

rho_ax = 1
u_ax = 1.5
v_ax = 0.5
p_ax = 2

rho_ay = 0.5
u_ay = 0.6
v_ay = 2.0 / 3
p_ay = 1

L = 1
sp_r = 1.4
R = 287
T = 350
x, y = symbols('x y')

init_printing(use_unicode=True)
rho = rho_0 + rho_x * sin(rho_ax * pi * x / L) + rho_y * cos(rho_ay * pi * y / L)
u = u_0 + u_x * sin(u_ax * pi * x / L) + u_y * cos(u_ay * pi * y / L)
v = v_0 + v_x * cos(v_ax * pi * x / L) + v_y * sin(v_ay * pi * y / L)
p = p_0 + p_x * cos(p_ax * pi * x / L) + p_y * sin(p_ay * pi * y / L)
e = R * T / (sp_r - 1) + (u * u + v * v) / 2

sm = diff(rho * u, x) + diff(rho * v, y)
sx = diff(rho * u * u + p, x) + diff(rho * v * u, y)
sy = diff(rho * v * v + p, y) + diff(rho * v * u, x)
sec = diff(rho * u * e + p * u, x) + diff(rho * v * e + p * v, y)

f_rho = lambdify([x, y], rho, 'numpy')
f_u = lambdify([x, y], u, 'numpy')
f_v = lambdify([x, y], v, 'numpy')
f_p = lambdify([x, y], p, 'numpy')
f_sm = lambdify([x, y], sm, 'numpy')
f_sx = lambdify([x, y], sx, 'numpy')
f_sy = lambdify([x, y], sy, 'numpy')
f_sec = lambdify([x, y], sec, 'numpy')

step = 100

sol_rho = np.zeros((step + 1, step + 1))
sol_u = np.zeros((step + 1, step + 1))
sol_v = np.zeros((step + 1, step + 1))
sol_p = np.zeros((step + 1, step + 1))
sol_sm = np.zeros((step + 1, step + 1))
sol_sx = np.zeros((step + 1, step + 1))
sol_sy = np.zeros((step + 1, step + 1))
sol_sec = np.zeros((step + 1, step + 1))

evalX = np.mgrid[0: 1 + 1 / step: 1 / step]
evalY = np.mgrid[0: 1 + 1 / step: 1 / step]

i = 0
j = 0
for i in range(0, step + 1):
    for j in range(0, step + 1):
        sol_rho[i, j] = f_rho(evalX[i], evalY[j])
        sol_u[i, j] = f_u(evalX[i], evalY[j])
        sol_v[i, j] = f_v(evalX[i], evalY[j])
        sol_p[i, j] = f_p(evalX[i], evalY[j])
        sol_sm[i, j] = f_sm(evalX[i], evalY[j])
        sol_sx[i, j] = f_sx(evalX[i], evalY[j])
        sol_sy[i, j] = f_sy(evalX[i], evalY[j])
        sol_sec[i, j] = f_sec(evalX[i], evalY[j])

# print(sol_sec)

fig_rho = plt.figure()
sol_rho_plot = plt.contourf(evalX, evalY, np.transpose(sol_rho), 10, cmap=plt.cm.plasma, origin='lower')
sol_rho_plot_2 = plt.contour(sol_rho_plot, levels=sol_rho_plot.levels[::1], colors='k', origin='lower')
cbar = plt.colorbar(sol_rho_plot)
cbar.set_label('kg/m3')
cbar.add_lines(sol_rho_plot_2)
plt.xlabel("x/l")
plt.ylabel("y/l")
plt.title("Manufactured solution for density for the Euler equations")
plt.show()

fig_u = plt.figure()
sol_u_plot = plt.contourf(evalX, evalY, np.transpose(sol_u), 10, cmap=plt.cm.plasma, origin='lower')
sol_u_plot_2 = plt.contour(sol_u_plot, levels=sol_u_plot.levels[::1], colors='k', origin='lower')
cbar = plt.colorbar(sol_u_plot)
cbar.set_label('m/s')
cbar.add_lines(sol_rho_plot_2)
plt.xlabel("x/l")
plt.ylabel("y/l")
plt.title("Manufactured solution for the x direction velocity for the Euler equations")
plt.show()

fig_v = plt.figure()
sol_v_plot = plt.contourf(evalX, evalY, np.transpose(sol_v), 10, cmap=plt.cm.plasma, origin='lower')
sol_v_plot_2 = plt.contour(sol_v_plot, levels=sol_v_plot.levels[::1], colors='k', origin='lower')
cbar = plt.colorbar(sol_v_plot)
cbar.set_label('m/s')
cbar.add_lines(sol_v_plot_2)
plt.xlabel("x/l")
plt.ylabel("y/l")
plt.title("Manufactured solution for the y direction velocity for the Euler equations")
plt.show()

fig_sec = plt.figure()
sol_sec_plot = plt.contourf(evalX, evalY, np.transpose(sol_sec), 10, cmap=plt.cm.plasma, origin='lower')
sol_sec_plot_2 = plt.contour(sol_sec_plot, levels=sol_sec_plot.levels[::1], colors='k', origin='lower')
cbar = plt.colorbar(sol_sec_plot)
cbar.set_label('sec')
cbar.add_lines(sol_sec_plot_2)
plt.xlabel("x/l")
plt.ylabel("y/l")
plt.title("analytics solution sec")
plt.show()
