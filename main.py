#main code for labo 3 equipe 1 Lucka Barbeau Matthew Coffey

import numpy as np
from sympy import *




#define varaibles

rho_0=1
u_0=800
v_0=800
p_0=1*10**5

rho_x=0.15
u_x=50
v_x=-75
p_x=0.2*10**5

rho_y=-0.1
u_y=-30
v_y=40
p_y=0.5*10**5

rho_ax=1
u_ax=1.5
v_ax=0.5
p_ax=2

rho_ay=0.5
u_ay=0.6
v_ay=2.0/3
p_ay=1

L=1

x , y = symbols('x y')
init_printing(use_unicode=True)
rho= rho_0+rho_x*np.sin(rho_ax*np.pi()*x/L)+rho_y*np.cos(rho_ay*np.pi()*y/L)
u= u_0+u_x*np.sin(u_ax*np.pi()*x/L)+u_y*np.cos(u_ay*np.pi()*y/L)
v= v_0+v_x*np.sin(v_ax*np.pi()*x/L)+v_y*np.cos(v_ay*np.pi()*y/L)
p= p_0+p_x*np.sin(p_ax*np.pi()*x/L)+p_y*np.cos(p_ay*np.pi()*y/L)

sm=diff(rho*u,x)+diff(rho*v,y)

