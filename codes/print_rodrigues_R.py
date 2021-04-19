from sympy import *
from rodrigues_R_utils import *

T_x, T_y, T_z = symbols('T_x T_y T_z')
s_x, s_y, s_z = symbols('s_x s_y s_z')

RT_wc = matrix44FromRodrigues(T_x, T_y, T_z, s_x, s_y, s_z)

print(RT_wc)
print(latex(RT_wc))
