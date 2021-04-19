from sympy import *
from rodrigues_R_utils import *
from rodrigues_Rutheta_utils import *

x, y, z = symbols('x y z')
px, py, pz = symbols('px py pz')
sx, sy, sz = symbols('sx sy sz')
ux, uy, uz, theta = symbols('ux uy uz theta')

point = Matrix([x, y, z, 1]).vec()
position_symbols = [px, py, pz]
rodriguez_symbols = [sx, sy, sz]
rodriguez_utheta_symbols = [ux, uy, uz, theta]
all_symbols = position_symbols + rodriguez_symbols

transformation_matrix = matrix44FromRodrigues(px, py, pz, sx, sy, sz)
transformedPoint = transformation_matrix * point

all_symbols_utheta = position_symbols + rodriguez_utheta_symbols
transformation_matrix_utheta = matrix44FromRodrigues_utheta(px, py, pz, ux, uy, uz, theta)

print(transformation_matrix)
print(transformedPoint)
print(latex(transformation_matrix_utheta))
