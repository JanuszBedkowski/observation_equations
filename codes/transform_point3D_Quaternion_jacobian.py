from sympy import *
from quaternion_R_utils import *

x, y, z = symbols('x y z')
px, py, pz = symbols('px py pz')
q0, q1, q2, q3 = symbols('q0 q1 q2 q3')

point = Matrix([x, y, z, 1]).vec()
position_symbols = [px, py, pz]
quaternion_symbols = [q0, q1, q2, q3]
all_symbols = position_symbols + quaternion_symbols

transformation_matrix = matrix44FromQuaternion(px, py, pz, q0, q1, q2, q3)

transformedPoint = transformation_matrix * point

print(transformedPoint)





