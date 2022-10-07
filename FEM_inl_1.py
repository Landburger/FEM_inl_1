import numpy as np
from matplotlib import pyplot as plt
import calfem.core as cfc
import calfem.vis_mpl as cfv
import math

# Materialparametrar
L = [math.sqrt(0.9**2 + 0.1**2),  math.sqrt(3.25), math.sqrt(0.5), math.sqrt(1.25), math.sqrt(2*0.4**2), math.sqrt(1.1**2+0.9**2), math.sqrt(0.5), math.sqrt(0.5), 0.5, 0.5, 0.5, 0.5]
E = 70 * 10 ** 9  # Young's, Pascal
F_ext = 1000  # Belastande kraft, Newton
sigma_y = 70 * 10 ** 6  # Sträckgräns, Pascal
a_y = 0.012  # Yttermått, Meter
a_i = 0.009  # Innermått, Meter
rho = 2800 #kg/m^-3

I = (a_y ** 4 - a_i ** 4) / 12  # Tröghetsmoment,Meter
A = a_y ** 2 - a_i ** 2

# Styvhetsmatris och kraftvektor
K = np.matrix(np.zeros((16, 16)))
f = np.matrix(np.zeros((16, 1)))

# Elementegenskaper
ep = [E, A]

# Konnektivitetsmatrisen och nodkoordinater
Edof = np.array([[1, 2, 3, 4],
                 [1,2,13,14],
                 [1,2,5,6],
                 [1, 2, 9, 10],

                 [3,4,5,6],
                 [3,4,15,16],

                 [5,6,7,8],
                 [5,6,11,12],

                 [7,8,9,10],
                 [9,10,11,12],
                 [11,12,13,14],
                 [13,14,15,16]])
# Elementkoordinater
ex = np.array([[0., .9],
               [0., 1.5],
               [0., .5],
               [0., .5],

               [.9,0.5],
               [.9,2.],

               [0.5,0.],
               [.5,1.],

               [0.,0.5],
               [0.5,1.],
               [1.,1.5],
               [1.5,2.]])

ey = np.array([[1., .9],
               [1., 0.],
               [1., 0.5],
               [1.,0.],

               [.9,0.5],
               [.9,0.],

               [0.5,0.],
               [0.5,0.],

               [0.,0.],
               [0.,0.],
               [0.,0.],
               [0.,0.]])

for elx, ely, eltopo in zip(ex, ey, Edof):
    Ke = cfc.bar2e(elx, ely, ep)
    cfc.assem(eltopo, K, Ke)

bc = np.array([1,2,7,8])
f[15] = -1000

a, r = cfc.solveq(K, f, bc)

print("Displacement: ")
print(a)
#print("Support force: ")
#print(r)

ed = cfc.extractEldisp(Edof,a)

N = np.zeros([Edof.shape[0]])

print("Element forces:")

#Beräkna stångkrafter
i = 0
for elx, ely, eld in zip(ex, ey, ed):
    N[i] = cfc.bar2s(elx, ely, ep, eld)
    print("N%d = %g" % (i+1,N[i]))
    i+=1

#Beräkna säkerhetsfaktorer
j = 0
for value in N:
    if value > 0:
        print("SF Spänning för element " + str(j+1) + " är: ", round(sigma_y / (value / A), 2))
    elif value < 0:
        print("SF knäckning för element " + str(j+1) + " är: ", round((math.pi ** 2 * E * I / (L[j] ** 2 * abs(value))), 2))
    else:
        print("Element " + str(j+1) + " är i kraftjämnvikt")
    j+=1

#Beräkna vikten
total = 0
for length in L:
    total += length
vikt = total*A*rho
print("Fackverket väger: ", vikt, "kg")

cfv.eldraw2(ex,ey,[1,1,1])
cfv.eldisp2(ex,ey,ed,[1,4,1])
cfv.show()