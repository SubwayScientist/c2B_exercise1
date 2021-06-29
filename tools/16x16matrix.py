import sympy as sy

kx, ky = sy.symbols("k_x k_y")
Kx, Ky = sy.symbols("K_x K_y")
t, tp, tpp, mu = sy.symbols("t tp tpp mu")

# choose the right Rx, Ry
R1x = 0
R1y = 0
R2x = 2
R2y = 0
# Number of sites in cluster
N = 16
# Change K according to Brillouin zone
Kx = [0,sy.pi/2,sy.pi,sy.pi * 3/2 \
     ,0,sy.pi/2,sy.pi,sy.pi * 3/2 \
     ,0,sy.pi/2,sy.pi,sy.pi * 3/2 \
     ,0,sy.pi/2,sy.pi,sy.pi * 3/2]

Ky = [0,0,0,0 \
     ,sy.pi/2,sy.pi/2,sy.pi/2,sy.pi/2 \
     ,sy.pi,sy.pi,sy.pi,sy.pi \
     ,sy.pi * 3/2,sy.pi * 3/2,sy.pi * 3/2,sy.pi * 3/2]

term = 0
for i in range(N):
  exp1 = sy.simplify(sy.exp(1.j * (Kx[i]*(R2x-R1x) + Ky[i]*(R2y-R1y))))
  #exp1 = sy.sympify('exp(I*(Kx*Rx+Ky*Ry))',locals={'Rx':Rx,'Ry':Ry,'Kx':Kx[i],'Ky':Ky[i]}) # we need to evaluate separatly otherwise sympy 

  term += (1 / N * exp1 * sy.exp(1.j * (kx*(R2x-R1x)+ky*(R2y-R1y)) ))  * \
                              ( -2*t  *( sy.cos(kx + Kx[i]) + sy.cos(ky + Ky[i]) ) \
                                -4*tp *( sy.cos(kx + Kx[i]) * sy.cos(ky + Ky[i]) ) \
                                -2*tpp*( sy.cos(2*(kx + Kx[i])) + sy.cos(2*(ky + Ky[i])) ) \
                                - mu )

#print(t)
sy.pprint(sy.simplify(term))
