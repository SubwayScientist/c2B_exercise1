import sympy as sy

kx, ky = sy.symbols("k_x k_y")
Kx, Ky = sy.symbols("K_x K_y")
t, tp, tpp, mu = sy.symbols("t tp tpp mu")

# choose the right Rx, Ry
Rx = -1
Ry = -1


Kx = [0,sy.pi,0,sy.pi]
Ky = [0,0,sy.pi,sy.pi]

term = 0
for i in range(4):
  exp1 = sy.sympify('exp(I*(Kx*Rx+Ky*Ry))',locals={'Rx':Rx,'Ry':Ry,'Kx':Kx[i],'Ky':Ky[i]}) # we need to evaluate separatly otherwise sympy 

  term += (0.25* exp1 * sy.exp(1.j * (kx*Rx+ky*Ry) ))  * \
                              ( -2*t  *( sy.cos(kx + Kx[i]) + sy.cos(ky + Ky[i]) ) \
                                -4*tp *( sy.cos(kx + Kx[i]) * sy.cos(ky + Ky[i]) ) \
                                -2*tpp*( sy.cos(2*(kx + Kx[i])) + sy.cos(2*(ky + Ky[i])) ) \
                                - mu )

#print(t)
sy.pprint(sy.simplify(term))
