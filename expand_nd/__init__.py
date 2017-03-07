import sympy as sy
import numpy as np

def series_nd(f, *args, n=6):
  temp = sy.symbols('temp')
  for a in args:
    f = f.subs(a, a*temp)
  return f.series(temp,n=n).removeO().subs(temp, 1)
  
def coeffMatrix_2d(E, x, y, n=5):
  # n is such that the biggest power on either x or y is n
  # x**n*y**n is the largest possible entry
  out = sy.Matrix(n+1, n+1, [0 for a in range((n+1)**2)])
  e = sy.Poly(E.expand(), x, y, domain='EX')
  for m in range(2*n, -1, -1):
    if m-n >= 0:
      for k in range(m-n, n+1):
        p = x**k*y**(m-k)
        out[k, m-k] = e.coeff_monomial(p)
    else:
      for k in range(0, m+1):
        p = x**k*y**(m-k)
        out[k, m-k] = e.coeff_monomial(p)
  return out

def Cplate(alpha):
  return 1.4-np.cos(2*alpha), 1.2*np.sin(2*alpha)

def GliderAuto(v, theta=0, C=Cplate):
  #v = [vx, vz]
  psi = -np.arctan2(v[1], v[0])
  Yp = [(v[0]**2+v[1]**2)*(C(psi + theta)[1]*np.sin(psi)-C(psi + theta)[0]*np.cos(psi)), (v[0]**2+v[1]**2)*(C(psi + theta)[1]*np.cos(psi)+C(psi + theta)[0]*np.sin(psi))-1]
  return np.array(Yp)
  
def eigenbasis(system, a, b, n=6):
  f1 = series_nd(system[0], a, b, n=2)
  f2 = series_nd(system[1], a, b, n=2)
  A = sy.Matrix([[f1.coeff(a), f1.coeff(b)],[f2.coeff(a), f2.coeff(b)]])
  P, lamda = A.diagonalize()
  p = sy.Matrix([sy.Matrix([P[:, 0]/(-P[0, 0])]).T, sy.Matrix([P[:, 1]/(-P[0, 1])]).T]).T
  # There's a bad assumption here re: P[0, 0] and P[0, 1] nonzero
  X, Y, m, n = sy.symbols('X Y m n')
  G = sy.Matrix([series_ND(system[0], a, b, n=n), series_ND(system[1], a, b, n=n)]).subs(a, -X-Y).subs(b, -m*X-n*Y)
  g1 = sy.cancel(sy.together(sy.Poly(G[0], X, Y), deep=True))
  g2 = sy.cancel(sy.together(sy.Poly(G[1], X, Y), deep=True))
  f = sy.Matrix([[n, -1],[-m, 1]])*sy.Matrix([g1, g2])/(m-n)
  m = -p[1, 0]; n=-p[1, 1]
  return f, m, n

def manifoldExpand(system, a, b, n=6):
  # Function to find b = h(a) from the system [a',b']
  # a'=f1(a, b)
  # b'=f2(a, b)
  H = 0
  HOut = 0
  d = {}
  dOut = {}
  for k in range(1, n+1):
    d['c{0}'.format(k)] = sy.symbols('c{0}'.format(k))
    H += d['c{0}'.format(k)]*a**(k)
  F = sy.Poly((sy.diff(H, a)*system[0].subs(b, H) - system[1].subs(b, h)).simplify().expand(), a, domain='EX')
  for k in range(1, n+1):
    term = a**k
	coeff = F.coeff_monomial(term)
	temp = sy.solve(coeff, d['c{0}'.format(k)])
	dOut['c{0}'.format(k)] = temp[np.argmin(np.abs(np.array(temp)-1))]
    HOut += dOut['c{0}'.format(k)]*term
  return sy.lambdify(HOut)