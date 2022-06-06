from sqlite3 import dbapi2
import sympy as sy
import numpy as np

def ItosLemma(mu, sigma, X, h):
    """
    This function applies Ito's lemma given a twice-differentiable,
    invertible function h: Y = h(X) and drift (mu, vector) and diffusion (sigma, matrix) for a
    stochastic process (X) satisfying SDE
      dX = mu(X,t) dt + sigma(X,t) dB

    This routine returns the new drift and diffusion
      dY = new_mu(X,t) dt + new_sigma(X,t) dB
    
    See Bernt Øksendal, Stochastic Differential Equations, section 4.2 for details.
    """

    # drift
    mu0 = sy.diff(h, t)
    # print(mu0)

    dh_dX = h.jacobian(X)
    mu1 = dh_dX * mu
    # print(sy.simplify(mu1))

    mu2 = []
    for k in range(len(h)):
      hk = h[k]
      mu2_k = 0
      for i in range(len(X)):
        for j in range(len(X)):
          dhk_dXidXj = sy.diff(hk, X[i], X[j])
          sum_s = 0
          for m in range(len(sigma[0, :])):
            sum_s += sigma[i, m]*sigma[j, m]
          mu2_k += dhk_dXidXj*sum_s
      mu2.append(mu2_k/2)
    mu2 = np.expand_dims(mu2, axis=1)

    # print(mu2)

    new_mu = mu0 + mu1 + mu2

    # diffusion
    dh_dX = h.jacobian(X)
    new_sigma = dh_dX * sigma

    return sy.simplify(new_mu), sy.simplify(new_sigma)

# time
t = sy.symbols('t')

# X:
r, theta, vr, omega, sigma_r, sigma_theta, mu_grav = sy.symbols('r theta vr omega sigma_r sigma_theta mu_grav')
X = sy.Matrix([r, theta, vr, omega])

# Y = h(X)
Y = sy.Matrix([r * sy.cos(theta), r * sy.sin(theta), vr * sy.cos(theta) - r * omega * sy.sin(theta),
            vr * sy.sin(theta) + r * omega * sy.cos(theta)])

# Stochastic Two-body problem:
mu = sy.Matrix([vr, omega, r*omega**2 - mu_grav/(r**2), -2*vr*omega/r])
sigma = sy.Matrix([[0, 0], [0, 0], [sigma_r, 0], [0, sigma_theta]])

new_mu, new_sigma = ItosLemma(mu, sigma, X, Y)
print('Here are drift and diffusion for two-body problem in cartesian coordinates:')
print(new_mu)
print(new_sigma)
print(' ')

V = sy.Matrix([vr, r*omega])
V_2 = sy.Matrix([V.dot(V)])

H = V_2/2 - sy.Matrix([mu_grav/r])

H_mu, H_sigma = ItosLemma(mu, sigma, X, H)

print('Here is the drift and diffusion for the total energy of the two-body problem:')
print(H_mu)
print(H_sigma)
print(' ')

