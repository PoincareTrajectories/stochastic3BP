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

    dh_dX = h.jacobian(X)
    mu1 = dh_dX * mu

    mu2 = []
    for k in range(len(X)):
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

    new_mu = mu0 + mu1 + mu2

    # diffusion
    dh_dX = h.jacobian(X)
    new_sigma = dh_dX * sigma

    return new_mu, new_sigma

# time
t = sy.symbols('t')

# X:
r, theta, vr, v_theta, sigma_r, sigma_theta = sy.symbols('r theta vr v_theta sigma_r sigma_theta')
X = sy.Matrix([r, theta, vr, v_theta])

# Y = h(X)
Y = sy.Matrix([r * sy.cos(theta), r * sy.sin(theta), vr * sy.cos(theta) - r * v_theta * sy.sin(theta),
            vr * sy.sin(theta) + r * v_theta * sy.cos(theta)])

# dX = mu dt + sigma dB
mu = sy.Matrix([vr, v_theta, 0, 0])
sigma = sy.Matrix([[0, 0], [0, 0], [sigma_r, 0], [0, sigma_theta]])

new_mu, new_sigma = ItosLemma(mu, sigma, X, Y)

print(new_mu)
print(new_sigma)