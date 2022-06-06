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
# X_polar = sy.Matrix([r, theta, vr, omega])

# # Stochastic Two-body problem:
# mu = sy.Matrix([vr, omega, r*omega**2 - mu_grav/(r**2), -2*vr*omega/r])
# sigma = sy.Matrix([[0, 0], [0, 0], [sigma_r, 0], [0, sigma_theta]])

# V = sy.Matrix([vr, r*omega])
# V_2 = sy.Matrix([V.dot(V)])

# H = V_2/2 - sy.Matrix([mu_grav/r])

# H_mu, H_sigma = ItosLemma(mu, sigma, X_polar, H)

# print('Here is the drift and diffusion for the total energy of the two-body problem:')
# print(H_mu)
# print(H_sigma)
# print(' ')

# # Y = h(X)
# Y_polar = sy.Matrix([r * sy.cos(theta), r * sy.sin(theta), vr * sy.cos(theta) - r * omega * sy.sin(theta),
#             vr * sy.sin(theta) + r * omega * sy.cos(theta)])

# mu_one, sigma_one = ItosLemma(sy.Matrix([0, 0, 0, 0]), sigma, X_polar, Y_polar)
# print('Here are drift and diffusion in cartesian coordinates from two-body problem sigma only:')
# print(mu_one)
# #Matrix([[0], [0], [0], [0]])
# print(sigma_one)
# # Matrix([[0, 0], [0, 0], [sigma_r*cos(theta), -r*sigma_theta*sin(theta)], [sigma_r*sin(theta), r*sigma_theta*cos(theta)]])
# print(' ')

# -----------------------------
# Now I forget about polar coordinates and start from cartesian in inertial space:
X, Y, VX, VY, x, y, vx, vy = sy.symbols('X Y VX VY x y vx vy')
X_vec = sy.Matrix([X, Y, VX, VY])

At = sy.Matrix([[sy.cos(t), sy.sin(t)], [-sy.sin(t), sy.cos(t)]]) # Matrix such that X = At x

mu_2 = 9.537e-4
mu_1 = 1 - mu_2

x_sun = sy.Matrix([-mu_2, 0])
x_jup = sy.Matrix([mu_1, 0])

X_sun = At @ x_sun 
X_jup = At @ x_jup

# mu comes from the deterministic 3BP.
R1 = sy.Matrix([X, Y]) - X_sun
R2 = sy.Matrix([X, Y]) - X_jup

X1 = sy.Matrix([X - X_sun[0]])
X2 = sy.Matrix([X - X_jup[0]])
Y1 = sy.Matrix([Y - X_sun[1]])
Y2 = sy.Matrix([Y - X_jup[1]])

R1_3 = sy.Matrix([R1.dot(R1)**1.5])
R2_3 = sy.Matrix([R2.dot(R2)**1.5])

mu = sy.Matrix([VX, VY, -mu_1*X1/R1_3[0] - mu_2*X2/R2_3[0], -mu_1*Y1/R1_3[0] - mu_2*Y2/R2_3[0]])
# print(mu)

# sigma comes from sigma_one, but expressed with respect to X_vec:
sigma = sy.Matrix([[0, 0], [0, 0], [sigma_r*X/sy.sqrt(X**2 + Y**2), -Y*sigma_theta], [sigma_r*Y/sy.sqrt(X**2 + Y**2), X*sigma_theta]])

# We go from inertial to rotating with the inverse/transpose of A_t populating the map for the state, see Koon Lo Marsden Ross chapter 2.
J = sy.Matrix([[0, 1], [-1, 0]])
Atinv = sy.Matrix([[sy.cos(t), - sy.sin(t)], [sy.sin(t), sy.cos(t)]])
AtinvJ = Atinv @ J

Ct = sy.Matrix([[sy.cos(t), -sy.sin(t), 0, 0], [sy.sin(t), sy.cos(t), 0, 0], [sy.sin(t), sy.cos(t), sy.cos(t), -sy.sin(t)], [-sy.cos(t), sy.sin(t), sy.sin(t), sy.cos(t)]])
h = Ct @ X_vec
# print(h)
mu_s3bp, sigma_s3bp = ItosLemma(mu, sigma, X_vec, h)
print('Here are drift and diffusion in cartesian coordinates from two-body problem sigma only:')
print(mu_s3bp)
# print(sigma_s3bp)

# https://github.com/MattCocci/ReusableCode/blob/master/Python/ItosLemma.py Implement inverse like this to have mu and sigma expressed in terms of
# corotating coordinates, so that the Jacobi "constant" equations is trivial, and also the numerical implementation is.