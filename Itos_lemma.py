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


    new_sigma = 0
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

# # print(G)

# dh1_dX = dh_dX[0, :]
# # print(dh1_dX)
# d2h1_dX2 = dh1_dX.jacobian(X)

# # print(d2h1_dX2)

# dh2_dX = dh_dX[1, :]
# d2h2_dX2 = dh2_dX.jacobian(X)
# # print(d2h2_dX2)


# dh3_dX = dh_dX[2, :]
# d2h3_dX2 = dh3_dX.jacobian(X)
# # print(d2h3_dX2)

# dh4_dX = dh_dX[3, :]
# d2h4_dX2 = dh4_dX.jacobian(X)
# # print(d2h4_dX2)

# d2h_dX2 = [d2h1_dX2, d2h2_dX2, d2h3_dX2, d2h3_dX2]
# print(d2h_dX2[0].row(1))

# 1/0
# # --------------------- Let's hard code it.
# m_111 = 0
# m_112 = -sy.sin(theta)
# m_113 = 0
# m_114 = 0

# m_11 = [m_111, m_112, m_113, m_114]

# m_121 = -sy.sin(theta)
# m_122 = -r * sy.cos(theta)
# m_123 = 0
# m_124 = 0

# m_12 = [m_121, m_122, m_123, m_124]
# m_13 = [0, 0, 0, 0]
# m_14 = [0, 0, 0, 0]

# m_21 = [0, sy.cos(theta), 0, 0]
# # test = type(m_21[1])
# # print(test)
# # print(type(m_121))
# # test = sy.Add(m_21[1], m_121)
# # print(test)

# m_22 = [sy.cos(theta), - r*sy.sin(theta), 0, 0]
# m_23 = [0, 0, 0, 0]
# m_24 = [0, 0, 0, 0]

# m_31 = [0, -v_theta, 0, -sy.sin(theta)]

# m_321 = -v_theta * sy.cos(theta)
# m_322 = r*v_theta*sy.sin(theta) - vr * sy.cos(theta)
# m_323 = - sy.sin(theta)
# m_324 = - r*sy.cos(theta)

# m_32 = [m_321, m_322, m_323, m_324]
# m_33 = [0, -sy.sin(theta), 0, 0]
# m_34 = [-sy.sin(theta), -r*sy.cos(theta), 0, 0]

# m_41 = [0, -v_theta * sy.sin(theta), 0, sy.cos(theta)]

# m_421 = -v_theta * sy.sin(theta)
# m_422 = -r*v_theta*sy.cos(theta) - vr*sy.sin(theta)
# m_423 = sy.cos(theta)
# m_424 = -r*sy.sin(theta)

# m_42 = [m_421, m_422, m_423, m_424]
# m_43 = [0, sy.cos(theta), 0, 0]
# m_44 = [sy.cos(theta), -r*sy.sin(theta), 0, 0]

# m_1 = [m_11, m_12, m_13, m_14]
# m_2 = [m_21, m_22, m_23, m_24]
# m_3 = [m_31, m_32, m_33, m_34]
# m_4 = [m_41, m_42, m_43, m_44]

# # print(m_1)
# # print(m_2)
# # print(m_3)
# # print(m_4)

# G_1 = [0, 0]
# G_2 = [0, 0]
# G_3 = [sigma_r, 0]
# G_4 = [0, sigma_theta]

# G = [G_1, G_2, G_3, G_4]

# # For the "trace" here - beginning
# y_drift_1 = 0
# for i in range(4):
#     for j in range(4):
#         for k in range(2):
#             GG = G[i][k] * G[j][k]
#         m_ij = m_1[i][j]
#         y_drift_1 += m_ij*GG


# print(y_drift_1)

# y_drift_2 = 0
# for i in range(4):
#     for j in range(4):
#         for k in range(2):
#             GG = G[i][k] * G[j][k]
#         m_ij = m_2[i][j]
#         y_drift_2 += m_ij*GG


# print(y_drift_2)


# y_drift_3 = 0
# for i in range(4):
#     for j in range(4):
#         for k in range(2):
#             GG = G[i][k] * G[j][k]
#         m_ij = m_3[i][j]
#         y_drift_3 += m_ij*GG


# print(y_drift_3)


# y_drift_4 = 0
# for i in range(4):
#     for j in range(4):
#         for k in range(2):
#             GG = G[i][k] * G[j][k]
#         m_ij = m_4[i][j]
#         y_drift_4 += m_ij*GG


# print(y_drift_4)

# # For the "trace" here - end

# M =[m_1, m_2, m_3, m_4]
# # M = np.array(M)

# # print('Ciao')
# # print(type(M[1][0][1]))

# GT_1 = [0, 0, sigma_r, 0]
# GT_2 = [0, 0, 0, sigma_theta]

# GT = [GT_1, GT_2]

# GT_d2h_dX2 = np.zeros((2, 4, 4))
# # print(GT_d2h_dX2)

# m = 2
# n = 4

# # exit()

# GT_d2h_dX2 = []
# for i in range(m):
#     # print("i is", i)
#     GT_d2h_dX2_i = []
#     for j in range(n):
#         GT_d2h_dX2_ij = []
#         # print("j is", j)
#         for g in range(n):
#             # print("g is", g)
#             GT_d2h_dX2_ijg = 0.0
#             for k in range(n):

#                 # print("k is", k)
#                 GT_ik = GT[i][k]
#                 # tg = type(GT_ik)
#                 # print(type(GT_ik))
#                 M_kjg = M[k][j][g]
#                 # tm = type(M_kjg)
#                 # print(type(M_kjg))
#                 prod = GT_ik * M_kjg
#                 # tp = type(prod)
#                 # print("tp", type(prod))

#                 GT_d2h_dX2_ijg += prod

#                 # ERROR HERE: IT DOES NOT LIKE NUMPY AND sy TOGETHER, NOW HARD CODED WITH LISTS.
#                 # GT_d2h_dX2[i, j, k] = sy.Add(GT_d2h_dX2[i, j, k], prod)
#                 # local = GT_d2h_dX2[i, j, k]
#                 # print("local", local)
#                 # print(GT_d2h_dX2[i, j, k])

#             GT_d2h_dX2_ij.append(GT_d2h_dX2_ijg)
#         GT_d2h_dX2_i.append(GT_d2h_dX2_ij)
#     GT_d2h_dX2.append(GT_d2h_dX2_i)

# # This is a 2x4x4
# print(GT_d2h_dX2)
# print(GT_d2h_dX2[1])
# print(GT_d2h_dX2[0][3])
# print(GT_d2h_dX2[0][0][3])
# print('')

# l = 2
# m = 4
# n = 2
# p = 4
# GT_d2h_dX2_G = []
# for i in range(l):
#     # print("i is", i)
#     GT_d2h_dX2_G_i = []
#     for j in range(m):
#         GT_d2h_dX2_G_ij = []
#         # print("j is", j)
#         for g in range(n):
#             # print("g is", g)
#             GT_d2h_dX2_G_ijg = 0.0
#             for k in range(p):

#                 GT_d2h_dX2_ijk = GT_d2h_dX2[i][j][k]

#                 G_kg = G[k][g]
#                 prod = GT_d2h_dX2_ijk * G_kg
#                 # tp = type(prod)
#                 # print("tp", type(prod))

#                 GT_d2h_dX2_G_ijg += prod

#             GT_d2h_dX2_G_ij.append(GT_d2h_dX2_G_ijg)
#         GT_d2h_dX2_G_i.append(GT_d2h_dX2_G_ij)
#     GT_d2h_dX2_G.append(GT_d2h_dX2_G_i)

# # This is a 2x4x2
# print(GT_d2h_dX2_G)
# # print(GT_d2h_dX2_G[1])
# # print(GT_d2h_dX2_G[0][0])
# # print(GT_d2h_dX2_G[0][0][1])