import numpy as np
from tristanVis.tmath import DOT, CROSS

def GAMMAtoBETA(gamma):
  return np.sqrt(1.0 - gamma**-2)

def BETAtoGAMMA(beta):
  return 1.0 / np.sqrt(1.0 - beta**2)

def UtoGAMMA(ux, uy, uz):
  return np.sqrt(1.0 + DOT(ux, uy, uz, ux, uy, uz))

def UtoBETA(ux, uy, uz):
  gamma = UtoGAMMA(ux, uy, uz)
  return (ux / gamma, uy / gamma, uz / gamma)

def boost_E(ex, ey, ez, bx, by, bz, ux, uy, uz):
  vx, vy, vz = UtoBETA(ux, uy, uz)
  gamma = UtoGAMMA(ux, uy, uz)
  beta_cross_B_x, beta_cross_B_y, beta_cross_B_z = CROSS(vx, vy, vz, bx, by, bz)
  beta_dot_e = DOT(vx, vy, vz, ex, ey, ez)
  dummy = gamma**2 / (gamma + 1.0)

  ex1 = gamma * (ex + beta_cross_B_x) - dummy * vx * beta_dot_e
  ey1 = gamma * (ey + beta_cross_B_y) - dummy * vy * beta_dot_e
  ez1 = gamma * (ez + beta_cross_B_z) - dummy * vz * beta_dot_e

  return (ex1, ey1, ez1)

def boost_B(bx, by, bz, ex, ey, ez, ux, uy, uz):
  vx, vy, vz = UtoBETA(ux, uy, uz)
  gamma = UtoGAMMA(ux, uy, uz)
  beta_cross_E_x, beta_cross_E_y, beta_cross_E_z = CROSS(vx, vy, vz, ex, ey, ez)
  beta_dot_b = DOT(vx, vy, vz, bx, by, bz)
  dummy = gamma**2 / (gamma + 1.0)

  bx1 = gamma * (bx - beta_cross_E_x) - dummy * vx * beta_dot_b
  by1 = gamma * (by - beta_cross_E_y) - dummy * vy * beta_dot_b
  bz1 = gamma * (bz - beta_cross_E_z) - dummy * vz * beta_dot_b

  return (bx1, by1, bz1)
