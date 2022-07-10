import numpy as np

def ADD(ax, ay, az, bx, by, bz):
  return (ax + bx, ay + by, az + bz)

def DIVIDE(ax, ay, az, bx, by, bz):
  return (ax / bx, ay / by, az / bz)

def OVER(ax, ay, az, b):
  return (ax / b, ay / b, az / b)

def MULT(ax, ay, az, bx, by, bz):
  return (ax * bx, ay * by, az * bz)

def TIMES(ax, ay, az, b):
  return (ax * b, ay * b, az * b)

def CROSS(ax, ay, az, bx, by, bz, i=None):
  if i is None:
    return (-(az * by) + ay * bz, az * bx - ax * bz, -(ay * bx) + ax * by)
  elif i == 'x' or i == 0:
    return -(az * by) + ay * bz
  elif i == 'y' or i == 1:
    return az * bx - ax * bz
  elif i == 'z' or i == 2:
    return -(ay * bx) + ax * by
  else:
    raise ValueError('Unknown component, ' + str(i) + '.')

def DOT(ax, ay, az, bx, by, bz):
  return ax * bx + ay * by + az * bz

def NORM(ax, ay, az):
  return np.sqrt(ax**2 + ay**2 + az**2)

def PARALLEL(ax, ay, az, bx, by, bz):
  return DOT(ax, ay, az, bx, by, bz) / NORM(bx, by, bz)

def PERPENDICULAR(ax, ay, az, bx, by, bz):
  return np.sqrt(DOT(ax, ay, az, ax, ay, az) - PARALLEL(ax, ay, az, bx, by, bz)**2)

def VEC3(obj, name):
  return (obj[name + 'x'], obj[name + 'y'], obj[name + 'z'])

def PRTLSUM(obj, name, species):
  f = obj[name + species[0]]
  for s in list(species)[1:]:
    f += obj[name + s]
  return f

def cart2sph(ax, ay, az, x, y, z, xc, yc, zc):
  rx = x - xc; ry = y - yc; rz = z - zc
  theta = np.arctan2(np.hypot(rx, ry), rz)
  phi = np.arctan2(ry, rx)
  ar = np.sin(theta) * (np.cos(phi) * ax + np.sin(phi) * ay) + np.cos(theta) * az
  atheta = np.cos(theta) * (np.cos(phi) * ax + np.sin(phi) * ay) - np.sin(theta) * az
  aphi = -np.sin(phi) * ax + np.cos(phi) * ay
  return (ar, atheta, aphi)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

class Fields:
  __slots__ = ("dataset",)
  vectors = ['e', 'b', 'j', 'curlB', 'vel', 'ExB']
  specbased = ['dens', 'dgca', 'enrg', 'momX', 'momY', 'momZ', 'nprt']
  sph = False
  def __init__(self, ds, xyz0):
    self.__class__.dataset.__set__(self, ds)
    for vec in self.vectors:
      try:
        xyz = [getattr(self.dataset, i*2) for i in ['x', 'y', 'z']]
        fxyz = [getattr(self.dataset, vec + i) for i in ['x', 'y', 'z']]
        self.dataset[vec+'r'], self.dataset[vec+'theta'], self.dataset[vec+'phi'] = cart2sph(*fxyz, *xyz, *xyz0)
      except:
        continue
  def __setattr__(self, k, v):
    setattr(self.dataset, k, v)
  def __getattr__(self, attr):
    if attr in self.vectors:
      return self.getVectorField(attr)
    elif attr[:4] in self.specbased and len(attr[4:]) > 1:
      specs = list(attr[4:])
      f = getattr(self.dataset, attr[:4] + specs[0])
      for s in specs[1:]:
        f = f + getattr(self.dataset, attr[:4] + s)
      return f.rename(attr)
    else:
      return getattr(self.dataset, attr)
  def getVectorField(self, attr):
    vecfield_xyz = [getattr(self.dataset, attr + i) for i in ['x', 'y', 'z']]
    vecfield_sph = [getattr(self.dataset, attr + i) for i in ['r', 'theta', 'phi']]
    xyz = [getattr(self.dataset, i*2) for i in ['x', 'y', 'z']]
    return VectorField(f_x_y_z=vecfield_xyz, f_r_theta_phi=vecfield_sph)

class VectorField:
  def __init__(self, f_x_y_z = None, f_r_theta_phi = None):
    self._cart = (f_x_y_z is not None)
    self._sph = (f_r_theta_phi is not None)
    if (f_x_y_z is not None):
      self._fx, self._fy, self._fz = f_x_y_z
    if (f_r_theta_phi is not None):
      self._fr, self._ftheta, self._fphi = f_r_theta_phi
  @property
  def x(self):
    assert self._cart
    return self._fx
  @property
  def y(self):
    assert self._cart
    return self._fy
  @property
  def z(self):
    assert self._cart
    return self._fz
  @property
  def norm(self):
    assert ((self._cart) or (self._sph))
    if (self._cart):
      return np.sqrt(self._fx**2 + self._fy**2 + self._fz**2)
    if (self._sph):
      return np.sqrt(self._fr**2 + self._ftheta**2 + self._fphi**2)
  @property
  def r(self):
    assert self._sph
    return self._fr
  @property
  def theta(self):
    assert self._sph
    return self._ftheta
  @property
  def phi(self):
    assert self._sph
    return self._fphi
