def integrateFieldline2D(point, xx, yy, fx, fy,
                         direction, stop_condition = (lambda x: False),
                         maxnsteps = 5000, dx = 0.01):
  import numpy as np
  sy, sx = fx.shape
  xmax = xx.max(); xmin = xx.min()
  ymax = yy.max(); ymin = yy.min()
  fieldline = [point]
  nstep = 0
  x, y = point
  # simple eulerian integrator ...
  # ... just make sure dx is small enough
  for nstep in range(maxnsteps):
    if (stop_condition(np.array([x, y]))):
      break
    i = int(min(max((x - xmin) * sx / (xmax - xmin), 0), sx - 1))
    j = int(min(max((y - ymin) * sy / (ymax - ymin), 0), sy - 1))
    fx_ = direction * fx[j, i]
    fy_ = direction * fy[j, i]
    ff_ = np.sqrt(fx_**2 + fy_**2)
    x += dx * fx_ / ff_
    y += dx * fy_ / ff_
    fieldline.append([x, y])
  return np.array(fieldline)

def slice1D(field, points):
  from scipy.interpolate import RegularGridInterpolator
  import numpy as np
  coords = list(field.coords.keys())
  fld = RegularGridInterpolator([field.coords[c] for c in coords[::-1]], field.data)
  return np.array([fld(point[::-1]) for point in points])

def vectorPotential(bx, by):
  import numpy as np
  my, mx = bx.shape
  bydx = -np.cumsum(by[0,:])
  bxdy = np.cumsum(bx[:,:], axis = 0)
  Az = bydx[np.newaxis,:] + bxdy
  return Az
