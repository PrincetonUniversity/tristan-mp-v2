import numpy as np

def getFieldline(x0, y0, z0,
                 bx, by, bz,
                 istep=2, direction=+1, maxit=1e4):
  x0 += istep / 2.
  y0 += istep / 2.
  z0 += istep / 2.
  sz, sy, sx = bx.shape
  sx *= istep; sy *= istep; sz *= istep
  tsteps = 0
  fieldline = np.array([[x0, y0, z0]])
  while((x0 >= 0 and x0 < sx) and (y0 >= 0 and y0 < sy) and (z0 >= 0 and z0 < sz) and (tsteps < maxit)):
    bx0 = bx[int(z0/istep), int(y0/istep), int(x0/istep)];
    by0 = by[int(z0/istep), int(y0/istep), int(x0/istep)];
    bz0 = bz[int(z0/istep), int(y0/istep), int(x0/istep)]
    b0 = np.sqrt(bx0**2 + by0**2 + bz0**2)
    x0 += direction * bx0 / (4 * b0)
    y0 += direction * by0 / (4 * b0)
    z0 += direction * bz0 / (4 * b0)
    tsteps += 1
    fieldline = np.append(fieldline, [[x0, y0, z0]], axis = 0)
  return fieldline[:-1]

def getFootpointFieldlines(xc, yc, rad,
                           bx, by, bz,
                           istep=2, resolution=30,
                           maxit=1e4, direction=+1):
  phi_list = np.linspace(0, 2 * np.pi, resolution)
  fieldlines = []
  for phi in phi_list:
    x0, y0, z0 = (xc + rad * np.sin(phi), yc + rad * np.cos(phi), 0)
    fld = getFieldline(x0, y0, z0, bx, by, bz, istep=istep, direction=direction, maxit=maxit)
    fieldlines.append(fld)
  return fieldlines
