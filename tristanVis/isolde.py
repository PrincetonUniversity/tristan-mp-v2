import struct
import numpy as np
import h5py
import os

def getParticles(fname):
  with h5py.File(fname, 'r') as file:
    keys = list(file.keys())
    species = np.unique([int(key.split('_')[1]) for key in keys])
    nspec = len(species)
    variables = np.unique([key.split('_')[0] for key in keys])
    nvars = len(variables)
    data = {}
    for s in range(nspec):
      data[str(s + 1)] = {}
      for i in range(nvars):
        (data[str(s + 1)])[variables[i]] = file[variables[i] + '_' + str(s + 1)][:]
  return data

def getFields(fname, nodes = False):
  # hdf5 file
  with h5py.File(fname, 'r') as file:
    keys = list(file.keys())
    data = {}
    for key in keys:
      data[key] = file[key][:]
  return data

def getParameters(fname):
  with h5py.File(fname, 'r') as file:
    keys = list(file.keys())
    params = {}
    for key in keys:
      params[key] = file[key][:][0]
  return params

def getSlice(output, proj, step):
  if (output[-1] != '/'):
    output += '/'
  proj, shift = proj.split('=')
  return getFields(output + '/slice' + proj.upper() + '=' + '%05d' % int(shift) + '.%05d' % step)

def convertToXarray(fields,
                    coordinateTransformation = {'x': lambda f: f,
                                                'y': lambda f: f,
                                                'z': lambda f: f},
                    additionalVariables = {}):
  import xarray as xr
  import numpy as np
  np.seterr(divide='ignore', invalid='ignore')
  for k in fields.keys():
    fields[k] = np.squeeze(fields[k])
  xr_data = xr.Dataset()
  dimension = len(fields[list(fields.keys())[0]].shape)
  if dimension == 1:
    xr_axes = np.array(list('xyz'))[[(np.min(fields[p+p]) != np.max(fields[p+p])) for p in list('xyz')]]
    if (len(xr_axes) != 1):
      raise ValueError("Incorrect `xr_axes`.")
    x1 = xr_axes
    xr_data.coords[x1] = ((x1), coordinateTransformation[x1](fields[x1*2][:]))
  elif dimension == 2:
    xr_axes = np.array(list('xyz'))[[(np.min(fields[p+p]) != np.max(fields[p+p])) for p in list('xyz')]]
    if (len(xr_axes) != 2):
      raise ValueError("Incorrect `xr_axes`.")
    x1, x2 = xr_axes
    xr_data.coords[x1] = ((x1), coordinateTransformation[x1](fields[x1*2][0,:]))
    xr_data.coords[x2] = ((x2), coordinateTransformation[x2](fields[x2*2][:,0]))
  elif dimension == 3:
    xr_axes = list('xyz')
    x1, x2, x3 = xr_axes
    xr_data.coords[x1] = ((x1), coordinateTransformation[x1](fields[x1*2][0,0,:]))
    xr_data.coords[x2] = ((x2), coordinateTransformation[x2](fields[x2*2][0,:,0]))
    xr_data.coords[x3] = ((x3), coordinateTransformation[x3](fields[x3*2][:,0,0]))
  for k in fields.keys():
    xr_data[k] = (xr_axes, fields[k][:])
  for k in additionalVariables.keys():
    xr_data[k] = (xr_axes, additionalVariables[k](xr_data)[:])
  return xr_data

# usage example for 2D uniform grid:
# ```
#   field_data = isolde.getFields("flds.tot.00000")
#   x_ = field_data['x']
#   y_ = field_data['y']
#   x_, y_ = np.mgrid[x_[0]: x_[-1] + 1 : x_[1] - x_[0],
#                   y_[0]: y_[-1] + 1 : y_[1] - y_[0]]
#   ex_ = field_data['ex'][:,:,0]
#   plt.pcolor(x_, y_, ex_) # <- 2D plot
# ```

class Spectra:
  def __init__(self, raw, radiation = False, gca = False):
    self.__radiation = radiation
    self.__gca = gca
    self.initialize(raw)
  def findSpecname(self, s, onlyGCA, onlyBoris):
    if (not self.__gca and (onlyGCA or onlyBoris)):
      raise ValueError('GCA set to `False` in `Spectra` class.')
    elif onlyGCA and onlyBoris:
      raise ValueError('At least one of `onlyGCA` and `onlyBoris` have to be `False`.')
    if (not onlyGCA) and (not onlyBoris):
      specname = 'n'
    elif (onlyGCA):
      specname = 'ngca'
    elif (onlyBoris):
      specname = 'nbor'
    return specname + str(s)
  def getBin(self, i, j, k):
    x0 = self.xbins[i]; y0 = self.ybins[j]; z0 = self.zbins[k]
    return (x0, y0, z0)
  def getTotal(self, s):
    ss = 'n' + str(s)
    return np.sum(self.__data[ss], axis=(0, 1, 2))
  def getBySpatialBin(self, s, ijk, onlyGCA = False, onlyBoris = False):
    i, j, k = ijk
    specname = self.findSpecname(s, onlyGCA, onlyBoris)
    xyz0 = self.getBin(i, j, k)
    if (len(self.xbins) == 1):
      sx = 1
    else:
      sx = self.xbins[1] - self.xbins[0]
    if (len(self.ybins) == 1):
      sy = 1
    else:
      sy = self.ybins[1] - self.ybins[0]
    if (len(self.zbins) == 1):
      sz = 1
    else:
      sz = self.zbins[1] - self.zbins[0]
    sxyz = (sx, sy, sz)
    return (xyz0, sxyz, self.__data[specname][i, j, k])
  def getByCoordinate(self, s, xyz, onlyGCA = False, onlyBoris = False):
    x, y, z = xyz
    ijk = [-1, -1, -1]
    bins = [self.xbins, self.ybins, self.zbins]
    for ind, (bn, crd) in enumerate(zip(bins, xyz)):
      for o in range(len(bn)):
        if bn[o] > crd:
          ijk[ind] = o - 1
          break
    return self.getBySpatialBin(s, ijk, onlyGCA, onlyBoris)
  def initialize(self, raw):
    import re
    self.__data = {}
    self.xbins = raw['xbins'][:]
    self.ybins = raw['ybins'][:]
    self.zbins = raw['zbins'][:]
    self.ebins = raw['ebins'][:]
    keylist = list(raw.keys())
    allspecies = ([int(re.findall("^n(\d+)", key)[0]) for key in keylist if re.match("^n\d+", key)])
    for s in allspecies:
      self.__data['n' + str(s)] = np.transpose(raw['n' + str(s)])
    if self.__radiation:
      self.rbins = raw['rbins'][:]
      radspecies = ([int(re.findall("^nr(\d+)", key)[0]) for key in keylist if re.match("^nr\d+", key)])
      for rs in radspecies:
        self.__data['nr' + str(rs)] = raw['nr' + str(rs)][:]
    if self.__gca:
      for s in allspecies:
        self.__data['nbor' + str(s)] = np.transpose(raw['nbor' + str(s)])
        self.__data['ngca' + str(s)] = np.transpose(raw['ngca' + str(s)])

def getSpectra(fname, radiation = False, gca = False):
  with h5py.File(fname, 'r') as file:
    spec = Spectra(file, radiation, gca)
  return spec

def getDomains(fname):
  with h5py.File(fname, 'r') as file:
    data = {}
    for k in file.keys():
      data[k] = file[k][:]
  return data

def parseReport(fname, nsteps = None, skip = 1, skip_every = 1e6):
  if (not nsteps):
    nsteps = 1e6
  def parseBlock(block, data, isfirst = False):
    for line in block.split('\n')[2:]:
      try:
        routine = line.split()[0]
      except:
        continue
      writing_particles = False
      if (routine == 'species'):
        writing_particles = True
        routine = line[:15].strip()
      if routine[-1] == ':':
        routine = routine[:-1]
      if (routine != '' and routine[0] != '['):
        line1 = line.split()
        if (isfirst):
          data[routine] = {}
          if not writing_particles:
            data[routine]['dt'] = np.array([])
            data[routine]['min'] = np.array([])
            data[routine]['max'] = np.array([])
          else:
            data[routine]['average'] = np.array([])
            data[routine]['min'] = np.array([])
            data[routine]['max'] = np.array([])
            data[routine]['total'] = np.array([])
        if len(line1) < 5:
          line1 = line1[-3:]
        else:
          line1 = line1[-4:]
        nums = [float(x.strip()) for x in line1]
        if (len(nums) < 3):
          print (nums, line)
          raise ValueError('len(nums) < 3')
        else:
          if not writing_particles:
            data[routine]['dt'] = np.append(data[routine]['dt'], [nums[0]])
            data[routine]['min'] = np.append(data[routine]['min'], [nums[1]])
            data[routine]['max'] = np.append(data[routine]['max'], [nums[2]])
          else:
            data[routine]['average'] = np.append(data[routine]['average'], [nums[0]])
            data[routine]['min'] = np.append(data[routine]['min'], [nums[1]])
            data[routine]['max'] = np.append(data[routine]['max'], [nums[2]])
            data[routine]['total'] = np.append(data[routine]['total'], [nums[3]])
  data = {}
  data['t'] = np.array([])
  with open(fname, 'r') as file:
    line = file.readline()
    isfirst = True
    ni = 0
    while line and (ni < nsteps):
      while (line.strip()[0:10] != '-'*10) and line:
        line = file.readline()
      block = ""
      line = file.readline()
      while (line.strip()[0:10] != '.'*10) and line:
        block += line
        line = file.readline()
      if (ni % skip == 0) and (ni % skip_every != 0) and line:
        parseBlock(block, data, isfirst = isfirst)
        isfirst = False
        data['t'] = np.append(data['t'], [ni])
      ni += 1
  return data

def parseHistory(fname):
  from itertools import groupby
  import re
  def make_grouper():
    counter = 0
    def key(line):
      nonlocal counter
      if line.startswith('===='):
        counter += 1
      return counter
    return key
  with open(fname, 'r') as f:
    data = {}
    for k, group in groupby(f, key=make_grouper()):
      fasta_section = ''.join(group)
      block = fasta_section.split("\n", 1)[1]
      if (k == 1):
        template = block
        template_keys = np.array(re.findall('\[.+?\]', template))
        mask = (template_keys != '[% Etot]')
        template_keys = np.array(list(map(lambda x: x[1:-1], template_keys[mask])))
        data = {key: np.array([]) for key in template_keys}
      else:
        block_values = np.array(re.findall('-?\ *[0-9]+\.?[0-9]*(?:[Ee]\ *[-|+]?\ *[0-9]+)?', block))
        if (len(block) == 0):
          break
        block_values = np.array(list(map(np.float, block_values[mask])))
        pairs = {key: v for key, v in zip(template_keys, block_values)}
        for key in template_keys:
          data[key] = np.append(data[key], pairs[key])
  return data
