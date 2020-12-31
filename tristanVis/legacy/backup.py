
# from abc import ABC, abstractmethod
# import tristanVis.aux
#
# class Simulation(ABC):
#   def __init__(self, root):
#     self._root = root
#     self.data = None
#     self.axes = None
#   @property
#   def step(self):
#     return self._step
#   @property
#   def root(self):
#     return self._root
#   @step.setter
#   def step(self, step):
#     self._step = step
#     self.loadData()
#   @abstractmethod
#   def loadData(self):
#     pass
#   @abstractmethod
#   def drawData(self):
#     tristanVis.aux.loadCustomStyles()
#   def requestField(self, key):
#     import h5py
#     try:
#       with h5py.File(self._root + 'flds.tot.%05d' % self._step, 'r') as fields:
#         self.data[key] = (self.axes, fields[key][:])
#     except:
#       print (key, 'not found in ...')
#       print ('...', self._root + 'flds.tot.%05d' % self._step)
#   def requestParam(self, key):
#     import h5py
#     try:
#       with h5py.File(self._root + 'params.%05d' % self._step, 'r') as params:
#         self.data.attrs[key] = params[key][0]
#     except:
#       print (key, 'not found in ...')
#       print ('...', self._root + 'params.%05d' % self._step)
#   def saveFig(self, savefig=None):
#     import matplotlib.pyplot as plt
#     if savefig is not None:
#       plt.savefig(savefig)
#       plt.close()
#     else:
#       plt.show()
#
# class SliceSimulation(ABC):
#   def __init__(self, root):
#     self._root = root
#     self.data = None
#     self.axes = None
#   @property
#   def step(self):
#     return self._step
#   @property
#   def root(self):
#     return self._root
#   @step.setter
#   def step(self, step):
#     self._step = step
#     self.loadData()
#   @abstractmethod
#   def loadData(self):
#     pass
#   @abstractmethod
#   def drawData(self):
#     tristanVis.aux.loadCustomStyles()
#   # def requestField(self, key):
#   #   import h5py
#   #   try:
#   #     with h5py.File(self._root + 'flds.tot.%05d' % self._step, 'r') as fields:
#   #       self.data[key] = (self.axes, fields[key][:])
#   #   except:
#   #     print (key, 'not found in ...')
#   #     print ('...', self._root + 'flds.tot.%05d' % self._step)
#   def requestParam(self, key):
#     import h5py
#     try:
#       with h5py.File(self._root + 'params.%05d' % 0, 'r') as params:
#         self.data.attrs[key] = params[key][0]
#     except:
#       print (key, 'not found in ...')
#       print ('...', self._root + 'params.%05d' % self._step)
#   def saveFig(self, savefig=None):
#     import matplotlib.pyplot as plt
#     if savefig is not None:
#       plt.savefig(savefig)
#       plt.close()
#     else:
#       plt.show()
#
# class GenericSimulation(Simulation):
#   def loadData(self):
#     import h5py
#     import xarray as xr
#     import numpy as np
#     np.seterr(divide='ignore', invalid='ignore')
#     with h5py.File(self._root + 'flds.tot.%05d' % self._step, 'r') as fields:
#       with h5py.File(self._root + 'params.%05d' % self._step, 'r') as params:
#         self.axes = ('z', 'y', 'x')
#         self.data = xr.Dataset()
#         self.data.attrs['t'] = params['timestep'][:][0]
#         try:
#           self.data.attrs['sx'] = params['grd:mx0'][:][0]
#         except:
#           self.data.attrs['sx'] = 1
#         try:
#           self.data.attrs['sy'] = params['grd:my0'][:][0]
#         except:
#           self.data.attrs['sy'] = 1
#         try:
#           self.data.attrs['sz'] = params['grd:mz0'][:][0]
#         except:
#           self.data.attrs['sz'] = 1
#
#         for key in fields.keys():
#           self.data[key] = (self.axes, fields[key][:])
#
#         self.data.coords['x'] = (('x'), fields['xx'][:][0,0,:])
#         self.data.coords['y'] = (('y'), fields['yy'][:][0,:,0])
#         self.data.coords['z'] = (('z'), fields['zz'][:][:,0,0])
#   def drawData(self):
#     super().drawData()
#     pass
#
# # import numpy as np
# # np.seterr(divide='ignore', invalid='ignore')
# #
# # # Generic data containers
# #
# # class ScalarFieldCartesian:
# #   def __init__(self, data, grid=None, istep=1):
# #     if isinstance(data, ScalarFieldCartesian):
# #       self.__data = data.data
# #     else:
# #       self.__data = data
# #     self.__data[np.isnan(self.__data)] = 0
# #     self.__istep = istep
# #     if grid is not None:
# #       self.__grid = grid
# #     else:
# #       self.__grid = {'x': np.arange(0, self.__data.shape[2] * istep, istep),
# #                      'y': np.arange(0, self.__data.shape[1] * istep, istep),
# #                      'z': np.arange(0, self.__data.shape[0] * istep, istep)}
# #   def slice2D(self, slicetext):
# #     plane, amount = slicetext.split('=')
# #     sz, sy, sx = self.__data.shape
# #     if '%' in amount:
# #       amount = np.float(amount[:-1]) / 100
# #       if plane == 'x':
# #         return ScalarFieldCartesian(self.__data[:,:,int(amount * sx)], self.__grid, self.__istep)
# #       elif plane == 'y':
# #         return ScalarFieldCartesian(self.__data[:,int(amount * sy),:], self.__grid, self.__istep)
# #       elif plane == 'z':
# #         return ScalarFieldCartesian(self.__data[int(amount * sz),:,:], self.__grid, self.__istep)
# #     else:
# #       if plane == 'x':
# #         return ScalarFieldCartesian(self.__data[:,:,int(amount)], self.__grid, self.__istep)
# #       elif plane == 'y':
# #         return ScalarFieldCartesian(self.__data[:,int(amount),:], self.__grid, self.__istep)
# #       elif plane == 'z':
# #         return ScalarFieldCartesian(self.__data[int(amount),:,:], self.__grid, self.__istep)
# #   def slice1D(self, xs, ys, zs):
# #     from scipy.interpolate import RegularGridInterpolator
# #     data_int = RegularGridInterpolator((self.__grid['z'], self.__grid['y'], self.__grid['x']), self.__data)
# #     return data_int(np.array([zs, ys, xs]).T)
# #   def rotateAround(self, axis, angle, order=0):
# #     import scipy.ndimage
# #     while (angle >= 2 * np.pi):
# #       angle -= 2 * np.pi
# #     while (angle <= -2 * np.pi):
# #       angle += 2 * np.pi
# #     if (angle == 0):
# #       data = self.__data
# #     else:
# #       if axis == 'x':
# #         data = scipy.ndimage.interpolation.rotate(self.__data, angle=angle*180/np.pi, axes=(1, 0), order=order, mode='nearest')
# #       elif axis == 'y':
# #         data = scipy.ndimage.interpolation.rotate(self.__data, angle=angle*180/np.pi, axes=(0, 2), order=order, mode='nearest')
# #       elif axis == 'z':
# #         data = scipy.ndimage.interpolation.rotate(self.__data, angle=angle*180/np.pi, axes=(2, 1), order=order, mode='nearest')
# #     return ScalarFieldCartesian(data, self.__istep)
# #   @property
# #   def data(self):
# #     return self.__data
# #   @property
# #   def grid(self):
# #     return self.__grid
# #   @property
# #   def istep(self):
# #     return self.__istep
# #   def __mul__(self, other):
# #     if type(self) is type(other):
# #       return ScalarFieldCartesian(self.__data * other.data, self.__grid, self.__istep)
# #     elif isinstance(other, (int, float, np.float32, np.float64, np.int)):
# #       return ScalarFieldCartesian(self.__data * other, self.__grid, self.__istep)
# #     else:
# #       raise TypeError('no idea how to multiply like that')
# #   __rmul__ = __mul__
# #   def __sub__(self, other):
# #     if type(self) is type(other):
# #       return ScalarFieldCartesian(self.__data - other.data, self.__grid, self.__istep)
# #     elif isinstance(other, (int, float, np.float32, np.float64, np.int)):
# #       return ScalarFieldCartesian(self.__data - other, self.__grid, self.__istep)
# #     else:
# #       raise TypeError('no idea how to sub like that')
# #   def __rsub__(self, other):
# #     if type(self) is type(other):
# #       return ScalarFieldCartesian(other.data - self.__data, self.__grid, self.__istep)
# #     elif isinstance(other, (int, float, np.float32, np.float64, np.int)):
# #       return ScalarFieldCartesian(other - self.__data, self.__grid, self.__istep)
# #     else:
# #       raise TypeError('no idea how to sub like that')
# #   def __add__(self, other):
# #     if type(self) is type(other):
# #       return ScalarFieldCartesian(self.__data + other.data, self.__grid, self.__istep)
# #     elif isinstance(other, (int, float, np.float32, np.float64, np.int)):
# #       return ScalarFieldCartesian(self.__data + other, self.__grid, self.__istep)
# #     else:
# #       raise TypeError('no idea how to add like that')
# #   __radd__ = __add__
# #   def __pow__(self, other):
# #     if isinstance(other, (int, float, np.float32, np.float64, np.int)):
# #       return ScalarFieldCartesian(self.__data**other, self.__grid, self.__istep)
# #     else:
# #       raise TypeError('no idea how to raise to power like that')
# #   def __truediv__(self, other):
# #     try:
# #       return ScalarFieldCartesian(self.__data / other.data, self.__grid, self.__istep)
# #     except:
# #       if isinstance(other, (int, float, np.float32, np.float64, np.int)):
# #         return ScalarFieldCartesian(self.__data / other, self.__grid, self.__istep)
# #       else:
# #         raise TypeError('no idea how to divide like that')
# #   def __rtruediv__(self, other):
# #     if type(self) is type(other):
# #       return ScalarFieldCartesian(other.data / self.__data, self.__grid, self.__istep)
# #     if isinstance(other, (int, float, np.float32, np.float64, np.int)):
# #       return ScalarFieldCartesian(other / self.__data, self.__grid, self.__istep)
# #     else:
# #       raise TypeError('no idea how to divide like that')
# #   def __neg__(self):
# #     return ScalarFieldCartesian(-self.__data, self.__grid, self.__istep)
# #
# # class VectorFieldCartesian:
# #   def __init__(self, data_x, data_y, data_z, grid=None, istep=1):
# #     self.__x = ScalarFieldCartesian(data_x, grid, istep)
# #     self.__y = ScalarFieldCartesian(data_y, grid, istep)
# #     self.__z = ScalarFieldCartesian(data_z, grid, istep)
# #     self.__istep = istep
# #     if grid is not None:
# #       self.__grid = grid
# #     else:
# #       self.__grid = {'x': np.arange(0, self.__x.shape[2] * istep, istep),
# #                      'y': np.arange(0, self.__x.shape[1] * istep, istep),
# #                      'z': np.arange(0, self.__x.shape[0] * istep, istep)}
# #   def slice2D(self, slicetext):
# #     return VectorFieldCartesian(self.__x.slice2D(slicetext),
# #                                 self.__y.slice2D(slicetext),
# #                                 self.__z.slice2D(slicetext), self.__grid, self.__istep)
# #   @property
# #   def x(self):
# #     return self.__x
# #   @property
# #   def y(self):
# #     return self.__y
# #   @property
# #   def z(self):
# #     return self.__z
# #   @property
# #   def grid(self):
# #     return self.__grid
# #   @property
# #   def istep(self):
# #     return self.__istep
# #   def __add__(self, other):
# #     if type(self) is type(other):
# #       return VectorFieldCartesian(self.__x + other.x, self.__y + other.y, self.__z + other.z, self.__grid, self.__istep)
# #     elif isinstance(other, (int, float, np.float32, np.float64, np.int)):
# #       return VectorFieldCartesian(self.__x + other, self.__y + other, self.__z + other, self.__grid, self.__istep)
# #     elif isinstance(other, (list, type(np.array([])))) and len(other) == 3:
# #       return VectorFieldCartesian(self.__x + other[0], self.__y + other[1], self.__z + other[2], self.__grid, self.__istep)
# #     else:
# #       raise TypeError('no idea how to add like that')
# #   __radd__ = __add__
# #   def __sub__(self, other):
# #     if type(self) is type(other):
# #       return VectorFieldCartesian(self.__x - other.x, self.__y - other.y, self.__z - other.z, self.__grid, self.__istep)
# #     elif isinstance(other, (int, float, np.float32, np.float64, np.int)):
# #       return VectorFieldCartesian(self.__x - other, self.__y - other, self.__z - other, self.__grid, self.__istep)
# #     elif isinstance(other, (list, type(np.array([])))) and len(other) == 3:
# #       return VectorFieldCartesian(self.__x - other[0], self.__y - other[1], self.__z - other[2], self.__grid, self.__istep)
# #     else:
# #       raise TypeError('no idea how to sub like that')
# #   def __rsub__(self, other):
# #     if type(self) is type(other):
# #       return VectorFieldCartesian(other.x - self.__x, other.y - self.__y, other.z - self.__z, self.__grid, self.__istep)
# #     elif isinstance(other, (int, float, np.float32, np.float64, np.int)):
# #       return VectorFieldCartesian(other - self.__x, other - self.__y, other - self.__z, self.__grid, self.__istep)
# #     elif isinstance(other, (list, type(np.array([])))) and len(other) == 3:
# #       return VectorFieldCartesian(other[0] - self.__x, other[1] - self.__y, other[2] - self.__z, self.__grid, self.__istep)
# #     else:
# #       raise TypeError('no idea how to sub like that')
# #   def __mul__(self, other):
# #     if isinstance(other, ScalarFieldCartesian) or isinstance(other, (int, float, np.float32, np.float64, np.int)):
# #       return VectorFieldCartesian(self.__x * other, self.__y * other, self.__z * other, self.__grid, self.__istep)
# #     else:
# #       raise TypeError('no idea how to multiply like that')
# #   __rmul__ = __mul__
# #   def __truediv__(self, other):
# #     if isinstance(other, ScalarFieldCartesian) or isinstance(other, (int, float, np.float32, np.float64, np.int)):
# #       return VectorFieldCartesian(self.__x / other, self.__y / other, self.__z / other, self.__grid, self.__istep)
# #     else:
# #       raise TypeError('no idea how to divide like that')
# #   def __rtruediv__(self, other):
# #     if isinstance(other, ScalarFieldCartesian) or isinstance(other, (int, float, np.float32, np.float64, np.int)):
# #       return VectorFieldCartesian(other / self.__x, other / self.__y, other / self.__z, self.__grid, self.__istep)
# #     else:
# #       raise TypeError('no idea how to divide like that')
# #   def __neg__(self):
# #     return VectorFieldCartesian(-self.__x, -self.__y, -self.__z, self.__grid, self.__istep)
# #
# # # Arythmetic operations
# #
# # def DotProduct(vField1, vField2):
# #   if type(vField1) is type(vField2):
# #     try:
# #       return ScalarFieldCartesian(vField1.x * vField2.x + vField1.y * vField2.y + vField1.z * vField2.z, vField1.grid, vField1.istep)
# #     except:
# #       raise TypeError('unknown vector field type')
# #   else:
# #     raise TypeError('vector fields are incompatible')
# #
# # def CrossProduct(vField1, vField2):
# #   if type(vField1) is type(vField2):
# #     try:
# #       return VectorFieldCartesian(-vField1.z * vField2.y + vField1.y * vField2.z,
# #                                   vField1.z * vField2.x - vField1.x * vField2.z,
# #                                   -vField1.y * vField2.x + vField1.x * vField2.y, vField1.grid, vField1.istep)
# #     except:
# #       raise TypeError('unknown vector field type')
# #   else:
# #     raise TypeError('vector fields are incompatible')
# #
# # def Norm(vField):
# #   try:
# #     return DotProduct(vField, vField)**0.5
# #   except:
# #     raise TypeError('unknown vector field type')
# #
# # def Abs(Field):
# #   try:
# #     return VectorFieldCartesian(Abs(Field.x), Abs(Field.y), Abs(Field.z), Field.grid, Field.istep)
# #   except:
# #     return ScalarFieldCartesian(np.abs(Field.data), Field.grid, Field.istep)
# #
# # def Sqrt(Field):
# #   try:
# #     return VectorFieldCartesian(Sqrt(Field.x), Sqrt(Field.y), Sqrt(Field.z), Field.grid, Field.istep)
# #   except:
# #     return ScalarFieldCartesian(np.sqrt(Field.data), Field.grid, Field.istep)
# #
# # # Simulation classes
# #
# # class Simulation:
# #   def __init__(self, root, **kwargs):
# #     self.__root = root
# #     self.__step = 0
# #     self.__loaded_step = -1
# #     self._field_data = None
# #     self._kwargs = kwargs
# #   @property
# #   def root(self):
# #     return self.__root
# #   @property
# #   def step(self):
# #     return self.__step
# #   @step.setter
# #   def step(self, step):
# #     self.__step = step
# #     self.update()
# #   @property
# #   def isLoaded(self):
# #     return ((self._field_data is not None) and (self.__loaded_step == self.__step))
# #   def get(self, name):
# #     try:
# #       return self._parameters[name]
# #     except:
# #       return self._field_data[name]
# #   def update(self):
# #     self.__readAllData()
# #     self.__loaded_step = self.__step
# #   def problemDependentParams(self, dataspace, params):
# #     pass
# #   def problemDependentFields(self, dataspace, fields):
# #     pass
# #   def __readAllData(self):
# #     import h5py
# #     flds_fname = self.__root + 'flds.tot.%05d' % self.__step
# #     params_fname = self.__root + 'params.%05d' % self.__step
# #     self._parameters = {}
# #     self._field_data = {}
# #     with h5py.File(params_fname, 'r') as f:
# #       for key in f.keys():
# #         self._parameters[key] = f[key][0]
# #       self._parameters['xmin'] = 0
# #       self._parameters['xmax'] = self._parameters['grd:mx0']
# #       self._parameters['ymin'] = 0
# #       self._parameters['ymax'] = self._parameters['grd:my0']
# #       self._parameters['zmin'] = 0
# #       self._parameters['zmax'] = self._parameters['grd:mz0']
# #       self.problemDependentParams(f, self._parameters)
# #     with h5py.File(flds_fname, 'r') as f:
# #       grid = {'x': f['xx'][0,0,:], 'y': f['yy'][0,:,0], 'z': f['zz'][:,0,0]}
# #       self._field_data['xyz'] = VectorFieldCartesian(f['xx'][:], f['yy'][:], f['zz'][:], grid, self._parameters['out:istep'])
# #       self._field_data['e'] = VectorFieldCartesian(f['ex'][:], f['ey'][:], f['ez'][:], grid, self._parameters['out:istep'])
# #       self._field_data['b'] = VectorFieldCartesian(f['bx'][:], f['by'][:], f['bz'][:], grid, self._parameters['out:istep'])
# #       self.problemDependentFields(f, self._field_data)
# #
# # class PulsarSimulation(Simulation):
# #   def problemDependentParams(self, dataspace, params):
# #     params.update({'xc': params['grd:mx0'] / 2, 'yc': params['grd:my0'] / 2, 'zc': params['grd:mz0'] / 2})
# #     params.update({'qe': params['alg:c']**2 / (params['pls:ppc0'] * params['pls:c_omp']**2)})
# #     params.update({'B0': params['alg:c']**2 * np.sqrt(params['pls:sigma']) / params['pls:c_omp']})
# #     params.update({'me': params['qe']})
# #     params.update({'LC': params['alg:c'] * params['prb:psr_period'] / (2 * np.pi)})
# #     params = {k.replace('alg:', ''): v for k, v in params.items()}
# #     params = {k.replace('pls:', ''): v for k, v in params.items()}
# #     params = {k.replace('prb:', ''): v for k, v in params.items()}
# #     params = {k.replace('sigma', 'sigma0'): v for k, v in params.items()}
# #     self._parameters = params
# #   def problemDependentFields(self, dataspace, fields):
# #     self._kwargs['light'] = self._kwargs.get('light', False)
# #     if not self._kwargs['light']:
# #       fields['rho+'] = ScalarFieldCartesian(dataspace['dens2'][:], self._field_data['xyz'].grid, self._field_data['xyz'].istep)
# #       fields['rho-'] = ScalarFieldCartesian(dataspace['dens1'][:], self._field_data['xyz'].grid, self._field_data['xyz'].istep)
# #       fields['enrg+'] = ScalarFieldCartesian(dataspace['enrg2'][:], self._field_data['xyz'].grid, self._field_data['xyz'].istep)
# #       fields['enrg-'] = ScalarFieldCartesian(dataspace['enrg1'][:], self._field_data['xyz'].grid, self._field_data['xyz'].istep)
# #       fields['sigma'] = self._parameters['sigma0'] * (Norm(fields['b']) / self._parameters['B0'])**2 * (self._parameters['ppc0'] / (fields['enrg+'] + fields['enrg-']))
# #       fields['j'] = VectorFieldCartesian(dataspace['jx'][:], dataspace['jy'][:], dataspace['jz'][:], self._field_data['xyz'].grid, self._field_data['xyz'].istep)
# #       fields['e_DOT_b'] = DotProduct(fields['e'], fields['b']) / Norm(fields['b'])**2
# #       fields['gmean'] = (fields['enrg+'] + fields['enrg-']) / (fields['rho+'] + fields['rho-'])
# #       fields['de'] = self._parameters['c_omp'] * Sqrt(fields['gmean']) * (self._parameters['ppc0'] / (fields['rho+'] + fields['rho-']))
# #       fields['rL'] = fields['gmean'] * (self._parameters['c_omp'] / np.sqrt(self._parameters['sigma0'])) *\
# #                                        (self._parameters['B0'] / Norm(fields['b']))
# #       fields['rvec'] = fields['xyz'] - [self._parameters['xc'], self._parameters['yc'], self._parameters['zc']]
# #       fields['r'] = Norm(fields['rvec'])
# #       fields['R'] = (fields['rvec'].x**2 + fields['rvec'].y**2)**0.5
# #       fields['Omega'] = (CrossProduct(fields['e'], fields['b']) / Norm(fields['b'])**2) * self._parameters['LC'] / fields['R']
# #       fields['rhoGJ'] = -2 * (2 * np.pi / self._parameters['psr_period']) * fields['b'].z / (self._parameters['c'] * self._parameters['qe'])
# #
# # class LoopSimulation(Simulation):
# #   def problemDependentParams(self, dataspace, params):
# #     self._kwargs['debug'] = self._kwargs.get('debug', False)
# #     if self._kwargs['debug']:
# #       params['xmin'] -= 5; params['ymin'] -= 5; params['zmin'] -= 5
# #       params['xmax'] += 5; params['ymax'] += 5; params['zmax'] += 5
# #     params = {k.replace('prb:', ''): v for k, v in params.items()}
# #     self._parameters = params
# #   def problemDependentFields(self, dataspace, fields):
# #
# #     fields['rho'] = ScalarFieldCartesian(dataspace['var1'][:], self._field_data['xyz'].grid, self._field_data['xyz'].istep)
# #
# # # Visuals
# #
# # class Figure:
# #   def __init__(self):
# #     self.__subplots = []
# #     self.__subplot_locs = []
# #     self.__is3d = []
# #   def addPlot(self, kind, loc, data_dict=None, **kwargs):
# #     if kind is Plot3DLines:
# #       self.__is3d.append(True)
# #     else:
# #       self.__is3d.append(False)
# #     if (loc[0] <= 0 or loc[1] <= 0):
# #       raise ValueError('invalid subplot location')
# #     else:
# #       self.__subplot_locs.append(loc)
# #       if (issubclass(kind, Subplot)):
# #         self.__subplots.append(kind(data_dict, **kwargs))
# #       else:
# #         self.__subplots.append(kind(**kwargs))
# #   def render(self, figsize, title=None, titleargs=None, fontsize=15, style='dark_background', **kwargs):
# #     import matplotlib.pyplot as plt
# #     from matplotlib import rc
# #     import matplotlib as mpl
# #     mpl.rcParams.update(mpl.rcParamsDefault)
# #     rc('font',**{'family':'monospace','sans-serif':['Verdana'],'size':fontsize})
# #     rc('text', usetex=True)
# #     plt.style.use(style)
# #     loadCustomColormaps()
# #     fig = plt.figure(figsize=figsize)
# #     imax = 0; jmax = 0
# #     for sp, loc in zip(self.__subplots, self.__subplot_locs):
# #       i, j = loc
# #       if i > imax:
# #         imax = i
# #       if j > jmax:
# #         jmax = j
# #     fig_dict = {}
# #     for sp, loc, is3d in zip(self.__subplots, self.__subplot_locs, self.__is3d):
# #       if loc not in fig_dict.keys():
# #         i, j = loc
# #         n = i + imax * (j - 1)
# #         if (is3d):
# #           from mpl_toolkits.mplot3d import Axes3D
# #           fig_dict.update({loc: plt.subplot(jmax, imax, n, projection='3d')})
# #         else:
# #           fig_dict.update({loc: plt.subplot(jmax, imax, n)})
# #     for sp, loc in zip(self.__subplots, self.__subplot_locs):
# #       sp.render(fig=fig, ax=fig_dict[loc], subplots=self.__subplots, subplot_locs=self.__subplot_locs)
# #     if title is not None:
# #       fig.suptitle(title, **titleargs)
# #     fig.subplots_adjust(**kwargs)
# #     return fig, fig_dict
# #   def show(self, savefig=None):
# #     import matplotlib.pyplot as plt
# #     if savefig is not None:
# #       plt.savefig(savefig)
# #       plt.close()
# #     else:
# #       plt.show()
# #
# # class Annotation:
# #   def __init__(self, **kwargs):
# #     self._kwargs = kwargs
# #
# # class Circle(Annotation):
# #   def render(self, fig, ax, subplots=None, subplot_locs=None):
# #     from matplotlib.patches import Circle
# #     ax.add_patch(Circle((self._kwargs['x0'], self._kwargs['y0']), self._kwargs['r'], zorder=self._kwargs.get('zorder', 1),
# #                         edgecolor=self._kwargs.get('ec', 'red'), facecolor=self._kwargs.get('fc', 'none'),
# #                         ls=self._kwargs.get('ls', '-'), lw=self._kwargs.get('lw', 1)))
# #
# # class Subplot:
# #   def __init__(self, data_dict, **kwargs):
# #     self._data_dict = data_dict
# #     self._kwargs = kwargs
# #   def commonVisuals(self, ax, subplots=None, subplot_locs=None):
# #     if 'sharex' in self._kwargs.keys():
# #       sp_i = subplot_locs.index(self._kwargs['sharex'])
# #       sharex_subplot = subplots[sp_i]
# #       sharex_kwargs = sharex_subplot._kwargs
# #       sharex_datadict = sharex_subplot._data_dict
# #     else:
# #       sharex_kwargs = self._kwargs
# #       sharex_datadict = self._data_dict
# #     if 'sharey' in self._kwargs.keys():
# #       sp_i = subplot_locs.index(self._kwargs['sharey'])
# #       sharey_subplot = subplots[sp_i]
# #       sharey_kwargs = sharey_subplot._kwargs
# #       sharey_datadict = sharey_subplot._data_dict
# #     else:
# #       sharey_kwargs = self._kwargs
# #       sharey_datadict = self._data_dict
# #     if 'xlabel' in sharex_kwargs.keys():
# #       ax.set_xlabel(sharex_kwargs['xlabel'])
# #     if 'xlim' in sharex_kwargs.keys():
# #       ax.set_xlim(*sharex_kwargs['xlim'])
# #     elif 'x' in sharex_datadict.keys():
# #       ax.set_xlim(np.min(sharex_datadict['x']), np.max(sharex_datadict['x']))
# #     if 'ylabel' in sharey_kwargs.keys():
# #       ax.set_ylabel(sharey_kwargs['ylabel'])
# #     if 'ylim' in  sharey_kwargs.keys():
# #       ax.set_ylim(*sharey_kwargs['ylim'])
# #     elif 'y' in sharey_datadict.keys():
# #       ax.set_ylim(np.min(sharey_datadict['y']), np.max(sharey_datadict['y']))
# #     if 'title' in self._kwargs.keys():
# #       ax.set_title(self._kwargs['title'])
# #
# # class Plot3DLines(Subplot):
# #   def render(self, fig, ax, subplots=None, subplot_locs=None):
# #     for line in self._data_dict['lines']:
# #       ax.plot3D(*line.T, c=self._kwargs.get('c', 'blue'), lw=self._kwargs.get('lw', 0.2))
# #     self.commonVisuals(ax)
# #     if 'sharez' in self._kwargs.keys():
# #       sp_i = subplot_locs.index(self._kwargs['sharez'])
# #       sharez_subplot = subplots[sp_i]
# #       sharez_kwargs = sharez_subplot._kwargs
# #     else:
# #       sharez_kwargs = self._kwargs
# #     if 'zlabel' in sharez_kwargs.keys():
# #       ax.set_zlabel(sharez_kwargs['zlabel'])
# #     if 'zlim' in  sharez_kwargs.keys():
# #       ax.set_zlim(*sharez_kwargs['zlim'])
# #     import magic.aux as aux
# #     setaspect = self._kwargs.get('setaspect', (1,1,1))
# #     ax.get_proj = aux.make_get_proj(ax, *setaspect)
# #     if 'view' in self._kwargs.keys():
# #       ax.view_init(*self._kwargs['view'])
# #
# # class Plot2D(Subplot):
# #   def render(self, fig, ax, subplots=None, subplot_locs=None):
# #     import matplotlib.pyplot as plt
# #     import matplotlib as mpl
# #     self._kwargs['interpolation'] = self._kwargs.get('interpolation', 'gaussian')
# #     self._kwargs['cmap'] = self._kwargs.get('cmap', 'viridis')
# #     self._kwargs['cbar'] = self._kwargs.get('cbar', '2%')
# #     self._kwargs['norm'] = self._kwargs.get('norm', 'lin')
# #     self._kwargs['zmin'] = self._kwargs.get('zmin', self._data_dict['z'].min())
# #     self._kwargs['zmax'] = self._kwargs.get('zmax', self._data_dict['z'].max())
# #     self._kwargs['zoom'] = self._kwargs.get('zoom', None)
# #     # # # if shared Z axis
# #     if 'sharez' in self._kwargs.keys():
# #       sp_i = subplot_locs.index(self._kwargs['sharez'])
# #       sharez_subplot = subplots[sp_i]
# #       zdata_kwargs = sharez_subplot._kwargs
# #       zdata_dict = sharez_subplot._data_dict
# #     else:
# #       zdata_kwargs = self._kwargs
# #       zdata_dict = self._data_dict
# #     if zdata_kwargs['norm'] == 'lin':
# #       zmin = zdata_kwargs['zmin']
# #       zmax = zdata_kwargs['zmax']
# #       if (zmin * zmax < 0):
# #         zmax = max(abs(zmin), abs(zmax))
# #         zmin = -zmax; zmax = zmax
# #       norm = mpl.colors.Normalize(vmin=zmin, vmax=zmax)
# #     elif zdata_kwargs['norm'] == 'log':
# #       if (zdata_kwargs['zmin']) == 0:
# #         zdata_kwargs['zmin'] = zdata_kwargs['zmax'] / (zdata_kwargs['zmax'] / zdata_kwargs['z'].mean())**2
# #       norm = mpl.colors.LogNorm(vmin=zdata_kwargs['zmin'], vmax=zdata_kwargs['zmax'])
# #     elif zdata_kwargs['norm'] == 'symlog':
# #       zmax = max(abs(zdata_kwargs['zmin']), abs(zdata_kwargs['zmax']))
# #       zmean = np.abs(zdata_dict['z']).mean()
# #       norm = mpl.colors.SymLogNorm(vmin=-zmax, vmax=zmax,
# #                                    linthresh=zdata_kwargs.get('linthresh', zmax / (zmax / zmean)**2),
# #                                    linscale=zdata_kwargs.get('linscale', 1))
# #     cmap = zdata_kwargs['cmap']
# #     # # #
# #
# #     # # # if shared X/Y axes
# #     if 'sharex' in self._kwargs.keys():
# #       sp_i = subplot_locs.index(self._kwargs['sharex'])
# #       sharex_subplot = subplots[sp_i]
# #       xdata_kwargs = sharex_subplot._data_dict
# #     else:
# #       xdata_kwargs = self._data_dict
# #     if 'sharey' in self._kwargs.keys():
# #       sp_i = subplot_locs.index(self._kwargs['sharey'])
# #       sharey_subplot = subplots[sp_i]
# #       ydata_kwargs = sharey_subplot._data_dict
# #     else:
# #       ydata_kwargs = self._data_dict
# #
# #     self._data_dict['istep'] = xdata_kwargs.get('istep', 1)
# #     if 'x' in xdata_kwargs.keys():
# #       if len(xdata_kwargs['x']) == 2:
# #         xmin = xdata_kwargs['x'][0]
# #         xmax = xdata_kwargs['x'][1]
# #       elif xdata_kwargs['x'].shape == xdata_kwargs['z'].shape:
# #         xmin = xdata_kwargs['x'][0,0]
# #         xmax = xdata_kwargs['x'][0,-1]
# #         dx = xdata_kwargs['x'][0,1] - xmin
# #         xmax += dx
# #       elif len(xdata_kwargs['x'].shape) == 1:
# #         xmin = xdata_kwargs['x'][0]
# #         xmax = xdata_kwargs['x'][-1]
# #         dx = xdata_kwargs['x'][1] - xmin
# #         xmax += dx
# #     else:
# #       xmin = 0; xmax = xdata_kwargs['z'].shape[1] * self._data_dict['istep']
# #     if 'y' in ydata_kwargs.keys():
# #       if len(ydata_kwargs['y']) == 2:
# #         ymin = ydata_kwargs['y'][0]
# #         ymax = ydata_kwargs['y'][1]
# #       elif ydata_kwargs['y'].shape == ydata_kwargs['z'].shape:
# #         ymin = ydata_kwargs['y'][0,0]
# #         ymax = ydata_kwargs['y'][-1,0]
# #         dy = ydata_kwargs['y'][1,0] - ymin
# #         ymax += dy
# #       elif len(ydata_kwargs['y'].shape) == 1:
# #         ymin = ydata_kwargs['y'][0]
# #         ymax = ydata_kwargs['y'][-1]
# #         dy = ydata_kwargs['y'][1] - ymin
# #         ymax += dy
# #     else:
# #       ymin = 0; ymax = ydata_kwargs['z'].shape[0] * ydata_kwargs['istep']
# #     # # #
# #
# #     im = ax.imshow(self._data_dict['z'],
# #                    origin='lower', cmap=cmap,
# #                    extent=(xmin,xmax,ymin,ymax), norm=norm, interpolation=self._kwargs['interpolation'])
# #
# #     if zdata_kwargs['cbar'] is not None:
# #       from mpl_toolkits.axes_grid1 import make_axes_locatable
# #       divider = make_axes_locatable(ax)
# #       cax = divider.append_axes("right", size=self._kwargs['cbar'], pad=0.05)
# #       plt.colorbar(im, cax=cax)
# #     self.commonVisuals(ax, subplots, subplot_locs)
# #     if self._kwargs['zoom'] is not None:
# #       to, sx_sy = self._kwargs['zoom']
# #       sx, sy = sx_sy
# #       if to == 'center':
# #         x0 = (xmax + xmin) * 0.5
# #         y0 = (ymax + ymin) * 0.5
# #       else:
# #         x0, y0 = to
# #       if sx > 0:
# #         self._kwargs['xlim']=(x0 - sx / 2, x0 + sx / 2)
# #         ax.set_xlim(x0 - sx / 2, x0 + sx / 2)
# #       if sy > 0:
# #         self._kwargs['ylim']=(y0 - sy / 2, y0 + sy / 2)
# #         ax.set_ylim(y0 - sy / 2, y0 + sy / 2)
# #
# # class Streamplot(Subplot):
# #   def render(self, fig, ax, subplots=None, subplot_locs=None):
# #     # # # if shared X/Y axes
# #     if 'sharex' in self._kwargs.keys():
# #       sp_i = subplot_locs.index(self._kwargs['sharex'])
# #       sharex_subplot = subplots[sp_i]
# #       xdata_kwargs = sharex_subplot._data_dict
# #     else:
# #       xdata_kwargs = self._data_dict
# #     if 'sharey' in self._kwargs.keys():
# #       sp_i = subplot_locs.index(self._kwargs['sharey'])
# #       sharey_subplot = subplots[sp_i]
# #       ydata_kwargs = sharey_subplot._data_dict
# #     else:
# #       ydata_kwargs = self._data_dict
# #     self._data_dict['istep'] = xdata_kwargs.get('istep', 1)
# #     if 'x' in xdata_kwargs.keys():
# #       if len(xdata_kwargs['x']) == 2:
# #         xmin = xdata_kwargs['x'][0]
# #         xmax = xdata_kwargs['x'][1]
# #       elif xdata_kwargs['x'].shape == xdata_kwargs['z'].shape:
# #         xmin = xdata_kwargs['x'][0,0]
# #         xmax = xdata_kwargs['x'][0,-1]
# #         dx = xdata_kwargs['x'][0,1] - xmin
# #         xmax += dx
# #       elif len(xdata_kwargs['x'].shape) == 1:
# #         xmin = xdata_kwargs['x'][0]
# #         xmax = xdata_kwargs['x'][-1]
# #         dx = xdata_kwargs['x'][1] - xmin
# #         xmax += dx
# #     else:
# #       xmin = 0; xmax = xdata_kwargs['z'].shape[1] * self._data_dict['istep']
# #     if 'y' in ydata_kwargs.keys():
# #       if len(ydata_kwargs['y']) == 2:
# #         ymin = ydata_kwargs['y'][0]
# #         ymax = ydata_kwargs['y'][1]
# #       elif ydata_kwargs['y'].shape == ydata_kwargs['z'].shape:
# #         ymin = ydata_kwargs['y'][0,0]
# #         ymax = ydata_kwargs['y'][-1,0]
# #         dy = ydata_kwargs['y'][1,0] - ymin
# #         ymax += dy
# #       elif len(ydata_kwargs['y'].shape) == 1:
# #         ymin = ydata_kwargs['y'][0]
# #         ymax = ydata_kwargs['y'][-1]
# #         dy = ydata_kwargs['y'][1] - ymin
# #         ymax += dy
# #     else:
# #       ymin = 0; ymax = ydata_kwargs['z'].shape[0] * ydata_kwargs['istep']
# #     xs = np.linspace(xmin, xmax, self._data_dict['fx'].shape[1])
# #     ys = np.linspace(ymin, ymax, self._data_dict['fx'].shape[0])
# #     ax.streamplot(xs, ys, self._data_dict['fx'], self._data_dict['fy'],
# #                   **{k: self._kwargs[k] for k in self._kwargs.keys() if k in ['density', 'color', 'linewidth', 'arrowsize', 'start_points']})
# #     ax.set_aspect(1)
# #     self.commonVisuals(ax, subplots, subplot_locs)
# #
# # class Plot1D(Subplot):
# #   def render(self, fig, ax, subplots=None, subplot_locs=None):
# #     import matplotlib.pyplot as plt
# #     self._kwargs['c'] = self._kwargs.get('c', 'royalblue')
# #     self._kwargs['label'] = self._kwargs.get('label', None)
# #     self._kwargs['legend'] = self._kwargs.get('legend', None)
# #     self._kwargs['lw'] = self._kwargs.get('lw', 1)
# #     self._kwargs['ls'] = self._kwargs.get('ls', '-')
# #
# #     ax.plot(self._data_dict['x'], self._data_dict['y'], **{k: self._kwargs[k] for k in ['c', 'label', 'lw', 'ls']})
# #
# #     if self._kwargs['legend'] is not None:
# #       ax.legend(loc=self._kwargs['legend'])
# #     self.commonVisuals(ax, subplots, subplot_locs)
