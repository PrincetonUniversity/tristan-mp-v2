import tristanVis.aux as aux

from contextlib import contextmanager
import ipywidgets as ipyW
from IPython.display import display

class FieldData2D():
  def __init__(self, root, step, slices, isSlice, mask,
               extraVariables, coordinateTransformation):
    self._root = root
    self._step = step
    self.slices = slices
    self._isSlice = isSlice
    self._mask = mask
    self._extraVariables = extraVariables
    self._coordinateTransformation = coordinateTransformation

    self.data = {}
    self.axes = []

  def maskData(self, mask):
    self._mask = mask

  def addVariables(self, variables):
    if (self._extraVariables is None):
      self._extraVariables = {}
    for k in variables.keys():
      if (not k in self._extraVariables.keys()):
        self._extraVariables.update({k: variables[k]})

  def addCoordinateTransformation(self, transform):
    for k in transform.keys():
      self._coordinateTransformation.update({k: transform[k]})

  def loadSlice(self, slice):
    import xarray as xr
    import h5py

    coord, value = slice.split('=')
    xr_slice = 'xyz'.replace(coord, '')
    xr_axes = tuple(xr_slice)[::-1]
    self.axes.append(xr_axes)
    if (not self._isSlice):
      fname = self._root + 'flds.tot.%05d' % self._step
    else:
      slice_name = coord.upper() + '=%05d' % int(value) + '.%05d' % self._step
      fname = self._root + 'slices/slice' + slice_name
    with h5py.File(fname, 'r') as fields:
      xr_data = xr.Dataset()
      for k in fields.keys():
        if (not self._isSlice):
          xr_data[k] = (xr_axes, fields[k][:][0])
        else:
          xr_data[k] = (xr_axes, fields[k][:])
      x1, x2 = xr_axes
      if (not self._isSlice):
        xr_data.coords[x1] = ((x1), self._coordinateTransformation[x1](fields[x1*2][:][0,:,0]))
        xr_data.coords[x2] = ((x2), self._coordinateTransformation[x2](fields[x2*2][:][0,0,:]))
      else:
        xr_data.coords[x1] = ((x1), self._coordinateTransformation[x1](fields[x1*2][:][:,0]))
        xr_data.coords[x2] = ((x2), self._coordinateTransformation[x2](fields[x2*2][:][0,:]))
      if (self._extraVariables is not None):
        for k in self._extraVariables.keys():
          xr_data[k] = self._extraVariables[k](xr_data)
      for k in xr_data.keys():
        if not self._mask is None:
          xr_data[k] = xr_data[k].where(self._mask(xr_data))
      self.data[slice] = xr_data

  def loadData(self):
    if self._step is not None:
      for s in self.slices:
        self.loadSlice(s)

class Simulation():
  def __init__(self,
               root, fld_steps=None, useSlices=False, mask=None,
               coordinateTransformation={'x': lambda x: x, 'y': lambda y: y, 'z': lambda z: z},
               extraVariables=None
              ):
    if (root[-1] != '/'):
      root += '/'
    self._root = root

    self.fields = {}
    self.spectra = {}
    self.particles = {}

    self._mask = mask
    self._fld_steps = fld_steps
    self._useSlices = useSlices
    self._extraVariables = extraVariables
    self._coordinateTransformation = coordinateTransformation

  def __del__(self):
    del self.fields
    del self.spectra
    del self.particles

  def maskData(self, mask):
    self._mask = mask

  def addVariables(self, variables):
    if (self._extraVariables is None):
      self._extraVariables = {}
    for k in variables.keys():
      if (not k in self._extraVariables.keys()):
        self._extraVariables.update({k: variables[k]})
        [fld.addVariables(variables) for st, fld in self.fields.items()]

  def addCoordinateTransformation(self, transform):
    for k in transform.keys():
      self._coordinateTransformation.update({k: transform[k]})
      [fld.addCoordinateTransformation(transform) for st, fld in self.fields.items()]

  def readFiles(self):
    import numpy as np
    # parse directory and record files
    files = aux.listFiles(self._root)
    files = files[['spec' in file for file in files]]
    self._spec_steps = np.sort([file[-5:] for file in files])
    self._spec_steps = [int(step) for step in self._spec_steps]

    if (self._useSlices):
      files = aux.listFiles(self._root + 'slices/')
      self._slices = np.unique([file[:-6] for file in files if (not 'xdmf' in file)])
      self._slices = np.array(['='.join((lambda x: [x[0], str(int(x[1]))])(sl.lower()[5:].split('='))) for sl in self._slices])
    else:
      files = aux.listFiles(self._root)
      self._slices = np.array(['z=0'])
    if self._fld_steps is None:
      self._slice_steps = np.sort(np.unique([file[-5:] for file in files]))
      self._slice_steps = [int(step) for step in self._slice_steps]
    else:
      self._slice_steps = np.sort(self._fld_steps)
    # preload all the files
    # TODO: preload spectra

    for st in self._slice_steps:
      fld = FieldData2D(self._root, st, self._slices, self._useSlices, self._mask,
                        extraVariables=self._extraVariables,
                        coordinateTransformation=self._coordinateTransformation)
      self.fields.update({st: fld})

  def loadData(self):
    self.readFiles()
    [fld.loadData() for st, fld in self.fields.items()]

class FieldPlot2D(ipyW.VBox):
  def __init__(self, sim, **kwargs):
    import matplotlib.pyplot as plt
    super().__init__()
    self._kwargs = kwargs
    self.simulation = sim
    self._kwargs['timestep'] = self._kwargs.get('timestep', list(self.simulation.fields.keys())[0])

    projections = list(self.simulation.fields[self._kwargs['timestep']].data.keys())
    variables = list(self.simulation.fields[self._kwargs['timestep']].data[projections[0]].keys())
    self._kwargs['var'] = self._kwargs.get('var', variables[0])
    self._kwargs['proj'] = self._kwargs.get('proj', projections[0])

    self._kwargs['cmap'] = self._kwargs.get('cmap', 'viridis')
    self._kwargs['vmin'] = self._kwargs.get('vmin', None)
    self._kwargs['vmax'] = self._kwargs.get('vmax', None)
    self._kwargs['logplot'] = self._kwargs.get('logplot', False)
    self._kwargs['controls'] = self._kwargs.get('controls', True)
    self._kwargs['figsize'] = self._kwargs.get('figsize', (6, 4))
    self._kwargs['zoomQ'] = self._kwargs.get('zoomQ', False)
    self._kwargs['interpolation'] = self._kwargs.get('interpolation', None)

    if (self._kwargs['vmin'] is None):
      self._kwargs['vmin'], _ = self.findMinMax()
    if (self._kwargs['vmax'] is None):
      _, self._kwargs['vmax'] = self.findMinMax()

    if self._kwargs['controls']:
      self.obj_var = ipyW.Dropdown(options=variables, description='var:',
                                   value=self._kwargs['var'], layout={'width': 'max-content'});
      self.obj_proj = ipyW.Dropdown(options=projections, description='proj:',
                                    value=self._kwargs['proj'], layout={'width': 'max-content'});
      self.obj_cmap = ipyW.Dropdown(options=plt.colormaps(), description='cmap:',
                                    value=self._kwargs['cmap'], layout={'width': 'max-content'});

      self.obj_minval = ipyW.FloatText(description='min:', value=self._kwargs['vmin'], layout={'width': '200px'})
      self.obj_maxval = ipyW.FloatText(description='max:', value=self._kwargs['vmax'], layout={'width': '200px'})

      self.obj_logplot = ipyW.Checkbox(value=self._kwargs['logplot'], description='logplot')

      self.controls = ipyW.Box([
          ipyW.VBox([self.obj_var, self.obj_minval, self.obj_maxval]),
          ipyW.VBox([self.obj_proj, self.obj_cmap, self.obj_logplot]),
      ])
      self.obj_cmap.value = self._kwargs['cmap']
      self.obj_cmap.observe(self.update_cmap, 'value')
      self.obj_var.observe(self.update_var, 'value')
      self.obj_proj.observe(self.update_proj, 'value')
      self.obj_logplot.observe(self.update_logplot, 'value')
      self.obj_minval.observe(self.update_minval, 'value')
      self.obj_maxval.observe(self.update_maxval, 'value')
    output = ipyW.Output()
    with output:
      self.fig, self.ax = plt.subplots(figsize=self._kwargs['figsize'])
    self.fig.canvas.toolbar_visible = self._kwargs['zoomQ']
    self.fig.canvas.header_visible = False
    self.fig.canvas.footer_visible = False

    self.generatePlot()
    if self._kwargs['controls']:
      self.children = [self.controls, output]
    else:
      self.children = [output]

  def __del__(self):
    del self.fig
    del self.children

  def findMinMax(self):
    import numpy as np
    data_ = self.simulation.fields[self._kwargs['timestep']].data[self._kwargs['proj']][self._kwargs['var']].values
    vmin = np.nanmin(data_[(data_ != -np.inf) & (data_ != np.inf)])
    vmax = np.nanmax(data_[(data_ != -np.inf) & (data_ != np.inf)])
    if (vmin * vmax < 0):
      vv = max(np.abs(vmin), np.abs(vmax))
      vmin = -vv; vmax = vv
    elif self._kwargs['logplot']:
      vmin = vmax / 1e5
    return (vmin, vmax)

  def findNorm(self):
    import matplotlib as mpl
    if (self._kwargs['logplot']):
      if (self._kwargs['vmin'] * self._kwargs['vmax'] < 0):
        norm_ = mpl.colors.SymLogNorm(vmin=self._kwargs['vmin'], vmax=self._kwargs['vmax'],
                                     linthresh=self._kwargs['vmax'] / 1e3, linscale=1, base=10)
      else:
        if (self._kwargs['vmin'] * self._kwargs['vmax'] == 0):
          self.autoMinMax(maxval=self._kwargs['vmax'])
        norm_ = mpl.colors.LogNorm(vmin=self._kwargs['vmin'], vmax=self._kwargs['vmax'])
    else:
      norm_ = mpl.colors.Normalize(vmin=self._kwargs['vmin'], vmax=self._kwargs['vmax'])
    return norm_

  def generatePlot(self):
    import matplotlib.pyplot as plt
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    data_ = self.simulation.fields[self._kwargs['timestep']].data[self._kwargs['proj']][self._kwargs['var']].values
    coords_ = self.simulation.fields[self._kwargs['timestep']].data[self._kwargs['proj']].coords
    x2_, x1_ = list(coords_)
    x1min_ = coords_[x1_].values.min()
    x1max_ = coords_[x1_].values.max()
    x2min_ = coords_[x2_].values.min()
    x2max_ = coords_[x2_].values.max()
    norm_ = self.findNorm()
    self.im = self.ax.imshow(data_, origin='lower',
                        cmap=self._kwargs['cmap'],
                        norm=norm_, interpolation=self._kwargs['interpolation'],
                        extent=(x1min_, x1max_, x2min_, x2max_))
    try:
      self.fig.delaxes(self.fig.axes[1])
    except:
      ...
    divider = make_axes_locatable(self.ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    self.cbar = self.fig.colorbar(self.im, cax=cax)
    self.cbar.set_label(self._kwargs['var'].replace('_', '\_'))
    crds_ = list(coords_.keys())
    self.ax.set_xlabel(crds_[1])
    self.ax.set_ylabel(crds_[0])
    self.ax.set_aspect(1)
    plt.tight_layout()

  def update_var(self, change):
    import difflib
    import re
    oldvar = self._kwargs['var']
    newvar = change.new
    self._kwargs['var'] = change.new
    data_ = self.simulation.fields[self._kwargs['timestep']].data[self._kwargs['proj']][self._kwargs['var']].values
    self.im.set_data(data_)
    if re.compile("^[x,y,z]+$").match(''.join(sorted([li[-1] for li in difflib.ndiff(oldvar, newvar) if li[0] != ' ']))) is None:
      self.autoMinMax()
      self.im.set_clim(vmin=self._kwargs['vmin'], vmax=self._kwargs['vmax'])
    self.cbar.set_label(self._kwargs['var'].replace('_', '\_'))

  def update_logplot(self, change):
    self._kwargs['logplot'] = change.new
    self.im.set_norm(self.findNorm())

  def update_proj(self, change):
    self._kwargs['proj'] = change.new
    data_ = self.simulation.fields[self._kwargs['timestep']].data[self._kwargs['proj']][self._kwargs['var']].values
    coords_ = self.simulation.fields[self._kwargs['timestep']].data[self._kwargs['proj']].coords
    self.im.set_data(data_)
    crds_ = list(coords_.keys())
    self.ax.set_xlabel(crds_[1])
    self.ax.set_ylabel(crds_[0])

  def update_cmap(self, change):
    self._kwargs['cmap'] = change.new
    self.im.set_cmap(self._kwargs['cmap'])

  def update_minval(self, change):
    self._kwargs['vmin'] = change.new
    self.im.set_clim(vmin=self._kwargs['vmin'], vmax=self._kwargs['vmax'])
    if ((self._kwargs['vmin'] < 0) and self._kwargs['logplot']):
      norm_ = self.findNorm()
      self.im.set_norm(norm_)

  def update_maxval(self, change):
    self._kwargs['vmax'] = change.new
    self.im.set_clim(vmin=self._kwargs['vmin'], vmax=self._kwargs['vmax'])
    if ((self._kwargs['vmin'] < 0) and self._kwargs['logplot']):
      norm_ = self.findNorm()
      self.im.set_norm(norm_)

  def update_timestep(self, new_timestep):
    self._kwargs['timestep'] = new_timestep
    data_ = self.simulation.fields[self._kwargs['timestep']].data[self._kwargs['proj']][self._kwargs['var']].values
    self.im.set_data(data_)

  def autoMinMax(self, maxval=None, minval=None):
    mn, mx = self.findMinMax()
    if (maxval is None):
      self._kwargs['vmax'] = mx
    if (minval is None):
      self._kwargs['vmin'] = mn
    self.im.set_norm(self.findNorm())
    if self._kwargs['controls']:
      self.obj_minval.value = self._kwargs['vmin']
      self.obj_maxval.value = self._kwargs['vmax']

class PlotGrid():
  def __init__(self, simulation, maxncols=2, zoom=False, interpolation=None, timestep=None, init=[], controls=True, figsize=None):
    self.simulation = simulation
    self.parameters = init
    self.controls = controls
    self.figsize = figsize
    self.zoomQ = zoom
    self.interpolation = interpolation
    self.panels = []
    self.maxncols = maxncols
    if self.parameters != []:
      try:
        self.figsize = self.parameters[0]['figsize']
      except:
        ...

    self.plotgrid = ipyW.Box()

    self.addPlot_button = ipyW.Button(description='Add field')
    self.addPlot_button.on_click(self.addPanel)

    # self.button2 = ipyW.Button(description='Next [%d] >' % self.simulation.step)
    # self.button3 = ipyW.Button(description='Save .png')
    # self.button3.on_click(self.savePng)

    timesteps = list(self.simulation.fields.keys())

    try:
      newvalue = self.parameters[0]['timestep']
    except:
      newvalue = timestep if (not timestep is None) else timesteps[0]
    self.timestep = newvalue
    try:
      dtimestep = timesteps[-1] - timesteps[-2]
    except:
      dtimestep = 0

    # self.button2.on_click(self.nextTimestep)
    self.step_slider = ipyW.IntSlider(min=min(timesteps),
                                      max=max(timesteps),
                                      step=dtimestep, value=self.timestep, layout={'width': '100%'})
    self.step_slider.observe(self.changeTimestep, names="value")

    self.button_panel = ipyW.HBox([self.addPlot_button, self.step_slider], layout={'margin': '0px 0px 20px 0px'})
    self.generateGrid()

  def __del__(self):
    del self.simulation
    del self.plotgrid
    del self.panels

  def addPanel(self, b):
    if (self.controls):
      self.parameters.append({})
      self.generateGrid()

  def changeTimestep(self, change):
    self.timestep = change.new
    [panel.update_timestep(change.new) for panel in self.panels]
    # self.simulation.step = change.new
    # self.redraw()

  def savePng(self, filename):
    NN, ncols, nrows = self.getNxM()
    import os
    import numpy as np
    import shutil
    import PIL
    from PIL import Image
    temp = 'temp_'
    if not os.path.exists(temp):
      os.mkdir(temp)
    filenames = []
    for ii, subplot in enumerate(self.plotgrid.children[0].children[1].children):
      fname = temp + '/pic_%03d.png'%ii
      filenames.append(fname)
      subplot.fig.savefig(fname)
    imgs = [PIL.Image.open(i) for i in filenames]
    min_shape = sorted([(np.sum(i.size), i.size ) for i in imgs])[0][1]
    img_rows = []
    for nr in range(nrows):
      img_row = np.hstack([np.asarray(i.resize(min_shape)) for i in imgs[nr * ncols : nr * ncols + ncols]])
      img_rows.append(img_row)
    imgs_comb = PIL.Image.fromarray(np.vstack(img_rows))
    shutil.rmtree(temp)
    imgs_comb.save(filename)

  def getNxM(self):
    import numpy as np
    NN = len(self.parameters)
    if (NN > 1):
      ncols = self.maxncols
    elif NN == 0:
      ncols = 0
    else:
      ncols = 1
    if (ncols > 0):
      nrows = int(np.ceil(NN / ncols))
    else:
      nrows = 0
    return (NN, ncols, nrows)

  def generateGrid(self):
    import matplotlib.pyplot as plt
    plt.close('all')
    NN, ncols, nrows = self.getNxM()
    if (nrows * ncols > 0):
      grid = ipyW.GridspecLayout(nrows, ncols)
      self._oldpanels = []
      n = 0
      for i in range(grid.n_rows):
        for j in range(grid.n_columns):
          if (n < NN):
            try:
              self.parameters[n] = self.panels[n]._kwargs
            except:
              pass
            if not self.figsize is None:
              self.parameters[n]['figsize'] = self.figsize
            self.parameters[n]['controls'] = self.controls
            self.parameters[n]['timestep'] = self.timestep
            self.parameters[n]['zoomQ'] = self.zoomQ
            self.parameters[n]['interpolation'] = self.interpolation
            panel = FieldPlot2D(self.simulation, **self.parameters[n])
            self._oldpanels.append(panel)
            grid[i, j] = panel
          n += 1
      self.panels = self._oldpanels
      self.plotgrid.children = [ipyW.VBox([self.button_panel, grid])]
    else:
      self.plotgrid.children = [ipyW.VBox([self.button_panel])]

  def exportParameters(self, filename=None):
    params = []
    for panel in self.panels:
      params.append(panel._kwargs)
    if not (filename is None):
      import json
      import numpy as np
      class jsonEncoder(json.JSONEncoder):
        def default(self, obj):
          if isinstance(obj, np.integer):
            return int(obj)
          elif isinstance(obj, np.floating):
            return float(obj)
          elif isinstance(obj, np.ndarray):
            return obj.tolist()
          else:
            return super(MyEncoder, self).default(obj)
      if filename[-5:] != '.json':
        filename = filename + '.json'
      with open(filename, 'w') as outfile:
        json.dump(params, outfile, cls=jsonEncoder)
      print ('data written to', filename)
    else:
      return params

  def readParameters(self, params=None, filename=None):
    if (not (params is None)):
      self.parameters = params
    elif (not (filename is None)):
      import json
      if filename[-5:] != '.json':
        filename = filename + '.json'
      with open(filename) as infile:
        self.parameters = json.load(infile)
    try:
      newvalue = self.parameters[0]['timestep']
    except:
      newvalue = 0
    self.step_slider.value = newvalue
    self.generateGrid()
    #     ...
    #   # add spectra

  def show(self):
    display(self.plotgrid)
