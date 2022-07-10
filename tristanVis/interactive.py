from tristanVis.algorithms import slice1D

import xarray as xr
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

class TwoDPlotAndSlice:
  def __init__(self, field_data):
    coords = list(field_data.coords.keys())
    self.x1 = field_data.coords[coords[0]]; self.x2 = field_data.coords[coords[1]]
    self.field_data = field_data
    self.inside = False
    self.xy_1 = (None, None)
    self.xy_2 = (None, None)
    self.mouse_hold = False
  def in_axes(self, event):
    self.inside = (event.inaxes == self.fig.axes[0])
  def leave_axes(self, event):
    self.inside = False
  def on_press(self, event):
    if (self.inside):
      self.mouse_hold = True
      self.xy_1 = (event.xdata, event.ydata)
  def on_release(self, event):
    def interpolate(s, p1, p2):
      return (np.array(p1) + s * (np.array(p2) - np.array(p1)), np.linalg.norm(s * (np.array(p2) - np.array(p1))))
    if (self.inside):
      self.xy_2 = (event.xdata, event.ydata)
      ss = np.linspace(0, 1, 100)
      pnts = np.array([interpolate(s, self.xy_1, self.xy_2)[0] for s in ss])
      ys = slice1D(self.field_data, pnts)
      xs = np.array([interpolate(s, self.xy_1, self.xy_2)[1] for s in ss])
      self.slice.set_xdata(xs)
      self.slice.set_ydata(ys)
      try:
        (self.slice_kwargs['xlim'] != None)
      except:
        self.ax_slice.set_xlim(xs.min(), xs.max())
      try:
        (self.slice_kwargs['ylim'] != None)
      except:
        self.ax_slice.set_ylim(ys.min(), ys.max())
      self.mouse_hold = False
      self.fig.canvas.draw()
      self.fig.canvas.flush_events()
  def on_mousemove(self, event):
    if (self.inside and self.mouse_hold):
      start = self.xy_1
      end = (event.xdata, event.ydata)
      self.line.set_xdata([start[0], end[0]])
      self.line.set_ydata([start[1], end[1]])
  def plot(self, figsize=(8, 4), pcolor_kwargs={}, line_kwargs={}, slice_kwargs={}):
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    self.pcolor_kwargs = pcolor_kwargs
    self.line_kwargs = line_kwargs
    self.slice_kwargs = slice_kwargs
    self.fig = plt.figure(figsize=figsize)
    ax_fld = self.fig.add_subplot(121)
    ax_fld.set_aspect(1)
    self.pcolor = ax_fld.pcolormesh(self.x1, self.x2, self.field_data, **self.pcolor_kwargs)
    im = self.pcolor
    divider = make_axes_locatable(ax_fld)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    plt.colorbar(im, cax=cax)
    self.ax_slice = self.fig.add_subplot(122)
    self.slice, = self.ax_slice.plot([], [])
    try:
      self.ax_slice.set_xlim(*self.slice_kwargs['xlim'])
    except:
      self.ax_slice.set_xlim(0, 1)
    try:
      self.ax_slice.set_ylim(*self.slice_kwargs['ylim'])
    except:
      self.ax_slice.set_ylim(0, 1)
    try:
      self.ax_slice.set_xlabel(self.slice_kwargs['xlabel'])
    except:
      self.ax_slice.set_xlabel(r'$s$')
    try:
      self.ax_slice.set_ylabel(self.slice_kwargs['ylabel'])
    except:
      self.ax_slice.set_ylabel(r'field')
    self.line, = ax_fld.plot([], [], **self.line_kwargs)
    in1 = self.fig.canvas.mpl_connect('axes_enter_event', self.in_axes)
    in2 = self.fig.canvas.mpl_connect('axes_leave_event', self.leave_axes)
    mouse1 = self.fig.canvas.mpl_connect('button_press_event', self.on_press)
    mouse2 = self.fig.canvas.mpl_connect('button_release_event', self.on_release)
    mouse3 = self.fig.canvas.mpl_connect('motion_notify_event', self.on_mousemove)
    plt.tight_layout()
    plt.show()

class TwoDPlotAndSpectra:
  """
  Plot 2D field plot and particle distribution function side-by-side. Clicking on the field plot will select particular regions for the particle distribution plot. Holding `Shift` while clicking will add the region to an already selected one.

  ...

  Methods
  -------
  plot(figsize=(12, 4)):
    Create the interactive plot.

  Example
  -------
    myplot = TwoDPlotAndSlice(...)
    myplot.plot()

  """
  def __init__(self, 
      x, y, z, field_data, spec_data, 
      coordinates=None, 
      imshow_kwargs={}, rectangle_kwargs={}, spectra_kwargs={}):
    """
    Initializer for the `TwoDPlotAndSpectra` class.

    Args
    ----
      x, y, z : N-D arrays
        array of coordinates for the field

      field_data : 2D array
        field data to plot

      spec_data : `Spectra`
        particle distribution data of type `Spectra`

      coordinates : str, optional
        2D plane of the field data, can be either of the following: `'xy'`, `'yz'`, `'xz'` (default is `None`)

      imshow_kwargs : dict, optional
        keyword arguments that go into `plt.imshow` (default is {})

      rectangle_kwargs : dict, optional
        properties of rectangles drawn over the selected region (default is {})

      spectra_kwargs : dict, optional
        keyword arguments for the spectral plot (default is {})
    """
    self.x = x; self.y = y; self.z = z
    self.field_data = field_data
    self.spec_data = spec_data
    self.coordinates = coordinates
    self.imshow_kwargs = imshow_kwargs
    self.rectangle_kwargs = rectangle_kwargs
    self.spectra_kwargs = spectra_kwargs
    self.spectra_kwargs.setdefault('smooth')
    self.spectra_kwargs.setdefault('xlim')
    self.spectra_kwargs.setdefault('ylim')
    self.SHIFT = False
    self.specbins = []
  def on_key_press(self, event):
    if event.key == 'shift':
      self.SHIFT = True
  def on_key_release(self, event):
    if event.key == 'shift':
      self.SHIFT = False
  def on_click(self, event):
    import matplotlib.patches as mpatches
    if (self.coordinates == 'xy'):
      xi = event.xdata; yi = event.ydata; zi = self.z.mean()
    elif (self.coordinates == 'xz'):
      xi = event.xdata; yi = self.y.mean(); zi = event.ydata
    elif (self.coordinates == 'yz'):
      xi = self.x.mean(); yi = event.xdata; zi = event.ydata
    else:
      raise ValueError('Wrong coordinate system')
    # add species #1 and #2
    xyz, sxyz, cnt1 = self.spec_data.getByCoordinate(1, (xi, yi, zi))
    xyz, sxyz, cnt2 = self.spec_data.getByCoordinate(2, (xi, yi, zi))
    cnt = cnt1 + cnt2
    if (self.coordinates == 'xy'):
      x0, y0 = (xyz[0], xyz[1])
      sx, sy = (sxyz[0], sxyz[1])
    elif (self.coordinates == 'xz'):
      x0, y0 = (xyz[0], xyz[2])
      sx, sy = (sxyz[0], sxyz[2])
    elif (self.coordinates == 'yz'):
      x0, y0 = (xyz[1], xyz[2])
      sx, sy = (sxyz[1], sxyz[2])
    else:
      raise ValueError('Wrong coordinate system')
    bns = np.copy(self.spec_data.ebins)
    if not self.spectra_kwargs['smooth'] is None:
      from scipy.ndimage import gaussian_filter1d as smooth
      cnt1 = smooth(np.copy(cnt1), self.spectra_kwargs['smooth'])
      cnt2 = smooth(np.copy(cnt2), self.spectra_kwargs['smooth'])
      cnt = smooth(np.copy(cnt), self.spectra_kwargs['smooth'])
    else:
      cnt1 = np.copy(cnt1)
      cnt2 = np.copy(cnt2)
      cnt = np.copy(cnt)
    if ((not self.SHIFT) or (self.specbins == [])):
      self.specbins = [(x0, y0)]
      self.spec_x = bns
      self.spec_y = (cnt1, cnt2, cnt)
      for i in range(3):
        self.spec[i].set_xdata(self.spec_x)
        self.spec[i].set_ydata(self.spec_y[i])
      while len(self.fig.axes[0].patches) != 0:
        self.fig.axes[0].patches[0].remove()
      self.fig.axes[0].add_patch(mpatches.Rectangle((x0, y0), sx, sy, fc='None', zorder=100, **self.rectangle_kwargs))
    else:
      if (not (x0, y0) in self.specbins):
        self.specbins.append((x0, y0))
        self.spec_x = bns
        try:
          self.spec_y = (self.spec_y[0] + cnt1, self.spec_y[1] + cnt2, self.spec_y[2] + cnt) 
        except:
          self.spec_y = (cnt1, cnt2, cnt)
        for i in range(3):
          self.spec[i].set_xdata(self.spec_x)
          self.spec[i].set_ydata(self.spec_y[i])
        self.fig.axes[0].add_patch(mpatches.Rectangle((x0, y0), sx, sy, fc='None', zorder=100, **self.rectangle_kwargs))
      else:
        self.specbins.remove((x0, y0))
        self.spec_y = (self.spec_y[0] - cnt1, self.spec_y[1] - cnt2, self.spec_y[2] - cnt) 
        for i in range(3):
          self.spec[i].set_xdata(self.spec_x)
          self.spec[i].set_ydata(self.spec_y[i])
        for rect in self.fig.axes[0].patches:
          if (rect.get_xy() == (x0, y0)):
            rect.remove()
            break
    self.fig.canvas.draw()
    self.fig.canvas.flush_events()

  def plot(self, figsize=(12, 4)):
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    self.fig = plt.figure(figsize=figsize)
    if (self.coordinates is None):
      if self.x.min() == self.x.max():
        self.coordinates = 'yz'
      elif self.y.min() == self.y.max():
        self.coordinates = 'xz'
      elif self.z.min() == self.z.max():
        self.coordinates = 'xy'

    if self.coordinates == 'yz':
      xmin = self.y.min(); xmax = self.y.max()
      ymin = self.z.min(); ymax = self.z.max()
    elif self.coordinates == 'xz':
      xmin = self.x.min(); xmax = self.x.max()
      ymin = self.z.min(); ymax = self.z.max()
    elif self.coordinates == 'xy':
      xmin = self.x.min(); xmax = self.x.max()
      ymin = self.y.min(); ymax = self.y.max()
    else:
      raise

    self.extent = (xmin, xmax, ymin, ymax)
    ax_fld = self.fig.add_subplot(121)
    im = ax_fld.imshow(self.field_data, origin='lower', extent=self.extent, **self.imshow_kwargs)
    ax_fld.set_xlabel(self.coordinates[0])
    ax_fld.set_ylabel(self.coordinates[1])
    ax_spec = self.fig.add_subplot(122)
    sp1, = ax_spec.plot([], [], label='1')
    sp2, = ax_spec.plot([], [], label='2')
    sp, = ax_spec.plot([], [], label='tot')
    self.spec = (sp1, sp2, sp)
    if not self.spectra_kwargs['xlim'] is None:
      ax_spec.set_xlim(*self.spectra_kwargs['xlim'])
    if not self.spectra_kwargs['ylim'] is None:
      ax_spec.set_ylim(*self.spectra_kwargs['ylim'])
    ax_spec.set_xscale('log')
    ax_spec.set_yscale('log')
    ax_spec.set_xlabel(r'$\gamma - 1$')
    ax_spec.set_ylabel(r'$(\gamma - 1)f(\gamma)$')
    plt.legend();
    divider = make_axes_locatable(ax_fld)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    plt.colorbar(im, cax=cax)
    rect = mpatches.Rectangle((0, 0), 0, 0, fc='None', zorder=100, **self.rectangle_kwargs)
    ax_fld.add_patch(rect)
    cid = self.fig.canvas.mpl_connect('button_press_event', self.on_click)
    self.fig.canvas.mpl_connect('key_press_event', self.on_key_press)
    self.fig.canvas.mpl_connect('key_release_event', self.on_key_release)
    plt.tight_layout()
    plt.show()

class ThreeDPlotAndSlice:
  """
  Plot 2D field plot (from 3D data) and a custom 2D slice side-by-side. Clicking and dragging on the field plot will select particular line across which to take the slice.

  ...

  Methods
  -------
  plot(figsize=(12, 4)):
    Create the interactive plot.

  Example
  -------
    myplot = ThreeDPlotAndSlice(...)
    myplot.plot()

  """
  def __init__(self, field_data):
    self.field_data = field_data
    self.inside = False
    self.p1 = (None, None)
    self.p2 = (None, None)
    self.mouse_hold = False
  def in_axes(self, event):
    self.inside = (event.inaxes == self.fig.axes[0])
  def leave_axes(self, event):
    self.inside = False
  def on_press(self, event):
    if (self.inside):
      self.mouse_hold = True
      self.p1 = (event.xdata, event.ydata)
  def on_release(self, event):
    def interpolate(s, p1, p2):
      return (np.array(p1) + s * (np.array(p2) - np.array(p1)), np.linalg.norm(s * (np.array(p2) - np.array(p1))))
    if (self.inside):
      self.p2 = (event.xdata, event.ydata)
      res = self.resolution
      theta = np.arctan2((np.array(self.p2) - np.array(self.p1))[0], (np.array(self.p2) - np.array(self.p1))[1])
      ss = np.linspace(0, 1, res)
      x1x2s = np.array([interpolate(s, self.p1, self.p2)[0] for s in ss])
      x12s = np.array([interpolate(s, self.p1, self.p2)[1] for s in ss])
      ds = interpolate(1, self.p1, self.p2)[1] / res
      x3s = np.arange(*self.hidden_lims, ds)
      # if (self.proj == 'xy'):
        # sliced_bx = np.array([
          # [self.interpolated_bx((x3, *x1x2[::-1])) for x1x2 in x1x2s
          # ] for x3 in x3s
        # ])
        # sliced_by = np.array([
          # [self.interpolated_by((x3, *x1x2[::-1])) for x1x2 in x1x2s
          # ] for x3 in x3s
        # ])
        # sliced_bxy = sliced_bx * np.sin(theta) + sliced_by * np.cos(theta)
        # sliced_bz = np.array([
          # [self.interpolated_bz((x3, *x1x2[::-1])) for x1x2 in x1x2s
          # ] for x3 in x3s
        # ])
      # else:
        # raise ValueError('Unkonwn projection')
      sliced_field = np.array([
        [self.interpolated_field((x3, *x1x2[::-1])) for x1x2 in x1x2s
        ] for x3 in x3s
      ])
      while self.ax_slice.collections != []:
        self.ax_slice.collections[0].remove()
      while self.ax_slice.patches != []:
        self.ax_slice.patches[0].remove()
      # self.ax_slice.streamplot(x12s, x3s, sliced_bxy, sliced_bz, **self.streamplot_kwargs)
      self.im.set_data(sliced_field)
      self.im.set_extent((x12s.min(), x12s.max(), *self.hidden_lims))
      self.mouse_hold = False
      self.fig.canvas.draw()
      self.fig.canvas.flush_events()
  def on_mousemove(self, event):
    if (self.inside and self.mouse_hold):
      start = self.p1
      end = (event.xdata, event.ydata)
      self.line.set_xdata([start[0], end[0]])
      self.line.set_ydata([start[1], end[1]])
  def plot(self, proj, shift, field,
           hidden_lims=None,
           transform=None,
           resolution=100,
           figsize=(8, 3),
           pcolor_kwargs={},
           line_kwargs={},
           slice_kwargs={},
           streamplot_kwargs={}
          ):
    from scipy.interpolate import RegularGridInterpolator
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    self.field = field(self.field_data)
    self.pcolor_kwargs = pcolor_kwargs
    self.line_kwargs = line_kwargs
    self.slice_kwargs = slice_kwargs
    self.streamplot_kwargs = streamplot_kwargs
    self.fig = plt.figure(figsize=figsize)
    self.shift = shift
    self.resolution = resolution
    self.proj = proj
    self.hidden_lims = hidden_lims
    ax_fld = self.fig.add_subplot(121)
    ax_fld.set_aspect(1)
    if (self.proj == 'xy'):
      if (transform is not None):
        xyz = transform(self.field_data['xx'], self.field_data['yy'], self.field_data['zz'])
        if (self.hidden_lims is None):
          self.hidden_lims = (xyz[2].min(), xyz[2].max())
        x1 = xyz[0][self.shift,:,:]
        x2 = xyz[1][self.shift,:,:]
      else:
        if (self.hidden_lims is None):
          self.hidden_lims = (self.field_data['zz'].min(), self.field_data['zz'].max())
        x1 = self.field_data['xx'][self.shift,:,:]
        x2 = self.field_data['yy'][self.shift,:,:]
      fld = self.field[self.shift,:,:]
    else:
      raise ValueError('Unkonwn projection')
    self.pcolor = ax_fld.pcolormesh(x1, x2, fld, shading='auto',
                                    **self.pcolor_kwargs)
    if (transform is not None):
      xyz = transform(self.field_data['xx'], self.field_data['yy'], self.field_data['zz'])
      self.interpolated_field = RegularGridInterpolator((xyz[2][:,0,0],
                                                         xyz[1][0,:,0],
                                                         xyz[0][0,0,:]), self.field)
      self.interpolated_bx = RegularGridInterpolator((xyz[2][:,0,0],
                                                      xyz[1][0,:,0],
                                                      xyz[0][0,0,:]), self.field_data['bx'])
      self.interpolated_by = RegularGridInterpolator((xyz[2][:,0,0],
                                                      xyz[1][0,:,0],
                                                      xyz[0][0,0,:]), self.field_data['by'])
      self.interpolated_bz = RegularGridInterpolator((xyz[2][:,0,0],
                                                      xyz[1][0,:,0],
                                                      xyz[0][0,0,:]), self.field_data['bz'])
    else:
      self.interpolated_field = RegularGridInterpolator((self.field_data['zz'][:,0,0],
                                                         self.field_data['yy'][0,:,0],
                                                         self.field_data['xx'][0,0,:]), self.field)
    im = self.pcolor
    divider = make_axes_locatable(ax_fld)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    plt.colorbar(im, cax=cax)
    self.line, = ax_fld.plot([], [], **self.line_kwargs)

    self.ax_slice = self.fig.add_subplot(122)
    self.im = self.ax_slice.imshow([[0]], origin='lower',
                              **self.slice_kwargs)
    self.ax_slice.set_xlabel(r'$xy$')
    self.ax_slice.set_ylabel(r'$z$')
    divider = make_axes_locatable(self.ax_slice)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    plt.colorbar(self.im, cax=cax)

    in1 = self.fig.canvas.mpl_connect('axes_enter_event', self.in_axes)
    in2 = self.fig.canvas.mpl_connect('axes_leave_event', self.leave_axes)
    mouse1 = self.fig.canvas.mpl_connect('button_press_event', self.on_press)
    mouse2 = self.fig.canvas.mpl_connect('button_release_event', self.on_release)
    mouse3 = self.fig.canvas.mpl_connect('motion_notify_event', self.on_mousemove)
    plt.tight_layout()
    plt.show()
