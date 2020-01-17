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

def getSpectra(fname):
    with h5py.File(fname, 'r') as file:
        keys = list(file.keys())
        spectra = [key[1:] for key in keys if key.startswith("n")]
        data = {}
        for sp in spectra:
            data[sp] = {}
            (data[sp])['bn'] = np.exp(file['e' + sp][:])
            (data[sp])['cnt'] = file['n' + sp][:]
    return data

def getDomains(fname):
    with h5py.File(fname, 'r') as file:
        data = {}
        for k in file.keys():
            data[k] = file[k][:]
    return data

def parseReport(fname, nsteps = None, skip = 1):
    if (not nsteps):
        nsteps = 1e6
    import re
    def parseBlock(block, data, isfirst = False):
        for line in block.split('\n')[2:]:
            routine = line.split(':', 1)[0].strip()
            if (routine != ''):
                line1 = line.split(':', 1)[1]
                if (isfirst):
                    data[routine] = {}
                    data[routine]['dt'] = np.array([])
                    data[routine]['min'] = np.array([])
                    data[routine]['max'] = np.array([])
                nums = [float(x.strip()) for x in re.findall(re.compile('-?\.? *[0-9]+\.?[0-9]*(?:[Ee]\ *[-+]?\ *[0-9]+)?'), line1)]
                if (len(nums) < 3):
                    raise ValueError('len(nums) < 3')
                else:
                    data[routine]['dt'] = np.append(data[routine]['dt'], [nums[0]])
                    data[routine]['min'] = np.append(data[routine]['min'], [nums[1]])
                    data[routine]['max'] = np.append(data[routine]['max'], [nums[2]])
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
            if (ni % skip == 0) and line:
                parseBlock(block, data, isfirst = isfirst)
                isfirst = False
                data['t'] = np.append(data['t'], [ni])
            ni += 1
    return data

# easy plotting functions
def plot2DField(ax, x, y, field, rotate=False,
                title='field', cmap='jet',
                vmin=None, vmax=None,
                scale='lin', region=[-np.inf, np.inf, -np.inf, np.inf],
                cbar = '2%', cbar_pad=0.05,
                **kwargs):
    import matplotlib.pyplot as plt
    import matplotlib as mpl
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    if rotate:
      field = np.rot90(field)
      x_dummy = 1.0 * np.array(x)
      x = 1.0 * np.array(y)
      y = 1.0 * np.array(x_dummy)
    xmin = x[x > region[0]].min()
    xmax = x[x <= region[1]].max()
    ymin = y[y > region[2]].min()
    ymax = y[y <= region[3]].max()
    ax.set_aspect(1)
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    if not vmin:
        vmin = field.min()
    if not vmax:
        vmax = field.max()
    if scale == 'lin':
      norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
    elif scale == 'log':
      vmax = max(vmax, 1e-10)
      norm = mpl.colors.LogNorm(vmin=max(vmin, vmax/1e10), vmax=vmax)
    elif scale == 'sym':
      vmax = max(np.abs(vmin), vmax)
      norm=mpl.colors.SymLogNorm(vmin=-vmax, vmax=vmax,
                                 linthresh=kwargs['lth'],
                                 linscale=kwargs['lsc'])

    im = ax.imshow(field, norm=norm,
                   cmap=cmap, origin='lower',
                   extent=(x.min(), x.max(), y.min(), y.max()))

    if 'xlabel' in kwargs:
        ax.set_xlabel(kwargs['xlabel'])
    else:
        ax.set_xlabel('x' if not rotate else 'y')
    if 'ylabel' in kwargs:
        ax.set_ylabel(kwargs['ylabel'])
    else:
        ax.set_ylabel('y' if not rotate else 'x')
    if cbar is not None:
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size=cbar, pad=cbar_pad)
        plt.colorbar(im, cax=cax)
    ax.set_title(title)


def plot2DScatterParticles(ax, x_list, y_list,
                           label='particles', legend=True,
                           color='black', **kwargs):
    ax.scatter(x_list, y_list, c=color, label=label, **kwargs)
    ax.set_aspect(1)
    if legend:
        ax.legend()

def plot2DDomains(ax, domain_data,
                  color='red', **kwargs):
    from matplotlib.patches import Rectangle
    x0_list = domain_data['x0']
    y0_list = domain_data['y0']
    sx_list = domain_data['sx']
    sy_list = domain_data['sy']
    for x0, y0, sx, sy in zip(x0_list, y0_list, sx_list, sy_list):
        rect = Rectangle((x0, y0), sx, sy,
                         edgecolor=color, facecolor='none')
        ax.add_patch(rect)

def plotReport(ax, data, only_fullstep = True, **kwargs):
    labels = []
    y_list = []
    for k in list(data.keys())[1:]:
        y_list.append
        if (k.split()[0] != 'nprt') and (((k.split()[0] != 'Full_step') and (not only_fullstep)) or (only_fullstep and k.split()[0] == 'Full_step')):
            labels.append(k)
            y_list.append(data[k]['dt'])
    y = np.vstack(y_list)
    if (only_fullstep and 'label' in kwargs):
        labels = [kwargs['label']]
    if (only_fullstep):
        ax.plot(data['t'], y[0], label = labels[0])
    else:
        ax.stackplot(data['t'], y, labels = labels)
    ax.legend()
    ax.set_ylabel(r'$\Delta t$ [ms]')
    ax.set_xlabel(r'timestep')
