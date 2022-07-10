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
