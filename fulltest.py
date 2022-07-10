import sys
import os
import glob
import shutil
from abc import ABC, abstractmethod

from datetime import datetime
import argparse


def dir_path(string):
    if os.path.isdir(string):
        return string
    else:
        raise NotADirectoryError(string)


def clustername(string):
    if string in ['perseus', 'stellar']:
        return string
    else:
        raise ValueError(string)


parser = argparse.ArgumentParser()
parser.add_argument('--path', required=True, type=dir_path)
parser.add_argument('--cluster', type=clustername)
parser.add_argument('--dry', action='store_true',
                    default=False, help='dry run.')
parser.add_argument('-c', action='store_true',
                    default=False, help='only test compilation.')
parser.add_argument('-r', action='store_true',
                    default=False, help='only run simulations.')
parser.add_argument('-v', action='store_true', default=False,
                    help='verbose mode (full output).')
parser.add_argument('-d', '--diag', type=int, default=None,
                    help='diagnostic mode [choose test #].')
parser.add_argument('-t', '--test', type=str,
                    default=','.join(map(str, list(range(15)))), help='tests to run.')
options = parser.parse_args()
tests = [int(t) for t in options.test.split(',')]

suffix = '' if options.v else ' >/dev/null 2>&1'


def callCommand(command):
    if options.dry:
        print(f'$ {command}')
    else:
        os.system(command)


# modules for compilation
if (options.cluster == 'perseus'):
    modules = ['intel-mkl/2017.4/5/64',
               'intel/17.0/64/17.0.5.239',
               'intel-mpi/intel/2017.5/64',
               'hdf5/intel-17.0/intel-mpi/1.10.0']
elif (options.cluster == 'stellar'):
    modules = ['module load intel/2021.1.2',
               'module load intel-mpi/intel/2021.1.1',
               'module load hdf5/intel-2021.1/intel-mpi/1.10.6']

# global variables
outdir = options.path

codedir = os.getcwd()
testdir = 'test_0'
testdir_full = outdir + '/' + testdir

# simulation environment (tristan & slurm)


class Simulation(ABC):
    path = None
    exe = None
    exe_full = None
    submit = None
    submit_full = None
    input = None
    input_full = None

    def __init__(self, flags, walltime='00:30:00', nproc=28, params={}):
        self.walltime = walltime
        self.flags = flags
        self.params = params
        self.nproc = nproc

    @abstractmethod
    def diag(self, ax):
        pass


class TwoStream(Simulation):
    jobid = 'twostream'
    userfile = 'user_twostream'
    dimension = 1

    def diag(self):
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(8, 6))
        hist = isolde.parseHistory(self.path + '/output/history')
        omegap0 = self.params['algorithm']['c'] / \
            self.params['plasma']['c_omp']
        rate = 0.5 * (0.5)**0.5 / \
            (self.params['problem']['shift_gamma'])**(1.5)
        time = hist['time'] * omegap0
        E2 = hist['E^2'] / hist['Etot'][0]
        ax.plot(time, E2, label='sim')
        xs = np.linspace(10, 40, 10)
        ys = np.exp(2 * rate * xs)
        ys = (E2[time > 10][0]) * (ys / ys[0])
        ax.plot(xs, ys, label=r'$\omega/\omega_{b}^p=\gamma_b^{-3/2}$')
        ax.set_ylim(1e-4, 1e-1)
        ax.set_xlim(0, 200)
        ax.set_yscale('log')
        ax.set_xlabel(r'$t\omega_{\rm p0}$')
        ax.set_ylabel(r'$U_E / E_{\rm tot}$')
        ax.axvline(ax.get_xlim()[0], color='black')
        ax.axhline(ax.get_ylim()[0], color='black')
        ax.text(110, 0.6e-3, r"""
            energy conservation
            by $t\omega_{{\rm p 0}}={{{}}}: \Delta E/E={{{}}}\%$
            """.format(int(hist['time'][-1]*omegap0), int(hist['% dEtot'][-1]*10000) / 10000))
        ax.set_title(self.jobid)
        plt.legend()
        plt.tight_layout()
        plt.show()


class PlasmaOsc(Simulation):
    jobid = 'plasmaosc'
    userfile = 'user_langmuir'
    dimension = 1

    def diag(self):
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(8, 6))
        exs = []
        steps = np.arange(50)
        for step in steps:
            flds = isolde.getFields(self.path + '/output/flds.tot.%05d' % step)
            exs.append(flds['ex'][0, 0, int(
                self.params['grid']['mx0'] / 4)] / 1e-4)
        ax.plot(steps * self.params['output']['interval'] /
                (2 * np.pi * self.params['plasma']['c_omp'] / 0.45), exs)
        ax.set_xlabel(r'$t\omega_{\rm p0}$')
        ax.set_ylabel(r'$E_x$')
        ax.set_title(self.jobid)
        for i in range(4):
            ax.axvline(i, c='gray', lw=1, ls='--')
        hist = isolde.parseHistory(self.path + '/output/history')
        ax.text(2, -1, r"""
            energy conservation
            by $t\omega_{{\rm p 0}}={{{}}}: \Delta E/E={{{}}}%$
            """.format(int(hist['time'][-1]*0.45 / self.params['plasma']['c_omp']), int(hist['% dEtot'][-1]*10000) / 10000))
        plt.tight_layout()
        plt.show()


class Weibel(Simulation):
    jobid = 'weibel'
    userfile = 'user_weibel'
    dimension = 2

    def diag(self):
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(6, 6))
        ax.set_title(self.jobid)
        ax.grid(False)
        from matplotlib.animation import FuncAnimation
        flds = isolde.getFields(self.path + '/output/flds.tot.%05d' % 0)
        xmin = flds['xx'][0].min() / self.params['plasma']['c_omp']
        xmax = flds['xx'][0].max() / self.params['plasma']['c_omp']
        ymin = flds['yy'][0].min() / self.params['plasma']['c_omp']
        ymax = flds['yy'][0].max() / self.params['plasma']['c_omp']
        im = ax.imshow(flds['bz'][0], origin='lower', cmap='RdBu',
                       vmin=-0.01, vmax=0.01, extent=(xmin, xmax, ymin, ymax))
        ax.set_xlabel(r'$x/d_{e0}$')
        ax.set_ylabel(r'$y/d_{e0}$')
        txt1 = ax.text(10, 120, r'$B_z$', color='k')
        txt2 = ax.text(90, 120, r'', color='k', zorder=100)

        def init():
            im.set_data(flds['bz'][0])
            txt2.set_text(r'$t\omega_{{\rm p0}}=0$')
            return im, txt1, txt2,

        def animate(i):
            flds = isolde.getFields(self.path + '/output/flds.tot.%05d' % i)
            im.set_data(flds['bz'][0])
            txt2.set_text(r'$t\omega_{{\rm p0}}={{{}}}$'.format(
                i * self.params['output']['interval'] * 0.45 / self.params['plasma']['c_omp']))
            return im, txt1, txt2,
        anim = FuncAnimation(fig, animate, init_func=init,
                             frames=50, interval=500, blit=True, repeat=True)
        # anim.save('animation.mp4', fps=10)
        plt.tight_layout()
        plt.show()


class Merging(Simulation):
    jobid = 'merging'
    userfile = 'unit_chargedmerging'
    dimension = 2

    def diag(self):
        import matplotlib.pyplot as plt
        from matplotlib.animation import FuncAnimation
        import matplotlib.collections as mcoll
        from mpl_toolkits.axes_grid1 import make_axes_locatable
        divE_max = 1e-11
        fig, ax = plt.subplots(figsize=(6, 6))

        def getMomenta(prtls):
            en = np.sum(np.sqrt(1.0 + prtls['1']['u']**2 + prtls['1']['v']**2 + prtls['1']['w']**2) * prtls['1']['wei']) + np.sum(
                np.sqrt(1.0 + prtls['2']['u']**2 + prtls['2']['v']**2 + prtls['2']['w']**2) * prtls['2']['wei'])
            mx = np.sum(prtls['1']['u'] * prtls['1']['wei']) + \
                np.sum(prtls['2']['u'] * prtls['2']['wei'])
            my = np.sum(prtls['1']['v'] * prtls['1']['wei']) + \
                np.sum(prtls['2']['v'] * prtls['2']['wei'])
            mz = np.sum(prtls['1']['w'] * prtls['1']['wei']) + \
                np.sum(prtls['2']['w'] * prtls['2']['wei'])
            return (en, mx, my, mz)
        energy0, momx0, momy0, momz0 = getMomenta(
            isolde.getParticles(self.path + '/output/prtl.tot.%05d' % 0))

        def template(prtls):
            energy1, momx1, momy1, momz1 = getMomenta(prtls)
            return r'npart: {} ($\times$2)'.format(len(prtls['1']['x'])) + \
                '\ntotal energy [err%]: {:.3f} [{:.4f}%]'.format(energy1, np.abs((energy1 - energy0) * 100 / energy0)) + \
                '\nmomX [err%]: {:.3f} [{:.4f}%]'.format(momx1, np.abs((momx1 - momx0) * 100 / momx0)) + \
                '\nmomY [err%]: {:.3f} [{:.4f}%]'.format(momy1, np.abs((momy1 - momy0) * 100 / momy0)) + \
                '\nmomZ [err%]: {:.3f} [{:.4f}%]'.format(
                momz1, np.abs((momz1 - momz0) * 100 / (momz0 + 1e-10)))
            sc1 = ax.scatter([-100], [-100], fc='blue', label='lecs', zorder=2)
        sc2 = ax.scatter([-100], [-100], fc='red', label='ions', zorder=2)
        lgd = ax.legend(loc='lower right')
        major_ticks = np.arange(
            0, self.params['grid']['mx0'], self.params['grid']['tileX'])
        minor_ticks = np.arange(0, self.params['grid']['my0'], 1)
        ax.set_xticks(major_ticks)
        ax.set_xticks(minor_ticks, minor=True)
        ax.set_yticks(major_ticks)
        ax.set_yticks(minor_ticks, minor=True)
        ax.grid(False)
        lines1 = ([[(x, y) for y in (0, self.params['grid']['mx0'])] for x in major_ticks]
                  + [[(x, y) for x in (0, self.params['grid']['mx0'])] for y in major_ticks])
        grid1 = mcoll.LineCollection(
            lines1, linestyles='solid', linewidths=1, color='k', zorder=1, alpha=0.5)
        ax.add_collection(grid1)
        lines2 = ([[(x, y) for y in (0, self.params['grid']['mx0'])] for x in minor_ticks]
                  + [[(x, y) for x in (0, self.params['grid']['mx0'])] for y in minor_ticks])
        grid2 = mcoll.LineCollection(
            lines2, linestyles='solid', linewidths=0.2, color='k', zorder=1, alpha=0.3)
        ax.add_collection(grid2)
        im = ax.imshow(np.zeros((self.params['grid']['mx0'], self.params['grid']['my0'])), origin='lower', extent=(0, self.params['grid']['mx0'], 0,
                                                                                                                   self.params['grid']['my0']), cmap='seismic', vmin=-divE_max, vmax=divE_max, zorder=0)
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="2%", pad=0.02)
        plt.colorbar(im, cax=cax, label=r'$\nabla\cdot \bf{E}$')
        ttl = ax.text(.1, 0.1, '', transform=ax.transAxes, va='center')
        diag = ax.text(2, 47, '', multialignment='left', va='top', bbox=dict(
            facecolor='white', edgecolor='gray', boxstyle='round,pad=0.2'))

        def init():
            ax.set_xlim(0, self.params['grid']['mx0'])
            ax.set_ylim(0, self.params['grid']['my0'])
            ax.set_xlabel(r'$x$')
            ax.set_ylabel(r'$y$')
            ax.set_title(self.jobid)
            prtls = isolde.getParticles(
                self.path + '/output/prtl.tot.%05d' % 0)
            flds = isolde.getFields(self.path + '/output/flds.tot.%05d' % 0)
            sc1.set_offsets(np.array([prtls['1']['x'], prtls['1']['y']]).T)
            sc2.set_offsets(np.array([prtls['2']['x'], prtls['2']['y']]).T)
            im.set_data(flds['divE'][0])
            ttl.set_text('t=0')
            diag.set_text(template(prtls))
            return sc1, sc2, im, ttl, diag, lgd, grid1, grid2,

        def animate(i):
            prtls = isolde.getParticles(
                self.path + '/output/prtl.tot.%05d' % i)
            flds = isolde.getFields(self.path + '/output/flds.tot.%05d' % i)
            sc1.set_offsets(np.array([prtls['1']['x'], prtls['1']['y']]).T)
            sc2.set_offsets(np.array([prtls['2']['x'], prtls['2']['y']]).T)
            im.set_data(flds['divE'][0])
            ttl.set_text('t=' + str(i))
            diag.set_text(template(prtls))
            return sc1, sc2, im, ttl, diag, lgd, grid1, grid2,
        anim = FuncAnimation(fig, animate, init_func=init,
                             frames=40, interval=250, blit=True, repeat=True)
        plt.tight_layout()
        plt.show()


class AdaptiveLB(Simulation):
    jobid = 'alb'
    userfile = 'unit_alb'
    dimension = 3

    def diag(self):
        import matplotlib.pyplot as plt
        from matplotlib.animation import FuncAnimation
        fig = plt.figure(figsize=(8, 8))
        sx = self.params['grid']['mx0']
        sy = self.params['grid']['my0']
        sz = self.params['grid']['mz0']
        sizes = {
            'x': sx,
            'y': sy,
            'z': sz
        }

        ax1 = plt.subplot(221)
        ax2 = plt.subplot(222)
        ax3 = plt.subplot(223)
        ax4 = plt.subplot(224, projection='3d')

        sc4_1 = ax4.scatter([-100], [-100], fc='C0')
        sc4_2 = ax4.scatter([-100], [-100], fc='C0')
        ax4.set_xlim(0, sx)
        ax4.set_ylim(0, sy)
        ax4.set_zlim(0, sz)
        ax4.set_xlabel('x')
        ax4.set_ylabel('y')
        ax4.set_zlabel('z')

        def draw(i):
            prtls = isolde.getParticles(
                self.path + '/output/prtl.tot.%05d' % i)
            flds = isolde.getFields(self.path + '/output/flds.tot.%05d' % i)
            dom = isolde.getDomains(self.path + '/output/domain.%05d' % i)
            # xz
            ims = []
            for ax, coords, avr in zip([ax1, ax2, ax3], ['xz', 'xy', 'yz'], [1, 0, 2]):
                ax.cla()
                ax.grid(False)
                s1 = coords[0]
                s2 = coords[1]
                ax.set_xlim(0, sizes[s1])
                ax.set_ylim(0, sizes[s2])
                ax.set_xlabel(s1)
                ax.set_ylabel(s2)

                im = ax.imshow(np.zeros((sizes[s1], sizes[s2])), vmin=0, vmax=300,
                               origin='lower', extent=(0, sizes[s1], 0, sizes[s2]), cmap='hot')
                for si in np.unique(dom[s1 + '0']):
                    ax.axvline(si, c='C0', lw=0.5)
                for si in np.unique(dom[s2 + '0']):
                    ax.axhline(si, c='C0', lw=0.5)

                dens = np.sum(flds['dens1'] + flds['dens2'], axis=avr)
                im.set_data(dens)
                ims.append(im)
            # scatter
            sc4_1._offsets3d = (
                prtls['1']['x'], prtls['1']['y'], prtls['1']['z'])
            sc4_2._offsets3d = (
                prtls['2']['x'], prtls['2']['y'], prtls['2']['z'])
            return ax1, ax2, ax3, ax4, ims, sc4_1, sc4_2,

        def init():
            return draw(0)

        def animate(i):
            return draw(i)
        anim = FuncAnimation(fig, animate, init_func=init,
                             frames=40, interval=250, blit=False, repeat=True)
        # anim.save('animation.mp4', fps=10)
        plt.tight_layout()
        plt.show()


# Here specify the test simulations and give additional specs of the environment
if options.cluster:
    common_flags = ' --cluster={} -hdf5 --debug=1'.format(options.cluster)
else:
    common_flags = ' -hdf5'
simulations = [
    TwoStream(common_flags,
              params={
                  'node_configuration': {'sizex': 28},
                  'time': {'last': 50000},
                  'grid': {'mx0': 5600, 'tileX': 100},
                  'algorithm': {'c': 0.35, 'nfilter': 8},
                  'output': {'enable': 0, 'hst_enable': 1, 'hst_interval': 5},
                  'plasma': {'ppc0': 64, 'sigma': 10, 'c_omp': 40},
                  'particles': {'nspec': 2, 'maxptl1': 1e5, 'm1': 1, 'ch1': -1, 'maxptl2': 1e5, 'm2': 1, 'ch2': -1},
                  'problem': {'shift_gamma': 2.5}}
              ),
    PlasmaOsc(common_flags, nproc=4,
              params={
                  'node_configuration': {'sizex': 4},
                  'time': {'last': 1000},
                  'grid': {'mx0': 1120, 'tileX': 10},
                  'algorithm': {'nfilter': 8},
                  'output': {'enable': 1, 'prtl_enable': 0, 'spec_enable': 0, 'interval': 10, 'hst_enable': 1, 'hst_interval': 5},
                  'plasma': {'ppc0': 500, 'sigma': 1, 'c_omp': 10},
                  'particles': {'nspec': 2, 'maxptl1': 1e8, 'm1': 1, 'ch1': -1, 'maxptl2': 1e8, 'm2': 1, 'ch2': 1},
                  'problem': {'upstream_T': 1e-5, 'amplitude': 0.01, 'nwaves': 1}}
              ),
    Weibel(common_flags, nproc=28,
           params={
               'node_configuration': {'sizex': 7, 'sizey': 4},
               'time': {'last': 500},
               'grid': {'mx0': 1400, 'my0': 1400, 'tileX': 20, 'tileY': 20},
               'algorithm': {'nfilter': 8},
               'output': {'enable': 1, 'prtl_enable': 0, 'spec_enable': 0, 'interval': 10, 'hst_enable': 1, 'hst_interval': 10},
               'plasma': {'ppc0': 16, 'sigma': 10, 'c_omp': 10},
               'particles': {'nspec': 2, 'maxptl1': 1e8, 'm1': 1, 'ch1': -1, 'maxptl2': 1e8, 'm2': 1, 'ch2': 1},
               'problem': {'backgr_T': 1e-5, 'shift_beta': 0.5}}
           ),
    Merging(common_flags + ' -dwn', nproc=1,
            params={
                'node_configuration': {'sizex': 1, 'sizey': 1},
                'time': {'last': 500},
                'grid': {'mx0': 50, 'my0': 50, 'tileX': 10, 'tileY': 10},
                'algorithm': {'nfilter': 0},
                'output': {'interval': 10, 'stride': 1, 'smooth_window': 0, 'write_nablas': 1},
                'plasma': {'ppc0': 50, 'sigma': 5, 'c_omp': 50},
                'particles': {'nspec': 2, 'maxptl1': 1e8, 'm1': 1, 'ch1': -1, 'dwn1': 1, 'maxptl2': 1e8, 'm2': 1, 'ch2': 1, 'dwn2': 1},
                'downsampling': {'interval': 1, 'start': 1, 'max_weight': 1e5, 'cartesian_bins': 1, 'energy_min': 0, 'energy_max': 1e5, 'int_weights': 0, 'dynamic_bins': 1, 'mom_bins': 1, 'mom_spread': 1e5}
            }
            ),
    AdaptiveLB(common_flags + ' -alb', nproc=192,
               params={
                   'node_configuration': {'sizex': 8, 'sizey': 8, 'sizez': 3},
                   'time': {'last': 400},
                   'grid': {'mx0': 112, 'my0': 112, 'mz0': 114, 'resize_tiles': 1, 'tileX': 5, 'tileY': 5, 'tileZ': 5},
                   'algorithm': {'nfilter': 0, 'fieldsolver': 0, 'currdeposit': 0},
                   'output': {'enable': 1, 'interval': 10, 'stride': 10, 'istep': 1, 'diag_enable': 1},
                   'adaptive_load_balancing': {'in_x': 1, 'in_y': 1, 'in_z': 1,
                                               'sx_min': 8, 'sy_min': 8, 'sz_min': 8,
                                               'interval_x': 5, 'interval_y': 5, 'interval_z': 5,
                                               'slab_x': 2, 'slab_y': 2, 'slab_z': 2
                                               },
                   'plasma': {'ppc0': 16, 'sigma': 5, 'c_omp': 10},
                   'particles': {'nspec': 2, 'maxptl1': 1e6, 'm1': 1, 'ch1': 1, 'maxptl2': 1e6, 'm2': 1, 'ch2': -1},
                   'problem': {'radius': 5}
               }
               )
]

if (options.diag):
    # diagnostic mode where you analize the test results
    import matplotlib.pyplot as plt
    import numpy as np
    import tristanVis.isolde as isolde
    plt.style.use('fivethirtyeight')
    for ii, simulation in enumerate(simulations):
        if (ii + 1 == options.diag):
            if (options.dry):
                print(f'Plotting for `{simulation.jobid}`')
            else:
                simulation.path = f'{testdir_full}/{ii + 1:02}_{simulation.jobid}'
                simulation.diag()
                break
else:
    # regular mode where you compile and run tests
    if not os.path.exists(testdir_full):
        os.makedirs(testdir_full)
    with open(testdir_full + '/test.log', 'w+') as testlog:
        # load modules
        for ii, simulation in enumerate(simulations):
            if not (ii + 1 in tests):
                continue
            # create directory for simulation
            simulation.path = f'{testdir_full}/{ii + 1:02}_{simulation.jobid}'
            if os.path.exists(simulation.path) and not options.r:
                callCommand(f'rm -r {simulation.path}')
            if not options.r and not options.dry:
                os.makedirs(simulation.path)
            testlog.write(('TEST_#{}_'.format(ii+1) +
                           simulation.jobid).ljust(50, '.') + '\n')
            # configure
            config_command = 'python3 configure.py'
            config_command += simulation.flags
            config_command += ' -{}d'.format(simulation.dimension)
            config_command += ' --{}='.format(
                simulation.userfile[:4]) + simulation.userfile
            if (not options.r):
                # compile everything if not in `r` mode
                callCommand(config_command + suffix)
                # clean
                callCommand('make clean' + suffix)
                # compile
                callCommand('make all -j ' + suffix)
                # check if compilation successfull
                simulation.exe = f'tristan-mp{simulation.dimension}d'
                simulation.exe_full = f'{simulation.path}/{simulation.exe}'
                if (not os.path.isfile(codedir + '/bin/' + simulation.exe)) and not options.dry:
                    raise RuntimeError("Something went wrong")
                testlog.write('compilation'.ljust(46, '.') + '[OK]\n')
                print('Compilation of `{}` done.'.format(simulation.jobid))
                # move executable
                callCommand(
                    f'mv {codedir}/bin/{simulation.exe} {simulation.path}')
                # clean
                callCommand('make clean' + suffix)
                # write input
                simulation.input = 'input.' + simulation.jobid
                simulation.input_full = simulation.path + '/' + simulation.input
                if not options.dry:
                    with open(simulation.input_full, 'w+') as inp:
                        for block in simulation.params.keys():
                            inp.write('\n<{}>\n\n'.format(block))
                            for var in simulation.params[block].keys():
                                inp.write('  {}  =  {}\n'.format(
                                    var, simulation.params[block][var]))
                                testlog.write(
                                    'input file'.ljust(46, '.') + '[OK]\n')
            # write submit
            simulation.submit = 'submit_' + simulation.jobid
            simulation.submit_full = simulation.path + '/' + simulation.submit
            if (not options.r) and (not options.dry) and (options.cluster):
                with open(simulation.submit_full, 'w+') as sub:
                    sub.write('#!/bin/bash\n')
                    sub.write('#SBATCH -t {}\n'.format(simulation.walltime))
                    sub.write('#SBATCH -n {}\n'.format(simulation.nproc))
                    sub.write('#SBATCH -J {}\n'.format(simulation.jobid))
                    sub.write(
                        '#SBATCH --output={}/tristan-v2.out\n'.format(simulation.path))
                    sub.write(
                        '#SBATCH --error={}/tristan-v2.err\n\n'.format(simulation.path))
                    sub.write('DIR={}\n'.format(simulation.path))
                    sub.write('EXECUTABLE=$DIR/{}\n'.format(simulation.exe))
                    sub.write('INPUT=$DIR/{}\n'.format(simulation.input))
                    sub.write(
                        'OUTPUT_DIR=$DIR/output\n'.format(simulation.path))
                    sub.write('SLICE_DIR=$DIR/slices\n'.format(simulation.path))
                    sub.write(
                        'REPORT_FILE=$DIR/report\n'.format(simulation.path))
                    sub.write(
                        'ERROR_FILE=$DIR/error\n\n'.format(simulation.path))
                    for module in modules:
                        sub.write('module load {}\n'.format(module))
                    sub.write('\nmkdir $OUTPUT_DIR\n\n')
                    sub.write(
                        'srun $EXECUTABLE -i $INPUT -o $OUTPUT_DIR -s $SLICE_DIR -r $RESTART_DIR -R $RESTART > $REPORT_FILE 2> $ERROR_FILE')
                    testlog.write('submit file'.ljust(46, '.') + '[OK]\n')
            testlog.write('\n\n')
        if (not options.c) or (options.r):
            for ii, simulation in enumerate(simulations):
                if not (ii + 1 in tests):
                    continue
                callCommand('sbatch ' + simulation.submit_full)
                testlog.write(('`{}`'.format(simulation.jobid)
                               ).ljust(41, '.') + 'submitted\n')
