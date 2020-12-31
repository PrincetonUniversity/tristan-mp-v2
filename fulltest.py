#!/usr/bin/env python3
import sys
import os
import glob
import shutil
from abc import ABC, abstractmethod

from datetime import datetime
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-c', action='store_true', default=False, help='only test compilation.')
parser.add_argument('-v', action='store_true', default=False, help='verbose mode (full output).')
parser.add_argument('-d', action='store_true', default=False, help='diagnostic mode.')
parser.add_argument('-t','--test', type=str, default=','.join(map(str, list(range(15)))), help='tests to run.')
options = parser.parse_args()
tests = [int(t) for t in options.test.split(',')]

suffix = '' if options.v else ' >/dev/null 2>&1'

# modules for compilation
modules = ['intel-mkl/2017.4/5/64',
            'intel/17.0/64/17.0.5.239',
            'intel-mpi/intel/2017.5/64',
            'hdf5/intel-17.0/intel-mpi/1.10.0']

# global variables
outdir = "/scratch/gpfs/hakobyan/tristan_2/test"

codedir = os.getcwd()
testdir = 'test_1'
# testdir = 'test_' + datetime.now().strftime("%H.%M_%d.%m.%Y")
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
  def diag(self, ax, fig=None):
    hist = isolde.parseHistory(self.path + '/output/history')
    omegap0 = self.params['algorithm']['c'] / self.params['plasma']['c_omp']
    rate = 0.5 * (0.5)**0.5 / (self.params['problem']['shift_gamma'])**(1.5)
    time = hist['time'] * omegap0; E2 = hist['E^2'] / hist['Etot'][0]; ax.plot(time, E2, label='sim')
    xs = np.linspace(10, 40, 10); ys = np.exp(2 * rate * xs); ys = (E2[time>10][0]) * (ys / ys[0]); ax.plot(xs, ys, label=r'$\omega/\omega_{\textrm{p}b}=\gamma_b^{-3/2}$')
    ax.set_ylim(1e-4, 1e-1); ax.set_xlim(0, 200); ax.set_yscale('log');
    ax.set_xlabel(r'$t\omega_{\rm p0}$'); ax.set_ylabel(r'$U_E / E_{\rm tot}$')
    ax.axvline(ax.get_xlim()[0], color='black'); ax.axhline(ax.get_ylim()[0], color='black')
    ax.text(110, 0.6e-3, r'energy conservation\\by $t\omega_{{\rm p 0}}={{{}}}: \Delta E/E={{{}}}\%$'.format(int(hist['time'][-1]*omegap0), int(hist['% dEtot'][-1]*10000) / 10000), bbox=dict(facecolor='white', edgecolor='gray', boxstyle='round,pad=0.5'))
    ax.set_title(self.jobid); plt.legend()

class PlasmaOsc(Simulation):
  jobid = 'plasmaosc'
  userfile = 'user_langmuir'
  dimension = 1
  def diag(self, ax, fig=None):
    exs = []; steps = np.arange(50)
    for step in steps:
      flds = isolde.getFields(self.path + '/output/flds.tot.%05d' % step)
      exs.append(flds['ex'][0, 0, int(self.params['grid']['mx0'] / 4)] / 1e-4)
    ax.plot(steps * self.params['output']['interval'] / (2 * np.pi * self.params['plasma']['c_omp'] / 0.45), exs)
    ax.set_xlabel(r'$t\omega_{\rm p0}$'); ax.set_ylabel(r'$E_x$')
    ax.set_title(self.jobid)
    for i in range(4):
      ax.axvline(i, c='gray', lw=1, ls='--')
    hist = isolde.parseHistory(self.path + '/output/history')
    ax.text(2, -1, r'energy conservation\\by $t\omega_{{\rm p 0}}={{{}}}: \Delta E/E={{{}}}\%$'.format(int(hist['time'][-1]*0.45 / self.params['plasma']['c_omp']), int(hist['% dEtot'][-1]*10000) / 10000), bbox=dict(facecolor='white', edgecolor='gray', boxstyle='round,pad=0.5'))

class Weibel(Simulation):
  jobid = 'weibel'
  userfile = 'user_weibel'
  dimension = 2
  def diag(self, ax, fig=None):
    ax.set_title(self.jobid); ax.grid(False)
    from matplotlib.animation import FuncAnimation
    flds = isolde.getFields(self.path + '/output/flds.tot.%05d' % 0)
    xmin = flds['xx'][0].min() / self.params['plasma']['c_omp']; xmax = flds['xx'][0].max() / self.params['plasma']['c_omp']
    ymin = flds['yy'][0].min() / self.params['plasma']['c_omp']; ymax = flds['yy'][0].max() / self.params['plasma']['c_omp']
    im = ax.imshow(flds['bz'][0], origin='lower', cmap='bipolar', vmin=-0.01, vmax=0.01, extent=(xmin,xmax,ymin,ymax))
    ax.set_xlabel(r'$x/d_{e0}$'); ax.set_ylabel(r'$y/d_{e0}$')
    txt1 = ax.text(10, 120, r'$B_z$', color='white')
    txt2 = ax.text(90, 120, r'', color='white', zorder=100)
    def init():
      im.set_data(flds['bz'][0])
      txt2.set_text(r'$t\omega_{{\rm p0}}=0$')
      return im, txt1, txt2,
    def animate(i):
      flds = isolde.getFields(self.path + '/output/flds.tot.%05d' % i); im.set_data(flds['bz'][0]);
      txt2.set_text(r'$t\omega_{{\rm p0}}={{{}}}$'.format(i * self.params['output']['interval'] * 0.45 / self.params['plasma']['c_omp']))
      return im, txt1, txt2,
    anim = FuncAnimation(fig, animate, init_func=init, frames=50, interval=1000, blit=True, repeat=True)

class Merging(Simulation):
  jobid = 'merging'
  userfile = 'unit_chargedmerging'
  dimension = 2
  def diag(self, ax, fig=None):
    from matplotlib.animation import FuncAnimation
    import matplotlib.collections as mcoll
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    divE_max = 1e-11;
    def getMomenta(prtls):
      en = np.sum(np.sqrt(1.0 + prtls['1']['u']**2 + prtls['1']['v']**2 + prtls['1']['w']**2) * prtls['1']['wei']) + np.sum(np.sqrt(1.0 + prtls['2']['u']**2 + prtls['2']['v']**2 + prtls['2']['w']**2) * prtls['2']['wei'])
      mx = np.sum(prtls['1']['u'] * prtls['1']['wei']) + np.sum(prtls['2']['u'] * prtls['2']['wei'])
      my = np.sum(prtls['1']['v'] * prtls['1']['wei']) + np.sum(prtls['2']['v'] * prtls['2']['wei'])
      mz = np.sum(prtls['1']['w'] * prtls['1']['wei']) + np.sum(prtls['2']['w'] * prtls['2']['wei'])
      return (en, mx, my, mz)
    energy0, momx0, momy0, momz0 = getMomenta(isolde.getParticles(self.path + '/output/prtl.tot.%05d' % 0))
    def template(prtls):
      energy1, momx1, momy1, momz1 = getMomenta(prtls)
      return r'npart: {} ($\times$2)'.format(len(prtls['1']['x'])) + \
              '\ntotal energy [err\%]: {:.3f} [{:.4f}\%]'.format(energy1, np.abs((energy1 - energy0) * 100 / energy0)) + \
              '\nmomX [err\%]: {:.3f} [{:.4f}\%]'.format(momx1, np.abs((momx1 - momx0) * 100 / momx0)) + \
              '\nmomY [err\%]: {:.3f} [{:.4f}\%]'.format(momy1, np.abs((momy1 - momy0) * 100 / momy0)) + \
              '\nmomZ [err\%]: {:.3f} [{:.4f}\%]'.format(momz1, np.abs((momz1 - momz0) * 100 / (momz0 + 1e-10)))
    sc1 = ax.scatter([-100], [-100], fc='blue', label='lecs', zorder=2)
    sc2 = ax.scatter([-100], [-100], fc='red', label='ions', zorder=2)
    lgd = ax.legend(loc='lower right');
    major_ticks = np.arange(0, self.params['grid']['mx0'], self.params['grid']['tileX']); minor_ticks = np.arange(0, self.params['grid']['my0'], 1)
    ax.set_xticks(major_ticks); ax.set_xticks(minor_ticks, minor=True)
    ax.set_yticks(major_ticks); ax.set_yticks(minor_ticks, minor=True)
    ax.grid(False);
    lines1 = ([[(x, y) for y in (0, self.params['grid']['mx0'])] for x in major_ticks]
             + [[(x, y) for x in (0, self.params['grid']['mx0'])] for y in major_ticks])
    grid1 = mcoll.LineCollection(lines1, linestyles='solid', linewidths=1, color='k', zorder=1, alpha=0.5); ax.add_collection(grid1)
    lines2 = ([[(x, y) for y in (0, self.params['grid']['mx0'])] for x in minor_ticks]
             + [[(x, y) for x in (0, self.params['grid']['mx0'])] for y in minor_ticks])
    grid2 = mcoll.LineCollection(lines2, linestyles='solid', linewidths=0.2, color='k', zorder=1, alpha=0.3); ax.add_collection(grid2)
    im = ax.imshow(np.zeros((self.params['grid']['mx0'], self.params['grid']['my0'])), origin='lower', extent=(0, self.params['grid']['mx0'], 0,
     self.params['grid']['my0']), cmap='seismic', vmin=-divE_max, vmax=divE_max, zorder=0)
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="2%", pad=0.02)
    plt.colorbar(im, cax=cax, label=r'$\nabla\cdot \bf{E}$');
    ttl = ax.text(.1, 0.1, '', transform = ax.transAxes, va='center')
    diag = ax.text(2, 47, '', multialignment='left', va='top', bbox=dict(facecolor='white', edgecolor='gray', boxstyle='round,pad=0.2'))
    def init():
      ax.set_xlim(0, self.params['grid']['mx0']); ax.set_ylim(0, self.params['grid']['my0']);
      ax.set_xlabel(r'$x$'); ax.set_ylabel(r'$y$');
      ax.set_title(self.jobid);
      prtls = isolde.getParticles(self.path + '/output/prtl.tot.%05d' % 0)
      flds = isolde.getFields(self.path + '/output/flds.tot.%05d' % 0)
      sc1.set_offsets(np.array([prtls['1']['x'], prtls['1']['y']]).T)
      sc2.set_offsets(np.array([prtls['2']['x'], prtls['2']['y']]).T)
      im.set_data(flds['divE'][0])
      ttl.set_text('t=0')
      diag.set_text(template(prtls))
      return sc1, sc2, im, ttl, diag, lgd, grid1, grid2,
    def animate(i):
      prtls = isolde.getParticles(self.path + '/output/prtl.tot.%05d' % i)
      flds = isolde.getFields(self.path + '/output/flds.tot.%05d' % i)
      sc1.set_offsets(np.array([prtls['1']['x'], prtls['1']['y']]).T)
      sc2.set_offsets(np.array([prtls['2']['x'], prtls['2']['y']]).T)
      im.set_data(flds['divE'][0])
      ttl.set_text('t=' + str(i))
      diag.set_text(template(prtls))
      return sc1, sc2, im, ttl, diag, lgd, grid1, grid2,
    anim = FuncAnimation(fig, animate, init_func=init,
                         frames=40, interval=1000, blit=True, repeat=True)

# Here specify the test simulations and give additional specs of the environment
common_flags = ' -perseus -hdf5 -debug'
simulations = [
               TwoStream(common_flags,
                          params={
                            'node_configuration': {'sizex' : 28},
                            'time': {'last' : 50000},
                            'grid': {'mx0' : 5600, 'tileX' : 100},
                            'algorithm': {'c' : 0.35, 'nfilter': 8},
                            'output': {'enable' : 0, 'hst_enable' : 1, 'hst_interval' : 5},
                            'plasma': {'ppc0' : 64, 'sigma' : 10, 'c_omp' : 40},
                            'particles': {'nspec' : 2, 'maxptl1' : 1e5, 'm1' : 1, 'ch1' : -1, 'maxptl2' : 1e5, 'm2' : 1, 'ch2' : -1},
                            'problem': {'shift_gamma' : 2.5}}
                          ),
               PlasmaOsc(common_flags, nproc=4,
                          params={
                            'node_configuration': {'sizex' : 4},
                            'time': {'last' : 1000},
                            'grid': {'mx0' : 1120, 'tileX' : 10},
                            'algorithm': {'nfilter': 8},
                            'output': {'enable' : 1, 'prtl_enable' : 0, 'spec_enable' : 0, 'interval': 10, 'hst_enable' : 1, 'hst_interval' : 5},
                            'plasma': {'ppc0' : 500, 'sigma' : 1, 'c_omp' : 10},
                            'particles': {'nspec' : 2, 'maxptl1' : 1e8, 'm1' : 1, 'ch1' : -1, 'maxptl2' : 1e8, 'm2' : 1, 'ch2' : 1},
                            'problem': {'upstream_T' : 1e-5, 'amplitude' : 0.01, 'nwaves' : 1}}
                         ),
               Weibel(common_flags, nproc=28,
                          params={
                            'node_configuration': {'sizex' : 7, 'sizey' : 4},
                            'time': {'last' : 500},
                            'grid': {'mx0' : 1400, 'my0' : 1400, 'tileX' : 20, 'tileY' : 20},
                            'algorithm': {'nfilter': 8},
                            'output': {'enable' : 1, 'prtl_enable' : 0, 'spec_enable' : 0, 'interval': 10, 'hst_enable' : 1, 'hst_interval' : 10},
                            'plasma': {'ppc0' : 16, 'sigma' : 10, 'c_omp' : 10},
                            'particles': {'nspec' : 2, 'maxptl1' : 1e8, 'm1' : 1, 'ch1' : -1, 'maxptl2' : 1e8, 'm2' : 1, 'ch2' : 1},
                            'problem': {'backgr_T' : 1e-5, 'shift_beta' : 0.5}}
                         ),
              Merging(common_flags + ' -dwn', nproc=1,
                         params={
                           'node_configuration': {'sizex' : 1, 'sizey' : 1},
                           'time': {'last' : 500},
                           'grid': {'mx0' : 50, 'my0' : 50, 'tileX' : 10, 'tileY' : 10},
                           'algorithm': {'nfilter': 0},
                           'output': {'interval': 10, 'stride' : 1, 'smooth_window' : 0, 'write_nablas' : 1},
                           'plasma': {'ppc0' : 50, 'sigma' : 5, 'c_omp' : 50},
                           'particles': {'nspec' : 2, 'maxptl1' : 1e8, 'm1' : 1, 'ch1' : -1, 'dwn1' : 1, 'maxptl2' : 1e8, 'm2' : 1, 'ch2' : 1, 'dwn2' : 1},
                           'downsampling' : {'interval' : 1, 'start' : 1, 'max_weight' : 1e5, 'cartesian_bins' : 1, 'energy_min' : 0, 'energy_max' : 1e5, 'int_weights' : 0, 'dynamic_bins' : 1, 'mom_bins' : 1, 'mom_spread' : 1e5}
                           }
                        )
               ]

if (options.d):
  # diagnostic mode where you analize the test results
  import matplotlib.pyplot as plt
  import numpy as np
  import tristanVis.isolde as isolde
  import tristanVis.aux as aux
  aux.loadCustomStyles(style='fivethirtyeight', fs=10)
  fig = plt.figure(figsize=(12, 8))
  for ii, simulation in enumerate(simulations):
    if not (ii + 1 in tests):
      continue
    ax = plt.subplot(2, 2, ii + 1)
    simulation.path = testdir_full + '/%02d_' % (ii + 1) + simulation.jobid
    simulation.diag(ax, fig)
  plt.tight_layout()
  plt.show()
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
      simulation.path = testdir_full + '/%02d_' % (ii + 1) + simulation.jobid
      if os.path.exists(simulation.path):
        shutil.rmtree(simulation.path)
      os.makedirs(simulation.path)

      testlog.write(('TEST_#{}_'.format(ii+1) + simulation.jobid).ljust(50, '.') + '\n')

      # configure
      config_command = 'python configure.py '
      config_command += simulation.flags
      config_command += ' -{}d'.format(simulation.dimension)
      config_command += ' --{}='.format(simulation.userfile[:4]) + simulation.userfile
      os.system(config_command + suffix)

      # clean
      os.system('make clean' + suffix)
      # compile
      os.system('make all' + suffix)

      # check if compilation successfull
      simulation.exe = 'tristan-mp{}d'.format(simulation.dimension)
      simulation.exe_full = simulation.path + '/' + simulation.exe
      if os.path.isfile(codedir + '/exec/' + simulation.exe):
        testlog.write('compilation'.ljust(46, '.') + '[OK]\n')
        print ('Compilation of `{}` done.'.format(simulation.jobid))

        # move executable
        os.system('mv {} {}'.format(codedir + '/exec/' + simulation.exe, simulation.path) + suffix)

        # clean
        os.system('make clean' + suffix)

        # write input
        simulation.input = 'input.' + simulation.jobid
        simulation.input_full = simulation.path + '/' + simulation.input
        with open(simulation.input_full, 'w+') as inp:
          for block in simulation.params.keys():
            inp.write('\n<{}>\n\n'.format(block))
            for var in simulation.params[block].keys():
              inp.write('  {}  =  {}\n'.format(var, simulation.params[block][var]))
        testlog.write('input file'.ljust(46, '.') + '[OK]\n')

        # write submit
        simulation.submit = 'submit_' + simulation.jobid
        simulation.submit_full = simulation.path + '/' + simulation.submit
        with open(simulation.submit_full, 'w+') as sub:
          sub.write('#!/bin/bash\n')
          sub.write('#SBATCH -t {}\n'.format(simulation.walltime))
          sub.write('#SBATCH -n {}\n'.format(simulation.nproc))
          sub.write('#SBATCH -J {}\n'.format(simulation.jobid))
          sub.write('#SBATCH --output={}/tristan-v2.out\n'.format(simulation.path))
          sub.write('#SBATCH --error={}/tristan-v2.err\n\n'.format(simulation.path))

          sub.write('DIR={}\n'.format(simulation.path))
          sub.write('EXECUTABLE=$DIR/{}\n'.format(simulation.exe))
          sub.write('INPUT=$DIR/{}\n'.format(simulation.input))
          sub.write('OUTPUT_DIR=$DIR/output\n'.format(simulation.path))
          sub.write('SLICE_DIR=$DIR/slices\n'.format(simulation.path))
          sub.write('REPORT_FILE=$DIR/report\n'.format(simulation.path))
          sub.write('ERROR_FILE=$DIR/error\n\n'.format(simulation.path))

          for module in modules:
            sub.write('module load {}\n'.format(module))

          sub.write('\nmkdir $OUTPUT_DIR\n\n')

          sub.write('srun $EXECUTABLE -i $INPUT -o $OUTPUT_DIR -s $SLICE_DIR -r $RESTART_DIR -R $RESTART > $REPORT_FILE 2> $ERROR_FILE')
        testlog.write('submit file'.ljust(46, '.') + '[OK]\n')
      testlog.write('\n')

    testlog.write('\n\n')
    if not options.c:
      for ii, simulation in enumerate(simulations):
        if not (ii + 1 in tests):
          continue
        os.system('sbatch ' + simulation.submit_full)
        testlog.write(('`{}`'.format(simulation.jobid)).ljust(41, '.') + 'submitted\n')
