#! /usr/bin/env python
# -----------------------------------------------------------------------------------------
import argparse
import glob
import re

# Set template and output filenames
makefile_input = 'Makefile.in'
makefile_output = 'Makefile'

# Step 1. Prepare parser, add each of the arguments
parser = argparse.ArgumentParser()

user_directory = 'user/'
user_choices = glob.glob(user_directory + '*.F90')
user_choices = [choice[len(user_directory):-4] for choice in user_choices]

unit_directory = 'unit/'
unit_choices = glob.glob(unit_directory + '*.F90')
unit_choices = [choice[len(unit_directory):-4] for choice in unit_choices]

clusters = ['perseus', 'frontera', 'stellar']

# system
parser.add_argument('--cluster',
                    default=None,
                    choices=clusters,
                    help='choose cluster-specific configurations.')

parser.add_argument('-intel',
                    action='store_true',
                    default=True,
                    help='enable intel compiler')
parser.add_argument('-hdf5',
                    action='store_true',
                    default=False,
                    help='enable HDF5 & use h5pfc compiler')
parser.add_argument('-serial',
                    action='store_true',
                    default=False,
                    help='enable serial output')
parser.add_argument('-mpinonblock',
                    action='store_true',
                    default=False,
                    help='enable non-blocking mpi calls (faster simulation but non-reproducible)')

vec_group = parser.add_mutually_exclusive_group(required=False)
vec_group.add_argument('-avx2',
                       action='store_true',
                       default=False,
                       help='enable avx2 vectorization')
vec_group.add_argument('-avx512',
                       action='store_true',
                       default=False,
                       help='enable avx512 vectorization')

parser.add_argument('-lowmem',
                    action='store_true',
                    default=False,
                    help='enable low memory regime')
parser.add_argument('-ifport',
                    action='store_true',
                    default=False,
                    help='enable IFPORT library (`mkdir` etc)')

mpi_group = parser.add_mutually_exclusive_group()
mpi_group.add_argument('-mpi',
                       action='store_true',
                       default=False,
                       help='enable mpi')
mpi_group.add_argument('-mpi08',
                       action='store_true',
                       default=False,
                       help='enable mpi_f08')

mpi_group.add_argument('-test',
                       action='store_true',
                       default=False,
                       help='enable test mode')

# user file
user_group = parser.add_mutually_exclusive_group(required=True)
user_group.add_argument('--user',
                        default=None,
                        choices=user_choices,
                        help='select user file')
user_group.add_argument('--unit',
                        default=None,
                        choices=unit_choices,
                        help='select unit file')

# algorithms
parser.add_argument('--nghosts',
                    action='store',
                    default=3,
                    help='specify the # of ghost cells')

parser.add_argument('-extfields',
                    action='store_true',
                    default=False,
                    help='apply external fields')

parser.add_argument('-absorb',
                    action='store_true',
                    default=False,
                    help='enable absorbing boundaries')

parser.add_argument('--debug',
                    action='store',
                    default='OFF',
                    help='enable the debug mode at specific level')

parser.add_argument('-safe',
                    action='store_true',
                    default=False,
                    help='enable dynamic memory allocations (safe regime)')

parser.add_argument('--gca',
                    action='store',
                    default='OFF',
                    help='enable GCA mover with specific # of iterations')

parser.add_argument('-vay',
                    action='store_true',
                    default=False,
                    help='enable Vay pusher')

parser.add_argument('-payload',
                    action='store_true',
                    default=False,
                    help='enable particle payloads')

parser.add_argument('-usroutput',
                    action='store_true',
                    default=False,
                    help='enable user-specified output routines and condition')

dim_group = parser.add_mutually_exclusive_group(required=True)
dim_group.add_argument('-1d',
                       action='store_true',
                       default=False,
                       help='enable 1d')
dim_group.add_argument('-2d',
                       action='store_true',
                       default=False,
                       help='enable 2d')
dim_group.add_argument('-3d',
                       action='store_true',
                       default=False,
                       help='enable 3d')

parser.add_argument('-dwn',
                    action='store_true',
                    default=False,
                    help='enable particle downsampling')

parser.add_argument('-alb',
                    action='store_true',
                    default=False,
                    help='enable adaptive load balancing')

parser.add_argument('-slb',
                    action='store_true',
                    default=False,
                    help='enable static load balancing')

args = vars(parser.parse_args())

# Step 2. Set definitions and Makefile options based on above arguments

makefile_options = {}

if (args['user']):
    makefile_options['USER_FILE'] = args['user']
    makefile_options['USER_DIR'] = user_directory
else:
    makefile_options['USER_FILE'] = args['unit']
    makefile_options['USER_DIR'] = unit_directory

makefile_options['COMPILER_COMMAND'] = ''
makefile_options['COMPILER_FLAGS'] = ''
makefile_options['PREPROCESSOR_FLAGS'] = ''

# specific cluster:
specific_cluster = False
if (args['cluster'] is not None):
    specific_cluster = True
    clustername = args['cluster'].capitalize()
    args['intel'] = True
    args['ifport'] = True
    if args['cluster'] == 'perseus':
        args['mpi'] = True
        args['avx2'] = True
    elif args['cluster'] == 'frontera':
        args['mpi08'] = True
        args['avx512'] = True
        args['lowmem'] = True
    elif args['cluster'] == 'stellar':
        args['mpi08'] = True
        args['avx512'] = True

# compilation command
if args['hdf5']:
    makefile_options['COMPILER_COMMAND'] += 'h5pfc '
    makefile_options['PREPROCESSOR_FLAGS'] += '-DHDF5 '
else:
    if ((not args['mpi']) and (not args['mpi08'])):
        makefile_options['COMPILER_COMMAND'] += 'gfortran '
    else:
        makefile_options['COMPILER_COMMAND'] += 'mpif90 ' if args['intel'] else 'mpiifort '
if args['ifport']:
    makefile_options['PREPROCESSOR_FLAGS'] += '-DIFPORT '
if args['lowmem']:
    makefile_options['PREPROCESSOR_FLAGS'] += '-DLOWMEM '
if args['serial']:
    makefile_options['PREPROCESSOR_FLAGS'] += '-DSERIALOUTPUT '
if args['mpinonblock']:
    makefile_options['PREPROCESSOR_FLAGS'] += '-DMPINONBLOCK '

# mpi version
if args['mpi']:
    makefile_options['PREPROCESSOR_FLAGS'] += '-DMPI '
elif args['mpi08']:
    makefile_options['PREPROCESSOR_FLAGS'] += '-DMPI08 '

# debug
if args['debug'] != 'OFF':
    # non-intel compilers are not supported
    if int(args['debug']) >= 0:
        makefile_options['PREPROCESSOR_FLAGS'] += '-DDEBUG '
    if int(args['debug']) >= 1:
        makefile_options['COMPILER_FLAGS'] += '-traceback -fpe0 '
    if int(args['debug']) >= 2:
        makefile_options['COMPILER_FLAGS'] += '-check all -check noarg_temp_created '
else:
    makefile_options['COMPILER_FLAGS'] += '-Ofast '

if args['test']:
    makefile_options['PREPROCESSOR_FLAGS'] += '-DTESTMODE '

# compiler (+ vectorization etc)
if args['intel']:
    makefile_options['MODULE'] = '-module '
    makefile_options['COMPILER_FLAGS'] += '-O3 -DSoA -ipo -qopenmp-simd -qopt-report=5 -qopt-streaming-stores auto '
else:
    print ('Non-intel compilers may not work properly.')
    makefile_options['MODULE'] = '-J '
    makefile_options['COMPILER_FLAGS'] += '-O3 -DSoA -fwhole-program -mavx2 -fopt-info-vec -fopt-info-vec-missed -ftree-vectorizer-verbose=5 '

if args['avx2']:
    makefile_options['COMPILER_FLAGS'] += '-xCORE-AVX2 '
elif args['avx512']:
    makefile_options['COMPILER_FLAGS'] += '-xCORE-AVX512 -qopt-zmm-usage:high '

if args['1d']:
    makefile_options['EXE_NAME'] = 'tristan-mp1d'
    makefile_options['PREPROCESSOR_FLAGS'] += '-DoneD '
elif args['2d']:
    makefile_options['EXE_NAME'] = 'tristan-mp2d'
    makefile_options['PREPROCESSOR_FLAGS'] += '-DtwoD '
elif args['3d']:
    makefile_options['EXE_NAME'] = 'tristan-mp3d'
    makefile_options['PREPROCESSOR_FLAGS'] += '-DthreeD '

if args['absorb']:
    makefile_options['PREPROCESSOR_FLAGS'] += '-DABSORB '

# extra algorithms
if args['dwn']:
    makefile_options['PREPROCESSOR_FLAGS'] += '-DDOWNSAMPLING '
if args['alb']:
    makefile_options['PREPROCESSOR_FLAGS'] += '-DALB '
if args['slb']:
    makefile_options['PREPROCESSOR_FLAGS'] += '-DSLB '

if args['gca'] != 'OFF':
    makefile_options['PREPROCESSOR_FLAGS'] += '-DGCA -DGCAITER=' + \
        str(args['gca']) + ' '
if args['vay']:
    makefile_options['PREPROCESSOR_FLAGS'] += '-DVAY '
if args['payload']:
    makefile_options['PREPROCESSOR_FLAGS'] += '-DPRTLPAYLOADS '
if args['usroutput']:
    makefile_options['PREPROCESSOR_FLAGS'] += '-DUSROUTPUT '

# extra physics
if args['extfields']:
    makefile_options['PREPROCESSOR_FLAGS'] += '-DEXTERNALFIELDS '

makefile_options['PREPROCESSOR_FLAGS'] += '-DNGHOST=' + \
    str(args['nghosts']) + ' '

# Step 3. Create new files, finish up
with open(makefile_input, 'r') as current_file:
    makefile_template = current_file.read()
for key, val in makefile_options.items():
    makefile_template = re.sub(r'@{0}@'.format(key), val, makefile_template)
makefile_template = re.sub('# Template for ', '# ', makefile_template)
with open(makefile_output, 'w') as current_file:
    current_file.write(makefile_template)

# Finish with diagnostic output
print('==============================================================================')
print('Your TRISTAN distribution has now been configured with the following options:')
if (args['cluster'] is not None):
    print('  Cluster configurations:  `{}`'.format(clustername))

print('SETUP ........................................................................')
print('  Userfile:                ' + makefile_options['USER_FILE'])
print('  Dim:                     ' +
      ('1D' if args['1d'] else ('2D' if args['2d'] else ('3D' if args['3d'] else 'None'))))
print('  # of ghost zones:        ' + str(args['nghosts']))
print('  Load balancing:          ' + ('static/adaptive' if (args['alb'] and args['slb']) else (
    'static' if args['slb'] else ('adaptive' if args['alb'] else 'OFF'))))
print('  Particle downsampling:   ' + ('ON' if args['dwn'] else 'OFF'))
print('  Particle pusher:         ' + ('Vay' if args['vay'] else 'Boris') + (
    '/GCA ({} iterations)'.format(args['gca']) if args['gca'] != 'OFF' else ''))
print('  Particle payloads:       ' + ('ON' if args['payload'] else 'OFF'))

print('EXTRA .........................................................................')
print('  External fields:         ' + ('ON' if args['extfields'] else 'OFF'))
print('  Absorbing boundaries:    ' + ('ON' if args['absorb'] else 'OFF'))

print('TECHNICAL .....................................................................')

print('  Compiler:                ' + ('intel' if args['intel'] else 'gcc') +
                                      (' [avx2]' if args['avx2'] else
                                       (' [avx512]' if args['avx512'] else '')
                                       ))
print('  Debug mode:              ' +
      ('level ' if args['debug'] != 'OFF' else '') + args['debug'])
print('  Low memory mode:         ' + ('ON' if args['lowmem'] else 'OFF'))
print('  Output:                  ' + (('HDF5' +
      (' (serial)' if args['serial'] else ' (parallel)')) if args['hdf5'] else 'N/A'))
print('  User output:             ' + ('ON' if args['usroutput'] else 'OFF'))
print('  MPI version:             ' +
      ('old' if not args['mpi08'] else 'MPI_08'))
print('  `IFPORT` mkdir:          ' + ('ON' if args['ifport'] else 'OFF'))

print('==============================================================================')

print('  Compilation command:     ' + makefile_options['COMPILER_COMMAND']
      + makefile_options['PREPROCESSOR_FLAGS'] + makefile_options['COMPILER_FLAGS'])
