# Tristan v2.6

For detailed tutorials and code description please visit our [wiki](https://princetonuniversity.github.io/tristan-v2/). If you are a new user testing the code on a new computer cluster, please consider contributing to [this chapter](https://princetonuniversity.github.io/tristan-v2/tristanv2-configure.html#cluster-specific-customization) about cluster-specific configuration to make the life easier for future generations.

## Getting Started

### Prerequisites

* MPI (either OpenMPI or Intel-MPI)
* GCC or Intel Fortran compiler
* (optional) Parallel HDF5 (compiled with either OpenMPI or Intel-MPI)

On clusters typically all you need to do is to load the specific modules see [here](https://princetonuniversity.github.io/tristan-v2/tristanv2-configure.html#cluster-specific-customization). 

If you are, however, running on a local machine make sure to install the following prerequisites (assuming non-Intel compiler and `apt` package manager):

```shell
# gcc + openmpi
sudo apt install build-essential libopenmpi-dev
# hdf5
sudo apt libhdf5-openmpi-dev hdf5-tools
```

> Also make sure you have a working `python3` installation to be able to configure the code.

### Usage

```shell
# to view all configuration options
python3 configure.py --help
# to configure the code (example)
python3 configure.py -mpi08 -hdf5 --user=user_2d_rec -2d
# compile and link (-j compiles in parallel which is much faster)
make all -j
# run the code (on clusters need to do `srun`)
mpirun -np <NCORES> ./bin/tristan-mp2d -input ../inputs/input.2d_rec
```

### Docker

Another way to avoid the tedium of installing libraries (especially for local development) is to use [Docker containers](https://docs.docker.com/get-started/). This approach allows to quickly create an isolated linux environment with all the necessary packages already preinstalled (similar to a VM, but much lighter). The best thing is that VSCode can natively attach to a running container, and all the development can be done there (make sure to install the [appropriate extension](https://code.visualstudio.com/docs/remote/containers)).

To get started with this approach, make sure to install the Docker (as well as the Docker-compose) then simply follow these steps.

```shell
# from tristan root diretory
cd docker
# launch the container in the background (first-time run might take a few mins)
docker-compose up -d
# if you are attempting to also force rebuild the container, use the `--build` flag
docker-compose up -d --build
# ensure the container is running
docker ps
```

Then you can attach to the container via VSCode, or if you prefer the terminal, simply attach to the running container by doing:
```shell
docker exec -it trv2 zsh
# then the code will be in the `/home/$USER/tristan-v2` directory
cd /home/$USER/tristan-v2
```
To stop the container run `docker-compose stop`. To stop and delete the container simply run `docker-compose down` from the same `docker/` directory.

> Each container has in principle its own isolated filesystem, aside from the shared `/root/tristan-v2` directory. Any changes to the rest of the container's filesystem are discarded when the container is deleted (either via `docker-compose down` or directly `docker rm <CONTAINER>`).

## For developers/users

### Branching Policy

To prevent `v2` from growing to become the Lovecraftian monster it once was we highly encourage both users and developers to follow the guidelines on branching policies.

* For development:
    * all the new features shall be added to `dev/<feature>` branch;
    * as soon as the new features are tested they can be pushed to the main development branch: `dev/main`.
* For users:
    * user-specific branches are allowed (e.g. to test userfiles), but as soon as it works we highly encourage people to merge their new userfiles to `dev`;
    * the naming for user-specific branches shall be the following: `user/<problem>` or `user/<username>`;
    * if any of the core routines is modified, we highly encourage to use the `dev/...` branching instead of the `user/...`.

The `dev/main` branch is the most up-to-date **tested** version of the code, and it will be merged to `master` as soon as all the new features are documented.

### Development

Since `v2.4` code formatting policy is employed. To follow the proper formatting automatically we use the `fprettify` tool. To install `fprettify` in the local directory (via `pip`) one can use the `dev-requirements.txt` file, by running the following: 
```shell
# create a local pip environment
python3 -m virtualenv venv
# activate it
source venv/bin/activate
# install required modules
pip install -r dev-requirements.txt
```

After that one can either use the tool in a stand-alone manner:
```shell
./venv/bin/fprettify -i 2 -w 4 --whitespace-assignment true --enable-decl --whitespace-decl true --whitespace-relational true --whitespace-logical true --whitespace-plusminus true --whitespace-multdiv true --whitespace-print true --whitespace-type true --whitespace-intrinsics true --enable-replacements -l 1000 [FILENAME.F90]
```
or in the VSCode environment (see the extension list in the `.vscode/settings.json` of the current repo).

---

## Contributors (alphabetical order)

* Fabio Bacchini (University of Colorado Boulder)
* Alexander Chernoglazov (University of Maryland)
* Daniel Groselj (Columbia)
* Hayk Hakobyan (PPPL/Columbia)
* Jens Mahlmann (Princeton)
* Arno Vanthieghem (Princeton)

## Board of trustees

* Prof. Anatoly Spitkovsky (Princeton)
* Prof. Sasha Philippov (University of Maryland)

## Publications

__@TODO__

## Latest Releases
* `v2.6` __Jan 2023__
  * Minor cleanup + bugfixes
* `v2.6r1` __Nov 2022__
  * New absorption treatment with target B-field specified in the userfile
  * Fieldsolvers separated into different files
* `v2.5.2` __Nov 2022__
  * Compton module now has `nph_over_ne` to mimic high photon-to-lepton ratio
  * `absorb_x` now takes care of absorbing boundary conditions, while `boundary_x` is for MPI
  * New userfiles for compton-mediated turbulence and reconnection
  * Minor non-critical bugs
* `v2.5.1` __Aug 2022__
  * Minor configure bug when working with non-intel compilers
  * Minor bug with `usroutput` flag
* `v2.5` __Jul 2022__
  * Proper makefile/compilation command
  * Double precision option
  * Auto converter for public version releases
* `v2.4` __Jun 2022__
  * Stress-energy tensor output
  * `mpi` particle alignment issue
  * Formatting (see the "for developers" section)
  * Compilation linking with non-intel compilers
  * minor warnings on gcc + intel fixed
* `v2.3` __Feb 2022__
  * Reproducibility (blocking MPI comms)
  * New boundary/injection conditions in the reconnection userfile
  * Minor bugfixes
* `v2.2` __May 2021__
  * Adaptive load balancing
  * Dynamic reallocation of tiles (more stable on low memory machines, slightly slower)
  * Explicit low-memory mode for clusters such as frontera (-lowmem flag)
  * User-specific output added at runtime
  * Updated the fulltest.py facility (see wiki)
  * Debug levels (0, 1, 2)
  * A facility for warnings and diagnostic output (see wiki)
  * Coupling of GCA with Vay pusher
  * Cooling limiter
  * Fluid velocity output
  * Major bugfix in the momentum output
  * Major restructuring in the initializer
  * Minor improvements, restructurings and bugfixes
* `v2.1.2` __Mar 2021__
  * Dynamic reallocation of tiles
  * Coupling of GCA with Vay pusher
  * Major restructuring in the initializer
  * A facility for warnings and diagnostic output
  * Cooling limiter
  * Fluid velocity output
  * Lots of subroutines to be used in the future for the adaptive load balancing
  * Minor improvements, restructurings and bugfixes
* `v2.1.1` __Jan 2021__
  * Patched the Compton cross section normalization
  * Minor bug fixes
* `v2.1` __Dec 2020__
  * Pair annihilation module (+ advanced test for all QED modules coupled)
  * Merging of charged particles
  * Particle payloads
  * Vay pusher
  * Individual particle current deposition (necessary for reflecting walls and downsampling of charged particles)
  * Spatial binning for spectra
  * Automated testing framework with fulltest.py
  * Major restructuring of output modules
  * Major restructuring of QED modules
  * Minor issues fixed for GCA pusher
* `v2.0.1` __Jul 2020__
  * Compton scattering and GCA pusher added
  * Slice outputs for 3d added
  * Particle momentum binning improved for downsampling
  * Minor bugs fixed

---

We employ [semantic versioning](https://semver.org/) for this code. Given a version number `v<MAJOR>.<MINOR>.<PATCH>r<CANDIDATE>`, increment the:
- `MAJOR` version when you make incompatible changes,
- `MINOR` version when you add functionality in a backwards compatible manner, and
- `PATCH` version when you make backwards compatible bug fixes.
- `CANDIDATE`: release candidate for "nightly" builds
