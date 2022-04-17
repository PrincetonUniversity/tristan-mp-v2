# Tristan-MP v2.3

For detailed tutorials and descriptions please visit our [wiki](https://ntoles.github.io/tristan-wiki/).

### Developed by
[Hayk Hakobyan](http://github.com/haykh), [Anatoly Spitkovsky](https://github.com/ntoles)

### Key contributors
[Jens Mahlmann](https://github.com/jmahlmann), [Sasha Philippov](https://github.com/Sashaph), [Daniel Groselj](https://github.com/dgroselj), [Sasha Chernoglazov](https://github.com/SChernoglazov), [Fabio Bacchini](https://github.com/fabsilfab)

By using this `Tristan-MP v2` you agree to the terms of the License Agreement (see `LICENSE` file in the repository). Please refer to the `CITATION` file for recommendations on how to properly cite this code if you use them in your research.

> Note, that certain modules (synchrotron/inverse Compton radiation, Breit-Wheeler pair production, Compton scattering, pair annihilation and guiding center particle pusher) are hidden in public release as they are still untested and/or under construction.

### Releases

* `v2.3` __Apr 2022__
  * Adaptive load balancing: see [wiki](https://ntoles.github.io/tristan-wiki/tristanv2-loadbal.html#adaptive-load-balancing)
  * Hybrid gyrokinetic (GCA) pusher see [wiki](https://ntoles.github.io/tristan-wiki/tristanv2-algorithms.html#guiding-center-approximation)
  * Non-blocking MPI communications support (enables precise reproducibility)