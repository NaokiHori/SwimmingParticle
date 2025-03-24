# Swimming particle

[![License](https://img.shields.io/github/license/NaokiHori/SwimmingParticle)](https://opensource.org/license/MIT)
[![Last Commit](https://img.shields.io/github/last-commit/NaokiHori/SwimmingParticle/main)](https://github.com/NaokiHori/SwimmingParticle/commits/main)

![cover image](https://github.com/NaokiHori/SwimmingParticle/blob/main/cover.gif)

## Overview

A finite-difference-based numerical simulator for diffusophoretic, self-propelling particles.

## Quick Start

Get the source:

```console
git clone https://github.com/NaokiHori/SwimmingParticle
cd SwimmingParticle
```

Build and run:

```console
make output
make all
./a.out
```

Check the result:

```console
python visualize.py
```

## Parameters

There are four parameters which can be changed:

- Péclet number
- Radius of the outer circle
- Number of radial and azimuthal grid points 

They are all defined in `src/main.c`.
**Note that the number of azimuthal grids are to be a power of 2.**

Re-compile the source (`make all`) when modified.

## Documentation

Governing equation and its numerical treatment are briefly documented [here](https://naokihori.github.io/SwimmingParticle/index.html).

## Reference

- Williamson, *J. Comput. Phys.* (**35**), 1980
- Hu et al., *Phys. Rev. Lett.* (**123**), 2019
- Wood and Porter, *J. Appl. Eng. Math.* (**6**), 2019

