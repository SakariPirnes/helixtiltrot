# HelixTiltRot

A python implementation of the method constructed in BSc thesis:

[Analytical method for tilt and rotation of
α-helix and β-strand like secondary
structures](https://github.com/SakariPirnes/helixtiltrot/blob/main/documentation-BSc_pirnes.pdf).

## Installation

Create python3 environment:

    python3 -m venv env
    
Activate environment:

    source env/bin/activate

Install dependencies using requirements.txt:

    pip install -r requirements.txt

Install HelixTiltRot:

    pip install -e .
    
    
## Documentation

> **_NOTE:_** The HelixTiltRot package mainly uses NumPy arrays as inputs and outputs, see https://numpy.org/devdocs/user/quickstart.html and https://numpy.org/ for information about NumPy.


The following functions are availible in HelixTiltRot:
1. `helixtiltrot.local_axes` [\[source\]](https://github.com/SakariPirnes/helixtiltrot/blob/main/helixtiltrot/core.py) - equation [(3.15)](https://github.com/SakariPirnes/helixtiltrot/blob/main/documentation-BSc_pirnes.pdf)
2. `helixtiltrot.axis` [\[source\]](https://github.com/SakariPirnes/helixtiltrot/blob/main/helixtiltrot/core.py) - equation [(3.2)](https://github.com/SakariPirnes/helixtiltrot/blob/main/documentation-BSc_pirnes.pdf)
3. `helixtiltrot.tilt_angle` [\[source\]](https://github.com/SakariPirnes/helixtiltrot/blob/main/helixtiltrot/core.py) - equation [(3.1)](https://github.com/SakariPirnes/helixtiltrot/blob/main/documentation-BSc_pirnes.pdf)
4. `helixtiltrot.kink_angle` [\[source\]](https://github.com/SakariPirnes/helixtiltrot/blob/main/helixtiltrot/core.py) - equation [(3.3)](https://github.com/SakariPirnes/helixtiltrot/blob/main/documentation-BSc_pirnes.pdf)
5. `helixtiltrot.tilt_vectors` [\[source\]](https://github.com/SakariPirnes/helixtiltrot/blob/main/helixtiltrot/core.py) - equation [(3.4)](https://github.com/SakariPirnes/helixtiltrot/blob/main/documentation-BSc_pirnes.pdf)
6. `helixtiltrot.local_tilt_angle` [\[source\]](https://github.com/SakariPirnes/helixtiltrot/blob/main/helixtiltrot/core.py) - equation [(3.5)](https://github.com/SakariPirnes/helixtiltrot/blob/main/documentation-BSc_pirnes.pdf)
7. `helixtiltrot.center_points` [\[source\]](https://github.com/SakariPirnes/helixtiltrot/blob/main/helixtiltrot/core.py) - subsubsection [Center of the local segment](https://github.com/SakariPirnes/helixtiltrot/blob/main/documentation-BSc_pirnes.pdf)
8. `helixtiltrot.rotation_vectors` [\[source\]](https://github.com/SakariPirnes/helixtiltrot/blob/main/helixtiltrot/core.py) - equation [(3.18)](https://github.com/SakariPirnes/helixtiltrot/blob/main/documentation-BSc_pirnes.pdf)
9. `helixtiltrot.local_rotation_angle` [\[source\]](https://github.com/SakariPirnes/helixtiltrot/blob/main/helixtiltrot/core.py) - equation [(3.29)](https://github.com/SakariPirnes/helixtiltrot/blob/main/documentation-BSc_pirnes.pdf)
10. `helixtiltrot.single_phase` [\[source\]](https://github.com/SakariPirnes/helixtiltrot/blob/main/helixtiltrot/core.py) - equation [(3.19)](https://github.com/SakariPirnes/helixtiltrot/blob/main/documentation-BSc_pirnes.pdf)
11. `helixtiltrot.rotation_angle` [\[source\]](https://github.com/SakariPirnes/helixtiltrot/blob/main/helixtiltrot/core.py) - equation [(3.20)](https://github.com/SakariPirnes/helixtiltrot/blob/main/documentation-BSc_pirnes.pdf)
12. `helixtiltrot.angle_diff` [\[source\]](https://github.com/SakariPirnes/helixtiltrot/blob/main/helixtiltrot/core.py) - equation [(3.21)](https://github.com/SakariPirnes/helixtiltrot/blob/main/documentation-BSc_pirnes.pdf)
13. `helixtiltrot.load_ca_dssp` [\[source\]](https://github.com/SakariPirnes/helixtiltrot/blob/main/helixtiltrot/core.py)
14. `helixtiltrot.sse_mask` [\[source\]](https://github.com/SakariPirnes/helixtiltrot/blob/main/helixtiltrot/core.py)
15. `helixtiltrot.circular_mean` [\[source\]](https://github.com/SakariPirnes/helixtiltrot/blob/main/helixtiltrot/core.py)
16. `helixtiltrot.circular_var` [\[source\]](https://github.com/SakariPirnes/helixtiltrot/blob/main/helixtiltrot/core.py)
17. `helixtiltrot.circular_std` [\[source\]](https://github.com/SakariPirnes/helixtiltrot/blob/main/helixtiltrot/core.py)

Detailed definitions of the functions 1-12 are found in [the BSc thesis](https://github.com/SakariPirnes/helixtiltrot/blob/main/documentation-BSc_pirnes.pdf).

The function 13 `helixtiltrot.load_ca_dssp` uses [`mdtraj.load`](https://mdtraj.org/1.9.4/api/generated/mdtraj.load.html?highlight=load#mdtraj.load), [`mdtraj.Trajectory`](https://mdtraj.org/1.9.4/api/generated/mdtraj.Trajectory.html?highlight=trajectory#mdtraj.Trajectory) and [`mdtraj.compute_dssp`](https://mdtraj.org/1.9.4/api/generated/mdtraj.compute_dssp.html?highlight=dssp#mdtraj.compute_dssp) to obtain alpha-carbon coordinates and DSSP assigment codes from given coordinate files. See https://mdtraj.org/1.9.4/index.html for more information about the MDTraj python library and for supported file formats, including pdb, xtc, trr, dcd, binpos, netcdf, mdcrd, prmtop.

The rotation angles are 2pi-periodic, therefore the [circular mean](https://en.wikipedia.org/wiki/Circular_mean), [circular variance](https://en.wikipedia.org/wiki/Directional_statistics#Dispersion) and [circular standard deviation](https://en.wikipedia.org/wiki/Directional_statistics#Dispersion) are implemented in the functions 15-17, respectively.


Plotting how quantities have evolve in time in given coordinate file is a common starting point for an analysis. The HelixTiltRot contains a module `plot` containing three functions:\
  18. `helixtiltrot.plot.rotation` [\[source\]](https://github.com/SakariPirnes/helixtiltrot/blob/main/helixtiltrot/plot.py)\
  19. `helixtiltrot.plot.angle_map` [\[source\]](https://github.com/SakariPirnes/helixtiltrot/blob/main/helixtiltrot/plot.py)\
  20. `helixtiltrot.plot.angle_density` [\[source\]](https://github.com/SakariPirnes/helixtiltrot/blob/main/helixtiltrot/plot.py)


