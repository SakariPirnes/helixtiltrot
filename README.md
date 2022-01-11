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

> **_NOTE:_** The HelixTiltRot package mainly uses NumPy arrays as inputs and outputs, see https://numpy.org/devdocs/user/quickstart.html and https://numpy.org/ for more information about NumPy.


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


Plotting how quantities have evolve in time in given coordinate file is a common starting point for an analysis. To quickly plot arrays obtained from example `helixtiltrot.local_tilt_angle` and `helixtiltrot.local_rotation_angle` the HelixTiltRot has a module `plot` containing three functions:\
  18. `helixtiltrot.plot.rotation` [\[source\]](https://github.com/SakariPirnes/helixtiltrot/blob/main/helixtiltrot/plot.py)\
  19. `helixtiltrot.plot.angle_map` [\[source\]](https://github.com/SakariPirnes/helixtiltrot/blob/main/helixtiltrot/plot.py)\
  20. `helixtiltrot.plot.angle_density` [\[source\]](https://github.com/SakariPirnes/helixtiltrot/blob/main/helixtiltrot/plot.py)


The documentation of these 20 function can be found from their docstrings. For example to see the documentation of `helixtiltrot.local_rotation_angle`:
```
>>> 
>>> import helixtiltrot
>>> print(helixtiltrot.local_rotation_angle.__doc__)


    Compute the local rotation angles from alpha-carbon coordinates, relative
    to chosen reference vector.

    Returns a new 1-D array of the local rotation angles in phase of residue
    `phase`, in radians, in range ]-pi,pi].
    

    Parameters
    ----------
    ca : array_like
        Array containing the Cartesian coordinates of alpha-carbon atoms.
        The shape of `ca` is (nf,nr,3), where nf is the number of frames,
        nr is the number of residues and the 3 values in the last axis
        are the Cartesian coordinates. If `ca` is not an
        array of `numpy.float64`s, a conversion is attempted.
    ref_vec : Array_like
        1-D array containing the Cartesian coordinates of the reference
        vector. Hence, the shape of `ref_vec` is (3,). If `ref_vec` is not
        an array of `numpy.float64`s, a conversion is attempted.
    mask : array_like or None, optional
        numpy mask with shape (nf,nr). Mask determines which `ca` atoms
        will be used. Default is None, which is equivalent
        to numpy.full((nf,nr), True), i.e. all `ca` atoms will be used.

    Returns
    -------
    local_rot : ndarray of `numpy.float64`s
        Returns a new 1-D array of rotation angle in phase of residue `phase`,
        in radians, in range ]-pi,pi]. The shape of `local_rot` is (nf,nr).

    
>>>
```
### Examples
To get an idea how the 20 functions provided in HelixTiltRot can be used, two example analysis done with HelixTiltRot can be found from chapter 4 of [the BSc thesis](https://github.com/SakariPirnes/helixtiltrot/blob/main/documentation-BSc_pirnes.pdf). The python scripts to produce the figures in the two examples can be found from folders [Example1](https://github.com/SakariPirnes/helixtiltrot/tree/main/Examples/Example1) and [Example2](https://github.com/SakariPirnes/helixtiltrot/tree/main/Examples/Example1).
