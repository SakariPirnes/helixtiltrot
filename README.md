# HelixSide

A python implementation of the method constructed in BSc thesis:

[Analytical method for calculating
geometrical metrics for protein structures](https://github.com/SakariPirnes/helixside/blob/main/Pirnes_Sakari_BSc_thesis.pdf).

> **_NOTE:_** In the BSc thesis the quantities (local) side angle and side vector are called (local) rotation angle and rotation vector, respectively. Also the package is called with its old name, HelixTiltRot. The naming convention was updated on 10th August 2025.

## Installation

Create python3 environment:

    python3 -m venv env
    
Activate environment:

    source env/bin/activate

Install dependencies using requirements.txt:

    pip install -r requirements.txt

Install HelixSide:

    pip install -e .
    
    
## Documentation

> **_NOTE:_** The HelixSide package mainly uses NumPy arrays as inputs and outputs, see https://numpy.org/devdocs/user/quickstart.html for quickstart and https://numpy.org/ for more information about NumPy.


The following functions are availible in HelixSide:
1. `helixside.local_axes` [\[source\]](https://github.com/SakariPirnes/helixside/blob/main/helixside/core.py) - equation [(3.15)](https://github.com/SakariPirnes/helixside/blob/main/Pirnes_Sakari_BSc_thesis.pdf)
2. `helixside.axis` [\[source\]](https://github.com/SakariPirnes/helixside/blob/main/helixside/core.py) - equation [(3.2)](https://github.com/SakariPirnes/helixside/blob/main/Pirnes_Sakari_BSc_thesis.pdf)
3. `helixside.tilt_angle` [\[source\]](https://github.com/SakariPirnes/helixside/blob/main/helixside/core.py) - equation [(3.1)](https://github.com/SakariPirnes/helixside/blob/main/Pirnes_Sakari_BSc_thesis.pdf)
4. `helixside.kink_angle` [\[source\]](https://github.com/SakariPirnes/helixside/blob/main/helixside/core.py) - equation [(3.3)](https://github.com/SakariPirnes/helixside/blob/main/Pirnes_Sakari_BSc_thesis.pdf)
5. `helixside.tilt_vectors` [\[source\]](https://github.com/SakariPirnes/helixside/blob/main/helixside/core.py) - equation [(3.4)](https://github.com/SakariPirnes/helixside/blob/main/Pirnes_Sakari_BSc_thesis.pdf)
6. `helixside.local_tilt_angle` [\[source\]](https://github.com/SakariPirnes/helixside/blob/main/helixside/core.py) - equation [(3.5)](https://github.com/SakariPirnes/helixside/blob/main/Pirnes_Sakari_BSc_thesis.pdf)
7. `helixside.center_points` [\[source\]](https://github.com/SakariPirnes/helixside/blob/main/helixside/core.py) - subsubsection [Center of the local segment](https://github.com/SakariPirnes/helixside/blob/main/Pirnes_Sakari_BSc_thesis.pdf)
8. `helixside.side_vectors` [\[source\]](https://github.com/SakariPirnes/helixside/blob/main/helixside/core.py) - equation [(3.18)](https://github.com/SakariPirnes/helixside/blob/main/Pirnes_Sakari_BSc_thesis.pdf)
9. `helixside.local_side_angle` [\[source\]](https://github.com/SakariPirnes/helixside/blob/main/helixside/core.py) - equation [(3.29)](https://github.com/SakariPirnes/helixside/blob/main/Pirnes_Sakari_BSc_thesis.pdf)
10. `helixside.single_phase` [\[source\]](https://github.com/SakariPirnes/helixside/blob/main/helixside/core.py) - equation [(3.19)](https://github.com/SakariPirnes/helixside/blob/main/Pirnes_Sakari_BSc_thesis.pdf)
11. `helixside.side_angle` [\[source\]](https://github.com/SakariPirnes/helixside/blob/main/helixside/core.py) - equation [(3.20)](https://github.com/SakariPirnes/helixside/blob/main/Pirnes_Sakari_BSc_thesis.pdf)
12. `helixside.angle_diff` [\[source\]](https://github.com/SakariPirnes/helixside/blob/main/helixside/core.py) - equation [(3.21)](https://github.com/SakariPirnes/helixside/blob/main/Pirnes_Sakari_BSc_thesis.pdf)
13. `helixside.load_ca_dssp` [\[source\]](https://github.com/SakariPirnes/helixside/blob/main/helixside/core.py)
14. `helixside.sse_mask` [\[source\]](https://github.com/SakariPirnes/helixside/blob/main/helixside/core.py)
15. `helixside.circular_mean` [\[source\]](https://github.com/SakariPirnes/helixside/blob/main/helixside/core.py)
16. `helixside.circular_var` [\[source\]](https://github.com/SakariPirnes/helixside/blob/main/helixside/core.py)
17. `helixside.circular_std` [\[source\]](https://github.com/SakariPirnes/helixside/blob/main/helixside/core.py)

Detailed definitions of the functions 1-12 are in [the BSc thesis](https://github.com/SakariPirnes/helixside/blob/main/Pirnes_Sakari_BSc_thesis.pdf).

The function 13 `helixside.load_ca_dssp` uses [`mdtraj.load`](https://mdtraj.org/1.9.4/api/generated/mdtraj.load.html?highlight=load#mdtraj.load), [`mdtraj.Trajectory`](https://mdtraj.org/1.9.4/api/generated/mdtraj.Trajectory.html?highlight=trajectory#mdtraj.Trajectory) and [`mdtraj.compute_dssp`](https://mdtraj.org/1.9.4/api/generated/mdtraj.compute_dssp.html?highlight=dssp#mdtraj.compute_dssp) to obtain alpha-carbon coordinates and DSSP assigment codes from given coordinate files. See https://mdtraj.org/1.9.4/index.html for more information about the MDTraj python library and for supported file formats, including pdb, xtc, trr, dcd, binpos, netcdf, mdcrd, prmtop.

The side angles are 2pi-periodic, therefore the [circular mean](https://en.wikipedia.org/wiki/Circular_mean), [circular variance](https://en.wikipedia.org/wiki/Directional_statistics#Dispersion) and [circular standard deviation](https://en.wikipedia.org/wiki/Directional_statistics#Dispersion) are implemented in the functions 15-17, respectively.


Plotting how quantities evolve in time in given coordinate file is a common starting point for an analysis. To quickly plot arrays obtained from example `helixside.local_tilt_angle` and `helixside.local_side_angle` the HelixSide has a module `plot` containing three functions:\
  18. `helixside.plot.polar` [\[source\]](https://github.com/SakariPirnes/helixside/blob/main/helixside/plot.py)\
  19. `helixside.plot.angle_map` [\[source\]](https://github.com/SakariPirnes/helixside/blob/main/helixside/plot.py)\
  20. `helixside.plot.angle_density` [\[source\]](https://github.com/SakariPirnes/helixside/blob/main/helixside/plot.py)


The documentation of these 20 functions can be found from their docstrings. For example to see the documentation of `helixside.local_side_angle`:
```
>>> 
>>> import helixside
>>> print(helixside.local_side_angle.__doc__)


    Compute the local side angles from alpha-carbon coordinates, relative
    to chosen reference vector.

    Returns a new 1-D array of the local side angles in phase of residue
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
    local_side : ndarray of `numpy.float64`s
        Returns a new 1-D array of side angle in phase of residue `phase`,
        in radians, in range ]-pi,pi]. The shape of `local_side` is (nf,nr).

    
>>>
```
### Examples
To get an idea how the 20 functions provided in HelixSide can be used, two example analysis done with HelixSide can be found from chapter 4 of [the BSc thesis](https://github.com/SakariPirnes/helixside/blob/main/Pirnes_Sakari_BSc_thesis.pdf). The python scripts to produce the figures of these two examples can be found from folders [Example1](https://github.com/SakariPirnes/helixside/tree/main/Examples/Example1) and [Example2](https://github.com/SakariPirnes/helixside/tree/main/Examples/Example1).
