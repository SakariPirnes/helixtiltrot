import numpy as np
from mdtraj import load, compute_dssp
from helixtiltrot.heart import *
from helixtiltrot import types
from warnings import catch_warnings, simplefilter


 

__all__ = ['load_ca_dssp','local_axes', 'axis', 'kink_angle','center_points','tilt_vectors','rotation_vectors','tilt_angle', 'local_tilt_angle','local_rotation_angle','rotation_angle','single_phase','sse_mask','angle_diff','circular_mean','circular_var','circular_std']







def load_ca_dssp(fname,top=None):
    

    """

    Read the alpha-carbon coordinates and compute the DSSP
    assigment codes from coordinate file.

    Returns tuple containing a ndarray of the alpha-carbon coordinates
    and a ndarray of the DSSP assigment codes.


    Parameters
    ----------
    fname : str
        Filename containing trajectory file.
    top : str or None, optional
        Most trajectory formats do not contain topology information.
        Pass in either the path to a RCSB PDB file, a trajectory,
        or a topology to supply this information. This option is not
        required for the .h5, .lh5, and .pdb formats, which already 
        contain topology information. Default is None.

    Returns
    -------
    ca,dssp : tuple
        A tuple containing a ndarray of the alpha-carbon coordinates
        and a ndarray of the DSSP assigment codes. The shape `ca.shape`
        of the ndarray containing the alpha-carbon coordinates is `(nf,nr,3)`,
        where `nf` and `nr` are respectively the number of the frames and 
        residues in `fname`. The shape `dssp.shape` of the ndarray containing
        the DSSP assigment codes is `(nf,nr)`
        
    """


    
    # Load mdtraj.Trajectory object
    if top == None:
        traj = load(fname)

    else:
        traj = load(fname,top=top)
    
    
    # take only protein
    prot_indexes = traj.topology.select('protein')
    traj.atom_slice(prot_indexes, inplace=True)
    
        
    # Alpha-carbon atom positions
    ca_indexes = traj.topology.select('name CA')
    ca = traj.xyz[:,ca_indexes,:].astype(np.float64)

    
    # DSSP
    dssp = compute_dssp(traj, simplified=False)
    
    return ca, dssp




def sse_mask(dssp, sse = 'H'):

    """

    Compute numpy mask of chosen secondary structures from the DSSP
    assigment codes.

    Returns a new numpy mask.


    Parameters
    ----------
    dssp : array_like
        Array of `character`s. If `dssp` is not an
        array of `str`s, a conversion is attempted.
    sse : str, optional
        String used to determine whether an element of the `mask`
        is `True` or `False`. Default is 'H', which stands for alpha-helix
        in the DSSP assigment codes.

    Returns
    -------
    mask : ndarray of `bool`s
        A new ndarray of `bool`s with the same shape as `dssp`. An element of
        the `mask` is `True` if the corresponding character of `dssp` is in
        `sse` string.

    """



    dssp = types.dssp(dssp, 'dssp')
    sse = types.sse(sse,'sse')

    


    masks = []
    for s in sse:
        masks.append(dssp == s)



    return np.any(masks, axis=0)




def local_axes(ca, mask=None):

    """

    Compute the local axes from the alpha-carbon coordinates of
    a helical structure.

    Returns a new 3-D array of the local axes for each frame.
    The shape of `h` is (nf,nr-3,3).


    Parameters
    ----------
    ca : array_like
        Array containing the Cartesian coordinates of alpha-carbon atoms.
        The shape of `ca` is (nf,nr,3), where nf is the number of frames,
        nr is the number of residues and the 3 values in the last axis
        are the Cartesian coordinates. If `ca` is not an
        array of `numpy.float64`s, a conversion is attempted.
    mask : array_like or None, optional
        numpy mask with shape (nf,nr). Mask determines which `ca` atoms
        will be used. Default is None, which is equivalent
        to numpy.full((nf,nr), True), i.e. all `ca` atoms will be used.

    Returns
    -------
    h : ndarray of `numpy.float64`s
        Returns a new 3-D array of the local axes for each frame.
        The shape of `h` is (nf,nr-3,3).

    """


    ca  = types.ca(ca,'ca').copy()


    if mask is not None:

        mask = types.array(mask.copy(),'mask', dtype = bool)
        types.equal_shape(ca.shape,'ca', mask.shape, 'mask', axes_to_equal=[0,1] )

        ca[~mask] = np.nan



    # four consecutive ca atoms
    r1 = ca[:,:-3]
    r2 = ca[:,1:-2]
    r3 = ca[:,2:-1]
    r4 = ca[:,3:]
    
    # vectors pointing to consecutive ca atoms
    r21 = normalize(r2-r1)
    r32 = normalize(r3-r2)
    r43 = normalize(r4-r3)
    

    # Vectors of plane which normal is the local axis
    s1 = normalize(r21 - r32)
    s2 = normalize(r32 - r43)
        
        
    # (unnormalized) APPROXb
    h = np.cross( s2, s1 )


    # mask for APPROXc's condition
    mc = dot(s1,s1) < -0.9


    # update APPROXb to APPROXc
    h[mc] = (r4-r2+r3-r1)[mc]


    # return APPROXc of local axes
    return normalize( h )




def axis(ca, mask=None, m=0, n=None):
    
    """

    Compute the axis from alpha-carbon coordinates of residues m,...,n-1.

    Returns a new 2-D array of the axis `H` for each frame.
    The shape of `H` is (nf,3).


    Parameters
    ----------
    ca : array_like
        Array containing the Cartesian coordinates of alpha-carbon atoms.
        The shape of `ca` is (nf,nr,3), where nf is the number of frames,
        nr is the number of residues and the 3 values in the last axis
        are the Cartesian coordinates. If `ca` is not an
        array of `numpy.float64`s, a conversion is attempted.
    mask : array_like or None, optional
        numpy mask with shape (nf,nr). Mask determines which `ca` atoms
        will be used. Default is None, which is equivalent
        to numpy.full((nf,nr), True), i.e. all `ca` atoms will be used.
    m : int or None, optional
        Only residues >=m will be used. Default is 0, the first residue.
    n : int or None, optional
        Only residues <n will be used. Default is None, in which case
        all residues >=m will be used.

    Returns
    -------
    H : ndarray of `numpy.float64`s
        Returns a new 2-D array of the axis `H` for each frame computed
        from residues m,...,n-1. The shape of `H` is (nf,3).

    """

    ca  = types.ca(ca,'ca').copy()

    if mask is not None:

        mask = types.array(mask.copy(),'mask', dtype = bool)
        types.equal_shape(ca.shape,'ca', mask.shape, 'mask', axes_to_equal=[0,1] )

        h = local_axes(ca[:,m:n], mask=mask[:,m:n])
    else:
        h = local_axes(ca[:,m:n])


    H_mn = np.nansum([h], axis=2)

    return normalize( H_mn )[0]




def kink_angle(ca, m1, n1, m2, n2, mask=None):


    """

    Compute the kink angle from alpha-carbon coordinates between residues
    m1,...,n1-1 and m2,...n2-1.

    Returns a new 1-D array of the local tilt angles in phase of residue
    `phase`, in radians, in range ]-pi,pi].


    Parameters
    ----------
    ca : array_like
        Array containing the Cartesian coordinates of alpha-carbon atoms.
        The shape of `ca` is (nf,nr,3), where nf is the number of frames,
        nr is the number of residues and the 3 values in the last axis
        are the Cartesian coordinates. If `ca` is not an
        array of `numpy.float64`s, a conversion is attempted.
    mask : array_like or None, optional
        numpy mask with shape (nf,nr). Mask determines which `ca` atoms
        will be used. Default is None, which is equivalent
        to numpy.full((nf,nr), True), i.e. all `ca` atoms will be used.
    m1 : int or None, optional
        Only residues >=m1 will be used. Default is 0, the first residue.
    n1 : int or None, optional
        Only residues <n1 will be used. Default is None, in which case
        all residues >=m1 will be used.
    m2 : int or None, optional
        Only residues >=m2 will be used. Default is 0, the first residue.
    n2 : int or None, optional
        Only residues <n2 will be used. Default is None, in which case
        all residues >=m2 will be used.

    Returns
    -------
    kink : ndarray of `numpy.float64`s
        Returns a new 1-D array of the kink angle between residues m1,...,n1-1
        and m2,...n2-1, in radians, in range ]-pi,pi]. The shape of `kink` is (nf,).

    """




    H1 = axis(ca, mask=mask, m=m1, n=n1)
    H2 = axis(ca, mask=mask, m=m2, n=n2)

    cosine_of_kink_angle = np.clip( (H1*H2).sum(axis=1), -1.0, 1.0)

    return np.arccos(cosine_of_kink_angle)





def tilt_angle(ca, ref_vec, mask=None, m=0, n=None):


    """

    Compute the tilt angle from alpha-carbon coordinates of residues m,...,n-1, relative
    to a chosen reference vector.

    Returns a new 1-D array of the local tilt angles in phase of residue
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
    m : int or None, optional
        Only residues >=m will be used. Default is 0, the first residue.
    n : int or None, optional
        Only residues <n will be used. Default is None, in which case
        all residues >=m will be used.

    Returns
    -------
    tilt : ndarray of `numpy.float64`s
        Returns a new 1-D array of tilt angle of residues m,...,n-1,
        in radians, in range ]-pi,pi]. The shape of `tilt` is (nf,).

    """

    H_mn = axis(ca, mask=mask, m=m, n=n)
    ref_vec = types.k(ref_vec, H_mn.shape)



    cosine_of_tilt_angle = np.clip( (H_mn*ref_vec).sum(axis=1), -1.0, 1.0)

    return np.arccos( cosine_of_tilt_angle )








def tilt_vectors(ca, mask=None):


    """

    Compute the tilt vectors from alpha-carbon coordinates.

    Returns a new 3-D array of (normalized) tilt vectors.


    Parameters
    ----------
    ca : array_like
        Array containing the Cartesian coordinates of alpha-carbon atoms.
        The shape of `ca` is (nf,nr,3), where nf is the number of frames,
        nr is the number of residues and the 3 values in the last axis
        are the Cartesian coordinates. If `ca` is not an
        array of `numpy.float64`s, a conversion is attempted.
    mask : array_like or None, optional
        numpy mask with shape (nf,nr). Mask determines which `ca` atoms
        will be used. Default is None, which is equivalent
        to numpy.full((nf,nr), True), i.e. all `ca` atoms will be used.

    Returns
    -------
    tilt_vec : ndarray of `numpy.float64`s
        Returns a new 3-D array of (normalized) rotation vectors. 
        The `tilt_vec` has the same shape as `ca`.

    """

    h = local_axes(ca, mask=mask)


    return h_to_tiltvecs(h)





def local_tilt_angle(ca, ref_vec, mask=None):

    """

    Compute the local tilt angles from alpha-carbon coordinates, relative
    to a chosen reference vector.

    Returns a new 2-D array of local tilt angles in phase of residue `phase`, 
    in radians, in range ]-pi,pi].
    

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
    local_tilt : ndarray of `numpy.float64`s
        Returns a new 2-D array of the local tilt angles in phase of residue
        `phase`, in radians, in range ]-pi,pi]. The shape of `local_tilt` 
        is (nf,nr).

    """





    t_vec = tilt_vectors(ca, mask=mask)

    ref_vec = types.k(ref_vec, ca.shape)


    cosine_of_local_tilt_angle = np.clip( dot(t_vec,ref_vec), -1.0, 1.0)

    
    return np.arccos( cosine_of_local_tilt_angle )




def center_points(ca, mask=None):

    """

    Compute the center points of local segments from
    alpha-carbon coordinates.

    Returns a new 3-D array of the center point.


    Parameters
    ----------
    ca : array_like
        Array containing the Cartesian coordinates of alpha-carbon atoms.
        The shape of `ca` is (nf,nr,3), where nf is the number of frames,
        nr is the number of residues and the 3 values in the last axis
        are the Cartesian coordinates. If `ca` is not an
        array of `numpy.float64`s, a conversion is attempted.
    mask : array_like or None, optional
        numpy mask with shape (nf,nr). Mask determines which `ca` atoms
        will be used. Default is None, which is equivalent
        to numpy.full((nf,nr), True), i.e. all `ca` atoms will be used.

    Returns
    -------
    center_points : ndarray of `numpy.float64`s
        Returns a new 3-D array of the center points. 

    """

    h = local_axes(ca, mask=mask)


    return h_to_centerpoints(h,ca)





def rotation_vectors(ca, mask=None):

    """

    Compute the rotation vectors from alpha-carbon coordinates.

    Returns a new 3-D array of (normalized) rotation vectors.


    Parameters
    ----------
    ca : array_like
        Array containing the Cartesian coordinates of alpha-carbon atoms.
        The shape of `ca` is (nf,nr,3), where nf is the number of frames,
        nr is the number of residues and the 3 values in the last axis
        are the Cartesian coordinates. If `ca` is not an
        array of `numpy.float64`s, a conversion is attempted.
    mask : array_like or None, optional
        numpy mask with shape (nf,nr). Mask determines which `ca` atoms
        will be used. Default is None, which is equivalent
        to numpy.full((nf,nr), True), i.e. all `ca` atoms will be used.

    Returns
    -------
    rot_vec : ndarray of `numpy.float64`s
        Returns a new 3-D array of (normalized) rotation vectors. 
        The `rot_vec` has the same shape as `ca`.

    """

    h = local_axes(ca, mask=mask)


    return h_to_rotvecs(h,ca)





def local_rotation_angle( ca, ref_vec, mask=None ):


    """

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

    """




    h = local_axes(ca, mask=mask)

    k = types.k(ref_vec, ca.shape)



    t_vec = h_to_tiltvecs(h)
    r_vec = h_to_rotvecs(h,ca)




    t_cross_k = np.cross(t_vec, k)


    r_dot_t_cross_k = r_vec * t_cross_k
    r_dot_t_cross_k = r_dot_t_cross_k[:,:,0] + r_dot_t_cross_k[:,:,1] + r_dot_t_cross_k[:,:,2]


    return np.arctan2( r_dot_t_cross_k, dot(r_vec,k) )



def single_phase(local_rot, turn_angle_deg, phase):
    

    """

    Shift local rotation angles in radians into local rotation angles in phase
    of residue `phase`, in radians .

    Returns a new array of local rotation angles in phase of residue `phase`, 
    in radians, in range ]-pi,pi].
    

    Parameters
    ----------
    local_rot : array_like
        Array containing local rotation angles. If `local_rot` is not an
        array of `numpy.float64`s, a conversion is attempted.
    turn_angle_deg : float or int
        The turn angle in degrees.
    phase : int
        Desired phase for `local_rot`

    Returns
    -------
    local_rot_phase : ndarray of `numpy.float64`s
        Returns a new array of local rotation angles in phase of residue `phase`, 
        in radians, in range ]-pi,pi].


    """

    local_rot = types.array(local_rot, 'local_rot')



    # Numper of residues
    n_res = local_rot.shape[-1]



    if isinstance( phase, int) == False:


        raise TypeError('\"phase\" argument should be int. You supplied {}'.format(type(phase)))

        
        
        
        

    if isinstance( turn_angle_deg, (int, float) ):


        # Change turn_angle from degrees to radians
        turn_angle  = turn_angle_deg*np.pi/180
    
        
        
        #,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
        # Angles to rotate residues to the chosen phase


        # Angles with 1 radians step between them.
        phase_angles_step1 = np.arange( -phase, (n_res - phase) )


        # Conver to the correct turn_angle
        phase_angles = phase_angles_step1 * turn_angle
        #''''''''''''''''''''''''''''''''''''''''''''''''''''''''

        
    else:

        raise TypeError('\"turn_angle_deg\" argument should be int or float. You supplied {}'.format(type(turn_angle_deg)))
    
    #''''''''''''''''''''''''''''''''''''''''''''''''''''''''

    
    local_rot_in_chosen_phase = local_rot + phase_angles

    not_nan = ~np.isnan(local_rot_in_chosen_phase)
    local_rot_in_chosen_phase[not_nan] =  np.pi - (( np.pi - local_rot_in_chosen_phase[not_nan] ) % (2 * np.pi ))




    return local_rot_in_chosen_phase




def circular_mean(data,  axis=None):

    """

    Compute the circular mean along the specified axis, ignoring NaNs.
    
    Returns the average of the array elements.  The average is taken over
    the flattened array by default, otherwise over the specified axis.


    Parameters
    ----------
    data : array_like
        Array containing numbers whose circular mean is desired. If `data` is not an
        array of `numpy.float64`s, a conversion is attempted.
    axis : {int, tuple of int, None}, optional
        Axis or axes along which the means are computed. The default is to compute
        the mean of the flattened array.

    Returns
    -------
    circ_mean : ndarray of `numpy.float64`s
        Returns a new array containing the circular mean values.
        Nan is returned for slices that contain only NaNs.

    """


    data = types.array(data,'data')


        
    # We asume period of the data is 2pi, i.e. array of angles

    data_in_complex_plane = np.exp(1j*data)

    with catch_warnings():
        simplefilter("ignore", category=RuntimeWarning)

        mean_of_data_in_complex_plane = np.nanmean(data_in_complex_plane, axis=axis)


    return np.angle(mean_of_data_in_complex_plane)




def rotation_angle(ca, ref_vec, turn_angle_deg, phase, mask=None, m=0, n=None):

    """

    Compute the rotation angle from alpha-carbon coordinates of residues
    m,...,n-1, relative to a chosen reference vector.

    Returns a new 1-D array of rotation angles in phase of residue `phase`, 
    in radians, in range ]-pi,pi].
    

    Parameters
    ----------
    ca : array_like
        Array containing the Cartesian coordinates of alpha-carbon atoms.
        The shape of `ca` is (nf,nr,3), where nf is the number of frames,
        nr is the number of residues and the 3 values in the last axis
        are the Cartesian coordinates. If `ca` is not an
        array of `numpy.float64`s, a conversion is attempted.
    ref_vec : Array_like
        1-D Array containing the Cartesian coordinates of the reference
        vector. Hence, the shape of `ref_vec` is (3,). If `ref_vec` is not
        an array of `numpy.float64`s, a conversion is attempted.
    turn_angle_deg : float or int
        The turn angle in degrees.
    phase : int
        Desired phase for `local_rot`.
    mask : array_like or None, optional
        numpy mask with shape (nf,nr). Mask determines which `ca` atoms
        will be used. Default is None which is equivalent
        to numpy.full((nf,nr), True), i.e. all `ca` atoms will be used.
    m : int or None, optional
        Only residues >=m will be used. Default is 0, the first residue.
    n : int or None, optional
        Only residues <n will be used. Default is None, in which case
        all residues >=m will be used.

    Returns
    -------
    rot : ndarray of `numpy.float64`s
        Returns a new 1-D array of rotation angle of residues
        m,...,n-1, relative to a chosen reference vector in phase of residue `phase`,
        in radians, in range ]-pi,pi]. The shape of `rot` is (nf,).

    """


    local_rot = local_rotation_angle( ca, ref_vec, mask=mask )

    rot = single_phase(local_rot, turn_angle_deg, phase)


    return circular_mean(rot[:,m:n],  axis=1)


 


def angle_diff(angles1, angles2):


    """


    Compute substraction between two broadcastable arrays containing angles
    in radians.

    Returns array of resulting difference angles in radians,
    in the range ]-pi,pi]. 


    Parameters
    ----------
    angles1 : array_like
        Array containing angles in radians, should be broadcastable with `angles2`. 
        If `angles1` is not an array of `numpy.float64`s,
        a conversion is attempted.
    angles2 : array_like
        Array containing angles in radians, should be broadcastable with `angles1`. 
        If `angles2` is not an array of `numpy.float64`s,
        a conversion is attempted.

    Returns
    -------
    diff : ndarray of `numpy.float64`s
        Returns a new array containing the difference angles.

    """

    
    angles1 = types.array(angles1, 'angles1')
    angles2 = types.array(angles2, 'angles2')

    #types.equal_shape(angles1.shape,'angles1', angles2.shape, 'angles2', axes_to_equal='all' )



    differense = angles1- angles2

    # We want the difference to be in range (-pi,pi]
    not_nan = ~np.isnan(differense)
    differense[not_nan] = np.pi - (( np.pi-differense[not_nan] ) % (2 * np.pi ))
    
    return differense





def circular_var(data,  axis=None):

    """

    Compute the circular variance along the specified axis, ignoring NaNs.
    
    Returns the variance of the array elements.  The variance is taken over
    the flattened array by default, otherwise over the specified axis.
    `float64` intermediate and return values are used for integer inputs.


    Parameters
    ----------
    data : array_like
        Array containing numbers whose circular variance is desired. If `data` is not an
        array of `numpy.float64`s, a conversion is attempted.
    axis : {int, tuple of int, None}, optional
        Axis or axes along which the means are computed. The default is to compute
        the mean of the flattened array.

    Returns
    -------
    circ_var : ndarray of `numpy.float64`s
        Returns a new array containing the circular variance values.
        Nan is returned for slices that contain only NaNs.

    """


    data = types.array(data,'data')




    # We asume period of the data is 2pi, e.g. an array of rotation angles

    data_in_complex_plane = np.exp(1j*data)
    
    with catch_warnings():
        simplefilter("ignore", category=RuntimeWarning)

        mean_of_data_in_complex_plane = np.nanmean(data_in_complex_plane, axis=axis)


    # circular variance
    return 1 - np.abs(mean_of_data_in_complex_plane)





def circular_std(data,  axis=None):

    """

    Compute the circular standard deviation along the specified axis, ignoring NaNs.
    
    Returns the standard deviation of the array elements.  The standard deviation is taken over
    the flattened array by default, otherwise over the specified axis.


    Parameters
    ----------
    data : array_like
        Array containing numbers whose circular standard deviation is desired. If `data` is not an
        array of `numpy.float64`s, a conversion is attempted.
    axis : {int, tuple of int, None}, optional
        Axis or axes along which the means are computed. The default is to compute
        the mean of the flattened array.

    Returns
    -------
    circ_std : ndarray of `numpy.float64`s
        Returns a new array containing the circular standard deviation values.
        Nan is returned for slices that contain only NaNs.


    """

    data = types.array(data,'angles')


    # We asume period of the data is 2pi, e.g. an array of rotation angles

    data_in_complex_plane = np.exp(1j*data)

    with catch_warnings():
        simplefilter("ignore", category=RuntimeWarning)
        
        mean_of_data_in_complex_plane = np.nanmean(data_in_complex_plane, axis=axis)



    mean_lenght = np.abs(mean_of_data_in_complex_plane)



    # Circular stantard deviation
    return np.sqrt( -2*np.log(mean_lenght) )



