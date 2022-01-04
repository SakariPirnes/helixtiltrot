import numpy as np
from mdtraj import load, compute_dssp
from helixtiltrot.heart import *
from helixtiltrot import types
from warnings import catch_warnings, simplefilter


 

__all__ = ['load_ca_dssp','local_axes', 'axis', 'kink_angle','center_points','tilt_vectors','rotation_vectors','tilt_angle', 'local_tilt_angle','local_rotation_angle','rotation_angle','single_phase','sse_mask','angle_diff','circular_mean','circular_var','circular_std']




######################################################################################
#                            Functions for user
######################################################################################



def load_ca_dssp(fname,top=None):
    
    
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

    dssp = types.dssp(dssp, 'dssp')
    sse = types.sse(sse,'sse')

    


    masks = []
    for s in sse:
        masks.append(dssp == s)



    return np.any(masks, axis=0)




def local_axes(ca, mask=None):


    ca  = types.ca(ca,'ca').copy()


    if mask is not None:

        mask = types.array(mask.copy(),'mask', dtype = bool)
        types.equal_shape(ca,'ca', mask, 'mask', axes_to_equal=[0,1] )

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
    #s1 = 2 * r2 - r1 - r3
    #s2 = 2 * r3 - r2 - r4
    s1 = r21 - r32
    s2 = r32 - r43
        
        
    # Local axes
    h = np.cross( s2, s1 )
        
    # return normalized local axes
    return normalize( h )


def axis(ca, mask=None, m=0, n=None):

    ca  = types.ca(ca,'ca').copy()

    if mask is not None:

        mask = types.array(mask.copy(),'mask', dtype = bool)
        types.equal_shape(ca,'ca', mask, 'mask', axes_to_equal=[0,1] )

        h = local_axes(ca[:,m:n], mask=mask[:,m:n])
    else:
        h = local_axes(ca[:,m:n])


    H_mn = np.nansum([h], axis=2)

    return normalize( H_mn )[0]




def kink_angle(ca, m1, n1, n2, m2, mask=None):

    H1 = axis(ca, mask=mask, m=m1, n=n1)
    H2 = axis(ca, mask=mask, m=m2, n=n2)

    cosine_of_kink_angle = np.clip( (H1*H2).sum(axis=1), -1.0, 1.0)

    return np.arccos(cosine_of_kink_angle)



def tilt_angle(ca, ref_vec, mask=None, m=0, n=None):

    H_mn = axis(ca, mask=mask, m=m, n=n)
    ref_vec = types.k(ref_vec)


    # we change H_mn to [H_mn] since the dot function is dumm
    H_mn = H_mn.reshape(1,*H_mn.shape)

    cosine_of_tilt_angle = np.clip( dot(H_mn,ref_vec), -1.0, 1.0)[0]

    return np.arccos( cosine_of_tilt_angle )








def tilt_vectors(ca, mask=None):


    h = local_axes(ca, mask=mask)


    return h_to_tiltvecs(h)





def local_tilt_angle(ca, ref_vec, mask=None):


    ref_vec = types.k(ref_vec)


    t_vec = tilt_vectors(ca, mask=mask)


    cosine_of_local_tilt_angle = np.clip( dot(t_vec,ref_vec), -1.0, 1.0)

    
    return np.arccos( cosine_of_local_tilt_angle )




def center_points(ca, mask=None):


    h = local_axes(ca, mask=mask)


    return h_to_centerpoints(h,ca)





def rotation_vectors(ca, mask=None):


    h = local_axes(ca, mask=mask)


    return h_to_rotvecs(h,ca)





def local_rotation_angle( ca, ref_vec, mask=None ):


    k = types.k(ref_vec)

    h = local_axes(ca, mask=mask)



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
        Array containing local rotation angles. If `data` is not an
        array of `numpy.float64`s, a conversion is attempted.
    turn_angle_deg : float or int
        The turn angle in degrees.
    phase : int
        Desired phase for `local_rot`

    Returns
    -------
    m : ndarray of `numpy.float64`s
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
    m : ndarray of `numpy.float64`s
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

    Shift local rotation angles in radians into local rotation angles in phase
    of residue `phase`, in radians .

    Returns a new array of local rotation angles in phase of residue `phase`, 
    in radians, in range ]-pi,pi].
    

    Parameters
    ----------
    local_rot : array_like
        Array containing local rotation angles. If `data` is not an
        array of `numpy.float64`s, a conversion is attempted.
    turn_angle_deg : float or int
        The turn angle in degrees.
    phase : int
        Desired phase for `local_rot`

    Returns
    -------
    rot : ndarray of `numpy.float64`s
        Returns a new array of local rotation angles in phase of residue `phase`, 
        in radians, in range ]-pi,pi].


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
    m : ndarray of `numpy.float64`s
        Returns a new array containing the difference angles.

    """

    
    angles1 = types.array(angles1, 'angles1')
    angles2 = types.array(angles2, 'angles2')

    #types.equal_shape(angles1,'angles1', angles2, 'angles2', axes_to_equal='all' )



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
    m : ndarray of `numpy.float64`s
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
    `float64` intermediate and return values are used for integer inputs.


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
    m : ndarray of `numpy.float64`s
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



