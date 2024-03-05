import numpy as np


__all__ = ['array','state','boolean']






def array(arg, arg_name, dtype=np.float64):

    

    try:
        arg = np.asarray(arg, dtype = dtype )
    

    except ValueError as err:
        raise ValueError('\"{}\" argument could not be converted to a numpy array: {}'.format(arg_name,err)) from None

    
    except Exception as err:
        raise err

    else:
        return arg







def boolean(arg,arg_name):

    
    if isinstance(arg, bool ):
        return arg
    

    else:
        raise TypeError('Boolean value expected for argument \"{}\". You supplied {}'.format(arg_name,type(arg)))




def equal_shape(shape1,arg1_name, shape2,arg2_name, axes_to_equal,  words = ['frames', 'residues', 'cartesian coordinates']):

    #shape1 = arg1.shape
    #shape2 = arg2.shape

    if shape1 != shape2:

        if axes_to_equal == 'all':

            msg = 'arguments \"{}\" and \"{}\" should have same shape: {}!={}'.format(arg1_name,arg2_name,shape1,shape2)
            raise ValueError(msg)
            
        else:


            for axis in axes_to_equal:
                if shape1[axis] != shape2[axis]:

                    msg = 'arguments \"{}\" and \"{}\" should have same number of {}: {}!={}'.format(arg1_name,arg2_name,words[axis],shape1[axis],shape2[axis])
                    raise ValueError(msg)




def dssp(arg, arg_name, command_line = False):

    arg =  array(arg,arg_name, dtype = str)

    return arg



def sse(arg, arg_name, command_line = False):

    if not isinstance( arg, str ):
        msg = 'String expected for argument \"{}\". You supplied {}'.format(arg_name,type(arg))
        raise ValueError(msg)

    return arg



def ca(arg,arg_name, ndim_res_gt5 = True, command_line = False):
    
    
    # ca to numpy array
    arg = array(arg,arg_name)

    
    # number of dimensions for residue axis greater than five aka ndim_res_gt5
    ndim_res_gt5 = boolean(ndim_res_gt5,'ndim_res_gt5')
    
    
    # ca should have three dimensions, time, residues and cartesian coordinates.
    ndim = arg.ndim


    if ndim != 3:



        if command_line: 
            msg = 'In the file between each two \'|\' symbols there should be three values representing cartesian coordinates x,y and z. You have only one.'


        else:
            msg = '\"{}\" argument should have three axes. You supplied {} axes.'.format(arg_name,ndim)
        

        raise ValueError(msg)


    else:


        
        # Checking the shape is correct
        n_frame, n_res, n_xyz = arg.shape

        
        # Making sure we are in the boring 3d world.
        if n_xyz != 3:




            if command_line: 
                msg = 'In the file between each two \'|\' symbols there should be three values representing cartesian coordinates x,y and z. You have {} values.'.format(n_xyz)


            else:
                msg = '\"{}\" arguments third axis should have three elements representing cartesian coordinates x,y and z. You supplied {} elements.'.format(arg_name,n_xyz)


            raise ValueError(msg)
        



        if ndim_res_gt5 and n_res < 5:

            raise ValueError('The used method needs atleast 5 residues. You supplied {} residues in {} argument.'.format(n_res,arg_name))

        
    return arg





def vec3d(arg, arg_name, command_line = False ):

    arg = ca(arg, arg_name, ndim_res_gt5 = False, command_line = command_line )

    return arg






xyzbase = {
    'x'  : np.array( [1,  0,  0] ),
    'y'  : np.array( [0,  1,  0] ),
    'z'  : np.array( [0,  0,  1] ),
    '-z' : np.array( [0,  0, -1] ),
    '-y' : np.array( [0, -1,  0] ),
    '-x' : np.array( [-1, 0,  0] )
}






def k(arg, data_shape, arg_name='ref_vec'):
    


    if isinstance(arg,str):
        
        try:
            arg = xyzbase[arg]
            
        except KeyError:
            raise ValueError('Expected "x", "y", "z", "-z", "-y" or "-x" for argument \"{}\". You supplied \"{}\".'.format(arg_name,arg))
        
    else:
        
        # Convert to numpy array
        arg = array(arg, arg_name)
        
        
        # Making sure the shape is correct.
        
        ndim = arg.ndim
        
        if ndim != 1 and ndim != 3 and ndim != 2:
            
            raise ValueError('\"{}\" argument should have either 1, 2 or 3 axes. You supplied {} axises.'.format(arg_name,ndim))

        elif ndim == 3 and len(data_shape) == 2:
            raise ValueError('\"{}\" argument should have either 1 or 2  axes. You supplied {} axises.'.format(arg_name,ndim))
        elif ndim == 2 and len(data_shape) == 3:

            arg_reshaped = np.empty((arg.shape[0], data_shape[1], arg.shape[1]))

            for i in range(data_shape[1]):
                arg_reshaped[:,i] = arg

            arg = arg_reshaped
            ndim+=1
            
        if ndim == 1:
            
            # Making sure we are in the boring 3d world.
            
            n_xyz = arg.shape[0]
            
            if n_xyz != 3:
                raise ValueError('\"{}\" arguments should have three elements representing cartesian coordinates x,y and z. You supplied {} elements.'.format(arg_name,n_xyz))

            else:
                
                # Making sure we don't have a null vector
                
                arg_lenght = ( arg[0]**2+arg[1]**2+arg[2]**2)**0.5
                
                if arg_lenght == 0:
                    
                    raise ValueError('Vector with length 0 is not valid for \"{}\" argument.'.format(arg_name))
                
                else:
                    
                    # Normalize
                    arg = arg/arg_lenght


        elif ndim == 3:
            equal_shape(arg.shape, arg_name, data_shape,'ca', 'all')

            arg_lenght = (arg[:,:,0]**2+arg[:,:,1]**2+arg[:,:,2]**2)**0.5

            if np.any(arg_lenght==0):
                raise ValueError('Array containing a vector with length 0 is not valid for \"{}\" argument.'.format(arg_name))
            else:
                arg = (arg.T/arg_lenght.T).T
        elif ndim == 2:

            words = ['frames', 'cartesian coordinates']
            equal_shape(arg.shape, arg_name, data_shape,'axis', [0,1], words=words)

            arg_lenght = (arg[:,0]**2+arg[:,1]**2+arg[:,2]**2)**0.5

            if np.any(arg_lenght==0):
                raise ValueError('Array containing a vector with length 0 is not valid for \"{}\" argument.'.format(arg_name))
            else:
                arg = (arg.T/arg_lenght.T).T


        
    return arg







