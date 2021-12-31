import numpy as np


#from helixstate import types
 

__all__ = ['dot','normalize','h_to_centerpoints','h_to_tiltvecs','h_to_rotvecs']




def dot( arr3d, vec ):
    """ 
    arr3d is numpy array with shape (m,n,3) and
    vec is numpy array with shape (3,) or (m,n,3).
    
    Returns the dot product between vec and all elements on
    the last axis of arr3d

    Outputs array with shape (m,n)

    
    """
    
    d = arr3d*vec
    
    return d[:,:,0]+d[:,:,1]+d[:,:,2]




def normalize(arr):


    arr_norm = (arr[:,:,0]**2+arr[:,:,1]**2+arr[:,:,2]**2)**0.5


    arr_norm[ arr_norm == 0 ] = np.nan


    arr_normalized = ( arr.T / arr_norm.T ).T


    return arr_normalized





def circMid(nArr,Aarr,Barr,Carr):       
    # Aarr,Barr,Carr are arrays containing points in a cirles in 3d
    # nArr is array containing normals to the plane where A,B,C are
    
    A1, A2, A3 = Aarr.T
    B1, B2, B3 = Barr.T
    C1, C2, C3 = Carr.T
    n1,n2,n3 = nArr.T
    
    BA1,BA2,BA3 = 2*(B1-A1), 2*(B2-A2), 2*(B3-A3)
    CA1,CA2,CA3 = 2*(C1-A1), 2*(C2-A2), 2*(C3-A3)
    
    
    DETinv = 1/(n3*(CA1*BA2 - CA2*BA1) + n2*(CA3*BA1 - CA1*BA3) + n1*(CA2*BA3 - CA3*BA2))
    
    AA = A1**2 + A2**2 + A3**2
    BB = B1**2 + B2**2 + B3**2
    CC = C1**2 + C2**2 + C3**2
    
    
    
    MID = DETinv*np.array([(CA3*n2 - CA2*n3)*(BB-AA) + (BA2*n3 - BA3*n2)*(CC-AA) + (CA2*BA3 - CA3*BA2)*(n1*A1+n2*A2+n3*A3),
    (CA1*n3 - CA3*n1)*(BB-AA) + (BA3*n1 - BA1*n3)*(CC-AA) + (CA3*BA1 - CA1*BA3)*(n1*A1+n2*A2+n3*A3),
    (CA2*n1 - CA1*n2)*(BB-AA) + (BA1*n2 - BA2*n1)*(CC-AA) + (CA1*BA2 - CA2*BA1)*(n1*A1+n2*A2+n3*A3)])
    
    return MID.T








def h_to_centerpoints(h,ca):


    # four consecutive ca atoms
    r1 = ca[:,:-3]
    r2 = ca[:,1:-2]
    r3 = ca[:,2:-1]
    r4 = ca[:,3:]


    # Calculate the center point of the cirle
    # Lenth of circle
    
    
    L = (r1-r4) * h
    L = L[:,:,0]+L[:,:,1]+L[:,:,2] 
    
    # Project points, points in the cirle are ri_c
    h_dist_1 = -L/2
    h_dist_2 = h_dist_1 + np.sum( (r1 - r2)*h , axis=-1 )
    h_dist_3 = h_dist_1 + np.sum( (r1 - r3)*h , axis=-1 )
    h_dist_4 = h_dist_1 + np.sum( (r1 - r4)*h, axis=-1 )
    


    r1C = r1 + (h_dist_1.T*h.T).T
    r2C = r2 + (h_dist_2.T*h.T).T
    r3C = r3 + (h_dist_3.T*h.T).T
    r4C = r4 + (h_dist_4.T*h.T).T
    



    # Calculate mean of points for the center of the cirle
    m = ( circMid(h,r1C,r2C,r3C)+circMid(h,r1C,r3C,r4C)+circMid(h,r1C,r2C,r4C)+circMid(h,r2C,r3C,r4C) )/4
        
    return m




def h_to_tiltvecs(h):

    

    # Sum between consecutive h is the tilt vector
    t_vecs = h[:,:-1] + h[:,1:]

    t_vecs = normalize( t_vecs )



    end = t_vecs.shape[1]


    t_vecs = np.insert(t_vecs, (0,0,end,end), np.nan, axis=1)


    return t_vecs





def h_to_rotvecs(h,ca):



    m = h_to_centerpoints(h,ca)


    # Tilt vectors
    t_vecs = h_to_tiltvecs(h)


    #,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
    # Calculate rotation vectors

    M = (m[:,1:]+m[:,:-1])/2
    end = M.shape[1]
    M = np.insert( M, (0,0,end,end), np.nan, axis=1 )
    
    
    # It is in the same plane with the vector between alpha carbon atom and middle point
    r_vecs = ca - M
    
    
    # Now we just need to make it perpendicular to the tilt vector
    r_vecs = r_vecs - ( np.sum(t_vecs*r_vecs, axis=-1).T * t_vecs.T ).T
    

    r_vecs = normalize( r_vecs )

    #''''''''''''''''''''''''''''''''''''''''''''''''''''''''
    
    return r_vecs
    




    


#def statistics(data, np_function, over_time, over_residues):
#
#
#
#    # Atleast few first and last residues are nan and this would
#    # cause RuntimeWarning so we ignore it
#
#
#    with catch_warnings():
#
#        simplefilter("ignore", category=RuntimeWarning)
#        
#
#        axis = ()
#
#        if over_residues: axis += (1,)
#
#  
#        if over_time: axis += (0,)
#
#        
#        data = np_function(data, axis=axis)
#
#
#    return data
#

