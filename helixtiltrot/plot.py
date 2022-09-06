import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

from helixtiltrot.core import *
from helixtiltrot import types


__all__ = [ 'rotation', "angle_map","angle_density" ]


def rotation(rot_angle,residues='all',ncols=3):

    """

    nice plott


    """
    
    rot_angle = types.array(rot_angle,'rot_angle')
    

    if rot_angle.ndim ==1:
        rot_angle = rot_angle.reshape(rot_angle.size,1)
  
    n_frame, n_res = rot_angle.shape
    
    
    # Check residues argument
    
    
    if residues == 'all':
        
        # Create array with indexes for all residues
        
        residues = np.arange(n_res)
    
    else:
        # Convert rersidues to numpy array
        try:
            residues = np.asarray(residues,dtype=int)
        except Exception as err:
            raise ValueError('\"residues\" argument should be array like and elements should be integers.\n')
        
    
    #,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
    # Statistics of residues
    
    
    rot_angle_mt = circular_mean(rot_angle, axis=0)
    
    
    Rs = 1-circular_var(rot_angle, axis=0)
    
    
    circ_std = circular_std(rot_angle, axis=0)
    #''''''''''''''''''''''''''''''''''''''''''''''''''''''''
    

    rot_angle = rot_angle.T

    
    
    #,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
    # Stuff for plotting
     
    deg_sign = '$^\circ$'
    

    nrows = int( (residues.size+ncols-1)/ncols )
    # Create figure
    fig = plt.figure(figsize=(10*ncols,10*nrows))
    fig.tight_layout(pad=8, w_pad=4, h_pad=3.7)
    #''''''''''''''''''''''''''''''''''''''''''''''''''''''''
    
    axs = []
    for i_fig, i_res in enumerate(residues):
        
        # Values for residue i_res
        angles = rot_angle[i_res]
        R = Rs[i_res]
        SD = circ_std[i_res] * 180/np.pi
        mean_angle = rot_angle_mt[i_res] * 180/np.pi
        
        
        #,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
        # Rounded values in degrees
        
        mean_angle = float( "{0:.1f}".format( mean_angle ) )
        
        SD = float( "{0:.1f}".format( SD ) )
        #''''''''''''''''''''''''''''''''''''''''''''''''''''''''
        

        #,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
        # Plot angles
        
        # Figure axis
        ax = fig.add_subplot(nrows,ncols,i_fig+1,projection='polar')
        

        # Limits
        ax.set_ylim(0,n_frame)
        

        # Time array
        time = np.arange(n_frame)
        

        # Filtter nan values
        
        nan_mask = ~np.isnan(angles)
        
        angles = angles[nan_mask]
        
        time = time[nan_mask]
        
        ax.plot(angles,time,'o',markersize=3)
        #''''''''''''''''''''''''''''''''''''''''''''''''''''''''
        
        
        #,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
        # Title and labels
        
        ax.set_title('RES{}, rotation={}$^\circ\pm${}$^\circ$'.format(i_res,mean_angle,
                                                             SD,), fontsize = 30)
        
        ax.set_xticklabels(['0'+deg_sign, '45'+deg_sign, '90'+deg_sign, '135'+deg_sign
                                , '180'+deg_sign,'-135'+deg_sign, '-90'+deg_sign, '-45'+deg_sign],
                               fontsize=20)
        
        # Time/radial labels
        label_position=ax.get_rlabel_position()
        ax.text(np.radians(label_position-3),ax.get_rmax()*1.15,'time (frame)',
                    rotation=label_position,ha='center',va='center',fontsize=20)
        #''''''''''''''''''''''''''''''''''''''''''''''''''''''''
        
        
        #,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
        # Plot mean vector
        
        
        ax.arrow(rot_angle_mt[i_res], 0, 0, R*n_frame, alpha = 0.7,
                                    width = 0.09, length_includes_head=True,
                                    head_width = 0.25, head_length = 0.2*n_frame*R,
                                    edgecolor = 'red', facecolor = 'red', lw = 2, zorder = 5)        
        #''''''''''''''''''''''''''''''''''''''''''''''''''''''''
        
        axs.append(ax)

    return fig, axs 


                        



def angle_map(angles, cmap='hsv', clim=[-180,180], figsize=(20,20)):
    
    
    
    fig, ax = plt.subplots(figsize=figsize)



    n_res = angles.shape[1]
        
    

    angles = angles.T * 180/np.pi


    
    plot = ax.pcolormesh(angles,cmap=cmap,clim=clim)
    plot.set_clim(clim)
     
    

    ax.grid(True, linewidth=1.7,color='gray')
    


    ax.set_xlabel('Time (frame)',fontsize=30)
    ax.set_ylabel('Residues',fontsize=30)
    
    

    ax.set_yticks(np.arange(n_res), minor=False)
    res_labels = [str(i) for i in range(n_res)]
    ax.set_yticklabels(res_labels, fontsize=20, va="bottom")
    
    
    
    cb = plt.colorbar(plot)
    cb.set_label("Angle (deg)", rotation=270, labelpad=30,fontsize=30)
    cb.ax.tick_params(labelsize=20)    
    

    
    return fig, ax, cb





def probability_density(data, bins,clim, normalize=False):
    
    
    clim = np.array(clim)*np.pi/180
    
    
    histo, edges = np.histogram(data, bins, density=normalize)
    
    
    
    # Mean of edges
    edges=(edges[1:]+edges[0:-1])/2
    
    
    
    function = interp1d(edges, histo, kind='cubic', bounds_error=False, fill_value=0)
    
    
    
    delta_x = np.arange(clim[0], clim[1], 0.01)
    

    N = function(delta_x)
    

    return delta_x, N



def add_density(angles, bins, clim,ax,label=1):


    # Filter nan values
    nan_mask = np.isnan(angles)
    angles = angles[~nan_mask]


    x, N = probability_density(angles, bins, clim, normalize=True)
    ax.plot(x*180/np.pi, N,label=label)





 
def angle_density(angles, bins=30, clim=[-180,180], figsize=(15,7), plot_many=False):
    
    
    


    fig, ax = plt.subplots(figsize=figsize)
    

    if plot_many:
        for i, angles_i in enumerate(angles):

            angles_i = types.array(angles_i,"angles_i")
            add_density(angles_i, bins, clim,ax,label=i+1)

    else:
        angles = types.array(angles,"angles")
        add_density(angles, bins, clim,ax)

    
    
    ax.set_xlim(clim)
    ax.set_ylim(0)
    

    
    ax.set_xlabel("Angle (deg)",fontsize = 20)
    ax.set_ylabel("angle density (1/rad)",fontsize = 20)

    return fig, ax
    
    






"""


def angle_map(angles, cmap='hsv', clim=[-180,180], figsize=(20,20)):
    
    
    
    
    n_frame, n_res = angles.shape
        
    
    # Radians to degrees. Also transpose to rotate the image
    angles = angles.T * 180/np.pi

    
    
    #,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
    # Plot
    
    
    # Create figure and axis
    fig, ax = plt.subplots(figsize=figsize)
    
    
    
    # yticks
    ax.set_yticks(np.arange(n_res), minor=False)
    
    
    
    # Plot the data
    plot = plt.pcolormesh(angles,cmap=cmap,clim=clim)
    #''''''''''''''''''''''''''''''''''''''''''''''''''''''''
     
    
    
    
    
    #,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
    # Residue labels
    
    
    res_labels = [str(i) for i in range(1,n_res+1)]

    
    
    ax.set_yticklabels(res_labels, fontsize=20, va="bottom")
    #''''''''''''''''''''''''''''''''''''''''''''''''''''''''
    
    
    plt.clim(clim)
    
    
    #,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
    # Color bar settings
    
    
    cb = plt.colorbar(plot)
    
    
    
    cb.set_label("Angle (deg)", rotation=270, labelpad=30,fontsize=30)
    
    
    
    # Change tick labelsize
    cb.ax.tick_params(labelsize=20)

    
    #''''''''''''''''''''''''''''''''''''''''''''''''''''''''
    
    
    
    
    
    ax.grid(True, linewidth=1.7,color='gray')
    
    
    
    ax.set_xlabel('Time (frame)',fontsize=30)
    ax.set_ylabel('Residues',fontsize=30)
    
    return fig, ax, cb


def orient_mean(state,title='orientation mean', saveas=False, show=True):
    
    # Check state
    state = types.state(state,'state')
    
    
    # number of frames
    n_frame = state.shape[0]
    
    
    # Mean orientation angles over residues
    state_mr = mean(state,residues=True, time=False)
    orient_mr = state_mr[:,1]
    
    #,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
    # Statistics
    
    # Mean angle
    state_mrt =  mean(state,residues=True, time=True)
    
    mean_angle = state_mrt[1]    
    
    # R is 1 - (circular variance)
    R =  Rlenghts(state, residues=True, time=True)[1]
    
    #''''''''''''''''''''''''''''''''''''''''''''''''''''''''
    
    
    #,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
    # Create figure and plot angles
    
    # Figure
    fig = plt.figure(figsize=(10,10))
    
    # Axis
    ax = fig.add_subplot(111, projection='polar')
    
    # Limits
    ax.set_ylim(0,n_frame)
    
    # Time array
    time = np.arange(n_frame)
    
    
    ax.plot(orient_mr,time,'o',markersize=4)
    #''''''''''''''''''''''''''''''''''''''''''''''''''''''''
    
    
    #,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
    # Title and labels
    
    ax.set_title(title, fontsize = 30)
    
    # Angle labels
    deg_sign = '$^\circ$'
    ax.set_xticklabels(['0'+deg_sign, '45'+deg_sign, '90'+deg_sign, '135'+deg_sign
                        , '180'+deg_sign,'-135'+deg_sign, '-90'+deg_sign, '-45'+deg_sign],
                        fontsize=20)
    
    
    label_position=ax.get_rlabel_position()
    ax.text(np.radians(label_position-3),ax.get_rmax()*1.15,'time (frame)',
            rotation=label_position,ha='center',va='center',fontsize=20)
    #''''''''''''''''''''''''''''''''''''''''''''''''''''''''
    
    #,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
    # Plot mean vector 
    
    ax.arrow(mean_angle, 0, 0, R*n_frame, alpha = 0.7,
                            width = 0.09, length_includes_head=True,
                            head_width = 0.25, head_length = 0.2*n_frame*R,
                            edgecolor = 'red', facecolor = 'red', lw = 2, zorder = 5)
    
    #''''''''''''''''''''''''''''''''''''''''''''''''''''''''
    
    
    if saveas != False:
        fig.savefig(saveas)
    
    if show == True:
        plt.show()
    else:
        plt.close()
"""