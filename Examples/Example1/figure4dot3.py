import helixside
from matplotlib.pyplot import subplots
import numpy as np

# input coordinate files
pdb = 'transmembrane-alpha-helix.pdb'
xtc = 'transmembrane-alpha-helix.xtc'



ca, dssp = helixside.load_ca_dssp(xtc, top=pdb)


# create numpy mask called sse_mask, such that
# sse_mask[n,i] == True if the residue i at time step n
# is either alpha-helix (H), pi-helix (I) or 310-helix
sse_mask = helixside.sse_mask(dssp, sse='HIG')



# membrane normal
k = [0,0,1] 


# tilt of residue 3,4,...,14 and 16,...,26
tilt_3to14 = helixside.tilt_angle(ca, ref_vec=k, mask=sse_mask, m=3,n=15)
tilt_16to26 = helixside.tilt_angle(ca, ref_vec=k, mask=sse_mask, m=16,n=27)


# kink angle between residues 3,4,...,14 and 16,...,26
kink_angle = helixside.kink_angle(ca, mask=None, m1=3, n1=15, m2=16, n2=27)



# Change angles from radians to degrees
tilt_3to14  *= 180/np.pi
tilt_16to26 *= 180/np.pi
kink_angle  *= 180/np.pi



# time array for the plot
time = np.arange(kink_angle.shape[0])


# Plot
fig, ax = subplots(figsize=(20,10))


ax.plot(time,tilt_3to14, "blue", alpha=0.25)
ax.plot(time,tilt_16to26, "green",alpha=0.25)
ax.plot(time,kink_angle, "red", alpha=0.25)



# To save the precice image used in the thesis, uncomment 
# the following lines to edit the figure accordingly:

"""
def rolling_mean(data, window=50, center = True, fill=np.nan):
    # function to compute rolling mean
    
    from math import floor, ceil
    if(center):
        start_offset = floor(window/2)
        end_offset   = -ceil(window/2)+1
        if(end_offset==0): end_offset=None
    else:
        start_offset = window-1
        end_offset   = None

    window-=1
    if(window==0): window=None
    cumsum = np.nancumsum(data, axis=0, dtype=float)
    mean   = np.full_like(cumsum, fill)
    mean[start_offset:end_offset] = (cumsum[window:]-cumsum[:-window])/window

    return mean
"""

#label3_14 = r'$\Theta_{5,14}$'
#label15_26 = r'$\Theta_{16,26}$'
#label_kink = r'$\xi_{(5,14),(16,26)}$'


"""
ax.plot(time,rolling_mean(tilt_3to14), 'blue',label=label3_14)
ax.plot(time, rolling_mean(tilt_16to26), 'green', label=label15_26)
ax.plot(time, rolling_mean(kink_angle), 'red',label=label_kink)
ax.legend(loc='upper left', fontsize=30)

ax.set_xlim(0,time[-1])
ax.set_ylim(0)

ax.tick_params(labelsize=25)

ax.set_xlabel('Time (ns)', fontsize=30)
ax.set_ylabel('Angle (deg)', fontsize=30)
"""


fig.savefig('figure4dot3.png')



