import helixtiltrot as htr
import numpy as np



##########################################
########  Experimental dimer
##########################################


pdb_1afo = '1afo.pdb'

ca_1afo, dssp_1afo = htr.load_ca_dssp(pdb_1afo)

# take the model 1
ca_1afo = ca_1afo[:1]
dssp_1afo = dssp_1afo[:1]


sse_mask_1afo = htr.sse_mask(dssp_1afo, sse='H')


# number of residues in chains A and B
n_res = int(ca_1afo.shape[1]/2) # 40


# Split into chain A and B
caA_1afo = ca_1afo[:,:n_res]
caB_1afo = ca_1afo[:,n_res:]
sse_maskA_1afo = sse_mask_1afo[:,:n_res]
sse_maskB_1afo = sse_mask_1afo[:,n_res:]


# use c2 symmetry of 1afo to estimate the bilayer normal k_1afo
H_A = htr.axis(caA_1afo, mask = sse_maskA_1afo)
H_B = htr.axis(caB_1afo, mask = sse_maskB_1afo)
k_1afo = (H_A+H_B)[0] # no need to normalize


# compute local rotation of chain B
local_rot_B = htr.local_rotation_angle(caB_1afo, k_1afo, mask=sse_maskB_1afo)



##########################################
##### rotation difference between 1afo and monomer trjectories
##########################################


k_M = [0,0,1] # bilayer normal for simulations

# we use only last 2000 nanosends/frames for the analysis
use_ns = 2000


# looping over trajectories and calculating rotation difference
diff_rot_angles = []
rot_labels = []

for j in range(1, 5+1):
    pdb = '{}monomer.pdb'.format(j)
    xtc = '{}monomer.xtc'.format(j)
    
    ca, dssp = htr.load_ca_dssp(xtc, top=pdb)
    
    sse_mask = htr.sse_mask(dssp, sse='H')
    
    
    # we use only last 2000 nanosends/frames for the analysis
    ca = ca[-2000:]
    sse_mask = sse_mask[-2000:]
    
    
    # tilt of monomer
    local_rot_M = htr.local_rotation_angle(ca, k_M, mask=sse_mask)
    
    
    # difference between dimer and monomer
    diff_local_rot = htr.angle_diff(local_rot_M,local_rot_B)
    diff_rot = htr.circular_mean(diff_local_rot,axis=1)
    diff_rot_angles.append(diff_rot)
    
    
    #label for figure legend
    rot_mean = htr.circular_mean(diff_rot)*180/np.pi
    rot_mean = '{0:.1f}'.format( rot_mean )
    rot_std  = htr.circular_std(diff_rot) *180/np.pi
    rot_std  = '{0:.1f}'.format( rot_std )
    label    = 'j={}: $(\mu_{}, \sigma_{})=({}^\circ, \ {}^\circ)$'.format(j,j,j,rot_mean,rot_std)
    rot_labels.append(label)
    

diff_rot_angles = np.asarray(diff_rot_angles)


rot_mean = htr.circular_mean(diff_rot_angles)*180/np.pi
rot_mean  ='{0:.1f}'.format(rot_mean)



rot_std  = htr.circular_std(diff_rot_angles)*180/np.pi
rot_std  ='{0:.1f}'.format(rot_std)



diff_rot_angles[diff_rot_angles<0] = 2*np.pi + diff_rot_angles[diff_rot_angles<0]
fig, ax = htr.plot.angle_density(diff_rot_angles, clim=[0,360], plot_many=True, bins=20,figsize=(15,7))


ax.legend(loc="upper right", fontsize=20, labels=rot_labels)


ax.set_xticks([0,45,90,135,180,360-135,360-90,360-45,360])
ax.set_xticklabels([0,45,90,135,180,-135,-90,-45,0],fontsize=15)


ax.set_xlabel('$\Theta^{M_j}(t)-\Theta^B$ (deg)', fontsize= 30)
ax.set_title('$(\mu_{all}, \sigma_{all})=({%s}^\circ, \ {%s}^\circ)$' %(rot_mean,rot_std), 
            fontsize = 40)
ax.tick_params(labelsize=20)
fig.set_size_inches(15,8)


fig.savefig('figure4dot10.png')