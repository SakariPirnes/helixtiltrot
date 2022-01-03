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


# compute tilt of chain B
tilt_B = htr.tilt_angle(caB_1afo, k_1afo, mask=sse_maskB_1afo)



##########################################
##### Tilt difference between 1afo and monomer trjectories
##########################################


k_M = [0,0,1] # bilayer normal for simulations



# we use only last 2000 nanosends/frames for the analysis
use_ns = 2000


# looping over trajectories and calculating tilt difference
diff_tilt_angles = []
tilt_labels = []


for j in range(1, 5+1):
    pdb = '{}monomer.pdb'.format(j)
    xtc = '{}monomer.xtc'.format(j)
    
    ca, dssp = htr.load_ca_dssp(xtc, top=pdb)
    
    sse_mask = htr.sse_mask(dssp, sse='H')
    
    
    # we use only last 2000 nanosends/frames for the analysis
    ca = ca[-2000:]
    sse_mask = sse_mask[-2000:]
    
    
    # tilt of monomer
    tilt_M = htr.tilt_angle(ca, k_M, mask=sse_mask)
    
    
    # difference between dimer and monomer
    diff_tilt = htr.angle_diff(tilt_M,tilt_B)
    diff_tilt_angles.append(diff_tilt)
    
    #label for figure legend
    tilt_mean = htr.circular_mean(diff_tilt)*180/np.pi
    tilt_mean = '{0:.1f}'.format( tilt_mean )
    tilt_std  = htr.circular_std(diff_tilt) *180/np.pi
    tilt_std  = '{0:.1f}'.format( tilt_std )
    label    = 'j={}: $(\mu_{}, \sigma_{})=({}^\circ, \ {}^\circ)$'.format(j,j,j,tilt_mean,tilt_std)
    tilt_labels.append(label)
    


    
tilt_mean = htr.circular_mean(diff_tilt_angles)*180/np.pi
tilt_mean  ='{0:.1f}'.format(tilt_mean)



tilt_std  = htr.circular_std(diff_tilt_angles)*180/np.pi
tilt_std  ='{0:.1f}'.format(tilt_std)



fig, ax = htr.plot.angle_density(diff_tilt_angles, clim=[-30,90], plot_many=True, bins=20,figsize=(15,7))
ax.legend(loc="upper right", fontsize=20, labels=tilt_labels)



ax.set_xlabel('$\Theta^{M_j}(t)-\Theta^B$ (deg)', fontsize= 30)
ax.set_title('$(\mu_{all}, \sigma_{all})=({%s}^\circ, \ {%s}^\circ)$' %(tilt_mean,tilt_std), 
            fontsize = 40)
ax.tick_params(labelsize=20)


fig.set_size_inches(15,8)
fig.savefig('figure4dot9.png')