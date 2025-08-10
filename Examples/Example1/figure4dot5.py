import helixtiltrot as htr


# input coordinate files
pdb = 'transmembrane-alpha-helix.pdb'
xtc = 'transmembrane-alpha-helix.xtc'



ca, dssp = htr.load_ca_dssp(xtc, top=pdb)



# create numpy mask called sse_mask, such that
# sse_mask[n,i] == True if the residue i at time step n
# is either alpha-helix (H), pi-helix (I) or 310-helix
sse_mask = htr.sse_mask(dssp, sse='HIG')



# membrane normal
k = [0,0,1] 



local_side = htr.local_side_angle(ca, ref_vec=k, mask=sse_mask)


# Change local side angles into phase of residue 10
local_side_10 = htr.single_phase(local_side,phase=10, turn_angle_deg=100.1 )


fig, ax, cb = htr.plot.angle_map(local_side_10)




# To save the precice figure used in the thesis, uncomment 
# the following lines to edit the figure accordingly:

"""
ax.set_xlabel('Time (ns)')
ax.tick_params(labelsize=25)

cb.set_label('Local side angle (deg)')
cb.ax.tick_params(labelsize=25)
"""


fig.savefig('figure4dot5.png')
