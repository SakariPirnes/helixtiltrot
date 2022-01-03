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



local_tilt = htr.local_tilt_angle(ca, ref_vec=k, mask=sse_mask)


# plot the colormap
fig, ax, cb = htr.plot.angle_map(local_tilt, clim=[0,90], cmap="jet")




# To save the same labels used in the thesis, uncomment 
# the following lines to edit the figure accordingly:


# The objects fig, ax, cb are the maplotlib figure, axis and colorbar, 
# respectively. These can be edited further, e.g. adding or changing titles
# and labels

"""
ax.set_xlabel('Time (ns)')
ax.tick_params(labelsize=25)
ax.set_title('b)'+' '*63, fontsize=50)

cb.set_label('Local tilt angle (deg)')
cb.ax.tick_params(labelsize=25)
"""

fig.savefig('figure4dot2b.png')