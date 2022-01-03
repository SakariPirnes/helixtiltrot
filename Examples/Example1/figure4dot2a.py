import helixtiltrot as htr


# input coordinate files
pdb = 'transmembrane-alpha-helix.pdb'
xtc = 'transmembrane-alpha-helix.xtc'



ca, dssp = htr.load_ca_dssp(xtc, top=pdb)



# membrane normal
k = [0,0,1] 



local_tilt = htr.local_tilt_angle(ca, ref_vec=k)


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
ax.set_title('a)'+' '*63, fontsize=50)

cb.set_label('Local tilt angle (deg)')
cb.ax.tick_params(labelsize=25)
"""

fig.savefig('figure4dot2a.png')
