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



l = 10
rot_5to14 = htr.rotation_angle(ca, k, mask=sse_mask, phase=l, turn_angle_deg=100.1, m=5,n=15 )
rot_15to26= htr.rotation_angle(ca, k, mask=sse_mask, phase=l, turn_angle_deg=100.1, m=15,n=27 )



# Merge rotation angles rot_5to14 and rot_15to26
rotation_angles = np.array([rot_5to14,rot_15to26]).T



# Plot rotaion angles
fig, axs = htr.plot.rotation(rotation_angles, ncols=2)



# To save the precice figure used in the thesis, uncomment 
# the following lines to edit the figure accordingly:

"""
avg_and_std = axs[0].get_title().split("=")[1]
new_title = r"$\Phi_{5,14}^{%s} = $%s" % (l, avg_and_std)
axs[0].set_title(new_title, fontsize=35)

avg_and_std = axs[1].get_title().split("=")[1]
new_title = r"$\Phi_{15,26}^{%s} = $%s" % (l, avg_and_std)
axs[1].set_title(new_title, fontsize=35)
"""



fig.savefig('figure4dot6.png')



