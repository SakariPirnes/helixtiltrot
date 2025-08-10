import helixtiltrot as htr


# input coordinate files
pdb = 'transmembrane-alpha-helix.pdb'
xtc = 'transmembrane-alpha-helix.xtc'



ca, dssp = htr.load_ca_dssp(xtc, top=pdb)



# create numpy mask called sse_mask, such that
# sse_mask[n,i] == True if the residue i at time step n
# is either alpha-helix (H), pi-helix (I) or 310-helix
sse_mask = htr.sse_mask(dssp, sse='HIG')



kappa = [1,1,0] # orthogonal vector to the membrane normal [0,0,1]



l = 10
turn_angle_deg=100.1
side_kappa_15to26 = htr.side_angle(ca, kappa, turn_angle_deg, l, mask=sse_mask, m=15,n=26 )



fig, axs = htr.plot.polar(side_kappa_15to26, ncols=1)



# To save the precice figure used in the thesis, uncomment 
# the following lines to edit the figure accordingly:

"""
avg_and_std = axs[0].get_title().split("=")[1]
new_title = r"$\Phi_{15,26}^{%s}(\hat{\mathbf{\kappa}}) = $%s" % (l, avg_and_std)
axs[0].set_title(new_title, fontsize=35)
"""


fig.savefig('figure4dot7.png')
