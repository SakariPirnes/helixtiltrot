import helixtiltrot as htr


# input coordinate files
pdb = 'transmembrane-alpha-helix.pdb'
xtc = 'transmembrane-alpha-helix.xtc'



ca, dssp = htr.load_ca_dssp(xtc, top=pdb)



# membrane normal
k = [0,0,1] 



local_side = htr.local_side_angle(ca, ref_vec=k)



# Plotting residues 6,...11. 
# Note that first residue is indexed by 0, because python indexing starts from 0.
residues = range(6,11+1)



# Plot local side angles for each residue.
fig, axs = htr.plot.polar(local_side,residues=residues, ncols=3)



# To save the precice figure used in the thesis, uncomment 
# the following lines to edit the figure accordingly:

"""
for res, ax in zip(residues, axs):
    avg_and_std = ax.get_title().split('=')[1]
    new_title = r'$\phi_{%s} = $%s' % (res, avg_and_std)
    ax.set_title(new_title, fontsize=35)
"""



fig.savefig('figure4dot4.png')
