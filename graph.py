 #!/usr/bin python3

# --------------------------------------------------------------
# Graphs Module (5/5)
# --------------------------------------------------------------

# ------------------------- WARNINGS!---------------------------
# This program requires Python 3 to work——hence, the '.py' extension. If you don't have Python 3, 
#+ please crawl out of whichever rock you've been living under and install it.

# This program requires Biopython to work. If you don't have Biopython installed, 
#+what are you waiting for? A sign from God? This is it.

# This program requires MatPlotLib to work. If you don't have it installed, check out
#+how to here: https://pypi.org/project/matplotlib/ 


# -------------------------- MODULES ---------------------------
import matplotlib.pyplot as plt
import numpy as np
import csv
from Bio import Phylo


# ------------------------- FUNCTIONS --------------------------

def the_perks_of_being_a_graphic(file):

    with open(file, 'r') as f:

        cov=[]
        seq=[]
        ident=[]

        blast_results= csv.reader(f, delimiter='\t')
        next(blast_results)

        for row in blast_results:
            cov.append(float(row[2]))
            seq.append(str(row[1]))
            ident.append(float(row[3]))

    N=len(cov)
    ind = np.arange(N)
    fig = plt.figure()
    ax = fig.add_axes([0,0,1,1])
    for i in range(N):
        if ident[i]<43:
            ax.barh(seq[i], cov[i], color='y')
            ax.set_yticklabels(seq[i], minor=False)
        elif ident[i]<46:
            ax.barh(seq[i], cov[i], color='m')
            ax.set_yticklabels(seq[i], minor=False)
        else:
            ax.barh(seq[i], cov[i], color='r')
            ax.set_yticklabels(seq[i], minor=False)
    
    # Show sequences names
    

    # Create legend
    ax.legend(labels=['I>80', '40<I<80', 'I<40'])

    # Title axes and figure
    ax.set_xlabel('Sequences')
    ax.set_ylabel('Coverage')
    ax.set_title('Bar Plot of BLAST results')

    # Show
    plt.show()


# This one shows the tree in the terminal. Yay!
def you_cant_handle_the_tree(tree_file):

    # Open the tree and show it to the world!
    tree = Phylo.read(tree_file, "newick")
    Phylo.draw(tree)

