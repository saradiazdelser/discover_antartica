#!/usr/bin python3

# --------------------------------------------------------------
# Graphs Module (5/5)
# --------------------------------------------------------------

# ------------------------- WARNINGS!---------------------------
# This program requires Python 3 to work——hence, the '.py' extension. 
#+If you don't have Python 3, please crawl out of whichever rock 
#you've been living under and install it.

# This program requires Biopython to work. If you don't have Biopython  
#+installed, what are you waiting for? A sign from God? This is it.

# This program requires MatPlotLib to work. If you don't have it 
#+installed, check out how to here: 
# https://pypi.org/project/matplotlib/ 


# -------------------------- MODULES ---------------------------
import matplotlib.pyplot as plt
import csv
from Bio import Phylo


# ------------------------- FUNCTIONS --------------------------

# This one creates a graphic out of a blast_results file and 
#+displays it in matplotlib.
def the_perks_of_being_a_graphic(file):
	# Open the file
    with open(file, 'r') as f:

    	# Create lists
        cov=[]
        seq=[]
        ident=[]

        # Parse through the file, except header
        blast_results= csv.reader(f, delimiter='\t')
        next(blast_results)

        # Create a list with each value type
        for row in blast_results:
            if seq not in seq:
                cov.append(float(row[2]))
                seq.append(str(row[1]))
                ident.append(float(row[3]))

    N=len(cov)
    fig = plt.figure()
    # Adds axes
    ax = fig.add_axes([0.2,0.05,0.75,0.85])

    # Design legend
    custom_lines = [plt.Line2D([0], [0], color='r', lw=4),
                    plt.Line2D([0], [0], color='m', lw=4),
                    plt.Line2D([0], [0], color='y', lw=4)]
    # Depending on the identity, create a different colored
    #+horizontal bar. 
    # Bar length = coverage
    for i in range(N):
        if ident[i]<40:
            ax.barh(seq[i], cov[i], color='y')
        elif ident[i]<80:
            ax.barh(seq[i], cov[i], color='m')
        else:
            ax.barh(seq[i], cov[i], color='r')

    # Display legend
    ax.legend(custom_lines, ['ident>80', '40<ident<80', 'ident<40'])

    # Make y labels (sequence names) smaller
    plt.tick_params(axis='y', direction='in', labelsize='xx-small')
    plt.xlabel('Coverage')

    # Add title
    fig.suptitle('BLAST Results\n'+\
    	'This barplot offers a simple representation of coverage and '+\
    	'identity of every hit according to given cut-off values.',\
    	 fontsize=10)

    # I feel it's better if the user saves the graphic manually
    # But this command will save it automatically:
    #fig.savefig(file+'_grafic.jpg')

    # Show
    plt.show()


# This one shows the tree in matplotlib. Yay!
def you_cant_handle_the_tree(tree_file):

    # Open the tree and show it to the world!
    tree = Phylo.read(tree_file, "newick")
    Phylo.draw(tree)

