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

def the_perks_of_being_a_graphic(file='data/filtered_results_file.tsv'):
    # Fixing random state for reproducibility
    np.random.seed(19680801)

    coverage=[]
    sequences=[]
    identity=[]

    with open(file, 'r') as f:

        plots= csv.reader(f, delimiter='\t')
        next(plots)

        for row in plots:
            coverage.append(float(row[2]))
            sequences.append(str(row[1]))
            identity.append(float(row[3]))

    plt.rcdefaults()
    fig, ax = plt.subplots()

    mask1 = (identity < 0.5)
    mask2 = (identity >= 0.5)

    # y_value represents bar width
    y_pos = np.arange(len(sequences))

    ax.barh(y_pos[mask1], coverage[mask1], identity[mask1],  align='center', color='red')
    ax.barh(y_pos[mask2], coverage[mask2], identity[mask2], align='center', color='pink')

    ax.set_yticks(y_pos)
    ax.set_yticklabels(sequences)
    ax.invert_yaxis()  # labels read top-to-bottom
    ax.set_xlabel('coverage')
    ax.set_title('BLAST results')

    plt.show()


# This one shows the tree in the terminal. Yay!
def you_cant_handle_the_tree(tree_file):

    # Open the tree and show it to the world!
    tree = Phylo.read(tree_file, "newick")
    Phylo.draw(tree)

