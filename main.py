#!/usr/bin python3

# -*- coding: utf-8 -*-

# --------------------------------------------------------------
# Discover Antartica 2.0
# --------------------------------------------------------------
# This wonderful, looong program has the ability do lots of interesting
#+things, including giving the coder (me!) grey hair.
# Here's what it'll do:

#    1. Run BLAST on a number of chosen query proteins, given a certain
#+    evalue, coverage and identity cut-off values. This will be done 
#+    locally, so you'll need a database. 
#+    It will show you the results using a pretty graphic.

#    2. Use MUSCLE to multiple-align the proteins aligned by BLAST.
#+    Then make a phylogenetic tree using the NJ method.
#+    It also show you tree using graphic.

#    3. Run all your proteins against a prosite database in order to
#+    identify their domains. 
#+    There's no graphic for this, sorry!

# Hope it helps! Good luck with your research!

# ------------------------- WARNING! ---------------------------
# This program requires Python 3 to work--hence, the '.py' extension. If you don't have Python 3, 
#+ please crawl out of whichever rock you've been living under and install it.

# This program requires Biopython to work. If you don't have Biopython installed, 
#+what are you waiting for? A sign from God? This is it.

# This program requires BLAST to work. If you don't have BLAST installed locally, 
# please reconsider your life choices and then download it from the NCBI's website. 
# Since I'm nice, here's a link: https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/).

# This program requires MUSCLE to work. If you don't have MUSCLE installed, 
# you should check out how to do it here: https://www.drive5.com/muscle/downloads.htm  

# This program requires Matplotlib to work. If you don't have it installed, 
#+ go do that now, it'll come in handy.


# -------------------------- MODULES ---------------------------
import argparse
import os
import sys
import datetime
from pathlib import Path


# Format Verification Module (1/5)
import formatt as ftt

# Sequence Alignment Module (2/5)
import blast as al

# Tree-Making Module (3/5)
import muscle as tr

# Domain Finder Module (4/5)
import prosite as dom

# Graphs Module (5/5)
import graph as gph 


# -------------------- ARGUMENT CONTROL ------------------------

# Only execute is this script is being used as main_script
if __name__ == '__main__':
    pass
else:
    print('\033[1;97;101mERROR:\033[0m Cannot be executed as module.')
    sys.exit(2)

# Create argument parser
ap = argparse.ArgumentParser(description="Discover Antartica 2.0",
                            epilog="Hope that helped! Good luck with your research!")

ap.version = '\033[1;92mVersion 2.0\033[0m'

# Add the arguments to the parser
ap.add_argument("-q", "--query", 
                action="append",
                required=True,
                help="query proteins files or directory\
                (\033[1;97mSequences must be single fasta, if they are to be analysed independently.\
                 There will be one (1) phylogenetic analysis per query FILE.\033[0m)")

ap.add_argument("-s", "--subject", 
                action="append",
                required=True,
                help="subject proteins files or directory \
                (\033[1;97mSequences can be either multifastas or genbank assemblies.\033[0m)")

ap.add_argument("-c", "--coverage", 
                type=float, 
                default=float(30),
                action="store",
                help="coverage cut-off value")

ap.add_argument("-i", "--identity", 
                type=float,
                default=float(40),
                action="store",
                help="identity cut-off value")

ap.add_argument("-e", "--evalue", 
                type=float,
                default=float(0.0000001), 
                help="evalue threshold")

ap.add_argument("-g", "--graphics",
                action='store_true',
                help="adds graphics generated with matplotlib")

ap.add_argument("-d", "--documentation",
                action='store_true',
                help="adds prodoc documentation")

ap.add_argument("-v", "--version",
                action='version',
                help="prints out version")

# Execute the parse_args() method
args = ap.parse_args()

# --------------------- FILE MANAGEMENT ---------------------
# Create results directory
results_directory='results_'+str(datetime.datetime.now().strftime("%Y_%m_%d_%H-%M-%S"))

if not os.path.exists(results_directory):
        os.mkdir(results_directory)

# Manage queries
# Find absolute path
path1=os.path.abspath(str(args.query[0]))
query_directory=[]

# If the query given is a directory, get the files
if os.path.isdir(path1):
    for file in os.listdir(path1):
        query_directory.append(path1+'/'+file)
else:
    # If not, it's probably a file
    query_directory.append(path1)


# Manage subject databases
# Find absolute path
path2=os.path.abspath(str(args.subject[0]))
subject_directory=[]

# If the subject given is a directory, get the files
if os.path.isdir(path2):
    for file in os.listdir(path2):
        subject_directory.append(path2+'/'+file)
else:
    # If not, it's probably a file
    subject_directory.append(path2)

# ------------------- MAIN FUNCTION -----------------------

def i_am_iron_main():
    
    # Turn list into string
    global query
    query = str(query).strip('\'[]\'')
    print(f"\n------------------------------------\
            \nQUERY: \033[1;96m{Path(query).stem}\033[0m \
            \n------------------------------------ \n")

    # Make sure it's a fasta file
    if ftt.the_format_is_strong_in_this_one(str(query))=='FA':

        # 1. Set new project
        projectname=results_directory+'/'+Path(query).stem

        query_dictionary=ftt.the_perks_of_being_a_dictionary(
            query_file=query)

        # 2. Create a database:
        for subject in subject_directory:
            # Turn list into string
            subject = str(subject).strip('\'[]\'')

            # Check the format
            format_is=ftt.the_format_is_strong_in_this_one(subject)

            if format_is=='GB':
                ftt.in_the_begining_there_was_genbank(
                    gb_file=subject, 
                    fasta_file=projectname+'.fasta'
                    )
            elif format_is=='GP':
                ftt.in_the_begining_there_was_genprot(
                    gp_file=subject, 
                    fasta_file=projectname+'.fasta'
                    )

            elif format_is=='FA':
                ftt.in_the_begining_there_was_fasta(
                    fa_file=subject, 
                    fasta_file=projectname+'.fasta'
                    )

        # 3. BLAST against that database:
        # Dictionary has all the results fasta sequences 
        print("\n------------------------------------\
            \nBLAST analysis...\
            \n------------------------------------ \n")
        try: 
            f_dictionary=al.to_blast_or_not_to_blast(
                query_file=query, 
                db_file=projectname+'.fasta', 
                evalue_co=args.evalue, 
                raw_results_file=projectname+'_raw_blast.tsv',
                filtered_results_file=projectname+'_blast.tsv', 
                fasta_results_file=projectname+'_blast.fasta', 
                iden_co=args.identity, 
                cov_co=args.coverage,
                query_dictionary=query_dictionary
                )

            # Show grafic if asked: 
            if args.graphics:
                print('Please save and/or close graphic to continue with the analysis.')
                gph.the_perks_of_being_a_graphic(
                    file=projectname+'_blast.tsv'
                    )
        
        # Make sure an exception is not raised just because there are no hits.
        except SystemExit as e:
            if e.code != 'ZERO':
                print ("\033[1;97;101mERROR:\033[0m Something wrong while BLASTin'!")
                sys.exit(2)
            else:
                sys.exit(ZERO)

        # 4. Make a tree with MUSCLE
        print("\n------------------------------------\
            \nPhylogenetic analysis...\
            \n------------------------------------ ")
        try:
            tr.run_muscle_run(
                fasta_file=projectname+'_blast.fasta',
                aligned_file=projectname+'_muscle.fasta', 
                tree_file=projectname+'_tree.nw'
                )
            
            # Show grafic if asked: 
            if args.graphics:
                print('Please save and/or close graphic to continue with the analysis.')
                gph.you_cant_handle_the_tree(
                    tree_file=projectname+'_tree.nw'
                    )
        except:
            print ("\033[1;97;101mERROR:\033[0m Something wrong while making phylogenetic tree'!")
            sys.exit(2)
    
        # 5. Look for domains in PROSITE database
        print("\n------------------------------------\
            \nDomain analysis...\
            \n------------------------------------ \n")
        #try:
        dom.its_a_wonderful_database(
            dictionary=f_dictionary,
            text_file=projectname+'_domains.txt',
            DOC=args.documentation,
            prosite_dat_file='discover_antartica/prosite.dat',
            prosite_doc_file='discover_antartica/prosite.doc'
            )
        #except:
            print ("\033[1;97;101mERROR:\033[0m Something wrong while analyzing domains!")
            sys.exit(2)

    else:
        print("\033[1;97;101mERROR:\033[0m Can't run the program with a non-valid query file!")
        sys.exit(2)


# --------------------- MAIN SCRIPT -----------------------

# One 'project' per query
try: 
    for query in query_directory:
        i_am_iron_main()
        
except SystemExit as e:
    if e.code != 'ZERO':
        print('\033[1;97;101mERROR:\033[0m Something went terribly wrong.')
        sys.exit(2)
    else:
        os._exit(0)
    

