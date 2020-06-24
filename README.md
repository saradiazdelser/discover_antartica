# Discover Antartica 2.0

This wonderful, looong program has the ability do lots of interesting things, including giving the coder (me!) grey hair.
Here's what it'll do:

  1. Run BLAST on a number of chosen query proteins, given a certain evalue, coverage and identity cut-off values. This will be done locally, so you'll need a database. It will show you the results using a pretty graphic.
  
  2. Use MUSCLE to multiple-align the proteins aligned by BLAST. Then make a phylogenetic tree using the NJ method. It also show you tree using graphic.
  
  3. Run all your proteins against a prosite database in order to identify their domains. There's no graphic for this, sorry!

Hope it helps! Good luck with your research!

# USAGE
In order to use this script, the user must download the following databases: prosite.dat and prosite.doc, and save them in the discover_antartica package.

Databases can be downloaded here: https://mega.nz/file/5iJygYQC#EuUTz_SInscP0aga9-u853sUQhV3usdE4rH0kgbkrbk

    discover_antartica/main.py [-h] -q QUERY_FILE -s SUBJECT_FILE [-c COVERAGE] [-i IDENTITY] [-e EVALUE] [-g] [-d] [-v]
  
    -h, --help              show this help message and exit
  
    -q --query QUERY_FILE   query proteins file or directory (sequences must be single fasta)
    
    -s --subject SUBJ_FILE  subject proteins file or directory (sequences can be multifasta or genbank)
    
    -c --coverage COVERAGE  coverage cut-off value
    
    -i --identity IDENTITY  identity cut-off value
    
    -e --evalue EVALUE      evalue threshold
    
    -g, --graphics          add graphics generated with matplotlib
    -d, --documentation     adds prodoc documentation
    -v, --version           prints out version

The program will run one analysis per query file. If there are multiple query files, please input them as a single directory Thus, if one of the query files is a multifasta, all the sequences will be analyzed together. 

#  REQUIREMENTS 
This program requires *Python 3* to work——hence, the '.py' extension. If you don't have Python 3, please crawl out of whichever rock you've been living under and install it.

This program requires *Biopython* to work. If you don't have Biopython installed, what are you waiting for? A sign from God? This is it.

This program requires *BLAST* to work. If you don't have BLAST installed locally, please reconsider your life choices and then download it from the NCBI's website. Since I'm nice, here's a link: https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/).

This program requires *MUSCLE* to work. If you don't have MUSCLE installed, you should check out how to do it here: https://www.drive5.com/muscle/downloads.htm  

This program requires *Matplotlib* to work. If it's not installed, go ahead and do it now, it'll come in handy!

