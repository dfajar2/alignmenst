Important Note:
This is a modified script. This is for sets of samples that 
do not contain '_S' in the sample name. 
I used this script in sample names that contain '_' in the name
and used it to split first string to the left to get the sample name.  

===============

USAGE:

    $ /path/to/script.sh OPTIONS
        Required:
        [ -g Define genome. Available organisms: 'Human', 'Ecoli', 'St43300', 'Staph', 'Rhodo',
                                                 'Rhodo241', 'Taq', 'Salmonella', 'Pseudomonas',
                                                 'Enterobacter' and 'Serratia' ]
        [ -f Full path to reads directory ]
        [ -o Output Directory ]
        Optional:
        [ -s Use BWA instead of Bowtie2 (Default) ]
        [ -e string to help filter files]
        [ -c Expected depth to subsample reads. i.e. '10', '5', '1', etc. Default: 10 (for 10x) ]
        [ -a Use all reads. No downsampling ]
