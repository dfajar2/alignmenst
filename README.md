Important Note:
I used this script in sample names that contain '_S' in the name
and used it to split first string to the left to get the sample name.  

===============

USAGE:

    $ /path/to/script.sh OPTIONS

    Required:
    [ -g Define genome. Available organisms: 'Human', 'Ecoli', 'St43300', 'Staph', 'Rhodo',
                                             'Rhodo241', 'Taq', 'Salmonella', 'Pseudomonas',
                                             'Enterobacter' and 'Serratia', 'Maize', '8gen',
					         '9gen', 'Seqwell_all_genomes' and 'Seqwell_all_genomes_renamed' ]
    [ -f Full path to reads directory ]
    [ -o Output Directory ]

    Optional:
    [ -s Use BWA instead of Bowtie2 (Default) ]
    [ -e string to help filter files]
    [ -c Expected depth to subsample reads. i.e. '10', '5', '1', etc. Default: 10 (for 10x) ]
    [ -a Use all reads. No downsampling ]
    [ -l Minimum length of reads after adapter trimming. Default: 90  ]
    [ -v Estimate full coverage stats ]
