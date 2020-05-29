Align and generate metrics. 

USAGE:

    $ /path/to/script.sh OPTIONS

        Required:
        [ -g Define genome. Available organisms: 'Human', 'Ecoli', 'St43300', 'Staph', 'Rhodo',
                                                 'Rhodo241', 'Taq', 'Salmonella', 'Pseudomonas',
                                                 'Enterobacter', 'Serratia', 'Maize', '8gen', 
                                                 '9gen', 'Seqwell_all_genomes' and 'Seqwell_all_genomes_renamed']
        [ -f Full path to reads directory ]
        [ -o Output Directory ]

        Optional:
        [ -s Use BWA instead of Bowtie2 (Default) ]
        [ -r If only want to use R1 from pair end reads. ]
        [ -e string to help filter files]
        [ -c Expected depth to subsample reads. i.e. '10', '5', '1', etc. Default: 10 (for 10x) ]
        [ -a Use all reads. No downsampling ]
        [ -l Minimum length of reads after adapter trimming. Default: 90  ]
        [ -v Estimate FULL coverage per base location  ]
