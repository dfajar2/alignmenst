#!/bin/bash
#set -eo pipefail
set -e
export PATH=/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:${PATH}
export PATH=/media/software/cutadapt/1.18/bin:${PATH}
export PATH=/media/software/bowtie2/2.3.4.3/bin:${PATH}
export PATH=/media/software/samtools/1.9/bin:${PATH}
export PATH=/media/software/ts/1.0/bin:${PATH}
export PATH=/media/software/freebayes/1.2.0-2-g29c4002/bin:${PATH}
export PATH=/media/software/bedtools/2.27.1/bin:${PATH}
export PATH=/media/software/vcftools/0.1.17/bin:/media/software/vcflib/bin:${PATH}
export PATH=/media/software/vcflib/bin:${PATH}
export PATH=/media/software/R/3.5.0/bin:${PATH}
export PATH=/media/software/picard/2.18.21/bin:${PATH}

# Export PERL5LIB for vcftools
export PERL5LIB=/media/software/vcftools/0.1.17/share/perl/5.26.1:${PERL5LIB}

# clear


echo Running pipeline; echo ==================
start=`date +%s`
echo $start

usage="
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

"
e="."
gen_size=""
c=10
s=0

while getopts "g:f:o:e:c:as" options; do
        case "${options}" in
                g)
                        g=${OPTARG} ;;
                f)
                        f=${OPTARG} ;;
                o)
                        o=${OPTARG} ;;
                s)
                        s=1 ;;
                e)
                        e=${OPTARG};;
                c)
                        c=${OPTARG};;
                a)
                        c=99999999 ;;
                *)
                        echo ${usage} ;;
        esac
done

shift $((OPTIND-1))
if [ -z "${g}" ] || [ -z "${f}" ] || [ -z "${o}" ] ; then
        echo ; echo "ERROR - Missing arguments"; echo "$usage"; exit 1

fi

if [ "$g" != "Human" ] && [ "$g" != "St43300" ] && [ "$g" != "Ecoli" ] && [ "$g" != "Staph" ] && [ "$g" != "Rhodo" ] && [ "$g" != "Rhodo241" ] &&[ "$g" != "Taq" ] && [ "$g" != "Salmonella" ] && [ "$g" != "Pseudomonas" ] && [ "$g" != "Enterobacter" ] && [ "$g" != "Serratia" ]
then
	echo ERROR - Invalid Genome \(-g\). ${g} genome not available.
	echo Valid options: 'Human', 'Ecoli', 'St43300', 'Staph', 'Rhodo', 'Rhodo241', 'Taq', 'Salmonella', 'Pseudomonas', 'Enterobacter' and 'Serratia'
	echo ; echo --------; echo "$usage"; exit 1; 
fi

genome=${g}
reads=${f}
output_dir=${o}
exp=${e}
st=1
en=100

# Define Reference Genome file
if [ ${genome} = "Human" ];
then
    reference="/media/ngs/ReferenceSequences/Homo_sapiens/hg38.fa";
    gen_size=3209286105 ;
    fastq_lines=$(echo 44000000*${c} | bc) ;
    st=22
    en=62
    echo ${genome} genome - $gen_size bases. $fastq_lines lines will be selected. ;
elif [ ${genome} = "Ecoli" ];
then
    reference="/media/ngs/ReferenceSequences/EcoliK12-MG1655-NC_000913.3.fasta" ;
    gen_size=4641652 ;
    fastq_lines=$(echo 61892*${c} | bc) ;
    st=35 ;
    en=65 ;
    echo ${genome} genome - $gen_size bases. $fastq_lines lines will be selected. ;
elif [ ${genome} = "Staph" ]; 
then
    reference="/media/ngs/ReferenceSequences/Staphylococcus_ATCC25923/Staph_ATCC25923.fasta" ;
    gen_size=2806345 ;
    fastq_lines=$(echo 42000*${c} | bc) ;
    st=21 ;
    en=44 ;
    echo ${genome} genome - $gen_size bases. $fastq_lines lines will be selected. ;
elif [ ${genome} = "St43300" ];
then
    reference="/media/ngs/ReferenceSequences/Staphylococcus_ATCC43300/GCF_003052445.1_ASM305244v1_genomic.fasta" ; 
    gen_size=2833835 ;
    fastq_lines=$(echo 42000*${c} | bc) ;
    st=21 ;
    en=44 ;
    echo ${genome} genome - $gen_size bases. $fastq_lines lines will be selected. ;
elif [ ${genome} = "Rhodo" ];
then
    reference="/media/ngs/ReferenceSequences/RhodobacterATCC17025/Rhodo_ATCC17025_All.fasta" ;
    gen_size=4557127 ;
    fastq_lines=$(echo 68000*${c} | bc) ;
    st=57
    en=81
    echo ${genome} genome - $gen_size bases. $fastq_lines lines will be selected. ;
#New_from_here
elif [ ${genome} = "Rhodo241" ];
then
    reference="/media/ngs/ReferenceSequences/Rhodobacter_sphaeroides_2.4.1/Rhodobacter_all.fasta" ;
    gen_size=4557127 ;
    fastq_lines=$(echo 68000*${c} | bc) ;
    st=57
    en=81
    echo ${genome} genome - $gen_size bases. $fastq_lines lines will be selected. ;
elif [ ${genome} = "Serratia" ];
then
    reference="/media/ngs/ReferenceSequences/Serratia_liquegaciens_GCF_000422085/Serratia_liquegaciens_GCF_000422085.1_ASM42208v1_genomic.fna" ;
    gen_size=5282719 ;
    fastq_lines=$(echo 70440*${c} | bc) ;
    st=35 ; #same as ecoli
    en=65 ;
    echo ${genome} genome - $gen_size bases. $fastq_lines lines will be selected. ;
elif [ ${genome} = "Enterobacter" ];
then
    reference="/media/ngs/ReferenceSequences/Enterobacter_cloacae_GCF_000025565/Enterobacter_cloacae_GCF_000025565.1_ASM2556v1_genomic.fna" ;
    gen_size=5598796 ;
    fastq_lines=$(echo 74652*${c} | bc) ;
    st=35 ; #same as ecoli
    en=65 ;
    echo ${genome} genome - $gen_size bases. $fastq_lines lines will be selected. ;
elif [ ${genome} = "Pseudomonas" ];
then
    reference="/media/ngs/ReferenceSequences/Pseudomonas_aeruginosa_GCF_000006765.1_ASM676v1_genomic.fna";
    gen_size=6264404 ;
    fastq_lines=$(echo 83528*${c} | bc) ;
    st=57 ; #same as rhodo
    en=81 ;
    echo ${genome} genome - $gen_size bases. $fastq_lines lines will be selected. ;
elif [ ${genome} = "Salmonella" ];
then
    reference="/media/ngs/ReferenceSequences/Salmonella_enterica_GCF_000195995/Salmonella_enterica_GCF_000195995.1_ASM19599v1_genomic.fna" ;
    gen_size=5133713 ;
    fastq_lines=$(echo 68448*${c} | bc) ;
    st=35 ; #same as ecoli
    en=65 ;
    echo ${genome} genome - $gen_size bases. $fastq_lines lines will be selected. ;
elif [ ${genome} = "Taq" ];
then
    reference="/media/ngs/ReferenceSequences/Thermus_aquaticus_GCF_001399775/Thermus_aquaticus_GCF_001399775.1_ASM139977v1_genomic.fna" ;
    gen_size=2338641 ;
    fastq_lines=$(echo 31180*${c} | bc) ;
    st=57 ; #same as rhodo
    en=81 ;
    echo ${genome} genome - $gen_size bases. $fastq_lines lines will be selected. ;
else
	echo "No valid reference."; exit 1
fi


# Check reference file name
if [[ "$reference" == *fasta  ]]
then
	echo Reference is in fasta format
elif [[ "$reference" == *fa  ]]
then
        echo Reference is in fasta format
elif [[ "$reference" == *fna  ]]
then
	echo Reference is in fasta format
else
        echo What the hell??? File Enterobacter_cloacae_GCF_000025565does not have .fasta, .fa or .fna suffix
        deactivate
	exit 1
fi


base_dir=$(pwd) 
ref_path=$(echo $reference | awk -F"/" '{OFS="/"; $(NF--)="";print}')
ref_complete_name=$(echo $reference | awk -F"/" '{print $NF}')
ref_basename=$(echo $reference | awk -F"/" '{print $NF}'| awk -F".fa" '{print $1}')


echo Subsampling Depth: ${c}x
num_reads=$((${fastq_lines}/4))
echo Number of read pairs: ${num_reads}
echo base_dir = ${base_dir}
echo output_dir = ${output_dir}
echo reads = ${reads}
echo genome = ${genome}
echo genome size = ${gen_size} bases
echo reference = ${reference}
head -n1 ${reference}
echo Reference_path = ${ref_path}
echo Reference_base_name = ${ref_basename}
echo expression = ${exp}
echo Lines to filter = ${fastq_lines}
echo ${genome} genome - $gen_size bases. $fastq_lines lines will be selected.

#Check paths
if [ ! -f ${reference} ] || [ ! -d ${reads} ] ;
then
	echo -----; echo ERROR !!! Reference File or Reads Directory do not exist. Check files and paths. ; exit 1
fi

# Get sample names from data folder 

if [ "${exp}" == '.' ] ;
then
	samples=$(for f in `ls --color=never -1 "${reads}" | grep _R1 | grep -v Unde | grep  -v _I | grep -v config.xml |  grep -v  FastqSummaryF1L1.txt | awk -F"_S" '{print $1}' | sort -u` ; do echo $f; done)  
else
	samples=$(for f in `ls --color=never -1 "${reads}" | grep ${exp} | grep _R1 | grep -v Unde | grep  -v _I | grep -v config.xml |  grep -v  FastqSummaryF1L1.txt | awk -F"_S" '{print $1}' | sort -u` ; do echo $f; done)  
fi

numsamp=$(printf '%s\n' $samples:q | wc -w)
echo Number of Samples: ${numsamp}
echo Names of the first 10 samples: `for f in {1..10} ; do echo ${samples} |  awk -v v=${f} '{print $v}' ; done`

if [ $s -eq 0 ]
then
        echo Aligner to use: Bowtie2
elif [ $s -eq 1  ]
then
        echo Aligner to use: BWA
else
        echo No aligner defined ???
fi

#Check if SE or PE

if [ $s -eq 0 ]
then
        echo Using Bowtie2 for alignment
        if ls -1 ${reads} | grep -q _R2 ;
        then
                echo ---------------------
                echo "This study contains PE reads"
                samp1=$(echo ${samples} | awk '{print $1}')
                echo Reads for Sample ${samp1}:
                ls ${reads}${samp1}_*_R1*
                ls ${reads}${samp1}_*_R2*
                echo ---------------------
                #ls ${reads}${samp1}_*_R1*
                popo1=$(ls ${reads}${samp1}_*_R1*)
                popo2=$(ls ${reads}${samp1}_*_R2*)
                echo Bowtie2 command:
                echo bowtie2 --threads 16 --rg-id ${samp1} --rg SM:${samp1} --rg LB:${samp1} --rg CN:LGC_Genomics --rg PL:Illumina -x ${base_dir}/../reference/${ref_basename} -1 `ls ${popo1} | tr "\n" "," | sed 's/,$//g'` -2 `ls ${popo2} | tr "\n" "," | sed 's/,$//g'`  2\>${samp1}.log \| samtools view -Sbu - \| samtools sort -@16 -m 4G \> ${samp1}_sorted.bam
        else
                echo ---------------------
                echo "This study contains SE reads"
                samp1=$(echo ${samples} | awk '{print $1}')
                echo Read for Sample ${samp1}:
                ls ${reads}${samp1}_*_R1*
                popo=$(ls ${reads}${samp1}_*_R1* | grep ${exp})
                echo Bowtie2 command:
                echo bowtie2 --threads 16 --rg-id ${samp1} --rg SM:${samp1} --rg LB:${samp1} --rg CN:LGC_Genomics --rg PL:Illumina -x ${base_dir}/../reference/${ref_basename} -U `ls ${popo} | tr "\n" "," | sed 's/,$//g'` 2\>${samp1}.log \| samtools view -Sbu - \| samtools sort -@16 -m 4G \> ${samp1}_sorted.bam
        fi
elif [ $s -eq 1  ]
then
        echo Using BWA for alignment
        if ls -1 ${reads} | grep -q _R2 ;
        then
                echo ---------------------
                echo "This study contains PE reads"
                samp1=$(echo ${samples} | awk '{print $1}')
                echo Reads for Sample ${samp1}:
                ls ${reads}${samp1}_*_R1*
                ls ${reads}${samp1}_*_R2*
                echo ---------------------
                #ls ${reads}${samp1}_*_R1*
                popo1=$(ls ${reads}${samp1}_*_R1*)
                popo2=$(ls ${reads}${samp1}_*_R2*)
                echo BWA command:
                echo bwa mem -M -t 16 ${reference} -1 `ls ${popo1} | tr "\n" "," | sed 's/,$//g'` -2 `ls ${popo2} | tr "\n" "," | sed 's/,$//g'` 2\>${samp1}.log \| samtools view -Sbu - \| samtools sort -@16 -m 4G \> ${samp1}_sorted.bam
        else
                echo ---------------------
                echo "This study contains SE reads"
                samp1=$(echo ${samples} | awk '{print $1}')
                echo Read for Sample ${samp1}:
                ls ${reads}${samp1}_*_R1*
                popo=$(ls ${reads}${samp1}_*_R1* | grep ${exp})
                echo BWA command:
                echo bwa mem -M -t 16 ${reference} -U `ls ${popo1} | tr "\n" "," | sed 's/,$//g'` 2\>${samp1}.log \| samtools view -Sbu - \| samtools sort -@16 -m 4G \> ${samp1}_sorted.bam
        fi
else
        echo No aligner defined????
fi

# Check for ${output_dir} folder - to continue
if [ -d ${output_dir} ]
then
echo ; echo ---------------; echo WARNING; echo; echo ${output_dir} folder already exist!!
echo Do you want to delete the \"${output_dir}\" folder and run the analysis again? \<Type Y to continue\>
read input2
if [ "$input2" != "Y" ] && [ "$input2" != "y" ];
then
	echo Exiting...
	exit 1
else
	# echo Uncomment line below to remove Directory
	rm -r ${output_dir}
fi
fi

#Safe stop
echo ---------------
echo Check files and paths
echo Press Y to continue, any other key to exit.
echo ---------------
read input
if [ "$input" != "Y" ] && [ "$input" != "y" ]; 
then
	echo Exiting...  ; 
	echo
	exit 0
fi


#############
### START ###
#############

#Project folder
mkdir ${output_dir} ; cd ${output_dir}
mkdir -p deliverables/
base_dir=$(pwd)
echo ${base_dir}

#Data folder
mkdir data; cd data
for file in $(ls --color=never -1 ${reads}); do ln -s ${reads}$file ; done

##Fastqc
#echo ---------------
#echo Do you want to run Fastqc?
#echo Press Y to run it, any other key to skip it.
#echo ---------------
#read input
#echo ${reads}
#cd ${reads}
#if [ "$input" != "Y" ] && [ "$input" != "y" ];
#then
#        echo Skipping Fastqc analysis...  ;
#        echo
#else
#       mkdir -p fastqc
#       if ls -1 --color=never ${reads} | grep -q _R2 ;
#       then
#               for f in `ls ${samp1}_*_R1*| awk -F"_R" '{print $1}'` ; do echo $f ; fastqc -o ${reads}/fastqc/${f}_fastqc -t 16 ${reads}/${f}_R1_001.fastq.gz ../${f}_R2_001.fastq.gz; done
#       else
#               for f in `ls ${samp1}_*_R1*| awk -F"_R" '{print $1}'` ; do echo $f ; fastqc -o ${reads}/fastqc/${f}_fastqc -t 16 ${reads}/${f}_R1_001.fastq.gz ; done
#       fi
#fi
#pwd
#exit 1
#cd ..

#Reference folder
if [ -d ${base_dir}../reference ]
then
if [ -f $ref_name ] && [ -f ${ref_basename}.1.bt2 ]
then
	echo Reference and Bowtie index already exist.
elif [ -f $ref_name ]
then
	mkdir ../reference; cd ../reference
	ln -s ${reference}* .
	cd ${base_dir}
else
	echo Fasta file does not exist to build index ; exit 1
fi
fi

cd ${base_dir}

echo ===== FUNCTIONS ==============

function phix(){
cd ${base_dir}
ref_path="/media/data/diego/work/reference/" 
ref_basename="PhiXgenome" 
mkdir bowtie2_phix
cd bowtie2_phix/ 
for f in ${samples} ; do echo bowtie2 --threads 48 --rg-id ${f} --rg SM:${f} --rg LB:${f} --rg CN:LGC_Genomics --rg PL:Illumina --dovetail --minins 0 --maxins 1000 -x ${ref_path}${ref_basename} -1 ${reads}${f}*_R1_*gz -2 ${reads}${f}*_R2_*gz  2\>${f}.log \| samtools view -Sbu - \| samtools sort -@16 -m 4G \> ${f}_sorted.bam; done > Bowtie_bash.sh ; #echo "ls *.bam > orig_bam.lst ; samtools merge -b orig_bam.lst --threads 60 merged_orig_bam_sorted.bam ; samtools index merged_orig_bam_sorted.bam" \>\> Bowtie_bash.sh 
sh Bowtie_bash.sh 
mkdir ${base_dir}/unmapped_to_phix_reads 
for f in `ls *bam | awk -F. '{print $1}'` ; do echo Extracting unmapped reads: Sample ${f}; samtools view -b -f 12 ${f}.bam > ${f}_unmapped.bam ; bamToFastq -i ${f}_unmapped.bam -fq ${base_dir}/unmapped_to_phix_reads/${f}_unmapped_R1.fq -fq2 ${base_dir}/unmapped_to_phix_reads/${f}_unmapped_R2.fq ; done 
cd ${base_dir}
}
echo phix
phix


function downsampling() {
echo Base_dir: ${base_dir}
cd ${base_dir}/data
mkdir downsampled ; cd downsampled
echo Base dir : ${base_dir}
echo
echo Files in ${base_dir}/data
ls ${base_dir}/data
echo 
echo Files in ${base_dir}/unmapped_to_phix_reads/
ls ${base_dir}/unmapped_to_phix_reads/
echo
echo Files in ${base_dir}
ls ${base_dir}/
echo
echo Samples: ${samples}
for f in ${samples} ; \
do echo Downsampling sample: ${f}; \
paste ${base_dir}/unmapped_to_phix_reads/${f}_sorted_unmapped_R1.fq ${base_dir}/unmapped_to_phix_reads/${f}_sorted_unmapped_R2.fq | \
awk '{ printf("%s",$0); n++; if(n%4==0) { printf("\n");} else { printf("\t\t");} }' | \
shuf | \
head -n ${num_reads} | awk '{if (/^@/) gsub(/ /,"_",$0) ;  print $0 }' | \
sed 's/\t\t/\n/g' | \
awk -v samp=${f} '{print $1 > samp"_1.fastq"; print $2 > samp"_2.fastq"}' ; \
cat ${f}_1.fastq | \
awk '{if (/^@/) gsub(/_/," ",$0) ;  print $0 }' > ${f}_R1.fastq; \
cat ${f}_2.fastq | \
awk '{if (/^@/) gsub(/_/," ",$0) ;  print $0 }' > ${f}_R2.fastq ;  \
rm ${f}_1.fastq; \
rm ${f}_2.fastq;  \
done 
echo Quality trimming
mkdir qual_trim
cd qual_trim
for f in ${samples}; do echo Trimming $f; java -jar /media/software/Trimmomatic-0.36/bin/trimmomatic-0.36.jar PE -phred33 -threads 48 ../${f}_R1.fastq ../${f}_R2.fastq ${f}_R1_trimmed.fastq ${f}_R1_Unpaired.seqs ${f}_R2_trimmed.fastq ${f}_R2_Unpaired.seqs SLIDINGWINDOW:5:20 LEADING:5 TRAILING:5 ; rm ${f}_R*_Unpaired.seqs ;done       
cd ${base_dir}
}

# Comment downsampling
echo downsampling
echo 
echo genome ${genome}
echo Reference $reference
echo genome Size $gen_size
echo num of lines to extract $fastq_lines
echo num of reads $num_reads
downsampling
#exit 1
#echo remember to uncomment after debugging

echo
ref_path=$(echo $reference | awk -F"/" '{OFS="/"; $(NF--)="";print}')
ref_complete_name=$(echo $reference | awk -F"/" '{print $NF}')
ref_basename=$(echo $reference | awk -F"/" '{print $NF}'| awk -F".fa" '{print $1}')
echo Reference : ${reference}
ref_name=$(echo $reference | awk -F"/" '{print $NF}')
echo Reference name: ${ref_name}
echo Reference base_name: ${ref_basename}


function CheckReference() {
cd ${ref_path}
if [[ "$ref_name" == *fasta  ]]
then
	echo Reference is in fasta format
elif [[ "$ref_name" == *fa ]]
then
	echo Reference is in fasta format
elif [[ "$ref_name" == *fna ]]
then
	echo Reference is in fasta format
else
	echo What the hell??? File does not have .fasta, .fa nor .fna suffix
	exit 1
fi

#Check Reference Index
if [ $s -eq 0 ]
then
        # Check Reference and Bowtie2 index
        if [ -f $ref_name ] && [ -f ${ref_basename}.1.bt2 ]
        then
                echo Reference and Bowtie index already exist.
        elif [ -f $ref_name ]
        then
                echo Reference exist. Building Bowtie2 reference index.
                bowtie2-build --threads 16  ${ref_complete_name} ${ref_basename}
                cd -
        else
                echo Fasta file does not exist to build index ; exit 1
        fi
else
        # Check Reference and BWA index
        if [ -f $ref_name ] && [ -f ${ref_name}.bwt ]
        then
                echo Reference and BWA index already exist.
        elif [ -f $ref_name ]
        then
                echo Reference exist. Building BWA reference index.
                bwa index ${ref_name}
                cd -
        else
                echo Fasta file does not exist to build index ; exit 1
        fi
fi
cd ${base_dir}
}


# Uncomment to build reference index again
echo CheckReference
CheckReference


function Bowtie2Alignment() {
cd ${base_dir}
mkdir alignments; cd alignments
# Bowtie2 alignment
echo Files list
lst=`ls --color=never -1 ${data_dir}`
echo Checking if SE or PE
data_dir=${base_dir}/data/downsampled/qual_trim/
if ls -1 $data_dir | grep -q _R2 ;
then
	echo "Paired data"
	for f in $samples ; do echo bowtie2 --threads 48 --rg-id ${f} --rg SM:${f} --rg LB:${f} --rg CN:LGC_Genomics --rg PL:Illumina --minins 0 -x ${ref_path}/${ref_basename} -1 `ls ${data_dir}${f}_*R1*fastq* |  tr "\n" "," | sed 's/,$//g'` -2 `ls ${data_dir}${f}_*R2*fastq* |  tr "\n" "," | sed 's/,$//g'` 2\>${f}.log \| samtools view -Sbu - \| samtools sort -@16 -m 4G \> ${f}_sorted.bam; done > Bowtie_bash.sh ; echo "ls *.bam > orig_bam.lst ; samtools merge -b orig_bam.lst --threads 60 merged_orig_bam_sorted.bam ; samtools index merged_orig_bam_sorted.bam" >> Bowtie_bash.sh
# echo Debug; echo STOP!!!!; cat Bowtie_bash.sh; exit 1
else
	echo "Unpaired data"
	for f in $samples ; do echo bowtie2 --threads 48 --rg-id ${f} --rg SM:${f} --rg LB:${f} --rg CN:LGC_Genomics --rg PL:Illumina --minins 0 -x ${ref_path}/${ref_basename} -U `ls ${data_dir}${f}_*R1*fastq* |  tr "\n" "," | sed 's/,$//g'` 2\>${f}.log \| samtools view -Sbu - \| samtools sort -@16 -m 4G \> ${f}_sorted.bam ; done > Bowtie_bash.sh ; echo "ls *.bam > orig_bam.lst ; samtools merge -b orig_bam.lst --threads 60 merged_orig_bam_sorted.bam ; samtools index merged_orig_bam_sorted.bam" >> Bowtie_bash.sh

#echo Debug; echo STOP!!!!; cat Bowtie_bash.sh; exit 1
fi

echo Running Bowtie scripts
# Check Reference and Bowtie2 index                              
if [ -f Bowtie_bash.sh ]                
then
	echo Aligning samples with bowtie2.
	sh Bowtie_bash.sh
else
	echo No bowtie file \(Bowtie_bash.sh\) generated.; exit 1
fi
#Exit bowtie2
cd ${base_dir} 
}


function BWAAlignment() {
cd ${base_dir}
mkdir alignments; cd alignments
# Bowtie2 alignment
echo Files list
lst=`ls --color=never -1 ${data_dir}`
echo Checking if SE or PE
data_dir=${base_dir}/data/downsampled/qual_trim/
if ls -1 $data_dir | grep -q _R2 ;
then
	echo "Paired data"
	for f in $samples ; \
	do echo bwa mem -M -t 16 ${ref_path}/${ref_name} \
	-1 `ls ${data_dir}${f}_*R1*fastq* | tr "\n" "," | sed 's/,$//g'` \
	-2 `ls ${data_dir}${f}_*R2*fastq* | tr "\n" "," | sed 's/,$//g'` 2\>${f}.log \| samtools view -Sbu - \| samtools sort -@16 -m 4G \> ${f}_sorted.bam ; \
	done  > BWA_bash.sh
# echo Debug; echo STOP!!!!; cat Bowtie_bash.sh; exit 1
else
	echo "Unpaired data"
	for f in $samples ; \
	do echo bwa mem -M -t 16 ${ref_path}/${ref_name} \
	-U `ls ${data_dir}${f}_*R1*fastq* |  tr "\n" "," | sed 's/,$//g'` 2\>${f}.log \| samtools view -Sbu - \| samtools sort -@16 -m 4G \> ${f}_sorted.bam ; \
	done  > BWA_bash.sh;

#echo Debug; echo STOP!!!!; cat Bowtie_bash.sh; exit 1
fi

echo Running Bowtie scripts
                              
if [ -f BWA_bash.sh ]                
then
	echo Aligning samples with BWA
	sh BWA_bash.sh
else
	echo No BWA file \(BWA_bash.sh\) generated.; exit 1
fi
#Exit bwa
cd ${base_dir} 
}

echo Alignment
if [ $s -eq 0 ]
then
        echo Aligner to use: Bowtie2
	echo Aligning samples
	Bowtie2Alignment
elif [ $s -eq 1  ]
then
        echo Aligner to use: BWA
	echo Aligning samples
	BWAAlignment
else
        echo No aligner defined ???
fi


 
#echo Comment exit to fix and continue after debugging
#deactivate; exit 1


 
function Markdup() {
cd ${base_dir}
mkdir markdup; cd markdup
cd ${base_dir}/alignments ; for f in `ls *bam| awk -F. '{print $1}'`; do echo ; echo ${f};  java -jar /media/software/picard/2.18.21/bin/picard.jar MarkDuplicates INPUT=${f}.bam  OUTPUT=../markdup/${f}.markdup.bam METRICS_FILE=../markdup/${f}_metrics.txt REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=LENIENT ; done
cd ${base_dir}/markdup ; for d in `ls *markdup.bam`; do echo ${d};  samtools index ${d}; done
cd ${base_dir}
} 

echo Markdup
Markdup

function Stats() {
mkdir stats
cd ${base_dir}/markdup
for f in `ls *bam`; do echo ; echo Stats for sample ${f}; samtools stats ${f}> ../stats/${f}.stats ; done
cd ${base_dir}/stats; for f in `ls *stats`; do echo ${f} ; grep ^SN ${f} | grep -e "sequences:" -e "reads mapped:" -e "bases mapped (cigar):" -e "insert size average:" -e "insert size standard deviation:" -e "reads duplicated:"; echo lolo; done | tr "\n" "\t" | sed 's/lolo\t/lolo\n/g'| sed 's/\tlolo//g' | awk -F"\t" -v gensize=$gen_size 'BEGIN{printf "Sample\tGen_size\tCoverage \tTotalReads\tReads_mapped\tPct\tBases_mapped_(cigar)\tReads_duplicated\tInsert_size_average\tInsert_size_standard_deviation\n"}{printf"%s\t%d\t%.2fx\t%d\t%d\t%.2f\t%d\t%d\t%d\t%d\n",$1,gensize,((($4)*150)/gensize),$4,$13,($13/$4)*100,$20,$16,$24,$27}' > stats_table.txt


#awk -F"\t" 'BEGIN{printf "Sample\tReads_mapped\tBases_mapped_(cigar)\tInsert_size_average\tInsert_size_standard_deviation\tReads_duplicated\n"}{printf "%s\t%d\t%d\t%d\t%.2f\t%d\n",$1,$4,$11,$15,$NF,$7}' > stats_table.txt

cd ${base_dir}/markdup
mkdir ../insert_size
for f in `ls *.bam | awk -F"." '{print $1}'` ; do echo Collect Insert Size Metric for sample ${f} ; java -jar /media/software/picard/2.18.21/bin/picard.jar CollectInsertSizeMetrics I=${f}.markdup.bam O=../insert_size/${f}_size_metrics.txt H=../insert_size/${f}.pdf ; done


mkdir ../complexity 
for f in `ls *.bam | awk -F"." '{print $1}'` ; do echo Library Complexity for sample ${f} ; java -jar /media/software/picard/2.18.21/bin/picard.jar EstimateLibraryComplexity I=${f}.markdup.bam O=../complexity/${f}_complexity.txt ; done

cd ${base_dir}
}

Stats
echo stats
echo Debug ; cat ${base_dir}/stats/stats_table.txt

function GCbias(){
cd ${base_dir}
mkdir gc_bias

cd ${base_dir}/markdup ; for f in `ls *markdup.bam| awk -F. '{print $1}'` ; do echo; java -jar /media/software/picard/2.18.21/bin/picard.jar CollectGcBiasMetrics I=${f}.markdup.bam O=../gc_bias/${f}_gc_bias_metrics.txt CHART=../gc_bias/${f}_gc_bias_metrics.pdf S=../gc_bias/${f}_summary_metrics.txt R=${ref_path}/${ref_name} ; echo; echo Creating table for GC plot - Sample ${f}; cat ../gc_bias/${f}_gc_bias_metrics.txt | awk '{if (/All/) print $4,$8}' > ../gc_bias/${f}_table_for_plot.txt ; done

cd ${base_dir}/gc_bias

merg_samp=$(ls *_table_for_plot.txt | awk -F"-" '{print $1}'  | sort | uniq -c | awk '{if ($1>1) print $2}')
echo $merg_samp

#for f in $merg_samp ; do join ${f}*_table_for_plot.txt > ${f}_joined.txt; done
####
####
#### commented because it is breaking the script
####
####
####if [ `ls ${merg_samp}*_table_for_plot.txt | wc -l` -eq 1 ]
####then
####    echo No files to merge
####    for f in $merg_samp ; do cp ${merg_samp}*_table_for_plot.txt > ${f}_joined.txt; done  
####elif [ `ls ${merg_samp}*_table_for_plot.txt | wc -l` -eq 2 ]
####then
####	for f in $merg_samp ; do join ${f}*_table_for_plot.txt > ${f}_joined.txt; done
####elif [ `ls ${merg_samp}*_table_for_plot.txt | wc -l` -eq 3 ]
####then
####	for f in $merg_samp ; do join <(join ${f}-{1..2}*_table_for_plot.txt) ${f}-3*_table_for_plot.txt > ${f}_joined.txt ; done
####elif [ `ls ${merg_samp}*_table_for_plot.txt | wc -l` -eq 4 ]
####then
####	for f in $merg_samp ; do join <(join ${f}-{1..2}*_table_for_plot.txt) <(join ${f}-{3..4}*_table_for_plot.txt) > ${f}_joined.txt ; done
####fi
####
####

#for f in `ls *_joined.txt | awk -F_ '{print $1}'`; do echo GC stats for ${f}; cat ${f}*_joined.txt | awk -v princ=$st -v fin=$en '{if (NR>princ && NR<=fin+1) print $1,($2+$3+$4)/(NF-1)}' | awk 'function abs(v) {return v < 0 ? -v : v}{x[NR]=$1; y[NR]=$2; sx+=x[NR]; sy+=y[NR]; sxx+=x[NR]*x[NR]; sxy+=x[NR]*y[NR]; dev1+=1-y[NR]; meanabs+=abs(1-y[NR])}END{det=NR*sxx - sx*sx; a=(NR*sxy - sx*sy)/det;print "slope:"a, "Dev_from_1:"dev1, "AbsMean:"meanabs/NR}' > ${f}_GC_stats.txt ; done

####
####for f in `ls *_joined.txt | awk -F_ '{print $1}'`; do echo GC stats for ${f}; cat ${f}*_joined.txt | awk -v princ=$st -v fin=$en '{if (NR>princ && NR<=fin+1) print $1,($2+$3+$4)/(NF-1)}' | awk 'function abs(v) {return v < 0 ? -v : v} BEGIN { FS = "[ ,\t]+" } NF == 2 { x_sum += $1 ; y_sum += $2 ; xy_sum += $1*$2 ; x2_sum += $1*$1 ; num += 1 ; x[NR] = $1 ; y[NR] = $2 ; dev1+=1-y[NR]; meanabs+=abs(1-y[NR])} END { mean_x = x_sum / num ; mean_y = y_sum / num ; mean_xy = xy_sum / num ; mean_x2 = x2_sum / num ; slope = (mean_xy - (mean_x*mean_y)) / (mean_x2 - (mean_x*mean_x)) ; inter = mean_y - slope * mean_x ; for (i = num; i > 0; i--) {ss_total += (y[i] - mean_y)**2 ; ss_residual += (y[i] - (slope * x[i] + inter))**2} ; r2 = 1 - (ss_residual / ss_total) ; printf("Slope      :  %g\n", slope) ; printf("Intercept  :  %g\n", inter) ; printf("R-Squared  :  %g\n", r2) ; printf("Dev from 1 :  %g\n", dev1) ; printf("MeanAbs Val:  %g\n", meanabs/NR) }' > ${f}_GC_stats.txt; done
####

}

echo GCbias 
GCbias
# echo debug
#cat *_GC_stats.txt
# exit 1;


function Gen0Cov() {
cd ${base_dir}
mkdir genome_coverage
echo Estimating genome regions 0 coverage
cd ${base_dir}/markdup
for f in `ls *bam | awk -F. '{print $1"."$2}'`; do echo ${f};  genomeCoverageBed -ibam ${f}.bam -bga > ../genome_coverage/${f}_coverage.bed ; echo ;done
cd ${base_dir}/genome_coverage
for f in `ls *bed`; do echo $f; grep -w 0$ ${f} | awk '{sum += $3-$2 }END { print NR" regions with no coverage", sum" total bases"}' ; done | tr "\n" "\t" | sed 's/bases\t/bases\n/g' | awk -F"[\t ]" 'BEGIN{printf "Sample\tRegions_with_no_coverage\tTotal_bases\n"}{printf "%s\t%d\t%d\n",$1,$2,$7}' > genome_0_coverage_table.txt
cd ${base_dir}
mkdir coverage
cd ${base_dir}/markdup
echo Creating bedfile
for f in `ls *markdup.bam | awk -F. '{print $1"."$2}'`; do echo $f ; genomeCoverageBed -bg -trackline -ibam ${f}.bam > ../coverage/${f}_bedgraph.bed ;done
}

echo Gen0Cov
Gen0Cov
echo Debug ; cat ${base_dir}/genome_coverage/genome_0_coverage_table.txt
echo Debug ; head ${base_dir}/coverage/*_bedgraph.bed

cd  ${base_dir}
mkdir ${base_dir}/deliverables
cp ${base_dir}/stats/stats_table.txt ${base_dir}/deliverables
cp ${base_dir}/gc_bias/*_table_for_plot.txt ${base_dir}/deliverables
cp ${base_dir}/genome_coverage/genome_0_coverage_table.txt  ${base_dir}/deliverables
#cp ${base_dir}/gc_bias/*joined.txt ${base_dir}/deliverables
cp ${base_dir}/insert_size/*size_metrics.txt ${base_dir}/deliverables
cp ${base_dir}/insert_size/*pdf ${base_dir}/deliverables
