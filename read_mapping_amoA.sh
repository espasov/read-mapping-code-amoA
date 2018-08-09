#!/bin/bash
# code for performing read-mapping of metagenomic reads to specific amoA sequences to determine relative abundance of amoA genes for Nitrospira
# The code is modified from Jackson Tsuji's ATLAS co-assembly code available at https://github.com/jmtsuji/atlas-extensions/blob/master/co-assembly.sh

#read_mapping_amoA.sh	
#Written by Emilie Spasov, Neufeld Lab, Univeristy of Waterloo, May 2018


set -e
set -u
set -o pipefail 
#see chapter 12 pg 397 in Vince Buffalo's Bioinformatics Data skills for bash scripting background



if [ "$#" -lt 5 ] #check if less than 5 arguments
then	
	echo "error: too few arguments, you provided $#, 5 required"
	echo "useage: read_mapping_amoA.sh altas_directory guide_filepath  mapping_directory mapping_references_seq threads"
	exit 1
fi

echo "script name: $0"
echo "atlas directory: $1"
echo "guide filepath: $2"
echo "mapping directory: $3"
echo "mapping reference sequences: $4"
echo "threads: $5"




# Atlas_output_dir: Path to the base directory where ATLAS files were output. 
# guide_filepath: TSV file with two columns: mapping_names (names of coassembly runs); read_mapping_samples (comma-separated names of previously run ATLAS samples to read map for binning)
# mapping_dir: path to directory where want mapping files to go
# mapping_ref: fasta file with the reference sequences to map to
# threads: number of threads to run


# Set variables from user input:
ATLAS_DIR=$1
GUIDE_FILEPATH=$2
MAPPING_DIR=$3
MAPPING_REF=$4
THREADS=$5


# Get GUIDE_FILEPATH params
	mapping_names=($(tail -n +2 ${GUIDE_FILEPATH} | cut -d $'\t' -f 1))
	read_mapping_samples=($(tail -n +2 ${GUIDE_FILEPATH} | cut -d $'\t' -f 2))
	
	
#Emilie's first run inputs:
#ATLAS_DIR=/Hippodrome/espasov/WWTP_metagenomics2017/180129_atlas_espasov_r1.0.22_bash_WWTP_full2/
#GUIDE_FILEPATH=/Hippodrome/espasov/WWTP_metagenomics2017/downstream_analysis_atlas/read_mapping_amoA/RBC_sample_names.tsv
#MAPPING_DIR=/Hippodrome/espasov/WWTP_metagenomics2017/downstream_analysis_atlas/read_mapping_amoA
#MAPPING_REF=/Hippodrome/espasov/WWTP_metagenomics2017/downstream_analysis_atlas/read_mapping_amoA/comammox_amoA_DNA_cluster99_clusternames.fasta
#THREADS=12

function read_map_to_coassemblies {
	# Description: iteratively maps read_mapping_samples to curated amoA sequences (RBC groups) like done within ATLAS
	# GLOBAL Params: OUTPUT_DIR; THREADS; mapping_names (array); read_mapping_samples (array)
	# Return: writes files/directories to disk

	echo "Read mapping to supplied sequences..."


	# Manually add additional settings needed for scripts.
	# TODO - pull these settings (at least MEMORY) from the .yaml file!
	local MEMORY=65 #TODO
	
	
	for i in $(seq 1 ${#mapping_names[@]}); do
		# Set counter to be based on zero, not one
		local j=$((${i}-1))
		
		# Make relevant directories for storing output
		mkdir -p ${MAPPING_DIR}/contig_stats
		mkdir -p ${MAPPING_DIR}/logs
	
		# Temporarily change the internal fields separator (IFS) to parse comma separators. See Vince Buffalo's "Bioinformatics Data Skills" (1st Ed.) chapter 12, pg 407 and corresponding Github page README at https://github.com/vsbuffalo/bds-files/tree/master/chapter-12-pipelines (accessed Nov 19, 2017)
		local OFS="$IFS"
		IFS=,

		# Get names of individual samples provided for that coassembly name
		local mapping_sample_IDs=(${read_mapping_samples[${j}]}) #read_mapping_samples is a comma separated file with the names of the samples

		# Fix the IFS
		IFS="$OFS"
		
		echo "${mapping_sample_IDs[@]}"
		
		# Read map iteratively for each mapping ID
		for mapping in ${mapping_sample_IDs[@]}; do
			# TODO - pull more settings from .yaml file (these are FIXED right now)

			echo "rule align_reads_to_amoA_seq (${mapping}):"
			echo "bbwrap.sh nodisk=t ref=${MAPPING_REF} in1=${ATLAS_DIR}/${mapping}/sequence_quality_control/${mapping}_QC_R1.fastq.gz,${ATLAS_DIR}/${mapping}/sequence_quality_control/${mapping}_QC_se.fastq.gz in2=${ATLAS_DIR}/${mapping}/sequence_quality_control/${mapping}_QC_R2.fastq.gz,null  outm=${MAPPING_DIR}/${mapping}.sam threads=${THREADS} trimreaddescriptions=t pairlen=1000 pairedonly=f mdtag=t xstag=fs nmtag=t sam=1.3 local=t ambiguous=best secondary=t ssao=t maxsites=10 -Xmx${MEMORY}G 2> ${MAPPING_DIR}/logs/contig_coverage_stats_${mapping}.log"
			echo ""
			bbwrap.sh nodisk=t ref=${MAPPING_REF} in1=${ATLAS_DIR}/${mapping}/sequence_quality_control/${mapping}_QC_R1.fastq.gz,${ATLAS_DIR}/${mapping}/sequence_quality_control/${mapping}_QC_se.fastq.gz in2=${ATLAS_DIR}/${mapping}/sequence_quality_control/${mapping}_QC_R2.fastq.gz,null  outm=${MAPPING_DIR}/${mapping}.sam threads=${THREADS} trimreaddescriptions=t pairlen=1000 pairedonly=f mdtag=t xstag=fs nmtag=t sam=1.3 local=t ambiguous=best secondary=t ssao=t maxsites=10 -Xmx${MEMORY}G 2> ${MAPPING_DIR}/logs/contig_coverage_stats_${mapping}.log
			echo ""

			echo "rule pileup (${mapping}):"
			echo "pileup.sh in=${MAPPING_DIR}/${mapping}.sam threads=${THREADS} -Xmx${MEMORY}G covstats=${MAPPING_DIR}/contig_stats/postfilter_coverage_stats_${mapping}.txt hist=${MAPPING_DIR}/contig_stats/postfilter_coverage_histogram_${mapping}.txt basecov=${MAPPING_DIR}/contig_stats/postfilter_base_coverage_${mapping}.txt.gz concise=t physcov=t secondary=f 2>> ${MAPPING_DIR}/logs/contig_coverage_stats_${mapping}.log"
			echo ""
			pileup.sh in=${MAPPING_DIR}/${mapping}.sam threads=${THREADS} -Xmx${MEMORY}G covstats=${MAPPING_DIR}/contig_stats/postfilter_coverage_stats_${mapping}.txt hist=${MAPPING_DIR}/contig_stats/postfilter_coverage_histogram_${mapping}.txt basecov=${MAPPING_DIR}/contig_stats/postfilter_base_coverage_${mapping}.txt.gz concise=t physcov=t secondary=f 2>> ${MAPPING_DIR}/logs/contig_coverage_stats_${mapping}.log
			echo ""

			echo "rule convert_sam_to_bam (${mapping}):"
			echo "samtools view -@ ${THREADS} -u ${MAPPING_DIR}/${mapping}.sam | samtools sort -m 4G -@ ${THREADS} -T ${MAPPING_DIR}/${mapping}_tmp -o ${MAPPING_DIR}/${mapping}.bam -O bam"
			echo ""
			samtools view -@ ${THREADS} -u ${MAPPING_DIR}/${mapping}.sam | samtools sort -m 4G -@ ${THREADS} -T ${MAPPING_DIR}/${mapping}_tmp -o ${MAPPING_DIR}/${mapping}.bam -O bam
			echo ""

		done
	
	done

	echo ""

}

function map_new_feature_counts {
	# Description: runs featureCounts with multi-mapped files as input 
	# GLOBAL Params: OUTPUT_DIR; THREADS; mapping_names (array); read_mapping_samples (array)
	# Return: writes files/directories to disk
	
	echo "Getting feature counts for individual samples onto reference sequences..."
	
	for i in $(seq 1 ${#mapping_names[@]}); do
		# Set counter to be based on zero, not one
		local j=$((${i}-1))	

		# Make relevant directories for storing output
		local fc_output_dir="${MAPPING_DIR}/feature_counts"
		mkdir -p ${fc_output_dir}
	
		# Temporarily change the internal fields separator (IFS) to parse comma separators. See Vince Buffalo's "Bioinformatics Data Skills" (1st Ed.) chapter 12, pg 407 and corresponding Github page README at https://github.com/vsbuffalo/bds-files/tree/master/chapter-12-pipelines (accessed Nov 19, 2017)
		local OFS="$IFS"
		IFS=,
		
		# Get names of individual samples provided for that coassembly name
		local mapping_sample_IDs=(${read_mapping_samples[${j}]}) #read_mapping_samples is a comma separated file with the names of the samples
	
		# Fix the IFS
		IFS="$OFS"
		
		# Build array of BAM file locations - iteratively add each filepath
		for i in $(seq 1 ${#mapping_sample_IDs[@]}); do
			
			# Set counter to zero-order
			j=$((${i}-1))

			# Get mapping name
			mapping=${mapping_sample_IDs[${j}]}
			
			if [ $i == 1 ]; then
				bam_filepaths=("${MAPPING_DIR}/${mapping}.bam") 
			elif [ $i > 1 ]; then
				bam_filepaths+=("${MAPPING_DIR}/${mapping}.bam")
			fi
			
		done
		
	
		
		# Run featureCounts
		echo "rule run_feature_counts_multi_mapping:"
		echo "featureCounts --minOverlap 1 -p -F SAF -T ${THREADS} --primary -M --fraction -R CORE -a ${MAPPING_DIR}/*.saf -o ${fc_output_dir}/gene_counts.txt ${bam_filepaths[@]} 2> ${MAPPING_DIR}/logs/featurecounts.log"
		echo ""
		featureCounts --minOverlap 1 -p -F SAF -T ${THREADS} --primary -M --fraction -R CORE -a ${MAPPING_DIR}/*.saf -o ${fc_output_dir}/gene_counts.txt ${bam_filepaths[@]} 2> ${MAPPING_DIR}/logs/featurecounts.log
		echo ""

	done
	

}

function main {
	# Get date and time of start
	start_time=$(date)
	#run all functions
	read_map_to_coassemblies
	map_new_feature_counts
	end_time=$(date)
	echo "Started read mapping at ${start_time} and finished at ${end_time}."
	echo ""
}

main