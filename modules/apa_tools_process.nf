#!/usr/bin/env nextflow

// set dsl version 
nextflow.enable.dsl=2
// defines process to run APA tools

process Sierra{

	label "Sierra"
	publishDir "${sample_name}_Sierra_${job_id}", mode: 'copy'
	
	conda "/mnt/mr01-home01/m57549qz/.conda/envs/sierra"
	input:
		path prog_dir
		tuple val(file_name), path(bam_file), path(bai_file)
		val sample_name
		path barcodes_file
		path reference_file
		val num_cores
		val job_id

	output:
		path "Sierra_counts/barcodes.tsv.gz"
		path "Sierra_counts/matrix.mtx.gz"
		path "Sierra_counts/sitenames.tsv.gz"
		path "Sierra_peak.txt"
		path "${sample_name}_junctions.bed"

	script:
	"""
	regtools junctions extract -s RF ${bam_file} -o ${sample_name}_junctions.bed
	zcat ${barcodes_file} > barcodes.tsv
	Rscript ${prog_dir}/s01_Sierra.R --bam ${bam_file} --junctions ${sample_name}_junctions.bed --reference ${reference_file} --whitelist barcodes.tsv --nthreads ${num_cores}
	"""
}

process polyApipe{

	label "polyApipe"
	publishDir "${sample_name}_polyApipe_${job_id}", mode: 'copy'

	conda "/mnt/mr01-home01/m57549qz/.conda/envs/polyApipe_env"
	input:
		path prog_dir
		tuple val(file_name), path(bam_file), path(bai_file)
		val sample_name
		val num_cores
		val job_id

	output:
		path "${sample_name}_counts.tab.gz"
		path "${sample_name}_polyA.bam"
		path "${sample_name}_polyA.bam.bai"
		path "${sample_name}_polyA_peaks.gff"

	script:
	"""
	python3 ${prog_dir}/polyApipe.py -i ${bam_file} -o ${sample_name} -t ${num_cores}
	"""
}

process prepare_utr_region{

	label "prepare_utr_region"
	publishDir "${genome_version}_utr_region", mode: 'copy'

	conda "/mnt/mr01-home01/m57549qz/.conda/envs/scape"

	input:
		path prog_dir
		path reference_file
		val genome_version

	output:
		path "${genome_version}_utr.bed", emit: scape_bed_file
		path "${genome_version}_intron.bed"

	script:
	"""
	bedtools sort -i ${reference_file} | bgzip > ${genome_version}.genes.gtf.gz
	python3 ${prog_dir}/main.py prepare --gtf ${genome_version}.genes.gtf.gz --prefix ${genome_version}
	"""
}

process SCAPE{

	label "SCAPE"
	publishDir "${sample_name}_SCAPE_${job_id}", mode: 'copy'

	conda "/mnt/mr01-home01/m57549qz/.conda/envs/scape"

	input:
		path prog_dir
		tuple val(file_name), path(bam_file), path(bai_file)
		val sample_name
		path bed_file
		path barcode_file
		val num_cores
		val job_id

	output:
		path "${sample_name}_SCAPE/pasite.csv.gz"
		path "${sample_name}_SCAPE/huge"
		path "${sample_name}_SCAPE/TimeConsulting"
		path "${sample_name}_SCAPE/TooLongRegion"

	script:
	"""
	python3 ${prog_dir}/main.py apamix --bed ${bed_file} --bam ${bam_file} --out ${sample_name}_SCAPE --cores ${num_cores} --cb ${barcode_file}
	"""
}

process prepare_SCAPTURE_annotation{

	label "prepare_SCAPTURE_annotation"
	//publishDir "${sample_name}_SCAPTURE_${job_id}", mode: 'copy'

	conda "/mnt/mr01-home01/m57549qz/.conda/envs/SCAPTURE_env"

	input:
		path prog_dir
		val sample_name
		path reference_file
		path genome_file
		path chromsize_file		
		val genome_version

	output:
		path "${genome_version}_SCAPTURE_annotation", emit: scapture_annotation_dir

	script:
	"""
	# annotation module
	scapture -m annotation -o ${genome_version}_SCAPTURE_annotation -g ${genome_file} --gtf ${reference_file} --cs ${chromsize_file} --extend 2000 &> ${genome_version}.scapture.annotation.log
	"""
}

process SCAPTURE{

	label "SCAPTURE"
	publishDir "${sample_name}_SCAPTURE_${job_id}", mode: 'copy'

	conda "/mnt/mr01-home01/m57549qz/.conda/envs/SCAPTURE_env"

	input:
		path prog_dir
		tuple val(file_name), path(bam_file), path(bai_file)
		val sample_name
		path barcode_file
		path scapture_annotation_dir
		path genome_file
		path polyaDB_file
		val genome_version
		val read_length
		val num_cores
		val job_id

	output:
		path "${sample_name}.PASquant.KeepCell.UMIs.tsv.gz"
		path "${sample_name}.exonic.peaks.bed" // filtered exonic peaks
		path "${sample_name}.exonic.peaks.annotated.bed" // assigned exonic peaks
		path "${sample_name}.exonic.peaks.evaluated.bed" // evaluated exonic peaks
		path "${sample_name}.intronic.peaks.bed" // filtered intronic peaks
		path "${sample_name}.intronic.peaks.annotated.bed" // assigned intronic peaks
		path "${sample_name}.intronic.peaks.evaluated.bed" // evaluated intronic peaks
		path "${sample_name}.3primeExtended.peaks.annotated.bed" // assigned peaks in dowstream region beyond 3' end of genes
		path "${sample_name}.3primeExtended.peaks.evaluated.bed" // evaluated peaks in dowstream region beyond 3' end of genes

	script:
	"""
	# PAScall module
	if [ "${genome_version}" = "mm10" ]; then
	species="mouse"
	elif [ "${genome_version}" = "GRCm39" ]; then
	species="mouse"
	else
	echo "Error：invalid ${genome_version}"
	exit 1
	fi

	${prog_dir}/scapture -m PAScall -a ${scapture_annotation_dir} -g ${genome_file} -b ${bam_file} -l ${read_length} -o ${sample_name} -p ${num_cores} --species ${species} --polyaDB ${polyaDB_file} &> ${sample_name}.PAScall.log

	# Select PASs with positive prediction and konwn sites overlapped (recomanded)
	perl -alne '\$,="\\t"; print @F[0..11] if \$F[12] > 0 | \$F[13] eq "positive";' ${sample_name}.exonic.peaks.evaluated.bed ${sample_name}.intronic.peaks.evaluated.bed > ${sample_name}.PASquant.bed

	# PASquant module
	${prog_dir}/scapture -m PASquant -b ${bam_file} --celllist ${barcode_file} --pas ${sample_name}.PASquant.bed -o ${sample_name}.PASquant -p ${num_cores} &> ${sample_name}.PASquant.log
	"""
}

process DarPars2{

	label "DarPars2"
	publishDir "${sample_name}_DarPars2", mode: 'copy'

	conda "/mnt/mr01-home01/m57549qz/.conda/envs/flair"
	input:
		path prog_dir
		tuple val(sample_name), path(bam_file), path(bai_file)
		val num_cores

	output:
		path "${sample_name}_counts.tab.gz"
		path "${sample_name}_polyA.bam"
		path "${sample_name}_polyA.bam.bai"
		path "${sample_name}_polyA_peaks.gff"

	script:
	"""
	# step 1 run utr annotation script
	# step 2 generate mapped reads files for all samples
	## generate bedgraph
	bedtools genomecov -ibam ${bamfile} -bga -split -trackline > ${output_dir}/${samplename}.wig

	depth=`samtools view -F 0x904 -c ${bamfile}`
	echo -e ${output_dir}/${samplename}.wig"\t"${depth} > ${output_dir}/mapping_wig_location_with_depth.txt

	# step3 run DaPars2 (DaPars2_Multi_Sample_Multi_Chr.py)
	## generate configure file
	echo -e "# Specify the reference of 3'UTR region

	Annotated_3UTR=${utr_annotation}

	# A comma separated list of wig files of all samples

	Aligned_Wig_files=${output_dir}/${samplename}.wig

	Output_directory=${output_dir}

	Output_result_file=Dapars2

	# Specify Coverage threshold

	Coverage_threshold=10

	# Specify the number of threads to process the analysis

	Num_Threads=${corenum}

	# Provide sequencing depth file for normalization

	sequencing_depth_file=${output_dir}/mapping_wig_location_with_depth.txt
	" > ${output_dir}/Dapars2_configure_file

	## run
	python3 ${prog_dir}/DaPars2_Multi_Sample_Multi_Chr.py ${output_dir}/Dapars2_configure_file ${chr_list}

	mv ${output_dir}_* ${output_dir}

	"""
}
