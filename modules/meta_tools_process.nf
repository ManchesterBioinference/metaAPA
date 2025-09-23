#!/usr/bin/env nextflow

// set dsl version 
nextflow.enable.dsl=2
// defines process to preprocess outputs of APA tools

process merge_bedfiles {

	publishDir("${params.apa_input}/${sample_name}", mode: "copy", overwrite: true)
	input:
		val sample_name
		path apa_bed_folder

	output:
		path "${sample_name}_pA_site_sorted_uniq_k4.bed", emit: merged_bed_file

	script:
		"""
		cat ${apa_bed_folder}/${sample_name}_*_pA_site.bed | sort -k1,1 -k2,2n -k3,3n -k4,4 -u > ${sample_name}_pA_site_sorted_uniq_k4.bed
		"""
}

process position_based_method {

	publishDir("${params.metatools_output}/${sample_name}_position_based_method", mode: "copy", overwrite: true)
	conda "/mnt/mr01-home01/m57549qz/.conda/envs/scape_apa_env"
	input:
		path merged_bed_file
		val sample_name
		val distance

	output:
		path "${sample_name}_cluster_uniq_${distance}nt.bed", emit: position_based_output

	script:
		"""
		bedtools cluster -d ${distance} -s -i ${merged_bed_file} > ${sample_name}_cluster_uniq_${distance}nt.bed
		"""
}

process similarity_based_method {

	publishDir("${params.metatools_output}/${sample_name}_similarity_based_method", mode: "copy", overwrite: true)
	conda "/mnt/mr01-home01/m57549qz/.conda/envs/r_env"
	input:
		path prog_dir
		path apa_count_folder
		path gene_list
		path barcode_file
		val sample_name
		val core_num
		val apa_reads_cutoff
		tuple val(distance_method), val(cluster_method)

	output:
		path "${sample_name}_${distance_method}_${cluster_method}_${apa_reads_cutoff}filter_site_cluster.qs", emit: similarity_based_output

	script:
		"""
		Rscript ${prog_dir}/s05_meta_tools_v2.R --sample ${sample_name} --apa_count_folder ${apa_count_folder} --nthreads ${core_num} --apacutoff ${apa_reads_cutoff} --gene_list ${gene_list} --barcodes ${barcode_file} --distance_method ${distance_method} --cluster_method ${cluster_method} --pipelinedir ${prog_dir}
		"""
}

