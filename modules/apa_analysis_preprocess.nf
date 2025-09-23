#!/usr/bin/env nextflow

// set dsl version 
nextflow.enable.dsl=2
// defines process to preprocess outputs of APA tools
process apa_annotation {

        publishDir("${params.apa_input}/${sample_name}", mode: "copy", overwrite: true)
	conda "/mnt/mr01-home01/m57549qz/.conda/envs/r_env"
        input:
				path prog_dir
				path gtf_file
				path apa_output_dir
				path barcode_file
				val sample_name
				val apa_method
				val core_num

        output:
				path "${apa_method}_peak_annotations.qs", emit: apa_annotation_file

        script:
				"""
				Rscript ${prog_dir}/s02_annotation.R --method ${apa_method} --nthread ${core_num} --gtf ${gtf_file} --APAdir ${apa_output_dir} --barcodes ${barcode_file} --pipelinedir ${prog_dir}
				"""
}

process apa_normalization {

        publishDir("${params.apa_input}/${sample_name}", mode: "copy", overwrite: true)
	conda "/mnt/mr01-home01/m57549qz/.conda/envs/r_env"
        input:
                                path prog_dir
                                path apa_output_dir
                                path apa_annotation
                                path barcode_file
                                path gene_dir
                                val sample_name
                                val apa_method
                                val core_num

        output:
                                tuple path("${sample_name}_${apa_method}_count_matrix.qs"), path("${sample_name}_${apa_method}_normalized_by_ranger_matrix.qs"), emit: apa_matrix_file

        script:
                                """
                                Rscript ${prog_dir}/s03_apa_normalization_v2.R --sample ${sample_name} --method ${apa_method} --APAannotation ${apa_annotation} --APAdir ${apa_output_dir} --genedir ${gene_dir} --barcodes ${barcode_file} --nthread ${core_num} --pipelinedir ${prog_dir}
                                """
}

process count_to_bed {

	publishDir("${params.apa_input}/${sample_name}", mode: "copy", overwrite: true)
	conda "/mnt/mr01-home01/m57549qz/.conda/envs/r_env"
	input:
		path prog_dir
		path apa_output_dir
		path barcode_file
		val sample_name
		val apa_method
		val core_num

	output:
		path "${sample_name}_${apa_method}_pA_site.bed", emit: apa_bed_file

	script:
		"""
		Rscript ${prog_dir}/s04_apa_count2bed.R --sample ${sample_name} --method ${apa_method} --barcodes ${barcode_file} --apa_count ${apa_output_dir} --nthread ${core_num}
		"""
}
