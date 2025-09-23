#!/usr/bin/env nextflow

// defines process to execute APA tools

// set dsl version 
nextflow.enable.dsl=2

// include modules for workflow
include {
	apa_annotation;
	apa_normalization;
	count_to_bed
} from "$launchDir/modules/apa_analysis_preprocess"


// define a new workflow
workflow {
	// Define inputs
	def sample_name = params.sample_name
	def apa_method = params.method
	def gtf_file = file(params.reference_file)
	def apa_rscript = file(params.Sierra_dir)
	def barcode_file = file(params.barcode_file)
	apa_output = channel.fromPath("${params.apa_input}/${sample_name}_${apa_method}_*", type: "dir", checkIfExists: true)
	gene_dir = channel.fromPath("${params.gene_dir}", type: "dir", checkIfExists: true)
	//apa_output = channel.fromPath("${launchDir}/data/${sample_name}_illumina/*_${apa_method}_*", type: "dir", checkIfExists: true)
	//barcode_file = channel.fromPath("${launchDir}/data/${sample_name}_illumina/${sample_name}_barcodes.tsv*", type: "file", checkIfExists: true)

	// annotation 	
	def apa_annotation_file = file("${params.annotation_dir}/${apa_method}_peak_annotations.qs")
	// Check if the file exists
        if (apa_annotation_file.exists()) {
                // If file exists, create a channel from the file
                apa_annotation_channel = channel.fromPath(apa_annotation_file)
                apa_annotation_channel.view { "File exists, importing: $it" }
        } else {
                // If file does not exist, run the process
                apa_annotation_channel = apa_annotation(apa_rscript, gtf_file, apa_output, barcode_file, sample_name, apa_method, params.core_num)
                apa_annotation_channel.view { "File does not exist, running process to generate: $it" }
        }
        
	// normalization
        def apa_count_file = file("${params.apa_input}/${sample_name}/${sample_name}_${apa_method}_count_matrix.qs")
        def apa_normalization_file = file("${params.apa_input}/${sample_name}/${sample_name}_${apa_method}_normalized_by_ranger_matrix.qs")

        // Check if the file exists
        if (apa_count_file.exists() && apa_normalization_file.exists()) {
                // If file exists, create a channel from the file
                apa_count_channel = channel.fromPath(apa_count_file)
                apa_normalization_channel = channel.fromPath(apa_normalization_file)

                apa_count_channel.view { "File1: $it" }
                apa_normalization_channel.view { "File2: $it" }
        } else {
                // If file does not exist, run the process
                apa_matrix_channel = apa_normalization(apa_rscript, apa_output, apa_annotation_channel.first(), barcode_file, gene_dir, sample_name, apa_method, params.core_num)

                apa_count_channel = apa_matrix_channel.map { it[0] }
                apa_normalization_channel = apa_matrix_channel.map { it[1] }

                apa_matrix_channel.view { result ->
    "File does not exist, running process to generate: File1: ${result[0]}, File2: ${result[1]}" }
        }

        // Convert count files into bed files 	
	def apa_bed_file = file("${params.apa_input}/${sample_name}/${sample_name}_${apa_method}_pA_site.bed")
	// Check if the file exists
        if (apa_bed_file.exists()) {
                // If file exists, create a channel from the file
                apa_bed_channel = channel.fromPath(apa_bed_file)
                apa_bed_channel.view { "File exists, importing: $it" }
        } else {
                // If file does not exist, run the process
                apa_bed_channel = count_to_bed(apa_rscript, apa_count_channel, barcode_file, sample_name, apa_method, params.core_num)
                apa_bed_channel.view { "File does not exist, running process to generate: $it" }
        }

}

workflow.onComplete {
	println "Pipeline successfully completed at: $workflow.complete"
}

workflow.onError {
	println "Pipeline stopped with following error: $workflow.errorMessage"
}





