#!/usr/bin/env nextflow

// defines process to execute APA tools

// set dsl version 
nextflow.enable.dsl=2

// include modules for workflow
include {
	merge_bedfiles;
	position_based_method;
	similarity_based_method
} from "$launchDir/modules/meta_tools_process"


// define a new workflow
workflow {
	// Define inputs
	def sample_name = params.sample_name
	def apa_rscript = file(params.Sierra_dir)
	def barcode_file = file(params.barcode_file)
	//def distance = params.distance
	def apa_reads_cutoff = 0
	//def distance_method = params.distance_method
	//def cluster_method = params.cluster_method
        apa_bed_folder = channel.fromPath("${params.apa_input}/${sample_name}", type: "dir", checkIfExists: true)
        apa_count_folder = channel.fromPath("${params.apa_input}/${sample_name}", type: "dir", checkIfExists: true)
	// Merge bed files
        def merged_bed_file = file("${params.apa_input}/${sample_name}/${sample_name}_pA_site_sorted_uniq_k4.bed")
        // Check if the file exists
        if (merged_bed_file.exists()) {
                // If file exists, create a channel from the file
                merged_bed_channel = channel.fromPath(merged_bed_file)
                merged_bed_channel.view { "File exists, importing: $it" }
        } else {
                // If file does not exist, run the process
                merged_bed_channel = merge_bedfiles(sample_name, apa_bed_folder)
                merged_bed_channel.view { "File does not exist, running process to generate: $it" }
        }
        
        // Distance threshold for strategy 1 to channel
        if (params.distance instanceof String && params.distance.startsWith('[') && params.distance.endsWith(']')) {
                def parsed_distance = params.distance.replaceAll(/\[|\]/, '').split(',').collect { it.trim().toInteger() }
                def (start, end, step) = parsed_distance
                distance = (start..end).step(step.toInteger()).toList()
        } else {
        // Treat win_size as a single number
                distance = [params.distance]
        }
        distance_channel = channel.from(distance)
        distance_channel.view{} 
        // Distance method for strategy 2
 	if (params.distance_method instanceof String && params.distance_method.startsWith('[') && params.distance_method.endsWith(']')) { //
                def parsed_distance_method = params.distance_method.replaceAll(/\[|\]/, '').split(',').collect { it.trim().toString() }
                distance_method = parsed_distance_method.toList()
        } else {
        // Treat win_size as a single number
                distance_method = [params.distance_method]
        }

        // Emit distance_method to channel
        distance_method_channel = channel.from(distance_method)
        
        // Cluster method for strategy 2
 	if (params.cluster_method instanceof String && params.cluster_method.startsWith('[') && params.cluster_method.endsWith(']')) { //
                def parsed_cluster_method = params.cluster_method.replaceAll(/\[|\]/, '').split(',').collect { it.trim().toString() }
                cluster_method = parsed_cluster_method.toList()
        } else {
        // Treat win_size as a single number
                cluster_method = [params.cluster_method]
        }

        // Emit cluster_method to channel
        cluster_method_channel = channel.from(cluster_method)        
           
        combination_methods = distance_method_channel.flatten().combine(cluster_method_channel.flatten())
        combination_methods.view{} 
        // Strategy 1: Position-based method
        // position_based_output = position_based_method(merged_bed_channel.first(), sample_name, distance_channel)
        // position_based_output.view { "Running position-based method to generate: $it" }
        // Strategy 2: Similarity-based method
	similarity_based_output = similarity_based_method(apa_rscript, apa_count_folder.first(), merged_bed_channel.first(), barcode_file, sample_name, params.core_num, apa_reads_cutoff, combination_methods)
	similarity_based_output.view { "Running similarity-based method to generate: $it" }
}

workflow.onComplete {
	println "Pipeline successfully completed at: $workflow.complete"
}

workflow.onError {
	println "Pipeline stopped with following error: $workflow.errorMessage"
}





