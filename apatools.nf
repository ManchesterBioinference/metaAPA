#!/usr/bin/env nextflow

// defines process to execute APA tools

// set dsl version 
nextflow.enable.dsl=2

// include modules for workflow
include {
	polyApipe;
	Sierra;
	SCAPE;
	prepare_utr_region
} from "$launchDir/modules/apa_tools_process"


// define a new workflow
workflow {
	bamfile_channel = channel.fromFilePairs("${params.data_dir}/*.bam{,.bai}", checkIfExists: true, flat: true)
	if(params.method == "polyApipe") {
		polyApipe(params.polyApipe_dir, bamfile_channel, params.sample_name, params.num_cores, params.job_id)
	}
	else if(params.method == "Sierra") {
		Sierra(params.Sierra_dir, bamfile_channel, params.sample_name, params.barcode_file, params.reference_file, params.num_cores, params.job_id)
	}
	else if(params.method == "SCAPE") {
		def scape_bed_file = file("${params.genome_version}_utr_region/${params.genome_version}_utr.bed")
		if(scape_bed_file.exists()) {
			scape_bed_file = Channel.fromPath(scape_bed_file)
		} else {
			scape_bed_file = prepare_utr_region(params.SCAPE_dir, params.reference_file, params.genome_version)
		}
		SCAPE(params.SCAPE_dir, bamfile_channel, params.sample_name, scape_bed_file.first(), params.barcode_file, params.num_cores, params.job_id) // server-specific problem module load apps/gcc/htslib/1.13
		
	} 
	else if(params.method == "SCAPTURE") {
		def scapture_annotation_dir = file("${genome_version}_SCAPTURE_annotation")
		if(scapture_annotation_dir.exists()) {
			scapture_annotation_dir = Channel.fromPath(scapture_annotation_dir, type: 'dir')
		} else {
			scapture_annotation_dir = prepare_SCAPTURE_annotation(params.SCAPTURE_dir, params.sample_name, params.reference_file, params.genome_file, params.chromsize_file, params.genome_version)
		}
		SCAPTURE(params.SCAPTURE_dir, bamfile_channel, params.sample_name, params.barcode_file, scapture_annotation_dir, params.genome_file, params.polyaDB_file, params.genome_version, params.read_length, params.num_cores, params.job_id)
	}
	else {
		println 'Wrong method!'
	}

}

workflow.onComplete {
	println "Pipeline successfully completed at: $workflow.complete"
}

workflow.onError {
	println "Pipeline stopped with following error: $workflow.errorMessage"
}





