process BBMERGE {

    /* */

	tag "${sample_id}"

	errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
	maxRetries 2

	input:
	tuple val(sample_id), path(reads1), path(reads2)

	output:
	tuple val(sample_id), path("${sample_id}.merged.preclump.fastq.gz")

	script:
	"""
	bbmerge.sh \
	-Xmx3g \
	in1=`realpath ${reads1}` \
	in2=`realpath ${reads2}` \
	out=${sample_id}.merged.preclump.fastq.gz \
	outu=${sample_id}.unmerged.preclump.fastq.gz \
	strict k=93 extend2=80 rem ordered \
	ihist=${sample_id}_ihist_merge.txt \
	threads=${task.cpus}
	"""

}

process FIND_ADAPTER_SEQS {

	/* */

	tag "${sample_id}"

	errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
	maxRetries 2

	cpus 2

	input:
	tuple val(sample_id), path(reads)

	output:
	tuple val(sample_id), path(reads), path("${sample_id}_adapters.fasta")

	script:
	"""
    bbmerge.sh -Xmx3g in=`realpath ${reads}` outa="${sample_id}_adapters.fasta" # ow qin=33
	"""

}

process TRIM_ADAPTERS {

	/* */

	tag "${sample_id}"

	errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
	maxRetries 2

	cpus 2

	input:
	tuple val(sample_id), path(reads), path(adapters)

    output:
    tuple val(sample_id), path("${sample_id}_no_adapters.fastq.gz")

    script:
    """
	reformat.sh -Xmx3g in=`realpath ${reads}` \
	out=${sample_id}_no_adapters.fastq.gz \
	ref=`realpath ${adapters}` \
	uniquenames=t overwrite=true t=${task.cpus} -Xmx8g
    """

}

process REMOVE_OPTICAL_DUPLICATES {

	/*
	This process removes optical duplicates from the Illumina flow cell.
	*/

	tag "${sample_id}"

	errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
	maxRetries 2

	cpus 1
	memory 3.GB

	input:
	tuple val(sample_id), path(reads)

	output:
	tuple val(sample_id), path("${sample_id}_deduped.fastq.gz")

	script:
	if ( sample_id.toString().contains("SRR") )
		"""
		rename.sh \
		in=${reads} \
		out=${sample_id}_deduped.fastq.gz \
		addpairnum=t -Xmx3g
		"""
	else
		"""
		clumpify.sh -Xmx2g in=`realpath ${reads}` \
		out=${sample_id}_deduped.fastq.gz \
		threads=${task.cpus} \
		optical tossbrokenreads reorder -Xmx3g
		"""

}

process REMOVE_LOW_QUALITY_REGIONS {

	/*
	Low quality regions of each read are removed in this process.
	*/

	tag "${sample_id}"

	errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
	maxRetries 2

	cpus 1
	memory 3.GB

	input:
	tuple val(sample_id), path(reads)

	output:
	tuple val(sample_id), path("${sample_id}_filtered.fastq.gz")

	script:
	if ( sample_id.toString().contains("SRR") )
		"""
		clumpify.sh -Xmx2g in=`realpath ${reads}` \
		out=${sample_id}_filtered.fastq.gz \
		threads=${task.cpus} \
		reorder markduplicates -Xmx3g
		"""
	else
		"""
		filterbytile.sh -Xmx2g in=`realpath ${reads}` \
		out=${sample_id}_filtered.fastq.gz \
		overwrite=true threads=${task.cpus}
		"""

}

process REMOVE_ARTIFACTS {

	/*
	Here we remove various contantimants that may have ended up in the reads,
	such as PhiX sequences that are often used as a sequencing control.
	*/

	tag "${sample_id}"

	errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
	maxRetries 2

	cpus 1
	memory 3.GB

	input:
	tuple val(sample_id), path(reads)

	output:
	tuple val(sample_id), path("${sample_id}_remove_artifacts.fastq.gz")

	script:
	"""
	bbduk.sh -Xmx2g in=`realpath ${reads}` \
	out=${sample_id}_remove_artifacts.fastq.gz \
	k=31 ref=artifacts,phix ordered cardinality \
	threads=${task.cpus}
	"""

}

process ERROR_CORRECT_PHASE_ONE {

	/*
	Bbmap recommends three phases of read error correction, the first of which
	goes through BBMerge.
	*/

	tag "${sample_id}"

	errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
	maxRetries 2

	cpus 1
	memory 3.GB

	input:
	tuple val(sample_id), path(reads)

	output:
	tuple val(sample_id), path("${sample_id}_error_correct1.fastq.gz")

	script:
	"""
	bbmerge.sh -Xmx2g in=`realpath ${reads}` \
	out=${sample_id}_error_correct1.fastq.gz \
	ecco mix vstrict ordered \
	ihist=${sample_id}_ihist_merge1.txt \
	threads=${task.cpus}
	"""

}

process ERROR_CORRECT_PHASE_TWO {

	/*
	The second phase of error correction goes through clumpify.sh
	*/

	tag "${sample_id}"

	cpus 1
	memory 3.GB

	input:
	tuple val(sample_id), path(reads)

	output:
	tuple val(sample_id), path("${sample_id}_error_correct2.fastq.gz")

	script:
	"""
	clumpify.sh -Xmx2g in=`realpath ${reads}` \
	out=${sample_id}_error_correct2.fastq.gz \
	ecc passes=4 reorder \
	threads=${task.cpus} \
	-Xmx2g
	"""

}

process ERROR_CORRECT_PHASE_THREE {

	/*
	The third phase of error correction uses tadpole.sh.
	*/

	tag "${sample_id}"

	errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
	maxRetries 2

	cpus 1
	memory 3.GB

	input:
	tuple val(sample_id), path(reads)

	output:
	tuple val(sample_id), path("${sample_id}_error_correct3.fastq.gz")

	script:
	"""
	tadpole.sh -Xmx2g in=`realpath ${reads}` \
	out=${sample_id}_error_correct3.fastq.gz \
	ecc k=62 ordered \
	threads=${task.cpus}
	"""

}

process QUALITY_TRIM {

	/*
	Here we quality trim reads from both ends to a minimum Phred quality of 10,
	and enforce a minimum read length of 70 bases.
	*/

	tag "${sample_id}"

	errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
	maxRetries 2

	cpus 1
	memory 3.GB

	input:
	tuple val(sample_id), path(reads)

	output:
	tuple val(sample_id), path("${sample_id}_qtrimmed.fastq.gz")

	script:
	"""
	bbduk.sh -Xmx2g in=`realpath ${reads}` \
	out=${sample_id}_qtrimmed.fastq.gz \
	qtrim=rl minavgquality=${params.min_qual} \
	minlength=${params.min_len} \
	ordered threads=${task.cpus}
	"""

}

process CLUMP_READS {

    /* */

	tag "${sample_id}"
    publishDir params.merged, mode: 'copy', overwrite: true

	errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
	maxRetries 2

	cpus 4

	input:
	tuple val(sample_id), path(reads)

	output:
    tuple val(sample_id), path("${sample_id}.merged.fastq.gz")

	script:
	"""
	clumpify.sh in=${reads} out=${sample_id}.merged.fastq.gz t=${task.cpus} reorder -Xmx3g
	"""

}
