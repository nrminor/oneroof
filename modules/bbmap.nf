process BBMERGE {

    /* */

	tag "${sample_id}"

	errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
	maxRetries 2

	time { "${5 * task.attempt}minutes" }
	cpus 4

	input:
	tuple val(sample_id), path(reads1), path(reads2)

	output:
	tuple val(sample_id), path("${sample_id}.merged.preclump.fastq.gz")

	script:
	def Xmx = 4 * task.attempt
	"""
	bbmerge.sh \
	-Xmx${Xmx}g \
	in1=`realpath ${reads1}` \
	in2=`realpath ${reads2}` \
	out=${sample_id}.merged.preclump.fastq.gz \
	outu=${sample_id}.unmerged.preclump.fastq.gz \
	strict k=93 extend2=80 rem ordered \
	ihist=${sample_id}_ihist_merge.txt \
	threads=${task.cpus} \
	-eoom
	"""

}

process FIND_ADAPTER_SEQS {

	/* */

	tag "${sample_id}"

	errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
	maxRetries 2

	time { "${5 * task.attempt}minutes" }

	cpus 2

	input:
	tuple val(sample_id), path(reads)

	output:
	tuple val(sample_id), path(reads), path("${sample_id}_adapters.fasta")

	script:
	def Xmx = 3 * task.attempt
	"""
    bbmerge.sh -Xmx${Xmx}g in=`realpath ${reads}` outa="${sample_id}_adapters.fasta" ow -eoom
	"""

}

process TRIM_ADAPTERS {

	/* */

	tag "${sample_id}"

	errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
	maxRetries 2

	time { "${5 * task.attempt}minutes" }

	cpus 2

	input:
	tuple val(sample_id), path(reads), path(adapters)

    output:
    tuple val(sample_id), path("${sample_id}_no_adapters.fastq.gz")

    script:
	def Xmx = 3 * task.attempt
    """
	reformat.sh -Xmx${Xmx}g in=`realpath ${reads}` \
	out=${sample_id}_no_adapters.fastq.gz \
	ref=`realpath ${adapters}` \
	uniquenames=t overwrite=true t=${task.cpus} -eoom
    """

}

process REMOVE_OPTICAL_DUPLICATES {

	/*
	This process removes optical duplicates from the Illumina flow cell.
	*/

	tag "${sample_id}"

	errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
	maxRetries 2

	time { "${5 * task.attempt}minutes" }

	cpus 1
	memory 3.GB

	input:
	tuple val(sample_id), path(reads)

	output:
	tuple val(sample_id), path("${sample_id}_deduped.fastq.gz")

	script:
	def Xmx = 3 * task.attempt
	if ( sample_id.toString().contains("SRR") )
		"""
		rename.sh \
		in=${reads} \
		out=${sample_id}_deduped.fastq.gz \
		addpairnum=t -Xmx${Xmx}g -eoom
		"""
	else
		"""
		clumpify.sh -Xmx${Xmx}g in=`realpath ${reads}` \
		out=${sample_id}_deduped.fastq.gz \
		threads=${task.cpus} \
		-eoom \
		optical tossbrokenreads reorder
		"""

}

process REMOVE_LOW_QUALITY_REGIONS {

	/*
	Low quality regions of each read are removed in this process.
	*/

	tag "${sample_id}"

	errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
	maxRetries 2

	time { "${5 * task.attempt}minutes" }

	cpus 1
	memory 3.GB

	input:
	tuple val(sample_id), path(reads)

	output:
	tuple val(sample_id), path("${sample_id}_filtered.fastq.gz")

	script:
	def Xmx = 3 * task.attempt
	if ( sample_id.toString().contains("SRR") )
		"""
		clumpify.sh -Xmx${Xmx}g in=`realpath ${reads}` \
		out=${sample_id}_filtered.fastq.gz \
		threads=${task.cpus} -eoom \
		reorder markduplicates
		"""
	else
		"""
		filterbytile.sh -Xmx${Xmx}g in=`realpath ${reads}` \
		out=${sample_id}_filtered.fastq.gz \
		-da overwrite=true threads=${task.cpus} -eoom
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

	time { "${5 * task.attempt}minutes" }

	cpus 1
	memory 3.GB

	input:
	tuple val(sample_id), path(reads)

	output:
	tuple val(sample_id), path("${sample_id}_remove_artifacts.fastq.gz")

	script:
	def Xmx = 3 * task.attempt
	"""
	bbduk.sh -Xmx${Xmx}g in=`realpath ${reads}` \
	out=${sample_id}_remove_artifacts.fastq.gz \
	k=31 ref=artifacts,phix ordered cardinality \
	threads=${task.cpus} -eoom
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

	time { "${5 * task.attempt}minutes" }

	cpus 1
	memory 3.GB

	input:
	tuple val(sample_id), path(reads)

	output:
	tuple val(sample_id), path("${sample_id}_error_correct1.fastq.gz")

	script:
	def Xmx = 3 * task.attempt
	"""
	bbmerge.sh -Xmx${Xmx}g in=`realpath ${reads}` \
	out=${sample_id}_error_correct1.fastq.gz \
	ecco mix vstrict ordered \
	ihist=${sample_id}_ihist_merge1.txt \
	threads=${task.cpus} -eoom
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
	def Xmx = 3 * task.attempt
	"""
	clumpify.sh -Xmx${Xmx}g in=`realpath ${reads}` \
	out=${sample_id}_error_correct2.fastq.gz \
	ecc passes=4 reorder \
	threads=${task.cpus} -eoom
	"""

}

process ERROR_CORRECT_PHASE_THREE {

	/*
	The third phase of error correction uses tadpole.sh.
	*/

	tag "${sample_id}"

	errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
	maxRetries 2

	time { "${5 * task.attempt}minutes" }

	cpus 1
	memory 3.GB

	input:
	tuple val(sample_id), path(reads)

	output:
	tuple val(sample_id), path("${sample_id}_error_correct3.fastq.gz")

	script:
	def Xmx = 3 * task.attempt
	"""
	tadpole.sh -Xmx${Xmx}g in=`realpath ${reads}` \
	out=${sample_id}_error_correct3.fastq.gz \
	ecc k=62 ordered \
	threads=${task.cpus} -eoom
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

	time { "${5 * task.attempt}minutes" }

	cpus 1
	memory 3.GB

	input:
	tuple val(sample_id), path(reads)

	output:
	tuple val(sample_id), path("${sample_id}_qtrimmed.fastq.gz")

	script:
	def Xmx = 3 * task.attempt
	"""
	bbduk.sh -Xmx${Xmx}g in=`realpath ${reads}` \
	out=${sample_id}_qtrimmed.fastq.gz \
	qtrim=rl minavgquality=${params.min_qual} \
	minlength=${params.min_len} \
	ordered threads=${task.cpus} -eoom
	"""

}

process CLUMP_READS {

    /* */

	tag "${sample_id}"
    publishDir params.merged, mode: 'copy', overwrite: true

	errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
	maxRetries 2

	time { "${5 * task.attempt}minutes" }

	cpus 4

	input:
	tuple val(sample_id), path(reads)

	output:
    tuple val(sample_id), path("${sample_id}.merged.fastq.gz")

	script:
	def Xmx = 3 * task.attempt
	"""
	clumpify.sh \
	in=${reads} out=${sample_id}.merged.fastq.gz \
	reorder \
	threads=${task.cpus} -Xmx${Xmx}g -eoom
	"""

}
