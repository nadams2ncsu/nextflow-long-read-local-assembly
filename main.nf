#!/usr/bin/env nextflow

/*
########################################################################
 Nextflow Long-Read Local Assembly
########################################################################
*/

nextflow.enable.dsl = 2


/*
####################################
 Import process modules
####################################
*/

include { VALIDATE_INPUTS }	from './modules/validate_inputs'
include { COVERAGE_QC }         from './modules/coverage_qc'
include { COMBINE_QC }		from './modules/combine_qc'
include { READ_STATS } 		from './modules/read_stats'
include { COMBINE_READ_STATS }  from './modules/combine_read_stats'
include { EXTRACT_READS }       from './modules/extract_reads'
include { BAM2FASTQ }           from './modules/bam2fastq'
include { HIFIASM }             from './modules/hifiasm'
include { GFA2FASTA }           from './modules/gfa2fasta'
include { ALIGN_CONTIGS }       from './modules/align_contigs'
include { SORT_INDEX }          from './modules/sort_index'
include { SUMMARIZE }           from './modules/summarize'
include { SELECT_CONTIGUOUS_PAIRS } 	from './modules/select_contiguous_pairs'
include { RECONSTRUCT_SEQUENCE }        from './modules/reconstruct_sequence'
include { ACROSS_GROUPS_COMPARISON }	from './modules/across_groups_comparison'
include { CLASSIFY_HAPLOTYPE }  from './modules/classify_haplotypes'
include { COMBINE_RESULTS }     from './modules/combine_results'
include { NOTE_ASSEMBLY_DIFFERENCES } 	from './modules/note_assembly_differences'
include { MERGE_ASSEMBLY_DIFFERENCES } 	from './modules/merge_assembly_differences'
include { COMPARE_GROUP_DIFFERENCES } 	from './modules/compare_group_differences'
include { GENERATE_REPORT }	from './modules/generate_report'

/*
####################################
 Parameters
####################################
*/

params.samples = null
params.reference = null
params.min_median_read_length = 5000
//params.haplotype_config = null
params.outdir = 'results'
params.hifiasm_args = ''
params.minimap2_preset = 'asm5'
params.minimap2_args = ''
params.group_similarity_threshold = 90.0

params.default_flanks = [
    50,
    100,
    200,
    300,
    400,
    500,
    1000
]


/*
####################################
 Parameter validation
####################################
*/

if (!params.samples) {
    error "Missing required parameter: --samples"
}

if (!params.reference) {
    error "Missing required parameter: --reference"
}


/*
####################################
 Helper functions
####################################
*/

// Parse flanking region values in samples.tsv file
def parseFlanks(String flankString) {
    if (!flankString) {
        error "FLANKS value cannot be empty."
    }

    flankString = flankString.trim()
    if (flankString.toLowerCase() == 'all') {
        return params.default_flanks
    }

    try {
        return flankString
            .split(',')
            .collect { it.trim() as Integer }

    } catch (Exception e) {
        error """
        Invalid FLANKS value: '${flankString}'
        Expected:
            all
        or:
            50,100,200,400
        """
    }
}

// converts 1000 region to 1Mb;
def getRunName(Integer flankKb) {
    if (flankKb == 1000) {
        return '1Mb'
    }
    return "${flankKb}kb"
}

// parsing user-supplied coordintes
def parseCoordinates(String coordinates) {
    def matcher =
        coordinates =~ /^([^:]+):(\d+)-(\d+)$/

    if (!matcher.matches()) {
        error """
        Invalid coordinate format:
            ${coordinates}
        Expected:
            chromosome:start-end
        Example:
            chr1:161105457-162078654
        """
    }

    def chromosome =
        matcher[0][1]
    def start =
        matcher[0][2] as Long
    def end =
        matcher[0][3] as Long
    if (start >= end) {
        error """
        Invalid coordinates:
            ${coordinates}
        Start must be smaller than end.
        """
    }

    return [
        chromosome: chromosome,
        start:      start,
        end:        end
    ]
}

// calculate total regoin including flanks
def calculateRegion(
    String coordinates,
    Integer flankKb
) {

    def parsed =
        parseCoordinates(coordinates)

    def flankBp =
        flankKb * 1000L

    def regionStart =
        Math.max(
            1L,
            parsed.start - flankBp
        )

    def regionEnd =
        parsed.end + flankBp

    def region =
        "${parsed.chromosome}:${regionStart}-${regionEnd}"

    def hgSize =
        regionEnd - regionStart + 1

    return [
        chromosome: parsed.chromosome,
        start:      regionStart,
        end:        regionEnd,
        region:     region,
        hg_size:    hgSize
    ]
}


/*
####################################
 Input channels
####################################
*/

// reference genome
reference_ch = Channel
    .fromPath(
        params.reference,
        checkIfExists: true
    )
    .first()

// reference genome index
reference_fai_ch = Channel
    .fromPath(
        "${params.reference}.fai",
        checkIfExists: true
    )
    .first()

// samples from samples.tsv
samples_ch = Channel
    .fromPath(
        params.samples,
        checkIfExists: true
    )
    .splitCsv(
        header: true,
        sep: '\t',
        strip: true
    )

/*
 * Raw samples.tsv file.
 *
 * Used by VALIDATE_INPUTS before workflow execution.
 */

samples_file_ch = Channel
    .fromPath(
        params.samples,
        checkIfExists: true
    )
    .first()

/*
########################################################################
 Validate sample rows

 Validation occurs before branching so both QC and assembly receive
 validated input.
########################################################################
*/

validated_samples_ch = samples_ch.map { row ->

    /*
     * Required columns
     */

    def requiredColumns = [
        'sample',
        'group',
        'bam',
        'gene',
        'coordinates',
        'data_type',
        'flanks',
        'qc_config',
        'haplotype_config'
    ]

    requiredColumns.each { column ->

        if (
            !row.containsKey(column) ||
            !row[column]
        ) {

            error """
            Missing required value:
                ${column}

            Problematic row:
                ${row}
            """
        }
    }


    /*
     * Validate sequencing technology
     */

    def dataType =
        row.data_type
            .trim()
            .toUpperCase()

    if (!(dataType in ['ONT', 'PB'])) {

        error """
        Invalid data_type:
            ${row.data_type}

        Sample:
            ${row.sample}

        Supported values:
            ONT
            PB
        """
    }


    /*
     * Validate locus coordinates
     */

    parseCoordinates(
        row.coordinates
    )


    /*
     * Validate flank configuration
     */

    parseFlanks(
        row.flanks
    )


    /*
     * Validate locus-specific coverage QC configuration.
     *
     * Every samples.tsv row must provide a QC YAML.
     */

    def qcConfig =
        row.qc_config
            .trim()

    def qcConfigFile =
        file(qcConfig)

    if (!qcConfigFile.exists()) {

        error """
        QC configuration file not found:
            ${qcConfig}

        Sample:
            ${row.sample}

        Gene/locus:
            ${row.gene}
        """
    }

    row.qc_config = qcConfig


    /*
     * Validate locus-specific structural haplotype configuration.
     *
     * Use NA when structural haplotype classification should be skipped
     * for a locus.
     */

    def haplotypeConfig =
        row.haplotype_config
            .trim()

    if (haplotypeConfig.toUpperCase() != 'NA') {

        def haplotypeConfigFile =
            file(haplotypeConfig)

        if (!haplotypeConfigFile.exists()) {

            error """
            Haplotype configuration file not found:
                ${haplotypeConfig}

            Sample:
                ${row.sample}

            Gene/locus:
                ${row.gene}
            """
        }
    }

    row.haplotype_config =
        haplotypeConfig.toUpperCase() == 'NA'
            ? 'NA'
            : haplotypeConfig


    /*
     * Validate BAM and BAM index
     */

    def bam =
        file(row.bam)

    def bai =
        file("${row.bam}.bai")

    if (!bam.exists()) {

        error """
        BAM file not found:
            ${row.bam}
        """
    }

    if (!bai.exists()) {

        error """
        BAM index not found:
            ${row.bam}.bai
        """
    }


    /*
     * Store normalized sequencing technology
     */

    row.data_type = dataType

    return row
}


/*
########################################################################
 Branch samples

 coverage QC: One job per original samples.tsv row.
 read statistics QC: calculate read stats once per sample
 Assembly: Each sample is subsequently expanded across requested flank sizes.
#################################################################################
*/

sample_branches = validated_samples_ch.multiMap { row ->

    qc: row
    read_stats: row
    assembly: row
}

/*
##########################################################
 Coverage QC input

 QC occurs before flank expansion.

 Therefore: one samples.tsv row = one COVERAGE_QC job
##########################################################
*/

qc_input_ch = sample_branches.qc.map { row ->

    def meta = [

        sample:
            row.sample,
        group:
            row.group,
        gene:
            row.gene,
        data_type:
            row.data_type,
        qc_config:
            row.qc_config

    ]

    def bam =
        file(row.bam)

    def bai =
        file("${row.bam}.bai")

    def qcConfig =
        file(row.qc_config)

    tuple(
        meta,
        bam,
        bai,
        qcConfig
    )
}

/*
#############################
 READ STATISTICS INPUT
#############################
*/

read_stats_input_ch = sample_branches.read_stats.map { row ->

    def bam = file(row.bam)
    def bai = file("${row.bam}.bai")

    def meta = [
        sample:      row.sample,
        group:       row.group,
        gene:        row.gene,
        coordinates: row.coordinates,
        data_type:   row.data_type
    ]

    tuple(
        meta,
        bam,
        bai
    )
}

/*
#############################
 Expand assembly samples across flank sizes
#############################################
*/

sample_flank_ch = sample_branches.assembly
    .flatMap { row ->

        def flanks =
            parseFlanks(
                row.flanks
            )

        flanks.collect { flankKb ->

            def regionInfo =
                calculateRegion(
                    row.coordinates,
                    flankKb
                )

            def runName =
                getRunName(
                    flankKb
                )

            def meta = [

                sample:
                    row.sample,
                group:
                    row.group,
                gene:
                    row.gene,
		coordinates:
		    row.coordinates,
                data_type:
                    row.data_type,
                haplotype_config:
                    row.haplotype_config,
                flank_kb:
                    flankKb,
                run_name:
                    runName,
                chromosome:
                    regionInfo.chromosome,
                region_start:
                    regionInfo.start,
                region_end:
                    regionInfo.end,
                region:
                    regionInfo.region,
                hg_size:
                    regionInfo.hg_size
            ]

            def bam =
                file(row.bam)

            def bai =
                file("${row.bam}.bai")

            tuple(
                meta,
                bam,
                bai
            )
        }
    }


/*
################################################################################
Main long-read local assembly workflow
################################################################################
*/

workflow {
     /*
      * QUALITY CONTROL
      */

    /*
     * Sanity check user samples.tsv input files
     */ 

     VALIDATE_INPUTS(
    	samples_file_ch,
    	reference_ch,
	reference_fai_ch,
    	projectDir.toString()
     )
    
    /*
     * Gate analysis channels on successful input validation
     */
     validated_qc_input_ch = qc_input_ch
    	 .combine(VALIDATE_INPUTS.out.done)
    	 .map { meta, bam, bai, qc_config, validation_ok ->
        	 tuple(meta, bam, bai, qc_config)
      }

     validated_read_stats_input_ch = read_stats_input_ch
     	.combine(VALIDATE_INPUTS.out.done)
     	.map { meta, bam, bai, validation_ok ->
        	tuple(meta, bam, bai)
      }


      validated_sample_flank_ch = sample_flank_ch
     	.combine(VALIDATE_INPUTS.out.done)
    	.map { meta, bam, bai, validation_ok ->
        	tuple(meta, bam, bai)
     }

     /*
     * 1: Depth of coverage QC (branch 1)
     *
     * Calculate:
     *     - average whole-genome depth
     *     - gene interval depth
     *     - normalized gene depth
     *     - low coverage warnings
     *
     * Runs once per original samples.tsv row.
     *
     * QC warnings do not prevent assembly.
     */

    COVERAGE_QC(
        validated_qc_input_ch
    )


    /*
     * 1.1: Combine coverage QC into 1 summary text file
     */

     qc_results_ch = COVERAGE_QC.out.qc.map { meta, qc_tsv ->
    	qc_tsv
     }

     COMBINE_QC(
    	 qc_results_ch.collect()
     )

    /*
     * 1.2 Calculate read statistics for user-supplied coordinate region
     *
     *	- Calculates number of reads and lengths (min, median, max) for region
     *  - Only done once per sample
     *  - Any reads that overlap the region are included, they can extend outside of the region 
     */

     READ_STATS(
    	 validated_read_stats_input_ch
     )

     read_stats_results_ch = READ_STATS.out.stats.map {
    	 meta, stats_tsv -> stats_tsv
     }

     COMBINE_READ_STATS(
    	 read_stats_results_ch.collect()
     )

    /*
     * 2: Extract reads mapping to the requested genomic region.
     *
     * Input:
     *     tuple(meta, bam, bai)
     */

    EXTRACT_READS(
        validated_sample_flank_ch
    )


    /*
     * 3: Convert regional BAM to FASTQ.
     */

    BAM2FASTQ(
        EXTRACT_READS.out.bam
    )


    /*
     * 4: Perform localized assembly using hifiasm.
     *
     * ONT samples use ONT mode.
     * PacBio HiFi samples use standard hifiasm mode.
     */

    HIFIASM(
        BAM2FASTQ.out.fastq
    )


    /*
     * 5: Convert phased hifiasm GFA outputs to FASTA.
     */

    GFA2FASTA(
        HIFIASM.out.gfa
    )


    /*
     * Add assembly haplotype labels to metadata.
     */

    hap1_ch = GFA2FASTA.out.hap1.map {
        meta,
        fasta ->

        tuple(
            meta + [haplotype: 'hap1'],
            fasta
        )
    }


    hap2_ch = GFA2FASTA.out.hap2.map {
        meta,
        fasta ->

        tuple(
            meta + [haplotype: 'hap2'],
            fasta
        )
    }


    /*
     * Combine hap1 and hap2 FASTA channels.
     */

    assembly_fasta_ch =
        hap1_ch.mix(
            hap2_ch
        )


    /*
     * 6: Align assembled contigs to the workflow reference genome.
     */

    ALIGN_CONTIGS(
        assembly_fasta_ch,
        reference_ch
    )


    /*
     * 7: Sort and index assembled-contig alignments.
     */

    SORT_INDEX(
        ALIGN_CONTIGS.out.bam
    )


   /*
    * 8: Branch sorted BAMs.
    *
    * One branch is used by SUMMARIZE to determine if an assembly is contiguous + stats.
    *
    * The second branch is retained so paired-contiguous BAMs can
    * later be reconstructed and published.
    */

    sorted_bam_branches = SORT_INDEX.out.bam.multiMap {
    	 meta, bam, bai ->

    	 summarize:
        	tuple(meta, bam, bai, meta.haplotype_config)

    	reconstruct:
        	tuple(meta, bam, bai)
     }


    /*
     * 8.1: Generate alignment/assembly summary for each
     *      assembled haplotype.
     */

     SUMMARIZE(
     	 sorted_bam_branches.summarize
     )


   /*
    * 9: Identify sample + group + gene + flank combinations
    *    where BOTH hap1 and hap2 are contiguous.
    */

    contiguity_summary_ch = SUMMARIZE.out.summary.map {
    	 meta, summary_tsv ->

    	 summary_tsv
    }

    SELECT_CONTIGUOUS_PAIRS(
    	 contiguity_summary_ch.collect()
    )


/*
 * 9.1: Read paired-contiguous manifest and create a set
 *      of valid sample + group + gene + flank keys.
 */

paired_contiguous_keyset_ch = SELECT_CONTIGUOUS_PAIRS.out.pairs
    .splitCsv(
        header: true,
        sep: '\t',
        strip: true
    )
    .map { row ->
        "${row.sample}|${row.group}|${row.gene}|${row.flank_kb}"
    }
    .collect()
    .map { keys ->
        keys as Set
    }


/*
 * 9.2: Retain BOTH hap1 and hap2 BAMs for every
 *      paired-contiguous sample + group + gene + flank.
 */

contiguous_bams_ch = sorted_bam_branches.reconstruct
    .combine(paired_contiguous_keyset_ch)
    .filter {
        meta,
        bam,
        bai,
        valid_keys ->

        def key =
            "${meta.sample}|${meta.group}|${meta.gene}|${meta.flank_kb}"

        valid_keys.contains(key)
    }
    .map {
        meta,
        bam,
        bai,
        valid_keys ->

        tuple(
            meta,
            bam,
            bai
        )
    }


/*
 * 9.3: Separate paired-contiguous hap1 and hap2 BAMs.
 */

contiguous_hap1_ch = contiguous_bams_ch
    .filter {
        meta, bam, bai ->
        meta.haplotype == 'hap1'
    }
    .map {
        meta, bam, bai ->
        tuple(
            "${meta.sample}|${meta.group}|${meta.gene}|${meta.flank_kb}",
            meta,
            bam,
            bai
        )
    }


contiguous_hap2_ch = contiguous_bams_ch
    .filter {
        meta, bam, bai ->
        meta.haplotype == 'hap2'
    }
    .map {
        meta, bam, bai ->
        tuple(
            "${meta.sample}|${meta.group}|${meta.gene}|${meta.flank_kb}",
            meta,
            bam,
            bai
        )
    }


/*
 * 9.4: Join hap1 and hap2 into paired-contiguous assemblies.
 *
 * At this point each key occurs exactly once in the hap1 channel
 * and exactly once in the hap2 channel.
 */

paired_contiguous_bams_ch = contiguous_hap1_ch
    .join(contiguous_hap2_ch)
    .map {
        key,
        meta1,
        hap1_bam,
        hap1_bai,
        meta2,
        hap2_bam,
        hap2_bai ->

        tuple(
            meta1,
            hap1_bam,
            hap1_bai,
            hap2_bam,
            hap2_bai
        )
    }


	/*
 	* 9.6: Branch paired-contiguous assemblies.
 	*
 	* All paired-contiguous assemblies are reconstructed.
 	* The comparison branches are used for across-group comparison.
 	*/

	paired_contiguous_branches = paired_contiguous_bams_ch.multiMap {
    		meta,
    		hap1_bam,
    		hap1_bai,
    		hap2_bam,
   	 	hap2_bai ->

    	reconstruct:
        	tuple(
            	meta,
            	hap1_bam,
            	hap1_bai,
            	hap2_bam,
            	hap2_bai
        	)

	differences:

        	tuple(
            	meta,
            	hap1_bam,
            	hap1_bai,
            	hap2_bam,
            	hap2_bai
        	)

    	comparison_hap1:
        	tuple(
            	meta + [haplotype: 'hap1'],
            	hap1_bam,
            	hap1_bai
        	)

    	comparison_hap2:
        	tuple(
            	meta + [haplotype: 'hap2'],
            	hap2_bam,
            	hap2_bai
        	)
	}


	/*
 	* 9.7: Reconstruct and trim paired-contiguous assemblies.
 	*
	* Runs for ALL paired-contiguous assemblies, including samples
 	* occurring in only one group.
 	*/

	RECONSTRUCT_SEQUENCE(
    		paired_contiguous_branches.reconstruct
	)


	/*
 	* 9.8: Prepare paired-contiguous haplotypes for
 	* assembly-reference difference annotation.
 	*/

	difference_haplotypes_ch = paired_contiguous_branches.differences
    		.flatMap {
        		meta,
        		hap1_bam,
        		hap1_bai,
        		hap2_bam,
        		hap2_bai ->

        		[
            	tuple(
                	meta + [haplotype: 'hap1'],
                	hap1_bam,
                	hap1_bai
            	),

            	tuple(
                	meta + [haplotype: 'hap2'],
                	hap2_bam,
                	hap2_bai
            		)
        	]
    	}


	/*
 	* 9.9: Record differences between each haplotype assembly
 	* and the workflow reference.
 	*/

	NOTE_ASSEMBLY_DIFFERENCES(
    		difference_haplotypes_ch,
    		reference_ch
	)


	/*
 	* 9.10: Separate hap1 and hap2 difference files.
 	*/

	difference_hap1_ch = NOTE_ASSEMBLY_DIFFERENCES.out.differences
    		.filter {
        		meta, bed ->

        		meta.haplotype == 'hap1'
    		}
    		.map {
        		meta, bed ->

        		tuple(
            			"${meta.sample}|${meta.group}|${meta.gene}|${meta.flank_kb}",
            			meta,
            			bed
        		)
    		}

	difference_hap2_ch = NOTE_ASSEMBLY_DIFFERENCES.out.differences
    		.filter {
        		meta, bed ->

        		meta.haplotype == 'hap2'
    		}
    		.map {
        		meta, bed ->

        		tuple(
            			"${meta.sample}|${meta.group}|${meta.gene}|${meta.flank_kb}",
            			meta,
            			bed
        		)
    		}


	/*
 	* 9.11: Pair hap1 and hap2 difference files.
 	*/

	paired_differences_ch = difference_hap1_ch
    		.join(difference_hap2_ch)
    		.map {
        		key,
        		meta1,
        		hap1_bed,
        		meta2,
        		hap2_bed ->

        	def merged_meta = new LinkedHashMap(meta1)
        	merged_meta.remove('haplotype')

        	tuple(
            	merged_meta,
            	hap1_bed,
            	hap2_bed
        	)
    	}


	/*
 	* 9.12: Merge hap1 and hap2 assembly-reference differences.
 	*/

	MERGE_ASSEMBLY_DIFFERENCES(
    		paired_differences_ch
	)


	/*
 	* 9.13: Prepare merged assembly-reference difference BEDs
 	* for across-group comparison.
	 *
 	* Each merged BED represents one:
 	*
 	* sample + group + gene + flank
 	*
 	* COMPARE_GROUP_DIFFERENCES will compare every available flank
 	* combination between different groups for the same sample + gene.
 	*/

	difference_comparison_branches =
    		MERGE_ASSEMBLY_DIFFERENCES.out.differences.multiMap {
        	meta,
        	bed ->

        	metadata:
			"${meta.sample}\t${meta.group}\t${meta.gene}\t${meta.flank_kb}\t${bed.name}"

        	beds:
            		bed
    	}


	/*
 	* 9.14: Compare merged assembly-reference differences
	 * across groups.
 	*
 	* Samples occurring in only one group are ignored by the
 	* comparison script.
 	*
 	* Haplotype labels are not used to determine whether a
 	* difference is shared across groups.
 	*/

	COMPARE_GROUP_DIFFERENCES(
    		difference_comparison_branches.metadata.collect(),
    		difference_comparison_branches.beds.collect()
	)


	/*
 	* 9.15: Combine hap1 and hap2 BAMs for sequence-based
 	* across-group comparison.
 	*/

	comparison_bams_ch =
    		paired_contiguous_branches.comparison_hap1.mix(
        	paired_contiguous_branches.comparison_hap2
    	)



	/*
 	* 9.16: Prepare all paired-contiguous BAMs for across-group
 	* comparison.
	*
 	* ACROSS_GROUPS_COMPARISON determines whether a sample + gene
 	* occurs in at least two different groups.
 	*/

	comparison_branches = comparison_bams_ch.multiMap {
		meta,
		bam,
		bai ->
            	
		metadata:
		        "${meta.sample}\t${meta.group}\t${meta.gene}\t${meta.flank_kb}\t${meta.haplotype}\t${bam.name}"
		bams:
			bam
            	bais:
			bai
    	 }


	ACROSS_GROUPS_COMPARISON(
    		comparison_branches.metadata.collect(),
		comparison_branches.bams.collect(),
		comparison_branches.bais.collect()
	)

    /*
     * 10: Optional locus-specific structural haplotype classification.
     *
     * The haplotype_config column in samples.tsv controls classification
     * independently for each locus:
     *
     *     path/to/config.yaml  -> classify this locus
     *     NA                   -> skip structural classification
     */

     haplotype_summary_branches = SUMMARIZE.out.summary.branch {
    	 meta, summary_tsv ->

    	classify:
       		meta.haplotype_config &&
        	meta.haplotype_config.toUpperCase() != 'NA'

    	unclassified:
        	true
      }


classify_input_ch = haplotype_summary_branches.classify.map {
    meta, summary_tsv ->

    tuple(
        meta,
        summary_tsv,
        file(meta.haplotype_config)
    )
}

    CLASSIFY_HAPLOTYPE(
        classify_input_ch
    )


classified_summary_ch = CLASSIFY_HAPLOTYPE.out.classified.map {
    meta, classified_tsv ->
    classified_tsv
}


unclassified_summary_ch = haplotype_summary_branches.unclassified.map {
    meta, summary_tsv ->
    summary_tsv
}

    final_summary_input_ch =
        classified_summary_ch
            .mix(unclassified_summary_ch)
            .collect()


    /*
     * 11: Combine all sample/assembly results.
     *
     * This process ALWAYS runs.
     *
     * With haplotype classification:
     *     local_assembly_summary.tsv contains
     *     structural_haplotype.
     *
     * Without haplotype classification:
     *     local_assembly_summary.tsv contains
     *     the SUMMARIZE results without that column.
     *
     * COMBINE_RESULTS publishes the final file to:
     *
     *     results/summary/local_assembly_summary.tsv
     */

    COMBINE_RESULTS(
        final_summary_input_ch
    )

    /*
     * 12: Generate final interactive HTML report.
     *
     * Inputs:
     *   - final assembly summary
     *   - coverage QC summary
     *   - regional read statistics
     *   - across-group comparison summary
     */

    GENERATE_REPORT(
        COMBINE_RESULTS.out.summary,
        COMBINE_QC.out.summary,
        COMBINE_READ_STATS.out.summary,
        ACROSS_GROUPS_COMPARISON.out.comparison,
	COMPARE_GROUP_DIFFERENCES.out.comparison
    )

}
