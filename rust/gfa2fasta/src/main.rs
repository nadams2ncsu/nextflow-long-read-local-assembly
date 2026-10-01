use clap::Parser;
use std::fs::File;
use std::io::{
    self,
    BufRead,
    BufReader,
    BufWriter,
    Write,
};
use std::path::PathBuf;
use std::process;


/*
===============================================================================
 Command-line interface
===============================================================================
*/

#[derive(Parser, Debug)]
#[command(
    name = "gfa2fasta",
    version,
    about = "Convert GFA segment records to FASTA and calculate assembly statistics."
)]
struct Args {
    /// Input GFA file
    #[arg(short, long)]
    input: PathBuf,

    /// Output FASTA file
    #[arg(short, long)]
    output: PathBuf,

    /// Optional TSV file containing assembly statistics
    #[arg(short, long)]
    stats: Option<PathBuf>,
}


/*
===============================================================================
 Assembly statistics
===============================================================================
*/

#[derive(Debug)]
struct AssemblyStats {
    contigs: usize,
    total_length: usize,
    min_length: usize,
    max_length: usize,
    mean_length: f64,
    n50: usize,
}


/*
===============================================================================
 Main
===============================================================================
*/

fn main() {

    let args = Args::parse();

    if let Err(error) = run(&args) {

        eprintln!("ERROR: {error}");

        process::exit(1);
    }
}


/*
===============================================================================
 Workflow
===============================================================================
*/

fn run(args: &Args) -> Result<(), Box<dyn std::error::Error>> {

    /*
     * Open input GFA.
     */
    let input_file = File::open(&args.input)?;

    let reader = BufReader::new(input_file);


    /*
     * Create output FASTA.
     */
    let output_file = File::create(&args.output)?;

    let mut writer = BufWriter::new(output_file);


    /*
     * Store contig lengths for statistics.
     */
    let mut lengths: Vec<usize> = Vec::new();


    /*
     * Track how many segment records were found.
     */
    let mut segment_count: usize = 0;


    /*
     * Read GFA one line at a time.
     *
     * This avoids loading the complete GFA file into memory.
     */
    for line_result in reader.lines() {

        let line = line_result?;


        /*
         * Ignore empty lines.
         */
        if line.trim().is_empty() {
            continue;
        }


        /*
         * GFA fields are tab-delimited.
         */
        let mut fields = line.split('\t');


        let record_type = match fields.next() {
            Some(value) => value,
            None => continue,
        };


        /*
         * We only need segment records.
         *
         * GFA segment format:
         *
         * S    segment_name    sequence    optional_tags...
         */
        if record_type != "S" {
            continue;
        }


        let segment_name = fields.next().ok_or_else(|| {

            format!(
                "Malformed segment record: missing segment name:\n{}",
                line
            )

        })?;


        let sequence = fields.next().ok_or_else(|| {

            format!(
                "Malformed segment record '{}': missing sequence.",
                segment_name
            )

        })?;


        /*
         * Some GFA files can represent an unavailable sequence
         * using '*'.
         *
         * Such records cannot be converted into a meaningful
         * FASTA sequence, so fail rather than silently producing
         * incorrect output.
         */
        if sequence == "*" {

            return Err(
                format!(
                    "Segment '{}' does not contain an embedded sequence.",
                    segment_name
                )
                .into()
            );
        }


        /*
         * Prevent invalid empty FASTA sequences.
         */
        if sequence.is_empty() {

            return Err(
                format!(
                    "Segment '{}' contains an empty sequence.",
                    segment_name
                )
                .into()
            );
        }


        /*
         * Write FASTA record.
         */
        writeln!(
            writer,
            ">{}",
            segment_name
        )?;

        write_wrapped_sequence(
            &mut writer,
            sequence,
            80
        )?;


        /*
         * Record sequence length.
         */
        lengths.push(sequence.len());

        segment_count += 1;
    }


    /*
     * Ensure buffered output has been written.
     */
    writer.flush()?;


    /*
     * Fail if the input GFA contained no segment records.
     */
    if segment_count == 0 {

        return Err(
            format!(
                "No GFA segment ('S') records were found in '{}'.",
                args.input.display()
            )
            .into()
        );
    }


    /*
     * Calculate assembly statistics.
     */
    let stats = calculate_stats(&lengths);


    /*
     * Optionally write statistics.
     */
    if let Some(stats_path) = &args.stats {

        write_stats(
            stats_path,
            &stats
        )?;
    }


    /*
     * Print a concise execution summary to stderr.
     *
     * Nextflow will capture this in the task log.
     */
    eprintln!(
        "Converted {} segment(s) from '{}' to '{}'.",
        stats.contigs,
        args.input.display(),
        args.output.display()
    );

    eprintln!(
        "Total assembly length: {} bp | N50: {} bp",
        stats.total_length,
        stats.n50
    );


    Ok(())
}


/*
===============================================================================
 FASTA output
===============================================================================
*/

fn write_wrapped_sequence<W: Write>(
    writer: &mut W,
    sequence: &str,
    width: usize,
) -> io::Result<()> {

    /*
     * DNA sequences produced by hifiasm are ASCII.
     *
     * Operating on bytes therefore avoids problems with slicing
     * UTF-8 strings by byte offsets.
     */
    let bytes = sequence.as_bytes();

    for chunk in bytes.chunks(width) {

        writer.write_all(chunk)?;

        writer.write_all(b"\n")?;
    }

    Ok(())
}


/*
===============================================================================
 Statistics
===============================================================================
*/

fn calculate_stats(
    lengths: &[usize]
) -> AssemblyStats {

    let contigs = lengths.len();


    let total_length: usize =
        lengths.iter().sum();


    let min_length =
        *lengths
            .iter()
            .min()
            .unwrap_or(&0);


    let max_length =
        *lengths
            .iter()
            .max()
            .unwrap_or(&0);


    let mean_length = if contigs > 0 {

        total_length as f64 / contigs as f64

    } else {

        0.0
    };


    let n50 =
        calculate_n50(
            lengths,
            total_length
        );


    AssemblyStats {
        contigs,
        total_length,
        min_length,
        max_length,
        mean_length,
        n50,
    }
}


/*
===============================================================================
 N50 calculation
===============================================================================
*/

fn calculate_n50(
    lengths: &[usize],
    total_length: usize
) -> usize {

    if lengths.is_empty() {
        return 0;
    }


    /*
     * Copy lengths so the original vector is not modified.
     */
    let mut sorted_lengths =
        lengths.to_vec();


    /*
     * Largest contigs first.
     */
    sorted_lengths.sort_unstable_by(
        |a, b| b.cmp(a)
    );


    let half_length =
        (total_length + 1) / 2;


    let mut cumulative_length: usize = 0;


    for length in sorted_lengths {

        cumulative_length += length;

        if cumulative_length >= half_length {
            return length;
        }
    }


    0
}


/*
===============================================================================
 Statistics output
===============================================================================
*/

fn write_stats(
    path: &PathBuf,
    stats: &AssemblyStats
) -> io::Result<()> {

    let stats_file =
        File::create(path)?;

    let mut writer =
        BufWriter::new(stats_file);


    /*
     * Header
     */
    writeln!(
        writer,
        "contigs\ttotal_length\tmin_length\tmax_length\tmean_length\tn50"
    )?;


    /*
     * Values
     */
    writeln!(
        writer,
        "{}\t{}\t{}\t{}\t{:.2}\t{}",
        stats.contigs,
        stats.total_length,
        stats.min_length,
        stats.max_length,
        stats.mean_length,
        stats.n50
    )?;


    writer.flush()?;


    Ok(())
}
