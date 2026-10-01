#!/usr/bin/env python3

import argparse, csv, os, re, subprocess, sys
from pathlib import Path
import yaml

###########################################################################################
# Validates input files exist and are the correct structure prior to analysis
#
# - inputs: sample.tsv && reference genome
# - outputs: sanity check print messages

################
# parameters

## expected columns for samples.tsv file
EXPECTED_COLUMNS = [
    'sample','group','bam','gene','coordinates','data_type','flanks',
    'qc_config','haplotype_config'
]
COORD_RE = re.compile(r'^([^:]+):(\d+)-(\d+)$')

###############################
# Helper functions

## error message printer
def fail(msg):
    raise ValueError(msg)

## path finder
def resolve_path(value, base_dir):
    p = Path(value)
    return p if p.is_absolute() else Path(base_dir) / p

#####
# Coordinates & BAM file checks

## read in user-supplied coordinates && check that they exist with BAM file header
def parse_coord(value, label):
    m = COORD_RE.fullmatch(value.strip())
    if not m:
        fail(f"{label}: invalid coordinate '{value}'; expected chromosome:start-end")
    chrom, start, end = m.group(1), int(m.group(2)), int(m.group(3))
    if start < 1 or start >= end:
        fail(f"{label}: invalid coordinate '{value}'; require 1 <= start < end")
    return chrom, start, end


def bam_contigs(bam):
    result = subprocess.run(['samtools','view','-H',str(bam)], text=True,
                            capture_output=True, check=True)
    contigs = {}
    for line in result.stdout.splitlines():
        if not line.startswith('@SQ\t'):
            continue
        fields = dict(x.split(':',1) for x in line.split('\t')[1:] if ':' in x)
        if 'SN' in fields and 'LN' in fields:
            contigs[fields['SN']] = int(fields['LN'])
    return contigs


def check_coord_in_bam(value, contigs, label):
    chrom, start, end = parse_coord(value, label)
    if chrom not in contigs:
        fail(f"{label}: chromosome '{chrom}' is not present in BAM header")
    if end > contigs[chrom]:
        fail(f"{label}: end {end} exceeds {chrom} length {contigs[chrom]} in BAM ({chrom}:{contigs[chrom]})")

#####################
# YAML file checks

## does yaml file exist?
def load_yaml(path, label):
    if not path.is_file():
        fail(f"{label}: file not found: {path}")
    with path.open() as handle:
        data = yaml.safe_load(handle)
    if not isinstance(data, dict):
        fail(f"{label}: YAML root must be a mapping: {path}")
    return data

## is the QC yaml file the correct structure?
def validate_qc_yaml(path, row, contigs):
    cfg = load_yaml(path, 'QC config')
    for key in ('region','datasets','references'):
        if key not in cfg:
            fail(f"QC config {path}: missing required key '{key}'")
    if str(cfg['region']) != row['gene']:
        fail(f"QC config {path}: region '{cfg['region']}' does not match samples.tsv gene '{row['gene']}'")
    if not isinstance(cfg['datasets'], dict) or row['group'] not in cfg['datasets']:
        fail(f"QC config {path}: group '{row['group']}' is not defined under datasets")
    ds = cfg['datasets'][row['group']]
    if not isinstance(ds, dict) or not ds.get('reference'):
        fail(f"QC config {path}: datasets.{row['group']}.reference is missing")
    ref = ds['reference']
    refs = cfg['references']
    if not isinstance(refs, dict) or ref not in refs:
        fail(f"QC config {path}: reference '{ref}' is not defined under references")
    loci = refs[ref]
    if not isinstance(loci, dict) or not loci:
        fail(f"QC config {path}: references.{ref} must contain at least one QC interval set")
    for name, settings in loci.items():
        if not isinstance(settings, dict) or not isinstance(settings.get('coordinates'), list) or not settings['coordinates']:
            fail(f"QC config {path}: references.{ref}.{name}.coordinates must be a non-empty list")
        for coord in settings['coordinates']:
            check_coord_in_bam(str(coord), contigs, f"QC config {path} [{ref}.{name}]")

## is HAPLOTYPE yaml file the correct structure
def validate_haplotype_yaml(path, row):
    cfg = load_yaml(path, 'Haplotype config')
    for key in ('region','defaults','haplotypes'):
        if key not in cfg:
            fail(f"Haplotype config {path}: missing required key '{key}'")
    if str(cfg['region']) != row['gene']:
        fail(f"Haplotype config {path}: region '{cfg['region']}' does not match samples.tsv gene '{row['gene']}'")
    if not isinstance(cfg['defaults'], dict):
        fail(f"Haplotype config {path}: defaults must be a mapping")
    if not isinstance(cfg['haplotypes'], dict) or not cfg['haplotypes']:
        fail(f"Haplotype config {path}: haplotypes must be a non-empty mapping")
    for hap, settings in cfg['haplotypes'].items():
        if not isinstance(settings, dict) or 'structural_events' not in settings:
            fail(f"Haplotype config {path}: {hap}.structural_events is required")
        events = settings['structural_events']
        if not isinstance(events, list):
            fail(f"Haplotype config {path}: {hap}.structural_events must be a list")
        for i, event in enumerate(events, 1):
            if not isinstance(event, dict):
                fail(f"Haplotype config {path}: {hap} event {i} must be a mapping")
            if event.get('type') not in ('insertion','deletion'):
                fail(f"Haplotype config {path}: {hap} event {i} type must be insertion or deletion")
            parse_coord(str(event.get('region','')), f"Haplotype config {path} [{hap} event {i}]")
            if not isinstance(event.get('length'), int) or event['length'] <= 0:
                fail(f"Haplotype config {path}: {hap} event {i} length must be a positive integer")
            if not isinstance(event.get('count'), int) or event['count'] <= 0:
                fail(f"Haplotype config {path}: {hap} event {i} count must be a positive integer")

###########################################################################################
# Run sanity checks

def main():

    ## required user-arguments
    ap = argparse.ArgumentParser()
    ap.add_argument('--samples', required=True)
    ap.add_argument('--reference', required=True)
    ap.add_argument('--reference-fai', required=True)
    ap.add_argument('--base-dir', required=True)
    ap.add_argument('--output', required=True)
    args = ap.parse_args()

    ## user-supplied required args?
    samples = Path(args.samples)
    reference = Path(args.reference)
    reference_fai = Path(args.reference_fai) 

    if not reference.is_file(): fail(f"Reference FASTA not found: {reference}")
    if not reference_fai.is_file(): fail(f"Reference FASTA index not found: {reference_fai}")

    ## reading in input samples from user samples.tsv
    with samples.open(newline='') as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        if reader.fieldnames != EXPECTED_COLUMNS:
            fail('samples.tsv header must be exactly:\n  ' + '\t'.join(EXPECTED_COLUMNS) +
                 '\nFound:\n  ' + '\t'.join(reader.fieldnames or []))
        rows = list(reader)
    if not rows: fail('samples.tsv contains no data rows')

    print('Validating workflow inputs...')
    print('samples.tsv structure................ PASS')
    print('Reference FASTA...................... PASS')
    print('Reference FASTA index................ PASS')

    ## correct number of columns && non-empty
    for n, row in enumerate(rows, 2):
        for col in EXPECTED_COLUMNS:
            if not row[col].strip(): fail(f"samples.tsv line {n}: '{col}' is empty")

        ## values accepted for data_type & flank
        if row['data_type'].strip().upper() not in ('ONT','PB'):
            fail(f"samples.tsv line {n}: data_type must be ONT or PB")
        flank = row['flanks'].strip().lower()
        if flank != 'all':
            try:
                vals = [int(x.strip()) for x in flank.split(',')]
                if not vals or any(x <= 0 for x in vals): raise ValueError
            except ValueError:
                fail(f"samples.tsv line {n}: invalid flanks '{row['flanks']}'")
        
        ## path to WGS BAM file exists?
        bam = resolve_path(row['bam'].strip(), args.base_dir)
        bai = Path(str(bam) + '.bai')
        if not bam.is_file(): fail(f"samples.tsv line {n}: BAM not found: {bam}")
        if not bai.is_file(): fail(f"samples.tsv line {n}: BAM index not found: {bai}")
        subprocess.run(['samtools','quickcheck','-v',str(bam)], check=True)
        contigs = bam_contigs(bam)
        check_coord_in_bam(row['coordinates'], contigs, f"samples.tsv line {n} coordinates")

        qc = resolve_path(row['qc_config'].strip(), args.base_dir)
        validate_qc_yaml(qc, row, contigs)

        hap = row['haplotype_config'].strip()
        if hap.upper() != 'NA':
            validate_haplotype_yaml(resolve_path(hap, args.base_dir), row)

        print(f"{row['sample']} / {row['group']} / {row['gene']}........ PASS")

    Path(args.output).write_text('PASS\n')
    print('Input validation PASSED.')

if __name__ == '__main__':
    try:
        main()
    except (ValueError, yaml.YAMLError, subprocess.CalledProcessError) as e:
        print(f"ERROR: input validation failed: {e}", file=sys.stderr)
        sys.exit(1)

