# Troubleshooting

This document describes common setup and execution issues that may occur when running the long-read local assembly workflow on an HPC system.

## Workflow Does Not Launch

Confirm that **Nextflow** and **Apptainer or Singularity** are installed and available in your environment.

```bash
nextflow -version
```

For Apptainer:

```bash
apptainer --version
```

or Singularity:

```bash
singularity --version
```

The workflow container includes the software used by individual workflow processes, but Nextflow and Apptainer/Singularity must be installed and available on the host HPC system.

If your HPC system uses environment modules, the required software may need to be loaded before running the workflow.

For example:

```bash
module load nextflow
module load apptainer
```

Available module names will vary between HPC systems.

## Permission Denied

The workflow uses executable shell scripts and Python utilities. If a `Permission denied` error occurs, confirm that the required files are executable.

From the repository root:

```bash
chmod +x submit_local_asm_workflow.sh run_nextflow.sh
chmod +x bin/*
```
You can inspect file permissions with:

```bash
ls -l submit_local_asm_workflow.sh run_nextflow.sh
ls -l bin/
```

## Input Files Cannot Be Found

Use direct paths to input files and configuration files in `samples.tsv`.

For example:

```text
/path/to/sample.bam
/path/to/qc_config.yaml
/path/to/haplotype_config.yaml
```

Confirm that each path exists and is accessible from the HPC compute nodes.

For BAM inputs, also confirm that the corresponding BAM index is present:

```bash
ls -lh /path/to/sample.bam
ls -lh /path/to/sample.bam.bai
```

Paths in `samples.tsv` should point to the actual locations of the files on the HPC filesystem rather than assuming that files are located relative to the workflow directory.

## Jobs Are Killed or Too Many Jobs Are Submitted

The number of tasks that Nextflow can submit concurrently is controlled by `--queue_size`.

If the workflow is submitting more jobs than appropriate for the HPC environment or encountering scheduler submission limits, rerun the workflow with a smaller queue size.

For example:

```bash
./submit_local_asm_workflow.sh samples.tsv reference.fasta --queue_size 5
```

The appropriate queue size depends on the HPC system, available resources, and scheduler limits.

A job may also be killed because it exceeds its requested memory or walltime. Check the scheduler information and job logs before assuming that workflow concurrency is the cause.

For SLURM systems, active jobs can be inspected with:

```bash
squeue -u $USER
```

Completed or failed jobs can be examined with:

```bash
sacct -u $USER
```

If a process consistently exceeds its allocated resources, the corresponding CPU, memory, or execution-time settings in `nextflow.config` may need to be adjusted for the HPC environment.

## Confirm the HPC Scheduler

The provided workflow configuration is designed for the **SLURM** workload manager.

Confirm that SLURM is available on the HPC system:

```bash
sbatch --version
squeue --version
```

The submission scripts and `nextflow.config` use SLURM for workflow execution.

If the HPC system uses a different scheduler, the Nextflow executor and associated scheduler configuration must be modified before running the workflow.

## Workflow Starts but a Process Fails

If Nextflow launches successfully but an individual process fails, first inspect the error reported by Nextflow and the corresponding SLURM logs.

Useful files and locations include:

- `.nextflow.log`
- SLURM output and error logs
- the failed process directory under `work/`
- Nextflow execution information in `pipeline_info/`, when available

Nextflow normally reports the `work/` directory associated with a failed process. Inspecting the `.command.*` files within that directory can help identify the command that failed and its output.

For example:

```bash
ls -la work/<process_directory>/
```

Common files include:

```text
.command.sh
.command.out
.command.err
.command.log
```

## Restarting an Interrupted Workflow

After correcting the cause of a failure, the workflow can be resubmitted using the same command. Nextflow can reuse successfully completed processes rather than rerunning them when `-resume` is enabled by the workflow launcher.

Before restarting, avoid deleting the `work/` directory or Nextflow execution metadata if you want previously completed tasks to be reused.

## Starting a Completely Fresh Run

If a completely fresh execution is required, previous Nextflow working files and results can be removed before resubmitting the workflow.

From the repository root:

```bash
rm -rf work .nextflow results
rm -f .nextflow.log*
rm -rf logs/*
mkdir -p logs results
```

This removes cached workflow tasks and previous results. Do not use these commands if you intend to resume an interrupted workflow.

## Additional Checks

If the problem is not resolved by the steps above, confirm that:

- Nextflow is installed and available in `PATH`.
- Apptainer or Singularity is installed and available in `PATH`.
- The HPC system uses SLURM or the workflow configuration has been adapted for its scheduler.
- Required shell scripts and Python utilities are executable.
- Paths in `samples.tsv` point to accessible files.
- BAM files have corresponding index files.
- The reference FASTA and its index are accessible.
- The workflow container is accessible from the compute nodes.
- The selected `--queue_size` is appropriate for the HPC scheduler.
- Requested CPU, memory, and walltime resources are compatible with the HPC system.

