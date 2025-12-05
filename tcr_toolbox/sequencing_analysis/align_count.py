import collections
import json
import os
import sys
import re
import shutil
import subprocess
import threading
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from typing import Tuple

import pandas as pd
import pysam
from dotenv import load_dotenv

from tcr_toolbox.sequencing_analysis.constants import trim_adapter_config
from tcr_toolbox.utils.constants import cigar_string_to_pysam_cigar_number, p20_pMX_trim_seqs
from tcr_toolbox.utils.utils import stream_stderr

load_dotenv()
tcr_toolbox_data_path = os.getenv("tcr_toolbox_data_path")


def return_fastq_basename(fastq: str | os.PathLike) -> str:
    fastq = Path(fastq).expanduser()
    if len(fastq.suffixes) > 1:
        return Path(fastq.stem).stem
    else:
        return fastq.stem


def umi_tools_extraction(
    fastq: str | os.PathLike, bc_pattern: str, whitelists_dir: str | os.PathLike, logs_dir: str | os.PathLike, tmp_dir: str | os.PathLike, barcode_file: str | os.PathLike
) -> Path:
    fastq = Path(fastq).expanduser()
    whitelists_dir = Path(whitelists_dir).expanduser()
    logs_dir = Path(logs_dir).expanduser()
    tmp_dir = Path(tmp_dir).expanduser()
    barcode_file = Path(barcode_file).expanduser()

    fastq_basename = return_fastq_basename(fastq)

    # umi_tools barcode whitelist extration here can probably be removed
    # but maybe I can also give users the option to choose between umi_tools extraction
    # and a user supplied barcode whitelist. In addition, keeping umi_tools whitelist extration here
    # might also be relevant for checking what kind of barcodes are present in the seq data.
    print("Started whitelist extraction:", fastq.name, file=sys.stderr)
    whitelist_file = whitelists_dir / f"{fastq_basename}_whitelist.txt"
    log_file = logs_dir / f"{fastq_basename}_whitelist.err"
    with open(whitelist_file, "w") as fout, open(log_file, "w") as ferr:
        subprocess.run(["umi_tools", "whitelist", "--stdin", str(fastq), "--bc-pattern=" + bc_pattern, "--set-cell-number=384", "--log2stderr"], stderr=ferr, stdout=fout)

    print("Started umi_tools extraction:", fastq.name, file=sys.stderr)
    log_file = logs_dir / f"{fastq_basename}_extraction.err"
    extracted_fastq = tmp_dir / f"{fastq_basename}_extracted.fastq.gz"
    with open(log_file, "w") as ferr:
        subprocess.run(
            ["umi_tools", "extract", "--bc-pattern=" + bc_pattern, "--stdin", str(fastq), "--stdout", str(extracted_fastq), "--whitelist=" + str(barcode_file), "--log2stderr"],
            stderr=ferr,
        )
    return extracted_fastq


def cutadapt_trim_reads(
    fastq: str | os.PathLike,
    logs_dir: str | os.PathLike,
    trim_adapter_config: dict,
    p20_pMX_trim_seqs: dict,
    threads: int,
    min_quality_3_end: int = None,
    remove_fastq: bool = False,
) -> Path:
    fastq = Path(fastq).expanduser()
    logs_dir = Path(logs_dir).expanduser()

    fastq_basename = return_fastq_basename(fastq)
    trimmed_fastq = fastq.with_name(f"{fastq_basename}_trimmed.fastq.gz")
    log_file_out = logs_dir / f"{fastq_basename}_cutadapt.out"
    log_file_err = logs_dir / f"{fastq_basename}_cutadapt.err"

    print(f"Started cutadapt trimming: {fastq.name}", file=sys.stderr)
    with log_file_err.open("w") as ferr, log_file_out.open("w") as fout:
        cmd = ["cutadapt"]

        if min_quality_3_end is not None:
            print(f"Applying 3' quality trimming: {min_quality_3_end}", file=sys.stderr)
            cmd.extend(["-q", str(min_quality_3_end)])

        for flag, keys in trim_adapter_config.items():
            for key in keys:
                if flag == "-g":
                    cmd.extend([flag, "^" + p20_pMX_trim_seqs[key]])
                else:
                    cmd.extend([flag, p20_pMX_trim_seqs[key]])

        cmd.extend(["-n", "2", "-j", str(threads)])
        cmd.extend(["-o", str(trimmed_fastq), str(fastq)])

        subprocess.run(cmd, stdout=fout, stderr=ferr, check=True)

    if remove_fastq:
        fastq.unlink()

    return trimmed_fastq


def align_reads(
    trimmed_fastq: str | os.PathLike,
    reference_file: str | os.PathLike,
    bam_outs_dir: str | os.PathLike,
    logs_dir: str | os.PathLike,
    threads: int,
    use_minimap2: bool = False,
    minimap2_kmer_length: int = 19,
) -> Tuple[Path, Path, Path]:
    trimmed_fastq = Path(trimmed_fastq).expanduser()
    reference_file = Path(reference_file).expanduser()
    bam_outs_dir = Path(bam_outs_dir).expanduser()
    logs_dir = Path(logs_dir).expanduser()

    fastq_basename = return_fastq_basename(trimmed_fastq)
    bam_file = bam_outs_dir / f"{fastq_basename}.bam"
    sorted_bam_file = bam_file.with_name(bam_file.stem + "_sorted.bam")
    align_err_file = logs_dir / f"{fastq_basename}_align.err"
    flagstat_out_file = logs_dir / f"{fastq_basename}_samtools_flagstat.out"

    if not use_minimap2:
        print("Started bwa mem alignment:", trimmed_fastq.name, file=sys.stderr)
        align_cmd = ["bwa", "mem", "-t", str(threads), "-M", str(reference_file), str(trimmed_fastq)]
    else:
        print("Started minimap2 alignment:", trimmed_fastq.name, file=sys.stderr)
        align_cmd = ["minimap2", "-ax", "sr", f"-k{minimap2_kmer_length}", f"-w{5}", "--MD", "-t", str(threads), str(reference_file), str(trimmed_fastq)]

    with open(align_err_file, "w") as err_log, open(bam_file, "wb") as bam_out:
        align_proc = subprocess.Popen(align_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

        stderr_thread = threading.Thread(target=stream_stderr, args=(align_proc.stderr, err_log, True))
        stderr_thread.start()

        subprocess.run(["samtools", "view", "-bSh", "-"], stdin=align_proc.stdout, stdout=bam_out, stderr=err_log, check=True)

        align_proc.stdout.close()
        retcode = align_proc.wait()
        stderr_thread.join()

        if retcode != 0:
            raise subprocess.CalledProcessError(retcode, align_cmd)

    with open(align_err_file, "a") as err_log:
        subprocess.run(["samtools", "sort", "-@", str(threads), "-o", str(sorted_bam_file), str(bam_file)], stderr=err_log, check=True)

        subprocess.run(["samtools", "index", "-@", str(threads), str(sorted_bam_file)], stderr=err_log, check=True)

    bam_file.unlink()

    with open(flagstat_out_file, "w") as fout:
        subprocess.run(["samtools", "flagstat", str(sorted_bam_file)], stdout=fout, check=True)

    print(f"Alignment {trimmed_fastq.name} completed.", file=sys.stderr)
    return bam_file, sorted_bam_file, flagstat_out_file


def samtools_filter_secondary_alignments(flagstat_out_file: str | os.PathLike, sorted_bam_file: str | os.PathLike, logs_dir: str | os.PathLike, threads: int) -> Path:
    flagstat_out_file = Path(flagstat_out_file).expanduser()
    sorted_bam_file = Path(sorted_bam_file).expanduser()
    logs_dir = Path(logs_dir).expanduser()

    with open(flagstat_out_file, "r") as flag_in:
        lines = flag_in.readlines()

    log_file = logs_dir / f"{sorted_bam_file.stem}_samtools_filter_secondary.err"

    if not int(lines[2].split(" +")[0]) == 0:
        print(int(lines[2].split(" +")[0]), "secondary alignments detected!", "Removing secondary alignments!", file=sys.stderr)

        second_filtered_bam_file = sorted_bam_file.with_name(sorted_bam_file.stem + "_second_filtered")
        sorted_bam_file_index = sorted_bam_file.with_suffix(".bam.bai")

        with open(log_file, "w") as ferr:
            subprocess.run(["samtools", "view", "-b", "-F 256", "-@", str(threads), "-o", str(second_filtered_bam_file), str(sorted_bam_file)], stderr=ferr, check=True)

            subprocess.run(["samtools", "index", "-@", str(threads), str(second_filtered_bam_file)], stderr=ferr)

        bam_file_path = second_filtered_bam_file

        sorted_bam_file.unlink()
        sorted_bam_file_index.unlink(missing_ok=True)

    else:
        print("filter_secondary_alignments is set to True, but no secondary alignments detected!", file=sys.stderr)
        bam_file_path = sorted_bam_file

    return bam_file_path


def print_cigarstring_md_tag_counts(bam_file_path: str | os.PathLike, write_cigar_json_file: str | os.PathLike = None, write_md_json_file: str | os.PathLike = None):
    bam_file_path = Path(bam_file_path).expanduser()
    write_cigar_json_file = Path(write_cigar_json_file).expanduser()
    write_md_json_file = Path(write_md_json_file).expanduser()

    cigar_string_counter = collections.Counter()
    md_string_counter = collections.Counter()
    with pysam.AlignmentFile(bam_file_path, "rb") as bam_in:
        for entry in bam_in:
            if entry.cigarstring:
                cigar_string_counter[entry.cigarstring] += 1
                md_string_counter[entry.get_tag("MD")] += 1
    print("Top 10 most common cigar strings:", cigar_string_counter.most_common(10), file=sys.stderr)
    print("Top 10 most common MD tags:", md_string_counter.most_common(10), file=sys.stderr)

    if write_cigar_json_file:
        with open(write_cigar_json_file, "w") as json_out:
            json.dump(cigar_string_counter, json_out)

    if write_md_json_file:
        with open(write_md_json_file, "w") as json_out:
            json.dump(md_string_counter, json_out)


def filter_bam_entry(
    entry: pysam.libcalignedsegment.AlignedSegment, max_soft_5_end: int, max_soft_3_end: int, max_insertion: int, max_deletion: int, minimal_overlap: int, max_mismatches: int
) -> bool:
    """Whether a read alignment entry should be filtered from a BAM file based on filter arguments.

    Parameters
    ----------
    entry : pysam.libcalignedsegment.AlignedSegment
    max_soft_5_end : int
        Maximal number of 5' end read bases soft-clipped from the alignment.
    max_soft_3_end : int
        Maximal number of 3' end read bases soft-clipped from the alignment.
    max_insertion : int
        Maximal number of read insertions from the reference.
    max_deletion : int
        Maximal number of read deletions from the reference.
    minimal_overlap : int
        Minimal alignment match length (can be match or mismatch).
    max_mismatches : int
        Maximal number of mismatches within minimal alignment match length.

    Returns
    -------
    True if entry should be filtered from BAM file.

    See Also
    --------
    print_cigarstring_md_tag_counts: Functions that can be used in alignment pipelines to check how BAM files should be
    filtered using this function and whether this filtering was successful.
    """
    max_indel_filter_cigar_dict = {"I": max_insertion, "D": max_deletion}

    if not entry.cigartuples:  # this should also filter unmapped reads
        return True

    cigar_op_counter = collections.defaultdict(int)
    for op, length in entry.cigartuples:
        cigar_op_counter[op] += length

    for cigar_filter_string in max_indel_filter_cigar_dict.keys():
        cigar_filter_number = cigar_string_to_pysam_cigar_number[cigar_filter_string]
        if cigar_op_counter[cigar_filter_number] > max_indel_filter_cigar_dict[cigar_filter_string]:
            return True

    if entry.cigartuples[0][0] == cigar_string_to_pysam_cigar_number["S"]:
        if entry.cigartuples[0][1] > max_soft_5_end:
            return True

    if entry.cigartuples[-1][0] == cigar_string_to_pysam_cigar_number["S"]:
        if entry.cigartuples[-1][1] > max_soft_3_end:
            return True

    if (cigar_op_counter[cigar_string_to_pysam_cigar_number["M"]] < minimal_overlap) or (len(re.findall(r"[a-zA-Z]", entry.get_tag("MD"))) > max_mismatches):
        return True

    return False


def umi_tools_count(bam_file: str | os.PathLike, counts_dir: str | os.PathLike, logs_dir: str | os.PathLike):
    bam_file = Path(bam_file).expanduser()
    counts_dir = Path(counts_dir).expanduser()
    logs_dir = Path(logs_dir).expanduser()

    bam_file_basename = bam_file.stem
    counts_tsv_file = counts_dir / f"{bam_file_basename}_counts.tsv"
    umi_tools_log_file = logs_dir / f"{bam_file_basename}_umi_tools_count.err"

    print("Started umi_tools counting:", bam_file.name, file=sys.stderr)
    with open(umi_tools_log_file, "w") as ferr:
        umi_proc = subprocess.Popen(
            ["umi_tools", "count", "--per-contig", "--log2stderr", "--per-cell", "-I", str(bam_file), "-S", str(counts_tsv_file)],
            stderr=subprocess.PIPE,
            stdout=subprocess.PIPE,
            universal_newlines=True,
            bufsize=1,
        )

        stderr_thread = threading.Thread(target=stream_stderr, args=(umi_proc.stderr, ferr, False))
        stderr_thread.start()
        umi_proc.wait()
        stderr_thread.join()

    if umi_proc.returncode != 0:
        print(f"Error: umi_tools count failed for {bam_file.stem}. Check the log file.", file=sys.stderr)


def count_umi_cell(
    project_dir: str | os.PathLike,
    reference_file: str | os.PathLike,
    barcode_file: str | os.PathLike = Path(tcr_toolbox_data_path)
    / "tcr_toolbox_datasets"
    / "pair_scan_luna_plate_seq"
    / "barcodes"
    / "barcodes_for_counting.tsv",
    bc_pattern: str = "CCCCCCCCNNNNNNN",
    use_minimap2: bool = False,
    filter_secondary_alignments: bool = True,
    remove_tmp_dir: bool = True,
    threads: int = 8,
    bwa_index_algo: str = "is",
    epitope_barcode: bool = True,
    minimal_overlap_epitope: int = 36,
    minimal_overlap_tcr: int = 105,
    max_mismatches_epitope: int = 0,
    max_mismatches_tcr: int = 0,
    max_soft_5_end_epitope: int = 76,
    max_soft_5_end_tcr: int = 0,
    max_soft_3_end_epitope: int = 0,
    max_soft_3_end_tcr: int = 0,
    max_insertion_epitope: int = 0,
    max_insertion_tcr: int = 0,
    max_deletion_epitope: int = 0,
    max_deletion_tcr: int = 0,
    min_quality_3_end: int = None,
    minimap2_kmer_length: int = 19,
    **kwargs,
):
    """
    Count unique TCR and epitope molecules per single cell from single-cell genomics FASTQ files.
    Writes UMI count TSV files per FASTQ to a `counts` subdirectory within the project directory.

    Parameters
    ----------
    project_dir : str | os.PathLike
        Path to the plate-based TCR and epitope sequencing project directory.
    reference_file : str | os.PathLike
        Path to the reference FASTA file.
    barcode_file : str | os.PathLike, optional
        Path to a TSV file containing cell barcodes. Defaults to the barcode file provided in
        tcr_toolbox datasets:
        '/tcr_toolbox_datasets/pair_scan_luna_plate_seq/barcodes/barcodes_for_counting.tsv'.
    bc_pattern : str, default='CCCCCCCCNNNNNNN'
        UMI-tools pattern describing the barcode structure.
    use_minimap2 : bool, default=False
        Whether to use minimap2 for alignment. Necessary for some processor architectures (e.g., M-series Macs).
    filter_secondary_alignments : bool, default=True
        If True, secondary alignments are filtered from the BAM files.
    remove_tmp_dir : bool, default=True
        If True, removes temporary pipeline directories after processing.
    threads : int, default=8
        Number of threads used for alignment and processing.
    bwa_index_algo : str, default='is'
        Algorithm for BWA indexing. Should not be changed without consultation.
    epitope_barcode : bool, default=True
        Whether epitopes are DNA-barcoded after the stop codon.
    minimal_overlap_epitope : int, default=36
        Minimum alignment length (including mismatches) for epitope reference matches.
    minimal_overlap_tcr : int, default=105
        Minimum alignment length (including mismatches) for TCR reference matches.
    max_mismatches_epitope : int, default=0
        Maximum mismatches allowed in epitope alignment.
    max_mismatches_tcr : int, default=0
        Maximum mismatches allowed in TCR alignment.
    max_soft_5_end_epitope : int, default=76
        Maximum number of 5' soft-clipped bases for epitope alignments.
    max_soft_5_end_tcr : int, default=0
        Maximum number of 5' soft-clipped bases for TCR alignments.
    max_soft_3_end_epitope : int, default=0
        Maximum number of 3' soft-clipped bases for epitope alignments.
    max_soft_3_end_tcr : int, default=0
        Maximum number of 3' soft-clipped bases for TCR alignments.
    max_insertion_epitope : int, default=0
        Maximum insertions allowed for epitope alignments.
    max_insertion_tcr : int, default=0
        Maximum insertions allowed for TCR alignments.
    max_deletion_epitope : int, default=0
        Maximum deletions allowed for epitope alignments.
    max_deletion_tcr : int, default=0
        Maximum deletions allowed for TCR alignments.
    min_quality_3_end : int, optional
        Minimum Phred quality for 3' end bases; bases below this are trimmed. Defaults to None (no trimming).
    minimap2_kmer_length : int, default=19
        K-mer length used by minimap2 if `use_minimap2` is True.
    **kwargs : dict
        Additional keyword arguments passed to underlying functions.

    See Also
    --------
    generate_assembly_nt_refs : Function in `tcr_toolbox.tcr_assembly.sequencing_analysis.reference` that can be used
    to generate the reference FASTA file.

    Examples
    --------
    >>> project_dir = '/path/to/project'
    >>> init_umi_count_analysis_dir(project_dir)
    >>> project_dir = '/path/to/your/project'
    >>> count_umi_cell(
    ...     project_dir=project_dir,
    ...     reference_file='/path/to/your/reference.fa',
    ...     bc_pattern='CCCCCCCCNNNNNNN',
    ...     use_minimap2=False,
    ...     remove_tmp_dir=True,
    ...     filter_secondary_alignments=True,
    ...     threads=16,
    ...     bwa_index_algo='is',
    ...     epitope_barcode=True,
    ...     minimal_overlap_epitope=36,
    ...     minimal_overlap_tcr=103,  # set to 105 with bwa-mem; minimap2 may not yield full 105M
    ...     max_mismatches_epitope=0,
    ...     max_mismatches_tcr=0,
    ...     max_soft_5_end_epitope=76,
    ...     max_soft_5_end_tcr=0,
    ...     max_soft_3_end_epitope=0,
    ...     max_soft_3_end_tcr=0,
    ...     max_insertion_epitope=0,
    ...     max_insertion_tcr=0,
    ...     max_deletion_epitope=0,
    ...     max_deletion_tcr=0,
    ...     min_quality_3_end=31
    ... )
        Started cutadapt trimming: 8191_10_15_GGTACCTT-GACGTCTT_S10_R1_001_extracted.fastq.gz
        Applying 3' quality trimming: 31
        Started bwa mem alignment: 8191_10_15_GGTACCTT-GACGTCTT_S10_R1_001_extracted_trimmed.fastq.gz
        Alignment 8191_10_15_GGTACCTT-GACGTCTT_S10_R1_001_extracted_trimmed.fastq.gz completed.
        3890 secondary alignments detected! Removing secondary alignments!
        Top 10 most common cigar strings: [('76S36M', 1017933), ('75S36M', 914122), ('74S36M', 180539), ('72S36M', 163743), ('104M', 111075), ('105M', 103442), ('71S36M', 99430), ('73S36M', 83782), ('68S36M', 59730), ('70S36M', 56200)]
        Top 10 most common MD tags: [('36', 4052143), ('104', 92746), ('105', 80194), ('34', 48366), ('33', 35205), ('101', 22897), ('103', 22836), ('32', 19193), ('9G26', 16646), ('35', 16246)]
        Filtering BAM file...
        Top 10 most common cigar strings: [('76S36M', 988783), ('75S36M', 893524), ('74S36M', 177532), ('72S36M', 161477), ('71S36M', 97667), ('104M', 92728), ('73S36M', 82096), ('105M', 79707), ('68S36M', 58351), ('70S36M', 55354)]
        Top 10 most common MD tags: [('36', 4042881), ('104', 92728), ('105', 79707), ('103', 22825)]


        UMI counting bam files...
        Started umi_tools counting: 8191_10_15_GGTACCTT-GACGTCTT_S10_R1_001_extracted_trimmed_sorted_second_filtered_cigarmd_filtered.bam
        UMI counting done!
    """

    if use_minimap2 and minimal_overlap_epitope < minimap2_kmer_length:
        raise Exception(
            "An epitope/epitope barcode length shorter than current minimap2 kmer length will not work:",
            minimap2_kmer_length,
            "! If epitope/epitope barcode length is really small, also consider changing the w parameter of minimap2!",
        )
    project_dir = Path(project_dir).expanduser()
    reference_path = Path(reference_file).expanduser()
    reference_dir = reference_path.parent
    barcode_file = Path(barcode_file).expanduser()

    fastq_list = list(project_dir.glob("*.fastq.gz"))
    if not fastq_list:
        fastq_list = list(project_dir.glob("*.fastq"))
    if not fastq_list:
        raise Exception("No fastq files found in work directory!")

    if not reference_path.is_file():
        raise Exception("Reference fasta file is not found!", reference_file)

    tmp_dir = project_dir / "tmp"
    logs_dir = tmp_dir / "logs"
    whitelists_dir = tmp_dir / "whitelists"
    bam_outs_dir = tmp_dir / "bam_outs"
    counts_dir = project_dir / "counts"

    for directory in [tmp_dir, logs_dir, whitelists_dir, bam_outs_dir, counts_dir]:
        directory.mkdir(parents=True, exist_ok=True)

    if not use_minimap2:
        ref_tmp_dir = reference_dir / "tmp"
        ref_logs_dir = ref_tmp_dir / "logs"

        # Create reference index if the index does not exist yet:
        sa_file = reference_path.with_name(reference_path.name + ".sa")
        log_file = ref_logs_dir / "bwa-mem-index.err"

        if not sa_file.is_file():
            print("Indexing reference...", file=sys.stderr)
            log_file.parent.mkdir(parents=True, exist_ok=True)  # Ensure log directory exists
            with open(log_file, "w") as ferr:
                subprocess.run(["bwa", "index", "-a", bwa_index_algo, str(reference_path)], stderr=ferr)

    ## 1 umi-tools extraction
    with ProcessPoolExecutor(max_workers=threads) as executor:
        umi_tools_extract_futures = [executor.submit(umi_tools_extraction, fastq, bc_pattern, whitelists_dir, logs_dir, tmp_dir, barcode_file) for fastq in fastq_list]
    extracted_fastq_list = [future.result() for future in umi_tools_extract_futures]

    print("\n\n", file=sys.stderr)
    for extracted_fastq in extracted_fastq_list:
        ## 2 cutadapt
        if not epitope_barcode:
            print("epitope_barcode is set to False, trimming epitope seqs with stop codon...", file=sys.stderr)
            trimmed_fastq = cutadapt_trim_reads(
                fastq=extracted_fastq,
                logs_dir=logs_dir,
                trim_adapter_config=trim_adapter_config["no_epi_barcode_rt_pair_scan_plate"],
                p20_pMX_trim_seqs=p20_pMX_trim_seqs,
                threads=threads,
                min_quality_3_end=min_quality_3_end,
                remove_fastq=True,
            )
        else:
            trimmed_fastq = cutadapt_trim_reads(
                fastq=extracted_fastq,
                logs_dir=logs_dir,
                trim_adapter_config=trim_adapter_config["epi_barcode_rt_pair_scan_plate"],
                p20_pMX_trim_seqs=p20_pMX_trim_seqs,
                threads=threads,
                min_quality_3_end=min_quality_3_end,
                remove_fastq=True,
            )

        ## 3 alignment
        _, sorted_bam_file, flagstat_out_file = align_reads(
            trimmed_fastq=trimmed_fastq,
            reference_file=reference_file,
            bam_outs_dir=bam_outs_dir,
            logs_dir=logs_dir,
            threads=threads,
            use_minimap2=use_minimap2,
            minimap2_kmer_length=minimap2_kmer_length,
        )

        ## 4 filter alignments
        if filter_secondary_alignments:
            bam_file_path = samtools_filter_secondary_alignments(flagstat_out_file=flagstat_out_file, sorted_bam_file=sorted_bam_file, logs_dir=logs_dir, threads=threads)
        else:
            bam_file_path = sorted_bam_file

        print_cigarstring_md_tag_counts(
            bam_file_path=bam_file_path,
            write_cigar_json_file=logs_dir / f"{bam_file_path.stem}_cigar_counter_before_filtering.json",
            write_md_json_file=logs_dir / f"{bam_file_path.stem}_MD_counter_before_filtering.json",
        )
        print("Filtering BAM file...", file=sys.stderr)
        bam_file_count_path = bam_file_path.with_name(bam_file_path.stem + "_cigarmd_filtered.bam")
        with pysam.AlignmentFile(bam_file_path, "rb") as bam_in, pysam.AlignmentFile(bam_file_count_path, "wb", template=bam_in) as bam_out:
            for entry in bam_in:
                if not entry.reference_name:
                    continue

                current_ref_epitope_or_tcr = entry.reference_name[:3]
                if current_ref_epitope_or_tcr not in ["epi", "tcr"]:
                    raise Exception(
                        'Incorrect reference file supplied!, Reference fasta file headers need to start with either "tcr" or "epi"!',
                        "Did you use the generate_assembly_nt_refs to generate the reference file?",
                    )

                if current_ref_epitope_or_tcr == "epi":
                    if not filter_bam_entry(
                        entry,
                        max_soft_5_end=max_soft_5_end_epitope,
                        max_soft_3_end=max_soft_3_end_epitope,
                        max_insertion=max_insertion_epitope,
                        max_deletion=max_deletion_epitope,
                        minimal_overlap=minimal_overlap_epitope,
                        max_mismatches=max_mismatches_epitope,
                    ):
                        bam_out.write(entry)

                if current_ref_epitope_or_tcr == "tcr":
                    if not filter_bam_entry(
                        entry,
                        max_soft_5_end=max_soft_5_end_tcr,
                        max_soft_3_end=max_soft_3_end_tcr,
                        max_insertion=max_insertion_tcr,
                        max_deletion=max_deletion_tcr,
                        minimal_overlap=minimal_overlap_tcr,
                        max_mismatches=max_mismatches_tcr,
                    ):
                        bam_out.write(entry)

        bam_file_path.unlink()
        bam_file_path.with_suffix(bam_file_path.suffix + ".bai").unlink()

        log_file = logs_dir / f"{bam_file_count_path.stem}_samtools_index_cigarmd_filter.err"
        with open(log_file, "w") as ferr:
            subprocess.run(["samtools", "index", "-@", str(threads), str(bam_file_count_path)], stderr=ferr)

        print_cigarstring_md_tag_counts(
            bam_file_path=bam_file_count_path,
            write_cigar_json_file=logs_dir / f"{bam_file_count_path.stem}_cigar_counter_after_filtering.json",
            write_md_json_file=logs_dir / f"{bam_file_count_path.stem}_MD_counter_after_filtering.json",
        )
        print("\n", file=sys.stderr)

    ## 5 umi counting
    print("\nUMI counting bam files...", file=sys.stderr)
    bam_files = list(bam_outs_dir.glob("*_cigarmd_filtered.bam"))
    with ProcessPoolExecutor(max_workers=threads) as executor:
        umi_count_futures = [executor.submit(umi_tools_count, bam_file, counts_dir, logs_dir) for bam_file in bam_files]
        for umi_count_future in umi_count_futures:
            umi_count_future.result()

    # remove temporary directory
    if remove_tmp_dir:
        shutil.move(str(logs_dir), str(Path(project_dir) / "logs"))
        shutil.rmtree(str(tmp_dir))

    print("UMI counting done!", file=sys.stderr)


def count_aligned_gdna(bam_file: str | os.PathLike, counts_dir: str | os.PathLike):
    bam_file = Path(bam_file)
    counts_dir = Path(counts_dir)

    alignment_counter = collections.Counter()

    bam_basename = bam_file.stem
    print("Counting:", bam_file.name, file=sys.stderr)
    with pysam.AlignmentFile(bam_file, "r") as bam_in:
        for entry in bam_in:
            alignment_counter[entry.reference_name] += 1

    pd.Series(alignment_counter).sort_values(ascending=False).to_frame(name="read_count").rename_axis("reference_name").to_csv(counts_dir / f"{bam_basename}_counts.csv")


def count_reads_bulk(
    project_dir: str | os.PathLike,
    reference_file: str | os.PathLike,
    minimal_overlap: int,
    max_mismatches: int,
    max_soft_5_end: int,
    max_soft_3_end: int,
    use_minimap2: bool = False,
    filter_secondary_alignments: bool = True,
    remove_tmp_dir: bool = True,
    threads: int = 8,
    bwa_index_algo: str = "is",
    epitope_barcode: bool = True,
    max_insertion: int = 0,
    max_deletion: int = 0,
    min_quality_3_end: int = None,
    minimap2_kmer_length: int = 19,
    **kwargs,
):
    """
    Count TCR and barcoded-epitope reads in bulk sequencing FASTQ files.

    This function performs read trimming, alignment, BAM filtering, and counting of aligned reads.
    Results are written to a `counts` subdirectory within the project directory. Temporary files
    are stored in a `tmp` folder and removed if `remove_tmp_dir` is True.

    Parameters
    ----------
    project_dir : str | os.PathLike
        Path to project directory containing the FASTQ files.
    reference_file : str | os.PathLike
        Path to the reference FASTA file.
    minimal_overlap : int
        Minimum alignment match length (can include mismatches).
    max_mismatches : int
        Maximum mismatches allowed within minimal alignment match length.
    max_soft_5_end : int
        Maximum number of 5' end read bases soft-clipped from the alignment.
    max_soft_3_end : int
        Maximum number of 3' end read bases soft-clipped from the alignment.
    use_minimap2 : bool, default=False
        Whether to use minimap2 for alignment. Necessary for some processor architectures (e.g., M-series Macs).
    filter_secondary_alignments : bool, default=True
        If True, secondary alignments are filtered from BAM files.
    remove_tmp_dir : bool, default=True
        If True, removes temporary pipeline directories after processing.
    threads : int, default=8
        Number of threads used for trimming, alignment, and counting.
    bwa_index_algo : str, default='is'
        Algorithm for BWA indexing. Should not be changed without consultation.
    epitope_barcode : bool, default=True
        Whether epitopes are DNA-barcoded after the stop codon.
    max_insertion : int, default=0
        Maximum insertions allowed for alignments.
    max_deletion : int, default=0
        Maximum deletions allowed for alignments.
    min_quality_3_end : int, optional
        Minimum Phred quality for 3' end bases; bases below this are trimmed. Defaults to None (no trimming).
    minimap2_kmer_length : int, default=19
        K-mer length used by minimap2 if `use_minimap2` is True.
    **kwargs : dict
        Additional keyword arguments passed to underlying functions.

    See Also
    --------
    generate_assembly_nt_refs : Function in `tcr_toolbox.tcr_assembly.sequencing_analysis.reference` that can be used
    to generate the reference FASTA.

    Examples
    --------
    # Epitope
    >>> bulk_project_dir = Path('/path/to/your/bulk/epi/project')
    >>> count_reads_bulk(
    ...     project_dir=bulk_project_dir,
    ...     reference_file='/path/to/your/reference_epi.fa',
    ...     minimal_overlap=36,
    ...     max_mismatches=0,
    ...     max_soft_5_end=78,
    ...     max_soft_3_end=0,
    ...     use_minimap2=False,
    ...     filter_secondary_alignments=True,
    ...     remove_tmp_dir=True,
    ...     threads=16,
    ...     bwa_index_algo='is',
    ...     epitope_barcode=True,
    ...     max_insertion=0,
    ...     max_deletion=0,
    ...     min_quality_3_end=31
    ... )
    Indexing reference...
    Started cutadapt trimming: 8191_33_B-YWE-2_TGGATCGA-TATCGCAC_S33_R1_001.fastq.gz
    Applying 3' quality trimming: 31 Started bwa mem alignment: 8191_33_B-YWE-2_TGGATCGA-TATCGCAC_S33_R1_001_trimmed.fastq.gz
    Alignment 8191_33_B-YWE-2_TGGATCGA-TATCGCAC_S33_R1_001_trimmed.fastq.gz completed.
    1 secondary alignments detected! Removing secondary alignments!
    Top 10 most common cigar strings: [('78S36M', 672164), ('80S36M', 39352), ('84S36M', 19364), ('77S36M', 15368), ('29S36M', 14358), ('32S36M', 14071), ('31S36M', 13947), ('79S36M', 13524), ('76S36M', 12767), ('34S36M', 12654)]
    Top 10 most common MD tags: [('36', 1337393), ('32', 9211), ('33', 5548), ('34', 5423), ('35', 5084), ('31', 4845), ('5G30', 3703), ('30', 3121), ('8C27', 2433), ('31A4', 2183)]
    Filtering BAM file...
    Top 10 most common cigar strings: [('78S36M', 649697), ('77S36M', 15002), ('29S36M', 13787), ('32S36M', 13577), ('31S36M', 13392), ('76S36M', 12456), ('34S36M', 12144), ('35S36M', 11864), ('27S36M', 11446), ('24S36M', 11059)]
    Top 10 most common MD tags: [('36', 1264645)]

    Counting bam files...
    Counting: 8191_33_B-YWE-2_TGGATCGA-TATCGCAC_S33_R1_001_trimmed_sorted_cigarmd_filtered.bam
    gdna reads counting done!


    # TCR
    >>> bulk_project_dir = Path('/path/to/your/bulk/tcr/project')
    >>> count_reads_bulk(
    ...     project_dir=bulk_project_dir,
    ...     reference_file='/path/to/your/reference_beta.fa',
    ...     minimal_overlap=110,
    ...     max_mismatches=0,
    ...     max_soft_5_end=3,
    ...     max_soft_3_end=0,
    ...     use_minimap2=False,
    ...     filter_secondary_alignments=True,
    ...     remove_tmp_dir=True,
    ...     threads=16,
    ...     bwa_index_algo='is',
    ...     epitope_barcode=True,
    ...     max_insertion=0,
    ...     max_deletion=0,
    ...     min_quality_3_end=31
    ... )
    Indexing reference...
    Started cutadapt trimming: 8191_37_T-YWE-2_GATTCTGC-CTCTCGTC_S37_R1_001.fastq.gz
    Applying 3' quality trimming: 31
    Started bwa mem alignment: 8191_37_T-YWE-2_GATTCTGC-CTCTCGTC_S37_R1_001_trimmed.fastq.gz
    Alignment 8191_37_T-YWE-2_GATTCTGC-CTCTCGTC_S37_R1_001_trimmed.fastq.gz completed.
    5731 secondary alignments detected! Removing secondary alignments!
    Top 10 most common cigar strings: [('2S110M', 106142), ('3S110M', 48031), ('108M', 14896), ('1S110M', 14352), ('110M', 11741), ('109M', 9654), ('105M', 7900), ('107M', 7801), ('106M', 6523), ('63M', 5667)]
    Top 10 most common MD tags: [('110', 153853), ('108', 13381), ('109', 8991), ('107', 7781), ('105', 7614), ('106', 6559), ('63', 5551), ('72', 4982), ('69', 4904), ('66', 4790)]
    Filtering BAM file...
    Top 10 most common cigar strings: [('2S110M', 91632), ('3S110M', 36261), ('1S110M', 12406), ('110M', 10275)]
    Top 10 most common MD tags: [('110', 150574)]

    Counting bam files...
    Counting: 8191_37_T-YWE-2_GATTCTGC-CTCTCGTC_S37_R1_001_trimmed_sorted_second_filtered_cigarmd_filtered.bam
    gdna reads counting done!
    """
    project_dir = Path(project_dir).expanduser()
    reference_path = Path(reference_file).expanduser()
    reference_dir = reference_path.parent

    fastq_list = list(project_dir.glob("*.fastq.gz"))
    if not fastq_list:
        fastq_list = list(project_dir.glob("*.fastq"))
    if not fastq_list:
        fastq_list = list(project_dir.glob("*.fq.gz"))
    if not fastq_list:
        fastq_list = list(project_dir.glob("*.fq"))
    if not fastq_list:
        raise Exception("No fastq files found in work directory!", project_dir)

    if not reference_path.is_file():
        raise Exception("Reference fasta file is not found!", reference_file)

    tmp_dir = project_dir / "tmp"
    logs_dir = tmp_dir / "logs"
    bam_outs_dir = tmp_dir / "bam_file_outs"
    counts_dir = project_dir / "counts"

    for directory in [tmp_dir, logs_dir, bam_outs_dir, counts_dir]:
        directory.mkdir(parents=True, exist_ok=True)

    if not use_minimap2:
        ref_tmp_dir = reference_dir / "tmp"
        ref_logs_dir = ref_tmp_dir / "logs"

        # Create reference index if the index does not exist yet:
        sa_file = reference_path.with_name(reference_path.name + ".sa")
        log_file = ref_logs_dir / "bwa-mem-index.err"

        if not sa_file.is_file():
            print("Indexing reference...", file=sys.stderr)
            log_file.parent.mkdir(parents=True, exist_ok=True)  # Ensure log directory exists
            with open(log_file, "w") as ferr:
                subprocess.run(["bwa", "index", "-a", bwa_index_algo, str(reference_path)], stderr=ferr)

    for fastq in fastq_list:
        ## 1 cutadapt
        if not epitope_barcode:
            trimmed_fastq = cutadapt_trim_reads(
                fastq=fastq,
                logs_dir=logs_dir,
                trim_adapter_config=trim_adapter_config["no_epi_barcode_pcr_gdna"],
                p20_pMX_trim_seqs=p20_pMX_trim_seqs,
                threads=threads,
                min_quality_3_end=min_quality_3_end,
                remove_fastq=False,
            )

        else:
            trimmed_fastq = cutadapt_trim_reads(
                fastq=fastq,
                logs_dir=logs_dir,
                trim_adapter_config=trim_adapter_config["epi_barcode_pcr_gdna"],
                p20_pMX_trim_seqs=p20_pMX_trim_seqs,
                threads=threads,
                min_quality_3_end=min_quality_3_end,
                remove_fastq=False,
            )

        ## 2 alignment
        _, sorted_bam_file, flagstat_out_file = align_reads(
            trimmed_fastq=trimmed_fastq,
            reference_file=reference_file,
            bam_outs_dir=bam_outs_dir,
            logs_dir=logs_dir,
            threads=threads,
            use_minimap2=use_minimap2,
            minimap2_kmer_length=minimap2_kmer_length,
        )

        ## 3 filter alignments
        if filter_secondary_alignments:
            bam_file_path = samtools_filter_secondary_alignments(flagstat_out_file=flagstat_out_file, sorted_bam_file=sorted_bam_file, logs_dir=logs_dir, threads=threads)
        else:
            bam_file_path = sorted_bam_file

        print_cigarstring_md_tag_counts(
            bam_file_path=bam_file_path,
            write_cigar_json_file=logs_dir / f"{bam_file_path.stem}_cigar_counter_before_filtering.json",
            write_md_json_file=logs_dir / f"{bam_file_path.stem}_MD_counter_before_filtering.json",
        )
        print("Filtering BAM file...", file=sys.stderr)
        bam_file_count_path = bam_file_path.with_name(bam_file_path.stem + "_cigarmd_filtered.bam")
        with pysam.AlignmentFile(bam_file_path, "rb") as bam_in, pysam.AlignmentFile(bam_file_count_path, "wb", template=bam_in) as bam_out:
            for entry in bam_in:
                if not filter_bam_entry(
                    entry,
                    max_soft_5_end=max_soft_5_end,
                    max_soft_3_end=max_soft_3_end,
                    max_insertion=max_insertion,
                    max_deletion=max_deletion,
                    minimal_overlap=minimal_overlap,
                    max_mismatches=max_mismatches,
                ):
                    bam_out.write(entry)

        bam_file_path.unlink()
        bam_file_path.with_suffix(bam_file_path.suffix + ".bai").unlink()

        log_file = logs_dir / f"{bam_file_count_path.stem}_samtools_index_cigarmd_filter.err"
        with open(log_file, "w") as ferr:
            subprocess.run(["samtools", "index", "-@", str(threads), str(bam_file_count_path)], stderr=ferr, check=True)

        print_cigarstring_md_tag_counts(
            bam_file_path=bam_file_count_path,
            write_cigar_json_file=logs_dir / f"{bam_file_count_path.stem}_cigar_counter_after_filtering.json",
            write_md_json_file=logs_dir / f"{bam_file_count_path.stem}_MD_counter_after_filtering.json",
        )
        print("\n", file=sys.stderr)

    ## 4 count alignments
    print("\nCounting bam files...", file=sys.stderr)
    bam_files = list(bam_outs_dir.glob("*_cigarmd_filtered.bam"))
    with ProcessPoolExecutor(max_workers=threads) as executor:
        count_futures = [executor.submit(count_aligned_gdna, bam_file, counts_dir) for bam_file in bam_files]
        for count_future in count_futures:
            count_future.result()

    # remove temporary directory
    if remove_tmp_dir:
        shutil.move(str(logs_dir), str(Path(project_dir) / "logs"))
        shutil.rmtree(str(tmp_dir))

    print("gdna reads counting done!", file=sys.stderr)
