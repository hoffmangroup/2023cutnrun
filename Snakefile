import math
import os
import json
from helper_functions import *
import sys

configfile: "config.json"

print(sys.executable)   # for debugging
print("\n\n")

# below code is to remove the python path from mordor when executing, 
# when snakemake runs, it will add in mordor python path despite 
# being in virtual env, which then uses mordor's numpy later on 
# which crashes MACS and Picard. 
print("before removal:")
print("\n".join(sys.path))

i = 0
while i != len(sys.path):
    print(f"PATH: {sys.path[i]}")
    if "multiqc" in sys.path[i]:
        print(f"removing: {sys.path[i]}")
        sys.path.remove(sys.path[i])
    else:
        i += 1

print("after removal:")
print("\n".join(sys.path))
print("\n\n")


# printing out the configuration for easier reproducibility
print("config: ")
for x in config:
    # putting "False" "True" strings in the config file will provide python boolean variables
    print(f"\t {x} : {str(config[x])} : {str(type(config[x]))}")


# Paths
RAW_DIR = config["raw_dir"]
BAM_DIR = config["bam_dir"]
QC_DIR = config["qc_dir"]
MULTIQC_DIR = config["multiqc_dir"]
FASTQC_DIR = config["fastqc_dir"]
PRE_FASTP_FASTQC_DIR = os.path.join(FASTQC_DIR + "/pre_fastp")
POST_FASTP_FASTQC_DIR = os.path.join(FASTQC_DIR + "/post_fastp")
OUTPUT_DIR = config["output_dir"]
MACS_DIR = config["macs_dir"]
TEMP_DIR = config["temp_dir"]
FASTP_OUT_DIR = config["fastp_out_dir"]
LOG_DIR = config["log_dir"]
BED_DIR = config["bed_dir"]
BDG_DIR = config["bdg_dir"]
SEACR_DIR = config["seacr_dir"]

THREADS = max(config["threads"] // 2 - 2, 1)

# preprocessing
FASTP_PARAM = config["fastp_param"]
FASTQ_END = ".fastq.gz"
FQ_END = ".fastp.fastq.gz"
STRANDS = ["_R1", "_R2"]

# alignment
BOWTIE2_CORES = max(config["threads"] // 2 - 2, 1)
BT2_EXT = [".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2", ".rev.1.bt2", ".rev.2.bt2"] # note: ".1" ".2" etc might change
MAIN_BOWTIE2_INDEX = config["bowtie2_index_main"]
MAIN_BOWTIE2_PARAM = config["bowtie2_main_param"]
SPIKE_IN_BOWTIE2_INDEX = config["bowtie2_index_spike_in"]
SPIKE_IN_BOWTIE2_PARAM = config["bowtie2_spike_in_param"]

# fragment length filtering
if config["filter_option"] == "min":
    FILTER_OPTION = "--minFragmentLength"
else: # by default, "max"
    FILTER_OPTION = "--maxFragmentLength"
MAX_FRAG_LEN = config["max_frag_len"]
MIN_FRAG_LEN = config["min_frag_len"]

# spike-in calibration
CHR_LEN_FILE = config["chr_len_file"]
MAIN_GEMOME = config["genome_main"]
MAIN_GEMOME_SIZE = config["genome_size"]

ARBITRARY_SCALE = config["scale"]
SPIKE_IN_MIN_LEN = config["spike_in_min_len"]
SPIKE_IN_MAX_LEN = config["spike_in_max_len"]
SPIKE_IN_SCRIPT = config["spike_in_script"]

# peak calling
PEAK_CALLER = config["peak_caller"]
PEAK_LOCATION = config["peak_location"]

# the macs version to use (macs2, macs3)
MACS = config["macs_version"]
MIN_PEAK_LEN = config["min_peak_length"]
MACS_FDR = config["macs_fdr"]
USE_MACS_CALLPEAK = config["macs_callpeak"]

SEACR_MODE = config["seacr_mode"]
SEACR_SCRIPT = config["seacr_script"]

# motif analysis
FASTA_EXT_LEN = config["fasta_ext_len"]
JASPAR = config["jaspar"]
SKIP_DREME = config["skip_dreme"]

# QC
QUALIMAP_ASSEMBLY = config["qualimap_assembly"]
MULTIQC_CONFIG = config["multiqc_config"]
ONLY_MULTIQC = config["manual"]

# condition
FILTER_FRAG_LEN = config["filter"]
SPIKE_IN_CALIBRATE = config["spkin"]
NEED_FLAGSTAT = config["flagstat"]
NEED_FASTQC = config["fastqc"]
NEED_QC = config["qc"]
NEED_MOTIF = config["motif"]
NEED_CUTOFF_ANALYSIS = config["cutoff"]

treatment_name = config["treatment"]
control_name = config["control"]

# these are separate functions to better generate the input file name strings
def get_flagstat_and_bam_index_input(wildcards):
    suffix = ["main"]
    if FILTER_FRAG_LEN:
        suffix.append("main_filtered")
    if SPIKE_IN_CALIBRATE:
        suffix.append("spkin")
    ext = ["bam.bai"]
    if NEED_FLAGSTAT:
        ext.append("flagstat")
    return expand("{dir}/{name}_{suffix}.{ext}", dir = BAM_DIR, name = [treatment_name, control_name], suffix = suffix, ext = ext)


def get_fastqc_input(wildcards):
    if not NEED_FASTQC:
        return []
    pre_fastp = expand("{dir}/{name}{st}{ext}", dir = PRE_FASTP_FASTQC_DIR, name = [treatment_name, control_name], 
        st = STRANDS, ext = FASTQ_END.replace(".fastq.gz", "_fastqc.html"))
    post_fastp = expand("{dir}/{name}{st}{ext}", dir = POST_FASTP_FASTQC_DIR, name = [treatment_name, control_name], 
        st = STRANDS, ext = FQ_END.replace(".fastq.gz", "_fastqc.html"))
    return pre_fastp + post_fastp


def create_symlink(wildcards, target, dst):
    target = target.replace("{name}", wildcards.name) 
    dst = dst.replace("{name}", wildcards.name)
    if os.path.exists(dst):
        os.remove(dst)
    os.symlink(target, dst)


def get_qc_input(wildcards):
    if not NEED_QC:
        return []
    if NEED_FASTQC:
        fastqc = expand("{dir}/{name}{st}.fastp_fastqc.zip", dir = MULTIQC_DIR, name = [treatment_name, control_name], st = STRANDS)
    else:
        fastqc = []
    qualimap = expand("{dir}/{name}_qualimap_main_stats", dir = MULTIQC_DIR, name = [treatment_name, control_name])    # from rule qualimap
    picard = expand("{dir}/{name}_main_picard_insert_size_metrics.txt", dir = MULTIQC_DIR, name = [treatment_name, control_name])  # from rule picard
    if FILTER_FRAG_LEN:
        picard_filtered = expand("{dir}/{name}_main_filtered_picard_insert_size_metrics.txt", dir = MULTIQC_DIR, name = [treatment_name, control_name])
    else:
        picard_filtered = []
    return fastqc + qualimap + picard + picard_filtered


def get_motif_input(wildcards):
    if not NEED_MOTIF:
        return []
    peak_caller = PEAK_CALLER
    peak_location = PEAK_LOCATION
    if peak_caller not in ["MACS2", "SEACR", "all"]:
        peak_caller = "MACS2"
    elif peak_caller == "all":
        peak_caller = ["SEACR", "MACS2"]
    if peak_location not in ["center", "summit", "all"]:
        peak_caller = "summit"
    elif peak_location == "all":
        peak_location = ["center", "summit"]
    return expand("{dir}/{caller}_{t}_{c}_motifs_{location}", dir = OUTPUT_DIR, caller = peak_caller, t = treatment_name, c = control_name, location = peak_location) 
        

if ONLY_MULTIQC:
    rule all:
        input:
            f"{OUTPUT_DIR}/{treatment_name}_{control_name}_multiqc_report.html"
else:
    rule all:
        input:
            get_fastqc_input, 
            get_flagstat_and_bam_index_input,
            get_qc_input,
            get_motif_input,
            f"{MACS_DIR}/{treatment_name}_{control_name}_cutoff_analysis.txt" if NEED_CUTOFF_ANALYSIS else []

# ---------------------------------------- before peak calling ----------------------------------------------------------------------

rule preprocessing: 
    input:
         r12 = expand("{dir}/{{name}}{st}{end}", dir = RAW_DIR, st = STRANDS, end = FASTQ_END)
    output:
        expand("{dir}/{{name}}{st}{end}", dir = FASTP_OUT_DIR, st = STRANDS, end = FQ_END),
        f"{FASTP_OUT_DIR}/{{name}}.json"
    log:
        os.path.join(LOG_DIR, "{name}_preprocessing.log")
    params:
        report = expand("{tag} {dir}/{{name}}.{end}", zip, tag = ["-h", "-j"], dir = [FASTP_OUT_DIR] * 2, end = ["html", "json"])
    threads: THREADS
    shell:
        "fastp {FASTP_PARAM} --thread {threads} -i {input.r12[0]} -I {input.r12[1]} -o {output[0]} -O {output[1]} {params.report} |& tee {log}"

rule get_insert_size:
    input:
        expand("{dir}/{{name}}.json", dir = FASTP_OUT_DIR)
    output:
        f"{FASTP_OUT_DIR}/{{name}}_insert_size.txt"
    run:  # might be easier with bash 'awk'
        averages = ""
        for file in input:
            with open(file, "r") as f:
                data = json.load(f)
                after_filter = data["summary"]["after_filtering"]
                average12 = int((after_filter["read1_mean_length"] + after_filter["read2_mean_length"]) / 2)
                averages += (str(average12) + "\n")
        with open(output[0], "w") as out:
            out.write(str(averages))

rule get_alignment_mode:
    input:
        insert_size = f"{FASTP_OUT_DIR}/{treatment_name}_insert_size.txt"
    output:
        f"{BAM_DIR}/mode.txt"
    run:
        with open(input[0], "r") as f:
            length = int(f.read())
        if length >= 80:
            bowtie2_mode = "--local"
        else:
            bowtie2_mode = "--end-to-end"
        with open(output[0], "w") as out:
            out.write(bowtie2_mode)

rule alignment: # requires module load bowtie2/2.4.4
    input:
        r12 = expand("{dir}/{{name}}{st}{end}", dir = FASTP_OUT_DIR, st = STRANDS, end = FQ_END),
        genome = expand("{index}{end}", index = MAIN_BOWTIE2_INDEX, end = BT2_EXT),
        mode = f"{BAM_DIR}/mode.txt"
    output: 
        f"{BAM_DIR}/{{name}}_main.bam"
    log:
        os.path.join(LOG_DIR, "{name}_main_alignment.log")
    threads: BOWTIE2_CORES
    shell: 
        """
        (bowtie2 -x {MAIN_BOWTIE2_INDEX} -1 {input.r12[0]} -2 {input.r12[1]} {MAIN_BOWTIE2_PARAM} $(cat {input.mode}) -p {threads} | 
        sambamba view --nthreads {threads} -S -f bam -l 0 /dev/stdin | 
        sambamba sort --nthreads {threads} --tmpdir={TEMP_DIR} -l 9 -o {output} /dev/stdin) |& tee {log}
        ln -sf {log} {MULTIQC_DIR}/$(basename {log})
        """

rule alignment_spike_in: # requires module load bowtie2/2.4.4
    input:
        r12 = expand("{dir}/{{name}}{st}{end}", dir = FASTP_OUT_DIR, st = STRANDS, end = FQ_END),
        genome = expand("{index}{end}", index = SPIKE_IN_BOWTIE2_INDEX, end = BT2_EXT),
        mode = f"{BAM_DIR}/mode.txt"
    output: 
        f"{BAM_DIR}/{{name}}_spkin.bam"
    log:
        os.path.join(LOG_DIR, "{name}_spkin_alignment.log")
    threads: BOWTIE2_CORES
    shell: 
        """
        (bowtie2 -x {SPIKE_IN_BOWTIE2_INDEX} -1 {input.r12[0]} -2 {input.r12[1]} {SPIKE_IN_BOWTIE2_PARAM} $(cat {input.mode}) -p {threads} | 
        sambamba view --nthreads {threads} -S -f bam -l 0 /dev/stdin | 
        sambamba sort --nthreads {threads} --tmpdir={TEMP_DIR} -l 9 -o {output} /dev/stdin) |& tee {log}
        ln -sf {log} {MULTIQC_DIR}/$(basename {log})
        """
     
rule fragment_length_filtering:
    input:
        "{dir}/{{name}}.bam".format(dir = BAM_DIR),
        "{dir}/{{name}}.bam.bai".format(dir = BAM_DIR)
    output:
        "{dir}/{{name}}_filtered.bam".format(dir = BAM_DIR)
    log:
        os.path.join(LOG_DIR, "{name}_fragment_length_filtering.log")
    params:
        length = MIN_FRAG_LEN if FILTER_OPTION == "--minFragmentLength" else MAX_FRAG_LEN,
    threads: THREADS
    shell:
        "alignmentSieve -b {input[0]} -o {output} {FILTER_OPTION} {params.length} --numberOfProcessors {threads} |& tee {log}"
        # "alignmentSieve -b {input[0]} -o {output} --minFragmentLength {MIN_FRAG_LEN} --maxFragmentLength {MAX_FRAG_LEN} --numberOfProcessors {threads} |& tee {log}"

# computes the spike-in genome background and spike-in calibration coefficients
rule compute_coefficients:
    input:
        t_bam = expand("{dir}/{t}_main_filtered.bam", dir = BAM_DIR, t = treatment_name) if FILTER_FRAG_LEN else expand("{dir}/{t}_main.bam", dir = BAM_DIR, t = treatment_name),
        spkin_t_bam = expand("{dir}/{t}_spkin_filtered.bam", dir = BAM_DIR, t = treatment_name) if FILTER_FRAG_LEN else expand("{dir}/{t}_spkin.bam", dir = BAM_DIR, t = treatment_name),
        spkin_c_bam = expand("{dir}/{c}_spkin_filtered.bam", dir = BAM_DIR, c = control_name) if FILTER_FRAG_LEN else expand("{dir}/{c}_spkin.bam", dir = BAM_DIR, c = control_name),
        insert_size = "{dir}/{name}_insert_size.txt".format(dir = FASTP_OUT_DIR, name = treatment_name)
    output:
        gb = "{dir}/gb.txt".format(dir = OUTPUT_DIR),    # the spike-in genome background, an "average" of the treatment reads
        t = "{dir}/{t}_spkin.txt".format(dir = TEMP_DIR, t = treatment_name), # spike-in calibration coefficient for treatment
        c = "{dir}/{c}_spkin.txt".format(dir = TEMP_DIR, c = control_name)    # spike-in calibration coefficient for control
    params:
        N = float(MAIN_GEMOME_SIZE),   # the total genome size of the treatment
    run:
        t_reads = []
        for file in input.t_bam:
            t_reads.append(BamReads(file))   
        t_spkin_reads = []
        for file in input.spkin_t_bam:
            t_spkin_reads.append(BamReads(file))
        c_spkin_reads = []
        for file in input.spkin_c_bam:
            c_spkin_reads.append(BamReads(file))
        with open(input.insert_size, "r") as file:
            lines = file.readlines()
            total_length = 0
            for i in range(0, len(t_reads), 1):
                total_length += t_reads[i] * int(lines[i])
        printsep("insert_size: " + str(total_length / sum(t_reads)))   # insert size used later when calculating genome background
        if sum(t_spkin_reads) == 0 or sum(c_spkin_reads) == 0:
            # if any sample gets 0 spike-in read, then spike-in calibration can not be used
            # so set both scales to 1
            t_spkin_scale = 1
            c_spkin_scale = 1
        else: # scale down the larger sample to the smaller one
            if t_spkin_reads >= c_spkin_reads:
                t_spkin_scale = sum(c_spkin_reads) / sum(t_spkin_reads) 
                c_spkin_scale = 1
            else:
                t_spkin_scale = 1
                c_spkin_scale = sum(t_spkin_reads) / sum(c_spkin_reads) 
        gb = total_length * t_spkin_scale / params.N  # CALCULATES THE SPIKE IN GENOME BACKGROUND
        printsep(f"t_spkin_reads: {str(sum(t_spkin_reads))}; t_spkin_scale: {str(t_spkin_scale)}")
        printsep(f"c_spkin_reads: {str(sum(c_spkin_reads))}; c_spkin_scale: {str(c_spkin_scale)}")
        printsep("gb: " + str(gb))
        open(output.t, "w").write(str(t_spkin_scale))
        open(output.c, "w").write(str(c_spkin_scale))
        open(output.gb, "w").write(str(gb))

# determines flag statistics of sorted bam. Output not used elsewhere in the pipeline, useful when checking if alignment succeeded
rule flagstat:
    input:
        "{dir}/{{name}}.bam".format(dir = BAM_DIR)
    output:
        "{dir}/{{name}}.flagstat".format(dir = BAM_DIR)
    log:
        os.path.join(LOG_DIR, "{name}_flagstat.log")
    threads: 1
    shell:
        "(sambamba flagstat --nthreads {threads} {input} > {output}) |& tee {log}"

rule bam_index:
    input:
        "{dir}/{{name}}.bam".format(dir = BAM_DIR)
    output:
        "{dir}/{{name}}.bam.bai".format(dir = BAM_DIR)
    log:
        os.path.join(LOG_DIR, "{name}_bam_index.log")
    threads: 1
    shell:
        "(sambamba index --nthreads {threads} {input} > {output}) |& tee {log}"

# ------------------------------------------------------peak calling------------------------------------------------------------

def get_pileup_input(wildcards):
    if treatment_name in wildcards.name:
        names = treatment_name
    elif control_name in wildcards.name:
        names = control_name
    else:
        raise ValueError("Invalid name wildcard")
    if FILTER_FRAG_LEN:
        return expand("{dir}/{names}_main_filtered.bam", dir = BAM_DIR, names = names)
    return expand("{dir}/{names}_main.bam", dir = BAM_DIR, names = names)

rule macs_pileup:
    input: 
        get_pileup_input  # takes input from rule 
    output:
        temp("{dir}/{{name}}_pileup.bdg".format(dir = TEMP_DIR))
    log:
        os.path.join(LOG_DIR, "{name}_pileup.log")
    shell:
        "{MACS} pileup -f BAMPE -i {input} -o {output} |& tee {log}"

rule macs_pseudocount:
    input:
        main = "{dir}/{{name}}_pileup.bdg".format(dir = TEMP_DIR)
    output:
        temp("{dir}/{{name}}_pileup_pseudo.bdg".format(dir = TEMP_DIR))
    log:
        os.path.join(LOG_DIR, "{name}_pileup_pseudo.log")
    shell:
        "{MACS} bdgopt -i {input.main} -m add -p 0.1 -o {output} |& tee {log}"

rule macs_spike_in:
    input:
        main = "{dir}/{{name}}_pileup_pseudo.bdg".format(dir = TEMP_DIR),
        spkin = "{dir}/{{name}}_spkin.txt".format(dir = TEMP_DIR)
    output:
        temp("{dir}/{{name}}_pileup_spkin.bdg".format(dir = TEMP_DIR))
    log:
        os.path.join(LOG_DIR, "{name}_pileup_spkin.log")
    shell:
        "{MACS} bdgopt -i {input.main} -m multiply -p $(cat {input.spkin}) -o {output} |& tee {log}"

rule macs_max:
    input:
        bdg = "{dir}/{c}_pileup_spkin.bdg".format(dir = TEMP_DIR, c = control_name),
        gb = "{dir}/gb.txt".format(dir = OUTPUT_DIR)
    output:
        temp("{dir}/{c}_pileup_max.bdg".format(dir = TEMP_DIR, c = control_name))
    log:
        os.path.join(LOG_DIR, "{c}_pileup_spkin_max.log".format(c = control_name))
    shell:
        "{MACS} bdgopt -i {input.bdg} -m max -p $(cat {input.gb}) -o {output} |& tee {log}"

rule compute_local_scales:
    input:
        t_bam = expand("{dir}/{t}_main.bam", dir = BAM_DIR, t = treatment_name),
        c_bam = expand("{dir}/{c}_main.bam", dir = BAM_DIR, c = control_name),
        d = "{dir}/d.txt".format(dir = MACS_DIR)
    output:
        tc = "{dir}/t_c_scale.txt".format(dir = TEMP_DIR),
        slocal = "{dir}/slocal.txt".format(dir = TEMP_DIR),
        llocal = "{dir}/llocal.txt".format(dir = TEMP_DIR)
    run:
        t_reads = []
        for file in input.t_bam:
            t_reads.append(BamReads(file)) 
        c_reads = []
        for file in input.c_bam:
            c_reads.append(BamReads(file))
        tcScale = sum(t_reads) / sum(c_reads)
        with open(input.d) as f:
            d = float(f.readline())
        slocal = d / 1000
        llocal = d / 10000
        printsep("tc_scale: " + str(tcScale))
        printsep("slocal: " + str(slocal))
        printsep("llocal: " + str(llocal))
        open(output.tc, "w").write(str(tcScale))
        open(output.slocal, "w").write(str(slocal))
        open(output.llocal, "w").write(str(llocal))

rule macs_predictd:
    input: 
        "{dir}/{c}_main.bam".format(dir = BAM_DIR, c = control_name)
    output:
        "{dir}/d.txt".format(dir = MACS_DIR)
    log:
        os.path.join(LOG_DIR, "{c}_predictd.log".format(c = control_name))
    shell:
        """
        {MACS} predictd -i {input} -f BAM --outdir {MACS_DIR} |& tee {log}
        grep "predicted fragment length is" {log} | grep -o "[0-9]*" | tail -1 > {output}
        """

rule macs_slocal:
    input: 
        bam = "{dir}/{c}_main.bam".format(dir = BAM_DIR, c = control_name),
        scale = "{dir}/slocal.txt".format(dir = TEMP_DIR)
    output:
        bg = temp("{dir}/1k_bg.bdg".format(dir = TEMP_DIR)),
        bg_norm = temp("{dir}/1k_bg_norm.bdg".format(dir = TEMP_DIR))
    shell:
        """
        {MACS} pileup -f BAMPE -i {input.bam} -B --extsize 500 -o {output.bg}
        {MACS} bdgopt -i {output.bg} -m multiply -p $(cat {input.scale}) -o {output.bg_norm}
        """

rule macs_llocal:
    input: 
        bam = "{dir}/{c}_main.bam".format(dir = BAM_DIR, c = control_name),
        scale = "{dir}/llocal.txt".format(dir = TEMP_DIR)
    output:
        bg = temp("{dir}/10k_bg.bdg".format(dir = TEMP_DIR)),
        bg_norm = temp("{dir}/10k_bg_norm.bdg".format(dir = TEMP_DIR))
    shell:
        """
        {MACS} pileup -f BAMPE -i {input.bam} -B --extsize 5000 -o {output.bg}
        {MACS} bdgopt -i {output.bg} -m multiply -p $(cat {input.scale}) -o {output.bg_norm}
        """

rule macs_local:
    input: 
        bdg = "{dir}/{c}_pileup_max.bdg".format(dir = TEMP_DIR, c = control_name) if SPIKE_IN_CALIBRATE else "{dir}/{c}_pileup_pseudo.bdg".format(dir = TEMP_DIR, c = control_name),
        onek = "{dir}/1k_bg_norm.bdg".format(dir = TEMP_DIR),
        tenk = "{dir}/10k_bg_norm.bdg".format(dir = TEMP_DIR),
        t_c_scale = "{dir}/t_c_scale.txt".format(dir = TEMP_DIR)
    output:
        "{dir}/local_lambda.bdg".format(dir = TEMP_DIR)
    shell:
        """
        {MACS} bdgcmp -m max -t {input.onek} -c {input.tenk} -o {TEMP_DIR}/1k_10k_bg_norm.bdg
        {MACS} bdgcmp -m max -t {TEMP_DIR}/1k_10k_bg_norm.bdg -c {input.bdg} -o {TEMP_DIR}/d_1k_10k_bg_norm.bdg
        {MACS} bdgopt -i {TEMP_DIR}/d_1k_10k_bg_norm.bdg -m multiply -p $(cat {input.t_c_scale}) -o {output}
        """

rule macs_qscore:
    input:
        t = "{dir}/{name}_pileup_spkin.bdg".format(dir = TEMP_DIR, name = treatment_name) if SPIKE_IN_CALIBRATE else "{dir}/{name}_pileup_pseudo.bdg".format(dir = TEMP_DIR, name = treatment_name),
        c = "{dir}/local_lambda.bdg".format(dir = TEMP_DIR)
    output:
        temp("{dir}/{t}_{c}_qscores.bdg".format(dir = TEMP_DIR, t = treatment_name, c = control_name))
    log:
        os.path.join(LOG_DIR, "{t}_{c}_qscores.log".format(t = treatment_name, c = control_name))
    params:
        method = "-m qpois" # the Q value of poisson distribution
    shell:
        "{MACS} bdgcmp {params.method} -t {input.t} -c {input.c} -o {output} |& tee {log}"

if not USE_MACS_CALLPEAK:
    rule call_peaks_macs:
        input:
            "{dir}/{t}_{c}_qscores.bdg".format(dir = TEMP_DIR, t = treatment_name, c = control_name)
        output:
            "{dir}/{t}_{c}_peaks.narrowPeak".format(dir = MACS_DIR, t = treatment_name, c = control_name)
        log:
            os.path.join(LOG_DIR, "{t}_{c}_peaks.log".format(t = treatment_name, c = control_name))
        params:
            fdr = "-c {x}".format(x = str(-math.log(float(MACS_FDR), 10)))
        shell:
            "{MACS} bdgpeakcall -i {input} {params.fdr} -l {MIN_PEAK_LEN} -o {output} |& tee {log}"

rule macs_cutoff_analysis:
    input:
        "{dir}/{t}_{c}_qscores.bdg".format(dir = TEMP_DIR, t = treatment_name, c = control_name)
    output:
        "{dir}/{t}_{c}_cutoff_analysis.txt".format(dir = MACS_DIR, t = treatment_name, c = control_name)
    log:
        os.path.join(LOG_DIR, "{t}_{c}_cutoff_analysis.log".format(t = treatment_name, c = control_name))
    params:
        fdr = "-c {x}".format(x = str(-math.log(float(MACS_FDR), 10)))
    shell:
        "{MACS} bdgpeakcall --cutoff-analysis -i {input} {params.fdr} -l {MIN_PEAK_LEN} -o {output} |& tee {log}"

def get_macs_input(wildcards):
    inputs = []
    suffix = ""
    if FILTER_FRAG_LEN:
        suffix = "_filtered"
    inputs.append(f"{BAM_DIR}/{treatment_name}_main{suffix}.bam")
    inputs.append(f"{BAM_DIR}/{control_name}_main{suffix}.bam")
    return inputs

if USE_MACS_CALLPEAK:
    rule macs_no_spkin:
        input:
            treatment = get_macs_input[0],
            control = get_macs_input[1]
        output:
            "{dir}/{t}_{c}_peaks.narrowPeak".format(dir = MACS_DIR, t = treatment_name, c = control_name)
        log:
            os.path.join(LOG_DIR, "{t}_{c}_macs_no_spkin.log".format(t = treatment_name, c = control_name))
        params:
            name = "{t}_{c}".format(t = treatment_name, c = control_name)
        shell:
            "{MACS} callpeak -t {input.treatment} -c {input.control} -f BAMPE --keep-dup all --outdir {MACS_DIR} -n {params.name}"

# conversion from bam to bed, preparing for spike in callibration
rule bam_to_bed:
    input:
        "{dir}/{{name}}.bam".format(dir = BAM_DIR)
    output:
        "{dir}/{{name}}.bed".format(dir = BED_DIR)
    shell:
        "bedtools bamtobed -i {input} | " \
        "awk 'BEGIN{{FS=OFS=\"\\t\";}} {{len = $3 - $2; print $0,len;}}' | " \
        "sort -T {TEMP_DIR} -k1,1 -k2,2n > {output}"

rule bed_to_bdg:
    input:
        "{dir}/{{name}}.bed".format(dir = BED_DIR)
    output:
        "{dir}/{{name}}.bedgraph".format(dir = BDG_DIR)
    shell:
        "bedtools genomecov -bg -i {input} -g {CHR_LEN_FILE} > {output}"

def get_seacr_bdg(wildcards):
    suffix = ""
    if FILTER_FRAG_LEN:
        suffix = "_filtered"
    return expand("{dir}/{IN}_main{suffix}.bedgraph", dir = BDG_DIR, IN = [treatment_name, control_name], suffix = suffix)

if not SPIKE_IN_CALIBRATE:
    rule seacr_no_spkin:
        input:
            get_seacr_bdg
        output:
            "{dir}/SEACR_{t}_{c}_peaks.auc.threshold.merge.bed".format(dir = SEACR_DIR, t = treatment_name, c = control_name)
        params:
            norm = "non"
        shell:
            "cd {SEACR_DIR}; bash {SEACR_SCRIPT} {input[0]} {input[1]} " \
            "{params.norm} {SEACR_MODE} {treatment_name}_{control_name}_peaks ;" \
            "cd .. ; mv {SEACR_DIR}/{treatment_name}_{control_name}_peaks.{SEACR_MODE}.bed {SEACR_DIR}/SEACR_{treatment_name}_{control_name}_peaks.auc.threshold.merge.bed"

def get_seacr_bed(wildcards):
    suffix = ""
    if FILTER_FRAG_LEN:
        suffix = "_filtered"
    return f"{BED_DIR}/{wildcards.name}_main{suffix}.bed"

rule seacr_spike_in:
    input:
        main = get_seacr_bed,
        spkin = "{dir}/{{name}}_spkin.txt".format(dir = TEMP_DIR)
    output:
        f"{BDG_DIR}/{{name}}_main.{SPIKE_IN_MIN_LEN}-{SPIKE_IN_MAX_LEN}.spkin.bedgraph", # TODO: remove the _main in output file
        f"{BDG_DIR}/{{name}}_spike_in_out.spkin.txt"
    params:
        output = "bg"
    shell:
        "cd {BDG_DIR} ; {SPIKE_IN_SCRIPT} {input.main} $(cat {input.spkin}) " \
        "{params.output} {CHR_LEN_FILE} {SPIKE_IN_MIN_LEN} {SPIKE_IN_MAX_LEN} > {output[1]} ;" \ 
        "mv {wildcards.name}_main*.{SPIKE_IN_MIN_LEN}-{SPIKE_IN_MAX_LEN}.bedgraph {output[0]}"
       
if SPIKE_IN_CALIBRATE:
    rule call_peaks_seacr:
        input:
            expand("{dir}/{IN}_main.{min_len}-{max_len}.spkin.bedgraph", dir = BDG_DIR, min_len= SPIKE_IN_MIN_LEN, max_len= SPIKE_IN_MAX_LEN, IN = [treatment_name, control_name])
        output:
            f"{SEACR_DIR}/SEACR_{treatment_name}_{control_name}_peaks.auc.threshold.merge.bed"
        params:
            norm = "non"
        shell:
            "cd {SEACR_DIR}; bash {SEACR_SCRIPT} {input[0]} {input[1]} " \
            "{params.norm} {SEACR_MODE} {treatment_name}_{control_name}_peaks ;" \
            "cd .. ; mv {SEACR_DIR}/{treatment_name}_{control_name}_peaks.{SEACR_MODE}.bed {SEACR_DIR}/SEACR_{treatment_name}_{control_name}_peaks.auc.threshold.merge.bed"

def GenerateFastaInput(wildcards):
    if wildcards.caller == "SEACR":
        return f"{SEACR_DIR}/SEACR_{treatment_name}_{control_name}_peaks.auc.threshold.merge.bed"
    elif wildcards.caller == "MACS2":
        return "{dir}/{t}_{c}_peaks.narrowPeak".format(dir = MACS_DIR, t = treatment_name, c = control_name)
    else:
        print("FALSE FASTA INPUT PEAK CALLER") # TODO: throw exception
        return False

# --------------------------------------------------motif analysis-------------------------------------------------------------------------

rule generate_fasta_macs2_summit:
    input:
        "{dir}/{t}_{c}_peaks.narrowPeak".format(dir = MACS_DIR, t = treatment_name, c = control_name)
    output:
        "{dir}/MACS2_{t}_{c}_summit.fa".format(dir = OUTPUT_DIR, t = treatment_name, c = control_name)
    log:
        os.path.join(LOG_DIR, "MACS2_{t}_{c}_fa.log".format(t = treatment_name, c = control_name))
    shell:
        "(awk 'BEGIN{{FS=OFS=\"\\t\"}} /^chr/{{print $1,$2+$10,$2+$10+1}}' {input} | " \
        "bedtools slop -i /dev/stdin -g {CHR_LEN_FILE} -l {FASTA_EXT_LEN} -r $(({FASTA_EXT_LEN} - 1)) | " \
        "cut -f 1-3 | fastaFromBed -fi {MAIN_GEMOME} -bed /dev/stdin -fo {output}) |& tee {log}"

rule generate_fasta_macs2_center:
    input:
        "{dir}/{t}_{c}_peaks.narrowPeak".format(dir = MACS_DIR, t = treatment_name, c = control_name)
    output:
        "{dir}/MACS2_{t}_{c}_center.fa".format(dir = OUTPUT_DIR, t = treatment_name, c = control_name)
    log:
        os.path.join(LOG_DIR, "MACS2_{t}_{c}_fa.log".format(t = treatment_name, c = control_name))
    shell:
        "(awk 'BEGIN{{FS=OFS=\"\\t\"}} /^chr/{{centre=($2+$3)/2; printf \"%s\\t%d\\t%d\\n\", $1, centre, centre+1}}' {input} | "
        "bedtools slop -i /dev/stdin -g {CHR_LEN_FILE} -l {FASTA_EXT_LEN} -r $(({FASTA_EXT_LEN} - 1)) | " \
        "cut -f 1-3 | fastaFromBed -fi {MAIN_GEMOME} -bed /dev/stdin -fo {output}) |& tee {log}"

rule generate_fasta_seacr_summit:
    input:
        "{dir}/SEACR_{t}_{c}_peaks.auc.threshold.merge.bed".format(dir = SEACR_DIR, t = treatment_name, c = control_name)
    output:
        "{dir}/SEACR_{t}_{c}_summit.fa".format(dir = OUTPUT_DIR, t = treatment_name, c = control_name)
    log:
        os.path.join(LOG_DIR, "SEACR_{t}_{c}_fa.log".format(t = treatment_name, c = control_name))
    shell:
        "(awk 'BEGIN{{FS=\"[-\\t:]\"; OFS=\"\\t\"}} /^chr/{{centre=($7+$8)/2; printf \"%s\\t%d\\t%d\\n\", $1, centre, centre+1}}' {input} | " \
        "bedtools slop -i /dev/stdin -g {CHR_LEN_FILE} -l {FASTA_EXT_LEN} -r $(({FASTA_EXT_LEN} - 1)) | " \
        "cut -f 1-3 | fastaFromBed -fi {MAIN_GEMOME} -bed /dev/stdin -fo {output}) |& tee {log}"

rule generate_fasta_seacr_center:
    input:
        "{dir}/SEACR_{t}_{c}_peaks.auc.threshold.merge.bed".format(dir = SEACR_DIR, t = treatment_name, c = control_name)
    output:
        "{dir}/SEACR_{t}_{c}_center.fa".format(dir = OUTPUT_DIR, t = treatment_name, c = control_name)
    log:
        os.path.join(LOG_DIR, "SEACR_{t}_{c}_fa.log".format(t = treatment_name, c = control_name))
    shell:
        "(awk 'BEGIN{{FS=OFS=\"\\t\"}} /^chr/{{centre=($2+$3)/2; printf \"%s\\t%d\\t%d\\n\", $1, centre, centre+1}}' {input} | "
        "bedtools slop -i /dev/stdin -g {CHR_LEN_FILE} -l {FASTA_EXT_LEN} -r $(({FASTA_EXT_LEN} - 1)) | " \
        "cut -f 1-3 | fastaFromBed -fi {MAIN_GEMOME} -bed /dev/stdin -fo {output}) |& tee {log}"

rule generate_motifs: # requires module load meme/5.5.2
    input:
        "{dir}/{{caller}}_{t}_{c}_{{location}}.fa".format(dir = OUTPUT_DIR, t = treatment_name, c = control_name)
    output:
        directory("{dir}/{{caller}}_{t}_{c}_motifs_{{location}}".format(dir = OUTPUT_DIR, t = treatment_name, c = control_name))
    log:
        os.path.join(LOG_DIR, "{{caller}}_{t}_{c}_motifs_{{location}}.log".format(t = treatment_name, c = control_name))
    params:
        minw = 7,
        maxw = 12,
        db = JASPAR,
        fast = "-meme-nmotifs 0 -streme-nmotifs 0 -fimo-skip -spamo-skip -centrimo-noseq" if SKIP_DREME else "" # this step skips meme, dreme, fima, spamo
        # which speeds up the process dramatically. If dreme and meme is kept, this step could take up to 2 weeks to complete. 
        # "-centrimo-noseq" makes the output html file a lot smaller. If this parameter is removed, opening up the html could drain all your computer's memory. 
        # If we are doing an actual run (meaning the results might actually go into a paper), this set of parameters should be removed. 
    threads: THREADS
    shell:
        "meme-chip -dna -minw {params.minw} -maxw {params.maxw} {params.fast} -meme-p {threads} -db {params.db} -oc {output} {input} |& tee {log}"

# ---------------------------------------------quality control---------------------------------------------------------

def get_fastqc_input_file(wildcards):
    if wildcards.name.endswith(".fastp"):
        return "{dir}/{name}.fastq.gz".format(dir = FASTP_OUT_DIR, name = wildcards.name)
    else:
        return "{dir}/{name}.fastq.gz".format(dir = RAW_DIR, name = wildcards.name)

rule fastqc_pre_fastp:
    input:
        "{dir}/{{name}}.fastq.gz".format(dir = RAW_DIR)
    output:
        "{dir}/{{name}}_fastqc.html".format(dir = PRE_FASTP_FASTQC_DIR),
        "{dir}/{{name}}_fastqc.zip".format(dir = PRE_FASTP_FASTQC_DIR)  
    params:
        option = "--noextract -d {TEMP_DIR}".format(TEMP_DIR = TEMP_DIR)
    log:
        os.path.join(LOG_DIR, "{name}_fastqc.log")
    threads: 2
    shell:
        "fastqc {params.option} -t {threads} -o {PRE_FASTP_FASTQC_DIR} {input} |& tee {log}"

rule fastqc_post_fastp:
    input:
        "{dir}/{{name}}.fastq.gz".format(dir = FASTP_OUT_DIR)
    output:
        "{dir}/{{name}}_fastqc.html".format(dir = POST_FASTP_FASTQC_DIR),
        "{dir}/{{name}}_fastqc.zip".format(dir = POST_FASTP_FASTQC_DIR)
    params:
        option = "--noextract -d {TEMP_DIR}".format(TEMP_DIR = TEMP_DIR)
    log:
        os.path.join(LOG_DIR, "{name}_fastqc.log")
    threads: 2
    shell:
        "fastqc {params.option} -t {threads} -o {POST_FASTP_FASTQC_DIR} {input} |& tee {log}"

rule fastqc_post_fastp_finish:
    input:
        "{dir}/{{name}}_fastqc.zip".format(dir = POST_FASTP_FASTQC_DIR)
    output:
        "{dir}/{{name}}_fastqc.zip".format(dir = MULTIQC_DIR)  
    run:
        create_symlink(wildcards, rules.fastqc_post_fastp_finish.input[0], rules.fastqc_post_fastp_finish.output[0])


rule picard:    # requires module load java/11.0.7, picard/2.10.9, R/3.6.1
    input:
        "{dir}/{{name}}.bam".format(dir = BAM_DIR)
    output:
        "{dir}/{{name}}_picard_insert_size_metrics.txt".format(dir = QC_DIR),
        "{dir}/{{name}}_picard_insert_size_plot.pdf".format(dir = QC_DIR)
    log:
        os.path.join(LOG_DIR, "{name}_picard.log")
    shell:
        "(java -jar $picard_dir/picard.jar CollectInsertSizeMetrics " \
        "INPUT={input} OUTPUT={output[0]} " \
        "HISTOGRAM_FILE={output[1]} METRIC_ACCUMULATION_LEVEL=ALL_READS) |& tee {log}"

rule picard_finish:
    input:
        "{dir}/{{name}}_picard_insert_size_metrics.txt".format(dir = QC_DIR)
    output:
        "{dir}/{{name}}_picard_insert_size_metrics.txt".format(dir = MULTIQC_DIR)
    run:
        create_symlink(wildcards, rules.picard_finish.input[0], rules.picard_finish.output[0])

rule qualimap:  # requires module load java/11.0.7 qualimap/2.2
    input:
        "{dir}/{{name}}_main.bam".format(dir = BAM_DIR),  # from rule alignment
    output:
        directory("{dir}/{{name}}_qualimap_main_stats".format(dir = QC_DIR))
    log:
        os.path.join(LOG_DIR, "{name}_qualimap_bamqc.log")
    params:
        memSize = "4g"  # or else will run out of memory and job will fail
    threads: 2
    shell:
        "qualimap bamqc -nt {threads} --java-mem-size={params.memSize} -outdir {output} -gd {QUALIMAP_ASSEMBLY} -ip -sd -c -bam {input} |& tee {log}"

rule qualimap_finish:
    input:
        "{dir}/{{name}}_qualimap_main_stats".format(dir = QC_DIR)
    output:
        directory("{dir}/{{name}}_qualimap_main_stats".format(dir = MULTIQC_DIR))
    run:
        create_symlink(wildcards, rules.qualimap_finish.input[0], rules.qualimap_finish.output[0])


def run_multiqc():
    log = os.path.join(LOG_DIR, "{t}_{c}_multiqc.log".format(t = treatment_name, c = control_name))
    p = "-f -v"
    c = "-c {}".format(MULTIQC_CONFIG)
    filename = "-n {t}_{c}_multiqc_report".format(t = treatment_name, c = control_name)
    cmd = "multiqc {p} -o {OUTPUT_DIR} {filename} {MULTIQC_DIR} |& tee {log}".format(p = p, \
        OUTPUT_DIR = OUTPUT_DIR, filename = filename, MULTIQC_DIR = MULTIQC_DIR, log = log)
    shell(cmd)


rule multiqc_manual:
    output:
        "{dir}/{t}_{c}_multiqc_report.html".format(dir = OUTPUT_DIR, t = treatment_name, c = control_name)
    run:
        run_multiqc()


onsuccess:
    if NEED_QC and not ONLY_MULTIQC:
        run_multiqc()
