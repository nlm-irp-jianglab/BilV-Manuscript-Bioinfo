#!/bin/bash
set -euo pipefail

# Per-sample SRA processing: [download or local] → QC → [optional host removal] → gene alignment → counts
# Usage: process_sample.sh [-r HOST_BT2_INDEX] [-s SRA_DIR] [-k] [-a] [-c CACHE_DIR] <SRA_ACCESSION> <OUTDIR> <GENE_BT2_INDEX_PREFIX>
#   -r  Host bowtie2 index prefix for host read removal (any organism)
#   -s  Directory containing pre-downloaded .sra files (from download_sra.sh)
#   -k  Use kingfisher to download from ENA (outputs .fastq.gz, no fasterq-dump needed)
#   -a  Use aria2c to download FASTQ from ENA FTP directly
#   -c  Directory with pre-downloaded fastq.gz (structure: CACHE_DIR/ACC/ACC_1.fastq.gz)

SCRIPTDIR="$(cd "$(dirname "$0")" && pwd)"
source "${PIPELINE_CONF:-${SCRIPTDIR}/pipeline.conf}"

HOST_IDX=""
SRA_DIR=""
CACHE_DIR=""
USE_KF=false
USE_ARIA=false
while getopts "r:s:c:ka" opt; do
    case $opt in
        r) HOST_IDX="$OPTARG" ;;
        s) SRA_DIR="$OPTARG" ;;
        c) CACHE_DIR="$OPTARG" ;;
        k) USE_KF=true ;;
        a) USE_ARIA=true ;;
        *) echo "Usage: $0 [-r HOST_BT2_INDEX] [-s SRA_DIR] [-c CACHE_DIR] [-k] [-a] <SRA_ACCESSION> <OUTDIR> <GENE_BT2_INDEX_PREFIX>"; exit 1 ;;
    esac
done
shift $((OPTIND - 1))

if [[ $# -ne 3 ]]; then
    echo "Usage: $0 [-r HOST_BT2_INDEX] [-s SRA_DIR] [-k] <SRA_ACCESSION> <OUTDIR> <GENE_BT2_INDEX_PREFIX>"
    exit 1
fi

SRA=$1
OUTDIR=$2
GENE_IDX=$3
[[ "$GENE_IDX" != /* ]] && GENE_IDX="${SCRIPTDIR}/../${GENE_IDX}"
THREADS=${THREADS:-8}
WORK="${TMPDIR}/${SRA}"
SAMPLE_OUT="${OUTDIR}/${SRA}"

log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"; }

count_reads() {
    local fq=$1
    if [[ "$fq" == *.gz ]]; then
        echo $(( $(zcat "$fq" | wc -l) / 4 ))
    else
        echo $(( $(wc -l < "$fq") / 4 ))
    fi
}

parse_bt2_log() {
    local log=$1 out=$2
    awk -F'[[:space:]]+' '
    /reads; of these:/                        { total=$1 }
    /aligned concordantly exactly 1 time/     { conc1=$1 }
    /aligned concordantly >1 times/           { concM=$1 }
    /aligned concordantly 0 times$/           { conc0=$1 }
    /overall alignment rate/                  { rate=$1 }
    END {
        print "metric\tvalue"
        print "total_pairs\t" total
        print "concordant_once\t" conc1
        print "concordant_multi\t" concM
        print "concordant_none\t" conc0
        print "overall_alignment_rate\t" rate
    }' "$log" > "$out"
}

mkdir -p "$WORK" "$SAMPLE_OUT"
cd "$WORK"

# --- Get FASTQs ---
RAW_R1=""
RAW_R2=""
if [[ -n "$CACHE_DIR" ]]; then
    # Use pre-downloaded fastq.gz from cache
    CACHE_R1="${CACHE_DIR}/${SRA}/${SRA}_1.fastq.gz"
    CACHE_R2="${CACHE_DIR}/${SRA}/${SRA}_2.fastq.gz"
    if [[ ! -f "$CACHE_R1" || ! -f "$CACHE_R2" ]]; then
        log "ERROR: cached fastq.gz not found: ${CACHE_DIR}/${SRA}/"
        exit 1
    fi
    log "Linking cached FASTQs for ${SRA}"
    ln -sf "$CACHE_R1" "${WORK}/${SRA}_1.fastq.gz"
    ln -sf "$CACHE_R2" "${WORK}/${SRA}_2.fastq.gz"
    RAW_R1="${SRA}_1.fastq.gz"
    RAW_R2="${SRA}_2.fastq.gz"
elif [[ -n "$SRA_DIR" ]]; then
    # Phase 2: extract from pre-downloaded .sra file
    SRA_FILE="${SRA_DIR}/${SRA}/${SRA}.sra"
    if [[ ! -f "$SRA_FILE" ]]; then
        log "ERROR: .sra file not found: ${SRA_FILE}"
        exit 1
    fi
    log "Extracting FASTQs from ${SRA_FILE}"
    fasterq-dump "$SRA_FILE" --split-files --threads "$THREADS" -O "$WORK" -t "$WORK"
    RAW_R1="${SRA}_1.fastq"
    RAW_R2="${SRA}_2.fastq"
elif [[ "$USE_KF" == true ]]; then
    # Kingfisher: download compressed FASTQ from ENA/AWS
    log "Downloading ${SRA} via kingfisher"
    COLUMNS=120 LINES=40 "$KINGFISHER" get -r "$SRA" \
        -m ena-ftp ena-ascp aws-http prefetch \
        -f fastq.gz --download-threads 4 -t "$THREADS" \
        --output-directory "$WORK"
    RAW_R1="${SRA}_1.fastq.gz"
    RAW_R2="${SRA}_2.fastq.gz"
elif [[ "$USE_ARIA" == true ]]; then
    # aria2c: download compressed FASTQ directly from ENA FTP
    log "Downloading ${SRA} via aria2c from ENA"
    PREFIX="${SRA:0:6}"
    ACC_LEN=${#SRA}
    if [[ $ACC_LEN -gt 9 ]]; then
        EXTRA=$((ACC_LEN - 9))
        TAIL="${SRA: -${EXTRA}}"
        # Force base-10 to avoid octal interpretation of 08/09
        SUBDIR=$(printf "%03d" "$((10#$TAIL))")
        ENA_DIR="vol1/fastq/${PREFIX}/${SUBDIR}/${SRA}"
    else
        ENA_DIR="vol1/fastq/${PREFIX}/${SRA}"
    fi
    ENA_BASE="http://ftp.sra.ebi.ac.uk/${ENA_DIR}"
    for attempt in 1 2 3; do
        ~/bin/aria2c -x 1 -s 1 --retry-wait=15 --max-tries=5 \
            --connect-timeout=60 --timeout=600 \
            -d "$WORK" "${ENA_BASE}/${SRA}_1.fastq.gz"
        ~/bin/aria2c -x 1 -s 1 --retry-wait=15 --max-tries=5 \
            --connect-timeout=60 --timeout=600 \
            -d "$WORK" "${ENA_BASE}/${SRA}_2.fastq.gz"
        if [[ -f "${WORK}/${SRA}_1.fastq.gz" && -f "${WORK}/${SRA}_2.fastq.gz" ]]; then
            break
        fi
        log "aria2c attempt ${attempt}: missing files, retrying in 60s..."
        rm -f "${WORK}/${SRA}_1.fastq.gz" "${WORK}/${SRA}_2.fastq.gz"
        sleep 60
        if [[ $attempt -eq 3 ]]; then
            log "ERROR: aria2c download failed after 3 attempts"
            exit 1
        fi
    done
    RAW_R1="${SRA}_1.fastq.gz"
    RAW_R2="${SRA}_2.fastq.gz"
else
    # Direct fasterq-dump download
    log "Downloading ${SRA}"
    for attempt in 1 2 3; do
        if fasterq-dump "$SRA" --split-files --threads "$THREADS" -O "$WORK" -t "$WORK"; then
            break
        fi
        log "fasterq-dump attempt ${attempt} failed, retrying in 30s..."
        sleep 30
        if [[ $attempt -eq 3 ]]; then
            log "ERROR: fasterq-dump failed after 3 attempts"
            exit 1
        fi
    done
    RAW_R1="${SRA}_1.fastq"
    RAW_R2="${SRA}_2.fastq"
fi
raw_reads=$(count_reads "$RAW_R1")
log "Download complete: ${raw_reads} read pairs"

# --- QC with fastp ---
log "Running fastp"
fastp \
    -i "$RAW_R1" -I "$RAW_R2" \
    -o "${SRA}_clean_1.fastq" -O "${SRA}_clean_2.fastq" \
    --detect_adapter_for_pe \
    --qualified_quality_phred 20 \
    --length_required 50 \
    --thread 4 \
    --json "${SAMPLE_OUT}/${SRA}.fastp.json" \
    --html "${SAMPLE_OUT}/${SRA}.fastp.html"
rm "$RAW_R1" "$RAW_R2"
qc_reads=$(count_reads "${SRA}_clean_1.fastq")
log "QC complete: ${qc_reads} read pairs passed"

# Track which FASTQs go into gene alignment
ALIGN_R1="${SRA}_clean_1.fastq"
ALIGN_R2="${SRA}_clean_2.fastq"

# --- Host read removal (optional) ---
if [[ -n "$HOST_IDX" ]]; then
    log "Removing host reads using index: ${HOST_IDX}"
    bowtie2 --sensitive -p "$THREADS" -x "$HOST_IDX" \
        -1 "$ALIGN_R1" -2 "$ALIGN_R2" \
        --un-conc "${WORK}/${SRA}_nonhost_%.fastq" \
        > /dev/null 2> "${WORK}/${SRA}.host_align.tmp"
    parse_bt2_log "${WORK}/${SRA}.host_align.tmp" "${SAMPLE_OUT}/${SRA}.host_align.tsv"
    rm "$ALIGN_R1" "$ALIGN_R2" "${SRA}.host_align.tmp"
    ALIGN_R1="${SRA}_nonhost_1.fastq"
    ALIGN_R2="${SRA}_nonhost_2.fastq"
    nonhost_reads=$(count_reads "$ALIGN_R1")
    log "Host removal complete: ${nonhost_reads} read pairs retained"
fi

# --- Gene alignment ---
log "Aligning to gene reference"
bowtie2 --no-unal -p "$THREADS" -x "$GENE_IDX" \
    -1 "$ALIGN_R1" -2 "$ALIGN_R2" \
    2> "${WORK}/${SRA}.gene_align.tmp" \
    | samtools sort -@ 2 -o "${SAMPLE_OUT}/${SRA}.gene.bam" -
parse_bt2_log "${WORK}/${SRA}.gene_align.tmp" "${SAMPLE_OUT}/${SRA}.gene_align.tsv"
samtools index "${SAMPLE_OUT}/${SRA}.gene.bam"
rm "$ALIGN_R1" "$ALIGN_R2" "${SRA}.gene_align.tmp"
aligned_reads=$(samtools view -c -F 0x904 "${SAMPLE_OUT}/${SRA}.gene.bam")
log "Gene alignment complete: ${aligned_reads} reads aligned"

# --- Read counts ---
log "Generating read counts"
samtools idxstats "${SAMPLE_OUT}/${SRA}.gene.bam" | awk '$1 != "*"' > "${SAMPLE_OUT}/${SRA}.idxstats.tsv"

# --- Read count summary ---
{
    echo -e "step\tread_pairs"
    echo -e "raw\t${raw_reads}"
    echo -e "after_qc\t${qc_reads}"
    if [[ -n "$HOST_IDX" ]]; then
        echo -e "after_host_removal\t${nonhost_reads}"
    fi
    echo -e "aligned_to_genes\t${aligned_reads}"
} > "${SAMPLE_OUT}/${SRA}.read_stats.tsv"

log "Done — results in ${SAMPLE_OUT}"
