#!/bin/bash -e

## v3.2 of STARsolo wrappers is set up to guess the chemistry automatically
## newest version of the script uses STAR v2.7.10a with EM multimapper processing
## in STARsolo which on by default; the extra matrix can be found in /raw subdir

# --- Function Definitions ---

# Finds paired-end FASTQ files based on common naming conventions.
# Args:
#   $1: FQDIR - The directory containing FASTQ files.
#   $2: TAG - The sample identifier.
# Returns:
#   A pipe-separated string "R1_FILES|R2_FILES".
function find_fastq_files() {
  local FQDIR=$1
  local TAG=$2
  local R1
  local R2

  ## Four popular cases: ENA, regular, Cell Ranger, and HRA.
  ## The command below will generate a comma-separated list for each read
  ## if there are >1 file for each mate.
  if [[ $(find "$FQDIR"/* | grep -P "/$TAG[\\/\\._]" | grep "_1\\.f.*q") != "" ]]; then
    R1=$(find "$FQDIR"/* | grep -P "/$TAG[\\/\\._]" | grep "_1\\.f.*q" | sort | tr '\n' ',' | sed "s/,$//g")
    R2=$(find "$FQDIR"/* | grep -P "/$TAG[\\/\\._]" | grep "_2\\.f.*q" | sort | tr '\n' ',' | sed "s/,$//g")
  elif [[ $(find "$FQDIR"/* | grep -P "/$TAG[\\/\\._]" | grep "R1\\.f.*q") != "" ]]; then
    R1=$(find "$FQDIR"/* | grep -P "/$TAG[\\/\\._]" | grep "R1\\.f.*q" | sort | tr '\n' ',' | sed "s/,$//g")
    R2=$(find "$FQDIR"/* | grep -P "/$TAG[\\/\\._]" | grep "R2\\.f.*q" | sort | tr '\n' ',' | sed "s/,$//g")
  elif [[ $(find "$FQDIR"/* | grep -P "/$TAG[\\/\\._]" | grep "_R1_.*\\.f.*q") != "" ]]; then
    R1=$(find "$FQDIR"/* | grep -P "/$TAG[\\/\\._]" | grep "_R1_.*\\.f.*q" | sort | tr '\n' ',' | sed "s/,$//g")
    R2=$(find "$FQDIR"/* | grep -P "/$TAG[\\/\\._]" | grep "_R2_.*\\.f.*q" | sort | tr '\n' ',' | sed "s/,$//g")
  elif [[ $(find "$FQDIR"/* | grep -P "/$TAG[\\/\\._]" | grep "_f1\\.f.*q") != "" ]]; then
    R1=$(find "$FQDIR"/* | grep -P "/$TAG[\\/\\._]" | grep "_f1\\.f.*q" | sort | tr '\n' ',' | sed "s/,$//g")
    R2=$(find "$FQDIR"/* | grep -P "/$TAG[\\/\\._]" | grep "_r2\\.f.*q" | sort | tr '\n' ',' | sed "s/,$//g")
  else
    >&2 echo "ERROR: No appropriate fastq files were found! Please check file formatting, and check if you have set the right FQDIR."
    exit 1
  fi

  echo "$R1|$R2"
}

# Checks for GZIP or BZIP2 compression and returns the appropriate command.
# Args:
#   $1: FQDIR - The directory containing FASTQ files.
#   $2: TAG - The sample identifier.
# Returns:
#   A pipe-separated string "GZIP_COMMAND|Z_COMMAND".
function check_compression() {
  local FQDIR=$1
  local TAG=$2
  local GZIP=""
  local ZCMD="cat" # Default to no compression

  if [[ $(find "$FQDIR"/* | grep -P "$TAG[\\/\\._]" | grep "\\.gz$") != "" ]]; then
    GZIP="--readFilesCommand zcat"
    ZCMD="zcat"
  elif [[ $(find "$FQDIR"/* | grep -P "$TAG[\\/\\._]" | grep "\\.bz2$") != "" ]]; then
    GZIP="--readFilesCommand bzcat"
    ZCMD="bzcat"
  fi

  echo "$GZIP|$ZCMD"
}

# Extracts a random subsample of reads for testing chemistry.
# Args:
#   $1: R1 - Comma-separated list of Read 1 FASTQ files.
#   $2: R2 - Comma-separated list of Read 2 FASTQ files.
#   $3: ZCMD - The command to decompress files (e.g., zcat, bzcat, cat).
#   $4: CMD - Optional prefix command (e.g., for singularity).
function extract_test_reads() {
  local R1=$1
  local R2=$2
  local ZCMD=$3
  local CMD=${4:-""}
  local COUNT=0

  for i in $(echo "$R1" | tr ',' ' '); do
    $ZCMD "$i" | head -4000000 > "$COUNT.R1_head" &
    COUNT=$((COUNT+1))
  done
  wait

  COUNT=0
  for i in $(echo "$R2" | tr ',' ' '); do
    $ZCMD "$i" | head -4000000 > "$COUNT.R2_head" &
    COUNT=$((COUNT+1))
  done
  wait

  # Use the same random seed to select corresponding reads from R1 and R2.
  cat *.R1_head | $CMD seqtk sample -s100 - 200000 > test.R1.fastq &
  cat *.R2_head | $CMD seqtk sample -s100 - 200000 > test.R2.fastq &
  wait
  rm *.R1_head *.R2_head
}

# Determines the correct 10x barcode whitelist and read lengths.
# Args:
#   $1: WL - Path to the directory containing whitelist files.
# Returns:
#   A pipe-separated string: "BC|NBC1|NBC2|NBC3|NBC4|NBC5|NBCA|R1LEN|R2LEN|R1DIS"
function determine_barcode_whitelist() {
    local WL=$1
    local BC=""

    # Elucidate the right barcode whitelist to use.
    local NBC1=$(cat test.R1.fastq | awk 'NR%4==2' | cut -c-14 | grep -F -f "$WL/737K-april-2014_rc.txt" | wc -l)
    local NBC2=$(cat test.R1.fastq | awk 'NR%4==2' | cut -c-16 | grep -F -f "$WL/737K-august-2016.txt" | wc -l)
    local NBC3=$(cat test.R1.fastq | awk 'NR%4==2' | cut -c-16 | grep -F -f "$WL/3M-february-2018.txt" | wc -l)
    local NBC4=$(cat test.R1.fastq | awk 'NR%4==2' | cut -c-16 | grep -F -f "$WL/3M-3pgex-may-2023.txt" | wc -l)
    local NBC5=$(cat test.R1.fastq | awk 'NR%4==2' | cut -c-16 | grep -F -f "$WL/3M-5pgex-jan-2023.txt" | wc -l)
    local NBCA=$(cat test.R1.fastq | awk 'NR%4==2' | cut -c-16 | grep -F -f "$WL/737K-arc-v1.txt" | wc -l)

    # Get average read lengths and distribution.
    local R1LEN=$(cat test.R1.fastq | awk 'NR%4==2' | awk '{sum+=length($0)} END {printf "%d\n",sum/NR+0.5}')
    local R2LEN=$(cat test.R2.fastq | awk 'NR%4==2' | awk '{sum+=length($0)} END {printf "%d\n",sum/NR+0.5}')
    local R1DIS=$(cat test.R1.fastq | awk 'NR%4==2' | awk '{print length($0)}' | sort | uniq -c | wc -l)

    if (( NBC3 > 50000 )); then
        BC="$WL/3M-february-2018.txt"
    elif (( NBC2 > 50000 )); then
        BC="$WL/737K-august-2016.txt"
    elif (( NBCA > 50000 )); then
        BC="$WL/737K-arc-v1.txt"
    elif (( NBC1 > 50000 )); then
        BC="$WL/737K-april-2014_rc.txt"
    elif (( NBC4 > 50000 )); then
        BC="$WL/3M-3pgex-may-2023.txt"
    elif (( NBC5 > 50000 )); then
        BC="$WL/3M-5pgex-jan-2023.txt"
    else
        >&2 echo "ERROR: No whitelist has matched a random selection of 200,000 barcodes! Match counts: $NBC1 (v1), $NBC2 (v2), $NBC3 (v3), $NBC4 (v4-3p), $NBC5 (v4-5p), $NBCA (multiome)."
        exit 1
    fi

    echo "$BC|$NBC1|$NBC2|$NBC3|$NBC4|$NBC5|$NBCA|$R1LEN|$R2LEN|$R1DIS"
}

# Checks read lengths and determines barcode/UMI parameters.
# Args:
#   $1: R1DIS, $2: R1LEN, $3: R2LEN, $4: BC, $5: WL
# Returns:
#   A pipe-separated string: "PAIRED|CBLEN|UMILEN"
function check_read_lengths() {
    local R1DIS=$1
    local R1LEN=$2
    local R2LEN=$3
    local BC=$4
    local WL=$5
    local PAIRED=False
    local CBLEN
    local UMILEN

    # Check for potential issues with read lengths.
    if (( R1DIS > 1 && R1LEN <= 30 )); then
        >&2 echo "ERROR: Read 1 (barcode) has varying length; possibly quality-trimmed. Please check the fastq files."
        exit 1
    elif (( R1LEN < 24 )); then
        >&2 echo "ERROR: Read 1 (barcode) is less than 24 bp in length. Please check the fastq files."
        exit 1
    elif (( R2LEN < 40 )); then
        >&2 echo "ERROR: Read 2 (biological read) is less than 40 bp in length. Please check the fastq files."
        exit 1
    fi

    if (( R1LEN > 50 )); then
        PAIRED=True
    fi

    if [[ $BC == "$WL/3M-february-2018.txt" || $BC == "$WL/737K-arc-v1.txt" || $BC == "$WL/3M-3pgex-may-2023.txt" || $BC == "$WL/3M-5pgex-jan-2023.txt" ]]; then
        CBLEN=16
        UMILEN=12
    elif [[ $BC == "$WL/737K-august-2016.txt" ]]; then
        CBLEN=16
        UMILEN=10
    elif [[ $BC == "$WL/737K-april-2014_rc.txt" ]]; then
        CBLEN=14
        UMILEN=10
    fi

    # Failsafe for short R1 reads with v3 chemistry.
    if (( CBLEN + UMILEN > R1LEN )); then
        local NEWUMI=$((R1LEN - CBLEN))
        local BCUMI=$((UMILEN + CBLEN))
        >&2 echo "WARNING: Read 1 length ($R1LEN) is less than the sum of barcode and UMI ($BCUMI). Changing UMI from $UMILEN to $NEWUMI!"
        UMILEN=$NEWUMI
    elif (( CBLEN + UMILEN < R1LEN )); then
        local BCUMI=$((UMILEN + CBLEN))
        >&2 echo "WARNING: Read 1 length ($R1LEN) is more than the sum of barcode and UMI ($BCUMI)."
    fi

    echo "$PAIRED|$CBLEN|$UMILEN"
}


# Determines strand specificity by running test alignments.
# Args:
#   $1: REF, $2: CBLEN, $3: UMILEN, $4: CPUS, $5: BC, $6: PAIRED
#   $7: CMD - Optional prefix command.
# Returns:
#   A pipe-separated string: "STRAND|PCTFWD|PCTREV|PAIRED"
function determine_strand_specificity() {
    local REF=$1
    local CBLEN=$2
    local UMILEN=$3
    local CPUS=$4
    local BC=$5
    local PAIRED=$6
    local CMD=${7:-""}
    local STRAND=Forward

    $CMD STAR --runThreadN "$CPUS" --genomeDir "$REF" --readFilesIn test.R2.fastq test.R1.fastq --runDirPerm All_RWX --outSAMtype None \
        --soloType CB_UMI_Simple --soloCBwhitelist "$BC" --soloBarcodeReadLength 0 --soloCBlen "$CBLEN" --soloUMIstart $((CBLEN + 1)) \
        --soloUMIlen "$UMILEN" --soloStrand Forward --genomeLoad LoadAndKeep \
        --soloUMIdedup 1MM_CR --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts --soloUMIfiltering MultiGeneUMI_CR \
        --soloCellFilter EmptyDrops_CR --clipAdapterType CellRanger4 --outFilterScoreMin 30 \
        --soloFeatures Gene GeneFull --soloOutFileNames test_forward/ features.tsv barcodes.tsv matrix.mtx &> /dev/null

    $CMD STAR --runThreadN "$CPUS" --genomeDir "$REF" --readFilesIn test.R2.fastq test.R1.fastq --runDirPerm All_RWX --outSAMtype None \
        --soloType CB_UMI_Simple --soloCBwhitelist "$BC" --soloBarcodeReadLength 0 --soloCBlen "$CBLEN" --soloUMIstart $((CBLEN + 1)) \
        --soloUMIlen "$UMILEN" --soloStrand Reverse --genomeLoad LoadAndKeep \
        --soloUMIdedup 1MM_CR --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts --soloUMIfiltering MultiGeneUMI_CR \
        --soloCellFilter EmptyDrops_CR --clipAdapterType CellRanger4 --outFilterScoreMin 30 \
        --soloFeatures Gene GeneFull --soloOutFileNames test_reverse/ features.tsv barcodes.tsv matrix.mtx &> /dev/null

    local PCTFWD=$(grep "Reads Mapped to GeneFull: Unique GeneFull" test_forward/GeneFull/Summary.csv | awk -F "," '{printf "%d\n",$2*100+0.5}')
    local PCTREV=$(grep "Reads Mapped to GeneFull: Unique GeneFull" test_reverse/GeneFull/Summary.csv | awk -F "," '{printf "%d\n",$2*100+0.5}')

    if (( PCTREV > PCTFWD )); then
        STRAND=Reverse
    fi

    if (( PCTREV < 50 && PCTFWD < 50 )); then
        >&2 echo "WARNING: Low percentage of reads mapping to GeneFull: forward = $PCTFWD, reverse = $PCTREV"
    fi

    # If a paired-end experiment turns out to be 3', process it as single-end.
    if [[ $STRAND == "Forward" && $PAIRED == "True" ]]; then
        PAIRED=False
    fi

    echo "$STRAND|$PCTFWD|$PCTREV|$PAIRED"
}

# Writes the final run configuration to a log file.
# Args: Takes 15 arguments representing all configuration parameters.
function write_config_file() {
    local TAG=$1
    local PAIRED=$2
    local STRAND=$3
    local PCTFWD=$4
    local PCTREV=$5
    local BC=$6
    local NBC1=$7
    local NBC2=$8
    local NBC3=$9
    local NBCA=${10}
    local NBC4=${11}
    local NBC5=${12}
    local CBLEN=${13}
    local UMILEN=${14}
    local GZIP=${15}
    local R1=${16}
    local R2=${17}

    echo "Done setting up the STARsolo run; here are final processing options:"
    echo "============================================================================="
    (
        echo "Sample: $TAG"
        echo "Paired-end mode: $PAIRED"
        echo "Strand (Forward = 3', Reverse = 5'): $STRAND, %reads mapped to GeneFull: forward = $PCTFWD , reverse = $PCTREV"
        echo "CB whitelist: $BC, matches out of 200,000: $NBC3 (v3), $NBC2 (v2), $NBC1 (v1), $NBC4 (v4-3p), $NBC5 (v4-5p), $NBCA (multiome)"
        echo "CB length: $CBLEN"
        echo "UMI length: $UMILEN"
        echo "GZIP: $GZIP"
        echo "-----------------------------------------------------------------------------"
        echo "Read 1 files: $R1"
        echo "-----------------------------------------------------------------------------"
        echo "Read 2 files: $R2"
        echo "-----------------------------------------------------------------------------"
    ) | tee strand.txt
}

# Executes the main STARsolo alignment.
# Args: Takes 13 arguments for all STARsolo parameters.
function run_starsolo() {
    local PAIRED=$1
    local R1=$2
    local R2=$3
    local GZIP=$4
    local BAM=$5
    local BC=$6
    local CBLEN=$7
    local UMILEN=$8
    local STRAND=$9
    local REF=${10}
    local CPUS=${11}
    local SOLOFILENAMES=${12}
    local CMD=${13:-""}

    if [[ $PAIRED == "True" ]]; then
        # Note the R1/R2 order and --soloStrand Forward for 5' paired-end.
        $CMD STAR --runThreadN "$CPUS" --genomeDir "$REF" --readFilesIn "$R1" "$R2" --runDirPerm All_RWX $GZIP $BAM --soloBarcodeMate 1 --clip5pNbases 39 0 \
            --soloType CB_UMI_Simple --soloCBwhitelist "$BC" --soloCBstart 1 --soloCBlen "$CBLEN" --soloUMIstart $((CBLEN + 1)) --soloUMIlen "$UMILEN" --soloStrand Forward \
            --soloUMIdedup 1MM_CR --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts --soloUMIfiltering MultiGeneUMI_CR \
            --soloCellFilter EmptyDrops_CR --outFilterScoreMin 30 --genomeLoad LoadAndRemove \
            --soloFeatures Gene GeneFull Velocyto --soloOutFileNames "$SOLOFILENAMES" features.tsv barcodes.tsv matrix.mtx --soloMultiMappers EM
    else
        $CMD STAR --runThreadN "$CPUS" --genomeDir "$REF" --readFilesIn "$R2" "$R1" --runDirPerm All_RWX $GZIP $BAM \
            --soloType CB_UMI_Simple --soloCBwhitelist "$BC" --soloBarcodeReadLength 0 --soloCBlen "$CBLEN" --soloUMIstart $((CBLEN + 1)) --soloUMIlen "$UMILEN" --soloStrand "$STRAND" \
            --soloUMIdedup 1MM_CR --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts --soloUMIfiltering MultiGeneUMI_CR \
            --soloCellFilter EmptyDrops_CR --clipAdapterType CellRanger4 --outFilterScoreMin 30 --genomeLoad LoadAndRemove \
            --soloFeatures Gene GeneFull Velocyto --soloOutFileNames "$SOLOFILENAMES" features.tsv barcodes.tsv matrix.mtx --soloMultiMappers EM
    fi
}

# Finalizes the output files.
# Args:
#   $1: TAG - The sample identifier.
#   $2: SOLOFILENAMES - The output directory for solo features.
#   $3: CMD - Optional prefix command (e.g., for singularity).
function process_output_files() {
    local CMD=${1:-""}

    if [[ -s "Aligned.sortedByCoord.out.bam" ]]; then
        $CMD samtools index -@16 Aligned.sortedByCoord.out.bam
    fi

    $CMD pbzip2 -9 "Unmapped.out.mate1" &
    $CMD pbzip2 -9 "Unmapped.out.mate2" &
    wait

    rm -rf test.R?.fastq test_forward test_reverse

    find . -name "*.tsv" -exec gzip {} +
    find . -name "*.mtx" -exec gzip {} +

    wait
    echo "ALL DONE!"
}


# --- Main Workflow Orchestrator ---

function starsolo_10x() {
    # --- 1. Setup and Argument Parsing ---
    local FQDIR=$1
    local TAG=$2
    local CPUS=$3
    local REF=$4
    local WL=$5
    local BAM=$6
    local CMD=${7:-""}

    local FQDIR_ABS
    FQDIR_ABS=$(readlink -f "$FQDIR")
    mkdir -p "$TAG" && cd "$TAG" || exit

    # --- 2. Find FASTQ Files ---
    local fastq_files
    fastq_files=$(find_fastq_files "$FQDIR_ABS" "$TAG")
    local R1
    local R2
    IFS='|' read -r R1 R2 <<< "$fastq_files"

    # --- 3. Determine Compression and Chemistry ---
    local compression_info
    compression_info=$(check_compression "$FQDIR_ABS" "$TAG")
    local GZIP
    local ZCMD
    IFS='|' read -r GZIP ZCMD <<< "$compression_info"

    extract_test_reads "$R1" "$R2" "$ZCMD" "$CMD"

    local whitelist_info
    whitelist_info=$(determine_barcode_whitelist "$WL")
    local BC NBC1 NBC2 NBC3 NBC4 NBC5 NBCA R1LEN R2LEN R1DIS
    IFS='|' read -r BC NBC1 NBC2 NBC3 NBC4 NBC5 NBCA R1LEN R2LEN R1DIS <<< "$whitelist_info"

    # --- 4. Check Reads and Determine Parameters ---
    local read_len_info
    read_len_info=$(check_read_lengths "$R1DIS" "$R1LEN" "$R2LEN" "$BC" "$WL")
    local PAIRED CBLEN UMILEN
    IFS='|' read -r PAIRED CBLEN UMILEN <<< "$read_len_info"

    local strand_info
    strand_info=$(determine_strand_specificity "$REF" "$CBLEN" "$UMILEN" "$CPUS" "$BC" "$PAIRED" "$CMD")
    local STRAND PCTFWD PCTREV
    # PAIRED is returned again as the function can modify it
    IFS='|' read -r STRAND PCTFWD PCTREV PAIRED <<< "$strand_info"

    # --- 5. Log, Run, and Cleanup ---
    write_config_file "$TAG" "$PAIRED" "$STRAND" "$PCTFWD" "$PCTREV" "$BC" "$NBC1" "$NBC2" "$NBC3" "$NBCA" "$NBC4" "$NBC5" "$CBLEN" "$UMILEN" "$GZIP" "$R1" "$R2"

    local SOLOFILENAMES="output/"
    run_starsolo "$PAIRED" "$R1" "$R2" "$GZIP" "$BAM" "$BC" "$CBLEN" "$UMILEN" "$STRAND" "$REF" "$CPUS" "$SOLOFILENAMES" "$CMD"

    process_output_files "$CMD"
}

# --- Script Entrypoint ---

function main () {
    if (( $# != 3 )); then
        >&2 echo "Usage: ./starsolo_10x_auto.sh <fastq_dir> <sample_id> <specie>"
        >&2 echo "(make sure you set the correct REF, WL, and BAM variables below)"
        >&2 echo "Specie can be 'human' or 'mouse' (or any other directory in /nfs/cellgeni/STAR/)"
        >&2 echo -e "You need to have the following software installed to run this script: \n- star_version=2.7.10a_alpha_220818\n- samtools_version=1.15.1\n- bbmap_version=38.97\n- rsem_version=1.3.3"
        >&2 echo "Use /nfs/cellgeni/singularity/images/reprocess_10x.sif if you work on FARM cluster."
        exit 1
    fi

    local FQDIR=$1
    local TAG=$2
    local SPECIE=$3

    # --- Configuration ---
    # Using local variables to keep all state contained.
    local SIF="/nfs/cellgeni/singularity/images/reprocess_10x.sif"
    local CMD="singularity run --bind /nfs,/lustre,/software $SIF"
    local CPUS=16
    local REF="/nfs/cellgeni/STAR/${SPECIE}/2020A/index"
    local WL="/nfs/cellgeni/STAR/whitelists"
    # Choose one of the two options, depending on whether you need a BAM file
    # local BAM="--outSAMtype BAM SortedByCoordinate --outBAMsortingBinsN 500 --limitBAMsortRAM 60000000000 --outSAMunmapped Within --outMultimapperOrder Random --runRNGseed 1 --outSAMattributes NH HI AS nM CB UB CR CY UR UY GX GN"
    local BAM="--outSAMtype None --outReadsUnmapped Fastx"

    starsolo_10x "$FQDIR" "$TAG" "$CPUS" "$REF" "$WL" "$BAM" "$CMD"
}

# This construct ensures that main() is called only when the script is executed directly.
if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
    main "$@"
fi