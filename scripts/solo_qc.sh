#!/bin/bash
# =============================================================================
# solo_qc.sh — Aggregate QC statistics across STARsolo runs
# =============================================================================
# Run from the parent directory that contains one subdirectory per sample
# (each produced by a starsolo run).
#
# Usage:
#   starsolo qc                   # via the CLI
#   ./scripts/solo_qc.sh          # standalone
#
# Pipe through `column -t` for pretty-printed output:
#   starsolo qc | column -t
# =============================================================================

# --- Check completion status ------------------------------------------------

>&2 echo "Checking that all STARsolo jobs went to completion …"
for i in */; do
    i="${i%/}"
    if [[ -d "$i/output" && -s "$i/Log.final.out" ]]; then
        if [[ -d "$i/_STARtmp" ]]; then
            >&2 echo "WARNING: Sample $i did not run to completion: _STARtmp is still present!"
        fi
    fi
done

# --- Check for cross-contamination -----------------------------------------

>&2 echo "Checking potential sample cross-contamination …"
for i in */; do
    i="${i%/}"
    if [[ -s "$i/output/Gene/filtered/barcodes.tsv.gz" ]]; then
        for j in */; do
            j="${j%/}"
            if [[ -s "$j/output/Gene/filtered/barcodes.tsv.gz" && "$i" != "$j" ]]; then
                if [[ ! -s "sorted.$i.barcodes.tsv" ]]; then
                    zcat "$i/output/Gene/filtered/barcodes.tsv.gz" | sort > "sorted.$i.barcodes.tsv" &
                fi
                if [[ ! -s "sorted.$j.barcodes.tsv" ]]; then
                    zcat "$j/output/Gene/filtered/barcodes.tsv.gz" | sort > "sorted.$j.barcodes.tsv" &
                fi
                wait
                N1=$(wc -l < "sorted.$i.barcodes.tsv")
                N2=$(wc -l < "sorted.$j.barcodes.tsv")
                MIN=$(( N1 <= N2 ? N1 : N2 ))
                COMM=$(comm -12 "sorted.$i.barcodes.tsv" "sorted.$j.barcodes.tsv" | wc -l)
                PCT=$(echo "$COMM" | awk -v v="$MIN" '{printf "%d\n",100*$1/v+0.5}')
                if (( PCT >= 20 )); then
                    >&2 echo "WARNING: Samples $i ($N1 barcodes) and $j ($N2 barcodes) share $COMM ($PCT%) barcodes — higher than expected!"
                fi
            fi
        done
    elif [[ -d "$i/output" && -s "$i/Log.out" ]]; then
        >&2 echo "Sample $i: no filtered output detected (probably plate-based)."
    fi
done

rm -f sorted.*.barcodes.tsv

# --- Extract and output statistics ------------------------------------------

>&2 echo "Extracting STARsolo stats …"
>&2 echo

echo -e "Sample\tRd_all\tRd_in_cells\tFrc_in_cells\tUMI_in_cells\tCells\tMed_nFeature\tGood_BC\tWL\tSpecies\tPaired\tStrand\tall_u+m\tall_u\texon_u+m\texon_u\tfull_u+m\tfull_u"

for i in */; do
    i="${i%/}"

    if [[ -d "$i/output/Gene/filtered/" && -s "$i/Log.out" ]]; then
        # --- Droplet-based ---
        PAIRED="Single"
        if grep -q "clip5pNbases 39 0" "$i/Log.out"; then
            PAIRED="Paired"
        fi

        REF="Other"
        if grep "^genomeDir" "$i/Log.out" | tail -n1 | grep -q "/human/"; then
            REF="Human"
        elif grep "^genomeDir" "$i/Log.out" | tail -n1 | grep -q "/mouse/"; then
            REF="Mouse"
        fi

        WL="Undef"
        if grep "^soloCBwhitelist" "$i/Log.out" | tail -n1 | grep -q "3M-february-2018.txt"; then
            WL="v3"
        elif grep "^soloCBwhitelist" "$i/Log.out" | tail -n1 | grep -q "737K-august-2016.txt"; then
            WL="v2"
        elif grep "^soloCBwhitelist" "$i/Log.out" | tail -n1 | grep -q "737K-april-2014_rc.txt"; then
            WL="v1"
        elif grep "^soloCBwhitelist" "$i/Log.out" | tail -n1 | grep -q "737K-arc-v1.txt"; then
            WL="arc"
        fi

        R1=$(grep "Number of Reads,"             "$i/output/Gene/Summary.csv"     | awk -F "," '{print $2}')
        B=$( grep "Reads With Valid Barcodes,"    "$i/output/Gene/Summary.csv"     | awk -F "," '{print $2}')
        G1=$(grep "Reads Mapped to Genome: Unique+Multiple," "$i/output/Gene/Summary.csv" | awk -F "," '{print $2}')
        G2=$(grep "Reads Mapped to Genome: Unique,"          "$i/output/Gene/Summary.csv" | awk -F "," '{print $2}')
        E1=$(grep "Reads Mapped to Gene: Unique+Multip.*e Gene,"     "$i/output/Gene/Summary.csv"     | awk -F "," '{print $2}')
        E2=$(grep "Reads Mapped to Gene: Unique Gene,"               "$i/output/Gene/Summary.csv"     | awk -F "," '{print $2}')
        F1=$(grep "Reads Mapped to GeneFull: Unique+Multip.*e GeneFull," "$i/output/GeneFull/Summary.csv" | awk -F "," '{print $2}')
        F2=$(grep "Reads Mapped to GeneFull: Unique GeneFull,"           "$i/output/GeneFull/Summary.csv" | awk -F "," '{print $2}')
        C=$( grep "Estimated Number of Cells,"    "$i/output/GeneFull/Summary.csv" | awk -F "," '{print $2}')
        R2=$(grep "Unique Reads in Cells Mapped to GeneFull," "$i/output/GeneFull/Summary.csv" | awk -F "," '{print $2}')
        CF=$(echo "$R1" | awk -v v="$R2" '{printf "%.3f\n",v/$1}')
        R3=$(grep "UMIs in Cells,"                "$i/output/GeneFull/Summary.csv" | awk -F "," '{print $2}')
        GC=$(grep "Median GeneFull per Cell,"      "$i/output/GeneFull/Summary.csv" | awk -F "," '{print $2}')
        ST=$(grep "^soloStrand" "$i/Log.out" | grep RE-DEFINED | awk '{print $2}')

        if [[ -z "$ST" ]]; then
            ST="Undef"
        elif [[ "$PAIRED" == "Paired" ]]; then
            ST="Reverse"
        fi

        F2PCT=$(echo "$F2" | awk '{printf "%d\n",$1*100+0.5}')
        if (( F2PCT <= 20 )); then
            >&2 echo "WARNING: Sample $i : GeneFull percentage ($F2PCT) is too low!"
        fi

        echo -e "$i\t$R1\t$R2\t$CF\t$R3\t$C\t$GC\t$B\t$WL\t$REF\t$PAIRED\t$ST\t$G1\t$G2\t$E1\t$E2\t$F1\t$F2"

    elif [[ -d "$i/output" && -s "$i/Log.out" ]]; then
        # --- Plate-based (Smart-seq2 etc.) ---
        PAIRED="Single"
        if grep -q "mate 2" "$i/Log.out"; then
            PAIRED="Paired"
        fi

        REF="Other"
        if grep "^genomeDir" "$i/Log.out" | tail -n1 | grep -q "/human/"; then
            REF="Human"
        elif grep "^genomeDir" "$i/Log.out" | tail -n1 | grep -q "/mouse/"; then
            REF="Mouse"
        fi

        WL="Undef"
        if grep -q "manifest" "$i/Log.out"; then
            WL="None(Smart-seq)"
        fi

        R1=$(grep "Number of Reads,"             "$i/output/Gene/Summary.csv"     | awk -F "," '{print $2}')
        B=$( grep "Reads With Valid Barcodes,"    "$i/output/Gene/Summary.csv"     | awk -F "," '{print $2}')
        G1=$(grep "Reads Mapped to Genome: Unique+Multiple," "$i/output/Gene/Summary.csv" | awk -F "," '{print $2}')
        G2=$(grep "Reads Mapped to Genome: Unique,"          "$i/output/Gene/Summary.csv" | awk -F "," '{print $2}')
        E1=$(grep "Reads Mapped to Gene: Unique+Multip.*e Gene,"     "$i/output/Gene/Summary.csv"     | awk -F "," '{print $2}')
        E2=$(grep "Reads Mapped to Gene: Unique Gene,"               "$i/output/Gene/Summary.csv"     | awk -F "," '{print $2}')
        F1=$(grep "Reads Mapped to GeneFull: Unique+Multip.*e GeneFull," "$i/output/GeneFull/Summary.csv" | awk -F "," '{print $2}')
        F2=$(grep "Reads Mapped to GeneFull: Unique GeneFull,"           "$i/output/GeneFull/Summary.csv" | awk -F "," '{print $2}')
        C=$( zcat "$i/output/Gene/raw/barcodes.tsv.gz" | wc -l)
        R2=$(grep "yessubWLmatch_UniqueFeature" "$i/output/GeneFull/Features.stats" | awk '{print $2}')
        CF=$(echo "$R1" | awk -v v="$R2" '{printf "%.3f\n",v/$1}')
        R3=$(grep "yesUMIs" "$i/output/GeneFull/Features.stats" | awk '{print $2}')
        GC=$(zcat "$i/output/GeneFull/raw/matrix.mtx.gz" | awk 'NR>3 {print $2}' | uniq -c | awk '{print $1}' | sort -n | awk '{a[i++]=$1} END {print a[int(i/2)]}')
        ST=$(grep "^soloStrand" "$i/Log.out" | grep RE-DEFINED | awk '{print $2}')

        [[ -z "$ST" ]] && ST="Undef"

        F2PCT=$(echo "$F2" | awk '{printf "%d\n",$1*100+0.5}')
        if (( F2PCT <= 20 )); then
            >&2 echo "WARNING: Sample $i : GeneFull percentage ($F2PCT) is too low!"
        fi

        echo -e "$i\t$R1\t$R2\t$CF\t$R3\t$C\t$GC\t$B\t$WL\t$REF\t$PAIRED\t$ST\t$G1\t$G2\t$E1\t$E2\t$F1\t$F2"
    fi
done
