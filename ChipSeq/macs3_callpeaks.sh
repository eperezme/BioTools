#!/bin/bash
export REFGENIE=/zpool_1TB/data/cram/genome_config.yaml
export REF_PATH=/zpool_1TB/data/cram/cache/%2s/%2s/%s:http://www.ebi.ac.uk/ena/cram/md5/%s
export REF_CACHE=/zpool_1TB/data/cram/cache/%2s/%2s/%s
# Usage: ./script.sh -s samplesheet.csv [-o outdir] [-g genome_size] [-q qvalue] [-b] [-p cores]
# spell-checker: disable samplesheet qvalue gsize qval bams BAMPE bedpe callpeak bedtools samtools annotatepeaks

# Default parameters
samplesheet=""
OUTDIR=./MACS3
gsize=2e9
qval=0.05
broad_flag=""
parallel_cores=4

# Additional default parameters
output_format="csv"

# Parse command-line arguments
while getopts "s:o:g:q:bp:G:f:" opt; do
    case ${opt} in
        s ) samplesheet="$OPTARG" ;;
        o ) OUTDIR="$OPTARG" ;;
        g ) gsize="$OPTARG" ;;
        q ) qval="$OPTARG" ;;
        b ) broad_flag="--broad" ;;
        p ) parallel_cores="$OPTARG" ;;
        G ) genome="$OPTARG" ;;
        f ) output_format="$OPTARG" ;;
        * ) echo "Usage: $0 -s <samplesheet.csv> [-o outdir] [-g genome_size] [-q qval] [-b] [-p cores] [-G genome] [-f csv|tsv]"; exit 1 ;;
    esac
done

# Check if samplesheet is provided and exists
if [[ -z "$samplesheet" ]]; then
    echo "Error: Sample sheet is required."
    echo "Usage: $0 -s <samplesheet.csv> [-o outdir] [-g genome_size] [-q qval] [-b] [-p cores]"
    exit 1
fi

if [[ ! -f "$samplesheet" ]]; then
    echo "Error: Sample sheet file '$samplesheet' not found."
    exit 1
fi

# Activate the conda environment for macs3
echo "Activating conda environment for macs3..."
source /zpool_1TB/lib/miniconda3/etc/profile.d/conda.sh
conda activate macs3

# Search for genome
echo " - Searching for genome"
echo " - Detecting map aligner for the first BAM file in the samplesheet"
first_bam=$(awk -F',' 'NR==2 {print $4}' "$samplesheet") # Get the first BAM file from the samplesheet
map_aligner=$(samtools view -H "$first_bam" | grep -m 1 '@PG' | awk '{print $2}' | cut -d':' -f2)
echo " - Map aligner detected: $map_aligner"
if [[ "$map_aligner" == "bwa" ]]; then
    map_aligner="BWAIndex"
elif [[ "$map_aligner" == "bowtie2" ]]; then
    map_aligner="Bowtie2Index"
elif [[ "$map_aligner" == "hisat2" ]]; then
    map_aligner="hisat2_index"
elif [[ "$map_aligner" == "star" ]]; then
    map_aligner="star_index"
elif [[ "$map_aligner" == "bismark" ]]; then
    map_aligner="bismark_bt2_index"
else
    echo " - Unknown map aligner: $map_aligner"
    exit 1
fi
echo " - Checking index file for reference genome"
ref=$(refgenie seek ${genome}/${map_aligner}) # genome
echo $ref # debug
# Check if the reference genome is already downloaded
if [[ -z "$ref" ]]; then
    echo " - Genome not found locally, downloading..."
    refgenie pull $genome/${map_aligner}
    ref=$(refgenie seek $genome/${map_aligner})
    echo " - Genome downloaded: $ref"
else
    echo " - Found genome: $ref"
fi

# Set up output directories
echo "Creating output directories..."
MERGED_PEAKS_DIR="$OUTDIR/MergedPeaks"
MERGED_BAM_DIR="$OUTDIR/MergedBAM"
mkdir -p "$OUTDIR" "$MERGED_PEAKS_DIR" "$MERGED_BAM_DIR/Peaks"

# Initialize associative arrays
declare -A factor_bams
declare -A factor_controls
declare -A factor_peaks
declare -A sample_controls
declare -A sample_is_paired
declare -A sample_format
factors=()
sample_count=0

# Parse the samplesheet line by line
echo "Parsing samplesheet..."
while IFS="," read -r SampleID Factor Replicate bamReads ControlID bamControl Peaks PeakCaller Tissue Condition PairedEnd; do
    [[ "$SampleID" == "SampleID" ]] && continue
    ((sample_count++))
    factor_bams[$Factor]+="$bamReads "
    factor_peaks[$Factor]+="$Peaks "
    sample_controls[$bamReads]="$bamControl"
    factor_controls[$Factor]+="$bamControl "

    # Detect format by checking if the BAM file contains paired-end reads
    if [[ "$bamReads" == *.bam ]]; then
        if samtools view -c -f 1 "$bamReads" > /dev/null 2>&1 && [ $(samtools view -f 1 "$bamReads" | wc -l) -gt 0 ]; then
            sample_format[$bamReads]="BAMPE"
        else
            sample_format[$bamReads]="BAM"
        fi
    elif [[ "$bamReads" == *.bed ]]; then
        sample_format[$bamReads]="BED"
    elif [[ "$bamReads" == *.bedpe ]]; then
        sample_format[$bamReads]="BEDPE"
    else
        sample_format[$bamReads]="BAM"
    fi

    [[ " ${factors[*]} " =~ " $Factor " ]] || factors+=("$Factor")
done < "$samplesheet"

echo "Found $sample_count samples across ${#factors[@]} factors: ${factors[*]}"


# Build MACS3 commands
echo "Generating MACS3 commands..."
macs3_cmds=()
for factor in "${factors[@]}"; do
    for bam in ${factor_bams[$factor]}; do
        name=$(basename "$bam" .sorted.bam)
        control="${sample_controls[$bam]}"
        format="${sample_format[$bam]}"
        echo "Processing sample: $name (format: $format)"

        if [[ ! -f "$OUTDIR/${name}_peaks.narrowPeak" ]]; then
            if [ -n "$control" ]; then
                echo "Using control: $control"
                macs3_cmds+=("macs3 callpeak -t $bam -c $control -n $name -f $format --outdir $OUTDIR -g $gsize -q $qval $broad_flag")
            else
                macs3_cmds+=("macs3 callpeak -t $bam -n $name -f $format --outdir $OUTDIR -g $gsize -q $qval $broad_flag")
            fi
        else
            echo "Skipping $name - peak file already exists."
        fi
    done
done

# Run MACS3 in parallel
echo "Running MACS3 peak calling in parallel using $parallel_cores cores..."
printf "%s\n" "${macs3_cmds[@]}" | parallel -j "$parallel_cores"


# Merge peak files and call MACS3 on merged BAMs
echo "Merging peak files and calling MACS3 on merged BAMs..."
for factor in "${factors[@]}"; do
    echo "Processing factor: $factor"

    peaks_master="$MERGED_PEAKS_DIR/${factor}_peaks_master.bed"
    merged_bam="$MERGED_BAM_DIR/${factor}_merged.bam"

    # Check if merged peaks and BAM files already exist
    if [[-f "$merged_bam" ]]; then
        echo " - Skipping $factor: Merged BAM file already exist."
        continue
    fi

    echo "Merging peaks for factor: $factor"
    peaks_concat="$MERGED_PEAKS_DIR/${factor}_peaks_concat.tmp"
    peaks_sorted="$MERGED_PEAKS_DIR/${factor}_peaks_sorted.bed"
    peaks_merged="$MERGED_PEAKS_DIR/${factor}_peaks_merged.bed"

    cat ${factor_peaks[$factor]} > "$peaks_concat"
    sort -k1,1 -k2,2n "$peaks_concat" > "$peaks_sorted"
    bedtools merge -i "$peaks_sorted" -c 4,5,6,7,8,9,10 \
                   -o first,sum,distinct,mean,mean,mean,mean > "$peaks_merged"
    sort -k1,1 -k2,2n "$peaks_merged" > "$peaks_master"

    echo "Merging BAM files for factor: $factor"
    samtools merge -f -o "$merged_bam" ${factor_bams[$factor]}
    samtools sort -o "$merged_bam" "$merged_bam"
    samtools index "$merged_bam"

    # Merge control BAMs for factor if any
    if [[ -n "${factor_controls[$factor]// }" ]]; then
        echo "Merging control BAM files for factor: $factor"
        merged_control="$MERGED_BAM_DIR/${factor}_control_merged.bam"
        samtools merge -f -o "$merged_control" ${factor_controls[$factor]}
        samtools sort -o "$merged_control" "$merged_control"
        samtools index "$merged_control"
        control_arg="-c $merged_control"
    else
        control_arg=""
    fi
done

for factor in "${factors[@]}"; do
    merged_bam="$MERGED_BAM_DIR/${factor}_merged.bam"
    merged_control="$MERGED_BAM_DIR/${factor}_control_merged.bam"

    # Check if peak file already exists
    if [[ -f "$MERGED_BAM_DIR/Peaks/${factor}_peaks.narrowPeak" ]]; then
        echo " - Skipping $factor: Peak file already exists."
        continue
    fi
    # Detect paired-end or single-end format
    if samtools view -c -f 1 "$merged_bam" > /dev/null 2>&1 && [ $(samtools view -f 1 "$merged_bam" | wc -l) -gt 0 ]; then
        format="BAMPE"
    else
        format="BAM"
    fi

    echo "Calling peaks on merged BAM for factor: $factor"
    macs3 callpeak -t "$merged_bam" $control_arg -n "$factor" -f "$format" --outdir "$MERGED_BAM_DIR/Peaks" -g "$gsize" -q "$qval" $broad_flag
    echo " - Peak calling completed for factor: $factor"
done

for factor in "${factors[@]}"; do
    merged_bam="$MERGED_BAM_DIR/${factor}_merged.bam"
    # Convert merged BAM to CRAM for storage efficiency
    echo " - Converting merged BAM to CRAM for factor $factor"
    merged_cram="$MERGED_BAM_DIR/${factor}_merged.cram"
    samtools view -C -T $ref -o "$merged_cram" "$merged_bam"
    # echo " - Deleting merged BAM file for factor $factor"
    # rm "$merged_bam"
    # rm "${merged_bam}.bai"
done

# Intersect peaks between all factor combinations
for (( i=0; i<${#factors[@]}; i++ )); do
        for (( j=i+1; j<${#factors[@]}; j++ )); do
                A="${factors[i]}"
                B="${factors[j]}"
                echo "Intersecting peaks: $A vs $B"
# TODO change name of intersect output files to know wich ones are the interesting ones
                bedtools intersect -a "$MERGED_PEAKS_DIR/${A}_peaks_master.bed" \
                                                     -b "$MERGED_PEAKS_DIR/${B}_peaks_master.bed" \
                                                     -wa > "$MERGED_PEAKS_DIR/${A}_vs_${B}_unique_peaks.bed"
                sort -k1,1 -k2,2n "$MERGED_PEAKS_DIR/${A}_vs_${B}_unique_peaks.bed" > "$MERGED_PEAKS_DIR/${A}_vs_${B}_unique_peaks.bed"

                bedtools intersect -a "$MERGED_BAM_DIR/Peaks/${A}_peaks.narrowPeak" \
                                                     -b "$MERGED_BAM_DIR/Peaks/${B}_peaks.narrowPeak" \
                                                     -wa > "$MERGED_BAM_DIR/Peaks/${A}_vs_${B}_intersected_peaks.bed"
                sort -k1,1 -k2,2n "$MERGED_BAM_DIR/Peaks/${A}_vs_${B}_intersected_peaks.bed" > "$MERGED_BAM_DIR/Peaks/${A}_vs_${B}_unique_peaks.bed"

                echo "Annotating intersected peaks for $A vs $B"
                # Check if bed files are empty before annotating
                if [[ ! -s "$MERGED_PEAKS_DIR/${A}_vs_${B}_unique_peaks.bed" && ! -s "$MERGED_BAM_DIR/Peaks/${A}_vs_${B}_unique_peaks.bed" ]]; then
                    echo "Both BED files are empty for $A vs $B. Skipping annotation."
                    continue
                fi

                if [[ -s "$MERGED_PEAKS_DIR/${A}_vs_${B}_unique_peaks.bed" ]]; then
                    echo "Annotating non-empty intersected peaks for $A vs $B (MERGED_PEAKS_DIR)"
                    ./annotatepeaks.R -i "$MERGED_PEAKS_DIR/${A}_vs_${B}_unique_peaks.bed" -g ${genome} -f $output_format
                fi

                if [[ -s "$MERGED_BAM_DIR/Peaks/${A}_vs_${B}_unique_peaks.bed" ]]; then
                    echo "Annotating non-empty intersected peaks for $A vs $B (MERGED_BAM_DIR)"
                    ./annotatepeaks.R -i "$MERGED_BAM_DIR/Peaks/${A}_vs_${B}_unique_peaks.bed" -g ${genome} -f $output_format
                fi

        done
done

# TODO: Create the samplesheet_merged.csv file
# Create the samplesheet_merged.csv file
echo "Creating samplesheet_merged.csv..."
merged_samplesheet="$MERGED_BAM_DIR/samplesheet_merged.csv"

# Extract the header from the original samplesheet
header=$(head -n 1 "$samplesheet")
echo "$header" > "$merged_samplesheet"

# Add merged BAM and peak information for each factor
for factor in "${factors[@]}"; do
    merged_bam="$MERGED_BAM_DIR/${factor}_merged.bam"
    peaks_master="$MERGED_PEAKS_DIR/${factor}_peaks_master.bed"
## Check if merged_control exists
    merged_control="$MERGED_BAM_DIR/${factor}_control_merged.bam"

   if [[ ! -f "$merged_control" ]]; then
        merged_control="NA"
        control_id="NA"
    else
        control_id=${factor}_control
    fi
    if [[ -f "$merged_bam" && -f "$peaks_master" ]]; then
        echo "$factor,$factor,1,$merged_bam,$control_id,$merged_control,$peaks_master,narrowPeak,NA,NA" >> "$merged_samplesheet"
    else
        echo "Warning: Missing merged BAM or peaks file for factor $factor. Skipping entry."
    fi
done

echo "samplesheet_merged.csv created at $merged_samplesheet"



# Cleanup
echo "Cleaning up intermediate files..."
rm "$MERGED_PEAKS_DIR"/*.tmp
rm "$MERGED_PEAKS_DIR"/*_sorted.bed
rm "$MERGED_PEAKS_DIR"/*_merged.bed
rm "$MERGED_BAM_DIR/Peaks/"*_intersected_peaks.bed