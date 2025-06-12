#!/bin/bash

# ==============================================================================
# A script to download and process GENCODE/Ensembl data to create
# custom IGV tracks.
#
# It performs the following main steps:
#   1. Downloads Gencode GTF and MANE summary data.
#   2. Creates a map of gene symbols to transcript IDs.
#   3. Processes a source BED file to:
#      - Rename features using the gene symbol.
#      - Colour MANE transcripts purple.
#      - Add GFF-style attributes for rich IGV pop-ups.
#   4. Creates a final subset track based on a gene list.
# ==============================================================================

# --- Script Configuration ---
# Stop on any error (`-e`) and on failures in pipelines (`-o pipefail`).
set -e
set -o pipefail

# --- File and URL Variables ---
# Makes updating versions and filenames easy.
GENCODE_URL="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_48/gencode.v48.annotation.gtf.gz"
MANE_URL="https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/release_1.2/MANE.GRCh38.v1.2.summary.txt.gz"
AVENIO_GENE_LIST="avenio/avenio_genes.txt"
SOURCE_BED="beds/gencode_ucsc_knownGene_wg.bed" # from UCSC genome table 

# --- Output Files ---
# Organize output into directories.
mkdir -p data beds/final avenio/final

GENCODE_GTF_GZ="data/gencode.v48.annotation.gtf.gz"
MANE_SUMMARY_GZ="data/MANE.GRCh38.v1.2.summary.txt.gz"
GENE_TX_MAP="data/ensembl_gene_tx_map.tsv"
MANE_ENSEMBL_TX="data/mane_ensembl_tx.txt"

FINAL_ENSEMBL_BED="beds/final/ensembl_final_tagged.bed"
FINAL_AVENIO_BED="avenio/final/ensembl_avenio_subset.bed"

# --- Headers for BED Tracks ---
ENSEMBL_HEADER='track name="Ensembl v48 (MANE Purple)" description="Gencode v48 with MANE transcripts in purple" itemRgb="On" db=hg38'
AVENIO_HEADER='track name="Avenio Gene Subset (Ensembl)" description="Avenio panel genes from Ensembl v48" itemRgb="On" db=hg38'


################################################################################
### STEP 1: Download Data and Prepare Mapping Files
################################################################################

echo "--> Downloading Gencode annotations..."
wget -O "$GENCODE_GTF_GZ" "$GENCODE_URL"

echo "--> Generating Gene-to-Transcript map from GTF..."
# Use a single, efficient awk command to parse the gzipped GTF directly.
# This avoids creating a huge intermediate uncompressed file.
gunzip -c "$GENCODE_GTF_GZ" | gawk -F'\t' '
    # This pattern only runs on lines where the 3rd field is "transcript"
    $3 == "transcript" {
        
        # Now that we are correctly splitting on tabs, $9 is guaranteed
        # to be the full attribute string. This logic should now work.
        if (match($9, /gene_name "([^"]+)"/, gene) && match($9, /transcript_id "([^"]+)"/, tx)) {
            print gene[1] "\t" tx[1]
        }
    }
' > "$GENE_TX_MAP"

echo "--> Downloading MANE summary..."
wget -O "$MANE_SUMMARY_GZ" "$MANE_URL"

echo "--> Extracting MANE Ensembl transcript IDs..."
# Extract the 3rd column (Ensembl ID) from the MANE summary, skipping the header.
gunzip -c "$MANE_SUMMARY_GZ" | gawk -F'\t' -v OFS='\t' '!/^#/ {split($3, parts, ":"); print parts[1], $6, $8}'  > "$MANE_ENSEMBL_TX"


################################################################################
### STEP 2: Process Ensembl Track (Rename, Colour, and Tag in ONE pass)
################################################################################

echo "--> Processing Ensembl BED file (Rename, Colour, Tag)..."

# This single awk command replaces three separate file passes.
# It reads the mapping files first, then processes the main BED file.
# The final result, including the header, is created in one go.
{
    echo "$ENSEMBL_HEADER"
    gawk -v OFS='\t' '
        # Block 1: Load gene-to-transcript map
        FNR==NR {
            map[$2] = $1; # map[transcript_id] = gene_symbol
            next;
        }

        # Block 2: Load MANE transcript IDs (gawk-specific ARGIND)
        ARGIND==2 {
            mane[$3] = 1; # Store MANE IDs in an array
            next;
        }

        # Block 3: Process the main BED file
        {
            current_tx = $4;
            # Step A: Rename the transcript using the map
            if (current_tx in map) {
                gene_symbol = map[current_tx];
                renamed_field = gene_symbol "(" current_tx ")";
            } else {
                gene_symbol = "UNKNOWN";
                renamed_field = "UNKNOWN(" current_tx ")";
            }

            # Step B: Determine colour and note based on MANE status
            colour = "062,163,104"; # Default to green
            note = "";
            if (current_tx in mane) {
                colour = "128,0,128"; # Purple for MANE
                note = ";note=MANE_Select";
            }

            # Step C: Construct the final GFF-style attributes for column 4
            $4 = "Name=" renamed_field ";id=" current_tx ";alias=" gene_symbol note;

            # Step D: Set the colour in column 9
            $9 = colour;

            print $0;
        }
    ' "$GENE_TX_MAP" "$MANE_ENSEMBL_TX" "$SOURCE_BED"

} > "$FINAL_ENSEMBL_BED"


################################################################################
### STEP 3: Create Avenio Subset
################################################################################

echo "--> Creating Avenio gene subset..."

# This awk command is more robust, using match() to find the alias tag
# instead of splitting the string multiple times.
{
    echo "$AVENIO_HEADER"
    gawk '
        # Block 1: Load the gene list
        NR==FNR {
            genes[$1] = 1;
            next;
        }

        # Block 2: Filter the final tagged BED file
        {
            # Efficiently find the gene symbol in the alias tag
            if (match($4, /alias=([^;]+)/, arr)) {
                if (arr[1] in genes) {
                    print $0;
                }
            }
        }
    ' "$AVENIO_GENE_LIST" "$FINAL_ENSEMBL_BED"

} > "$FINAL_AVENIO_BED"

echo "---"
echo "✅ All tasks complete."
echo "Final Ensembl Track: $FINAL_ENSEMBL_BED"
echo "Final Avenio Subset: $FINAL_AVENIO_BED"