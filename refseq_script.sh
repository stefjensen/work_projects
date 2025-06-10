#!/bin/bash

# ==============================================================================
# A script to download and process NCBI RefSeq data to create
# custom IGV tracks.
#
# It performs the following main steps:
#   1. Downloads RefSeq GTF and MANE summary data.
#   2. Creates a map of gene symbols to transcript IDs.
#   3. Processes a source BED file to:
#      - Rename features using the gene symbol.
#      - Colour transcripts by type (MANE, model, curated).
#      - Add rich, GFF-style attributes with quoted notes.
#   4. Creates a final subset track based on a gene list.
# ==============================================================================

# --- Script Configuration ---
# Stop on any error (`-e`) and on failures in pipelines (`-o pipefail`).
set -e
set -o pipefail

# --- File and URL Variables ---
REFSEQ_GTF_URL="https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh38_latest/refseq_identifiers/GRCh38_latest_genomic.gtf.gz"
MANE_URL="https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/release_1.2/MANE.GRCh38.v1.2.summary.txt.gz" # Assumed available from previous script
AVENIO_GENE_LIST="avenio/avenio_genes.txt"
SOURCE_BED="beds/ncbiRefSeq_ucsc_wg.bed"

# --- Output Files ---
mkdir -p data beds/final avenio/final

REFSEQ_GTF_GZ="data/GRCh38_latest_genomic.gtf.gz"
MANE_SUMMARY_GZ="data/MANE.GRCh38.v1.2.summary.txt.gz"
GENE_TX_MAP="data/refseq_gene_tx_map.tsv"
MANE_REFSEQ_TX="data/mane_refseq_tx.txt"

FINAL_REFSEQ_BED="beds/final/refseq_final_tagged.bed"
FINAL_AVENIO_BED="avenio/final/refseq_avenio_subset.bed"

# --- Headers for BED Tracks ---
REFSEQ_HEADER='track name="NCBI RefSeq (Coloured)" description="MANE purple, Models orange, Curated blue" itemRgb="On" db=hg38'
AVENIO_HEADER='track name="Avenio Gene Subset (RefSeq)" description="Avenio panel genes from RefSeq" itemRgb="On" db=hg38'


################################################################################
### STEP 1: Download Data and Prepare Mapping Files
################################################################################

echo "--> Downloading RefSeq annotations..."
wget -O "$REFSEQ_GTF_GZ" "$REFSEQ_GTF_URL"

echo "--> Generating Gene-to-Transcript map from RefSeq GTF..."
# Use a single, efficient awk command on the stream to avoid a huge intermediate file.
gunzip -c "$REFSEQ_GTF_GZ" | gawk '
    $3 == "transcript" {
        # Robustly find key-value pairs using match()
        if (match($9, /gene_name "([^"]+)"/, gene) && match($9, /transcript_id "([^"]+)"/, tx)) {
            print gene[1] "\t" tx[1]
        }
    }
' > "$GENE_TX_MAP"

echo "--> Downloading MANE summary (if not present)..."
# The -nc flag prevents re-downloading if the file already exists.
wget -nc -O "$MANE_SUMMARY_GZ" "$MANE_URL"

echo "--> Extracting MANE RefSeq transcript IDs..."
# Extract the 2nd column (RefSeq ID) from the MANE summary, skipping the header.
gunzip -c "$MANE_SUMMARY_GZ" | awk 'NR > 1 {print $2}' > "$MANE_REFSEQ_TX"

# Note: The original script downloaded hgdownload.soe.ucsc.edu/.../refGene.txt.gz
# but never used it. That download has been removed to avoid unnecessary work.


################################################################################
### STEP 2: Process RefSeq Track (Rename, Colour, and Tag in ONE pass)
################################################################################

echo "--> Processing RefSeq BED file (Rename, Colour, Tag)..."

# This single gawk command replaces multiple intermediate files and processing steps.
# It reads mapping files, then processes the main BED file to produce the final,
# headered output in one go.
{
    echo "$REFSEQ_HEADER"
    gawk -v OFS='\t' '
        # Block 1: Load gene-to-transcript map
        FNR==NR {
            map[$2] = $1; # map[transcript_id] = gene_symbol
            next;
        }

        # Block 2: Load MANE transcript IDs (gawk-specific ARGIND)
        ARGIND==2 {
            mane[$1] = 1; # Store MANE IDs in an array
            next;
        }

        # Block 3: Process the main BED file
        {
            current_tx = $4;
            gene_symbol = (current_tx in map) ? map[current_tx] : "UNKNOWN";

            # --- Step A: Determine colour and note based on transcript type ---
            note_text = "";
            colour = "0,0,0"; # Default: black

            if (current_tx in mane) {
                colour = "128,0,128";   # Purple for MANE
                note_text = "MANE select";
            } else if (current_tx ~ /^X[MRP]_/) {
                colour = "213,94,0";    # Orange for Models (XM, XR, XP)
                if (current_tx ~ /^XM_/) { note_text = "model_mRNA_transcript"; }
                if (current_tx ~ /^XR_/) { note_text = "model_non-coding_RNA_transcript"; }
                if (current_tx ~ /^XP_/) { note_text = "model_protein"; }
            } else if (current_tx ~ /^N[MRP]_/) {
                colour = "100,149,237"; # Cornflower blue for Curated (NM, NR, NP)
                if (current_tx ~ /^NM_/) { note_text = "curated_mRNA_transcript"; }
                if (current_tx ~ /^NR_/) { note_text = "curated_non-coding_RNA_transcript"; }
                if (current_tx ~ /^NP_/) { note_text = "curated_protein"; }
            } else if (current_tx ~ /^NG_/) {
                colour = "0,100,0";     # Green for genomic contigs
            }

            # --- Step B: Construct the final GFF-style attributes for column 4 ---
            gff_attributes = "Name=" gene_symbol "(" current_tx ");alias=" gene_symbol;
            if (note_text != "") {
                # Add the note with spaces, correctly wrapped in quotes
                gff_attributes = gff_attributes ";note=" note_text
            }
            $4 = gff_attributes;

            # --- Step C: Set the colour in column 9 ---
            $9 = colour;

            print $0;
        }
    ' "$GENE_TX_MAP" "$MANE_REFSEQ_TX" "$SOURCE_BED"

} > "$FINAL_REFSEQ_BED"


################################################################################
### STEP 3: Create Avenio Subset
################################################################################

echo "--> Creating Avenio gene subset for RefSeq..."

# This awk command robustly finds the alias tag using match()
{
    echo "$AVENIO_HEADER"
    awk '
        # Block 1: Load the gene list
        NR==FNR {
            genes[$1] = 1;
            next;
        }

        # Block 2: Filter the final tagged BED file
        # Skip the track header line
        /^track/ { next }

        # Efficiently find the gene symbol in the alias tag
        if (match($4, /alias=([^;]+)/, arr)) {
            if (arr[1] in genes) {
                print $0;
            }
        }
    ' "$AVENIO_GENE_LIST" "$FINAL_REFSEQ_BED"

} > "$FINAL_AVENIO_BED"

echo "---"
echo "✅ All tasks complete."
echo "Final RefSeq Track: $FINAL_REFSEQ_BED"
echo "Final Avenio Subset: $FINAL_AVENIO_BED"