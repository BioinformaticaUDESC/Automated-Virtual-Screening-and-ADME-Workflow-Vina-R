#!/bin/bash
set -euo pipefail

# ============================================================================
# USAGE: Run this script from anywhere. It will automatically locate the
#        receptors directory at ~/Vina_Workspace/receptors/
# 
# Expected directory structure (flat structure):
#   ~/Vina_Workspace/
#     ├── receptors/
#     │   ├── WNV_E/              (or WNV/NS5/, YFV/E/, etc.)
#     │   │   ├── receptor_WNV_E.pdb    (each subdirectory should have its own receptor file)
#     │   │   └── *.pdbqt (ligand files, copied from compounds directory)
#     │   ├── YFV_NS5/
#     │   │   ├── receptor_YFV_NS5.pdb    (different receptor structure for this protein)
#     │   │   └── ...
#     │   └── ...
#     └── compounds/
#         └── (compound folders with *_antechamber.mol2 files, converted to *.pdbqt by script)
#
# The script processes each direct subdirectory in the receptors folder.
# Each protein subdirectory should contain its own receptor file named "receptor_{PROTEIN_NAME}.pdb"
# where {PROTEIN_NAME} matches the subdirectory name.
# ============================================================================

# Load aliases and functions
source ~/.bashrc

# Base directory: automatically locate ~/Vina_Workspace/receptors/
# This allows the script to be run from any location
VINA_WORKSPACE="$HOME/Vina_Workspace"
RECEPTORS_DIR="$VINA_WORKSPACE/receptors"

if [ ! -d "$RECEPTORS_DIR" ]; then
    echo "Error: Receptors directory not found at $RECEPTORS_DIR"
    echo "Please ensure ~/Vina_Workspace/receptors/ exists and contains your protein subdirectories."
    exit 1
fi

BASE_DIR="$RECEPTORS_DIR"
cd "$BASE_DIR"
BASE_DIR=$(realpath "$BASE_DIR")

# Path to MGLTools
MGLTOOLS_UTILS="/home/your_path/Downloads/mgltools_x86_64Linux2_1.5.7/MGLToolsPckgs/AutoDockTools/Utilities24"

# Compounds directory (default: ~/Vina_Workspace/compounds)
COMPOUNDS_SRC="$HOME/Vina_Workspace/compounds"

if [ ! -d "$COMPOUNDS_SRC" ]; then
    echo "Error: Compounds directory not found at $COMPOUNDS_SRC"
    echo "Please ensure ~/Vina_Workspace/compounds/ exists and contains your *_antechamber.mol2 ligand files."
    exit 1
fi

# Check for MOL2 files in compounds directory
if ! ls "$COMPOUNDS_SRC"/*_antechamber.mol2 1>/dev/null 2>&1; then
    echo "Error: No *_antechamber.mol2 ligand files found in $COMPOUNDS_SRC"
    echo "Please ensure ~/Vina_Workspace/compounds/ contains *_antechamber.mol2 ligand files."
    exit 1
fi

# General Vina parameters
DOCKING_DIR="dockings"
FINAL_REPORT="Final_report.txt"
THRESHOLD=0.1
SIZE_X=20
SIZE_Y=20
SIZE_Z=20
NUM_MODES=9
EXHAUSTIVENESS=8
ENERGY_RANGE=4

# Initialize main report
{
    echo "Docking Report - $(date +'%Y-%m-%d %H:%M:%S')"
    echo "----------------------------------------"
} > "$FINAL_REPORT"

# Convert all ligands to PDBQT format in compounds directory (once, before processing proteins)
echo "===== Converting ligands to PDBQT format ====="
cd "$COMPOUNDS_SRC"
converted_count=0
for LIGAND_MOL2 in *_antechamber.mol2; do
    [ -f "$LIGAND_MOL2" ] || continue
    
    LIGAND_NAME="${LIGAND_MOL2%_antechamber.mol2}"
    LIGAND_PDBQT="${LIGAND_NAME}.pdbqt"
    
    # Skip if PDBQT already exists and is newer than MOL2
    if [ -f "$LIGAND_PDBQT" ] && [ "$LIGAND_PDBQT" -nt "$LIGAND_MOL2" ]; then
        echo "  Skipping $LIGAND_NAME (PDBQT already exists and is up-to-date)"
        continue
    fi
    
    echo "  Converting $LIGAND_NAME to PDBQT..."
    if pythonsh "$MGLTOOLS_UTILS/prepare_ligand4.py" \
        -l "$LIGAND_MOL2" \
        -o "$LIGAND_PDBQT" \
        -A bonds_hydrogens \
        &> "Prepare_ligand_${LIGAND_NAME}_OUT.out"; then
        converted_count=$((converted_count + 1))
    else
        echo "  Warning: Failed to convert $LIGAND_NAME" >> "$FINAL_REPORT"
    fi
done
echo "Converted $converted_count ligand(s) to PDBQT format"
cd "$BASE_DIR"

# Loop through each protein subdirectory
for DIR in "$BASE_DIR"/*/; do
    [ -d "$DIR" ] || continue
    echo "===== Entering $DIR ====="
    cd "$DIR"
    
    # Create dockings directory in each protein subdirectory
    mkdir -p "$DOCKING_DIR"

    # Copy all ligand files (*.pdbqt) into each protein subfolder
    echo "Copying ligand files to $DIR"
    if ! cp "$COMPOUNDS_SRC"/*.pdbqt . 2>/dev/null; then
        echo "Error: No *.pdbqt ligand files found in $COMPOUNDS_SRC" >> "$FINAL_REPORT"
        echo "Error: No *.pdbqt ligand files found in $COMPOUNDS_SRC"
        echo "Please ensure ligands have been converted to PDBQT format in the compounds directory."
        cd "$BASE_DIR"
        continue
    fi

    # Determine receptor filename based on subdirectory name
    PROTEIN_NAME="${DIR##*/}"
    RECEPTOR_FILE="receptor_${PROTEIN_NAME}.pdb"
    
    # Check for receptor file in current subdirectory
    if [ ! -f "$RECEPTOR_FILE" ]; then
        echo "Error: $RECEPTOR_FILE not found in ${PROTEIN_NAME}" >> "$FINAL_REPORT"
        echo "Error: $RECEPTOR_FILE not found. Please ensure each protein subdirectory contains its own receptor file named 'receptor_{PROTEIN_NAME}.pdb' where {PROTEIN_NAME} matches the subdirectory name."
        cd "$BASE_DIR"
        continue
    else
        echo "Using $RECEPTOR_FILE from ${PROTEIN_NAME} subdirectory"
    fi
    
    # Receptor PDBQT preparation
    echo "Preparing receptor..."
    pythonsh "$MGLTOOLS_UTILS/prepare_receptor4.py" \
        -r "$RECEPTOR_FILE" \
        -o receptor.pdbqt \
        -A bonds_hydrogens \
        &> Prepare_receptor_OUT.out

    # Run Concavity
    echo "Running Concavity..."
    concavity "$RECEPTOR_FILE" concavity_output.txt.scores
    if [ ! -s concavity_output.txt.scores ]; then
        echo "Concavity failed or generated empty file in $DIR" >> "$FINAL_REPORT"
        cd "$BASE_DIR"
        continue
    fi

    # Identify pockets
    pockets=()
    current_group=()
    while read -r line; do
        [[ "$line" =~ ^# ]] && continue
        res_id=$(awk '{print $1}' <<< "$line")
        score=$(awk '{print $3}' <<< "$line")
        if awk "BEGIN {exit !($score > $THRESHOLD)}"; then
            if [ ${#current_group[@]} -eq 0 ]; then
                current_group=("$res_id")
            else
                last=${current_group[-1]}
                if [ "$res_id" -eq $((last + 1)) ]; then
                    current_group+=("$res_id")
                else
                    pockets+=("$(printf "%s " "${current_group[@]}")")
                    current_group=("$res_id")
                fi
            fi
        else
            if [ ${#current_group[@]} -gt 0 ]; then
                pockets+=("$(printf "%s " "${current_group[@]}")")
                current_group=()
            fi
        fi
    done < concavity_output.txt.scores
    [ ${#current_group[@]} -gt 0 ] && pockets+=("$(printf "%s " "${current_group[@]}")")

    echo "Identified ${#pockets[@]} pockets in $DIR"

    # Convert receptor to PDBQT (reinforcement)
    echo "Converting receptor PDB→PDBQT..."
    pythonsh "$MGLTOOLS_UTILS/prepare_receptor4.py" \
        -r "$RECEPTOR_FILE" \
        -o receptor.pdbqt

    # Calculate pocket centers once (before ligand loop)
    declare -a pocket_centers
    for idx in "${!pockets[@]}"; do
        read -r -a res_array <<< "${pockets[$idx]}"
        center=$(awk -v res_list="$(IFS=,; echo "${res_array[*]}")" '
            BEGIN {
                split(res_list,R,",");
                for(i in R) seen[R[i]]=1;
                sumx=sumy=sumz=count=0;
            }
            /^ATOM/ && substr($0,13,4)==" CA " {
                resn=int(substr($0,23,4));
                if(resn in seen){
                    sumx+=substr($0,31,8);
                    sumy+=substr($0,39,8);
                    sumz+=substr($0,47,8);
                    count++;
                }
            }
            END {
                if(count>0) printf "%.3f %.3f %.3f", sumx/count, sumy/count, sumz/count;
                else print "0 0 0";
            }
        ' "$RECEPTOR_FILE")
        pocket_centers[$idx]="$center"
    done

    # Loop through each ligand file individually
    for LIGAND_PDBQT in *.pdbqt; do
        [ -f "$LIGAND_PDBQT" ] || continue
        
        # Extract ligand name (remove .pdbqt suffix)
        LIGAND_NAME="${LIGAND_PDBQT%.pdbqt}"
        
        echo "Processing ligand: $LIGAND_NAME"

        # Docking for each pocket
        for idx in "${!pockets[@]}"; do
            read cx cy cz <<< "${pocket_centers[$idx]}"

            CONFIG="config_${LIGAND_NAME}_pocket_${idx}.txt"
            OUT_PDBQT="${DOCKING_DIR}/resultado_${DIR##*/}_${LIGAND_NAME}_pocket_${idx}.pdbqt"
            LOG="${DOCKING_DIR}/${DIR##*/}_${LIGAND_NAME}_pocket${idx}.log"

            {
                echo "Subfolder: ${DIR##*/} | Ligand: $LIGAND_NAME | Pocket #$idx"
                echo "  Center: ($cx,$cy,$cz)"
            } >> "$FINAL_REPORT"

            cat > "$CONFIG" <<EOF
receptor = receptor.pdbqt
ligand = ${LIGAND_PDBQT}
center_x = $cx
center_y = $cy
center_z = $cz
size_x = $SIZE_X
size_y = $SIZE_Y
size_z = $SIZE_Z
num_modes = $NUM_MODES
exhaustiveness = $EXHAUSTIVENESS
energy_range = $ENERGY_RANGE
out = $OUT_PDBQT
EOF

            echo "Running Vina for $LIGAND_NAME, pocket #$idx in ${DIR##*/}..."
            if vina --config "$CONFIG" > "$LOG" 2>&1; then
                echo "  ✔ Completed $LIGAND_NAME, pocket #$idx"
            else
                echo "  ✖ Failed $LIGAND_NAME, pocket #$idx (see $LOG)" >> "$FINAL_REPORT"
            fi
        done
    done

    cd "$BASE_DIR"
done

echo "===== All dockings completed ====="
echo "Complete report in: $FINAL_REPORT"
