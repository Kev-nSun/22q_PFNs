#!/bin/bash

# Function to display help message
show_help() {
cat <<'EOF'
Usage: unzip_files.sh [--add_ses_freesurfer] <input_dir> <output_dir> [sub_list|file_pattern] [file_pattern]

Extract all .zip files from the specified input directory to the output directory using 7z.
Options:
  --add_ses_freesurfer  After extracting, rename sourcedata/freesurfer/sub-<id> to
                        sourcedata/freesurfer/sub-<id>_ses-<ses> using the ses from the zip filename.
  sub_list              Path to a text file with subject IDs (one per line).
  file_pattern          A pattern to extract specific files from each archive (using 7z's -r flag).

Important notes:
  * When --add_ses_freesurfer is enabled, your subject_list must include both subject and session IDs.
    For example:
      sub-1_ses-1
      sub-1_ses-2
      sub-1_ses-3

Examples:
  unzip_files.sh /path/to/input /path/to/output
  unzip_files.sh --add_ses_freesurfer /path/to/input /path/to/output
  unzip_files.sh /path/to/input /path/to/output subject_list.txt
  unzip_files.sh /path/to/input /path/to/output "*/sub-*/anat/*_space-MNI152NLin6Asym_res-2_desc-preproc_T1w.nii.gz"
  unzip_files.sh /path/to/input /path/to/output subject_list.txt "*/sub-*/anat/*_space-MNI152NLin6Asym_res-2_desc-preproc_T1w.nii.gz"
EOF
exit 0
}

# Check if help is requested
if [ "$1" == "-h" ] || [ "$1" == "--help" ]; then
    show_help
fi

# Optional flag parsing
add_ses_freesurfer=false
while true; do
    case "$1" in
        --add_ses_freesurfer)
            add_ses_freesurfer=true
            shift
            ;;
        *)
            break
            ;;
    esac
done

# Check if the correct number of arguments is provided (minimum 2, maximum 4 after flags)
if [ "$#" -lt 2 ] || [ "$#" -gt 4 ]; then
    echo "Error: Invalid number of arguments."
    show_help
fi

# Assign input and output directory variables
input_dir="$1"
output_dir="$2"
sub_list=""
file_pattern=""

# Determine if the third argument is a subject list file or a file pattern.
if [ "$#" -ge 3 ]; then
    if [ -f "$3" ]; then
        sub_list="$3"
    else
        file_pattern="$3"
    fi
fi

# If four arguments are provided, assume the third is sub_list and the fourth is file_pattern.
if [ "$#" -eq 4 ]; then
    sub_list="$3"
    file_pattern="$4"
fi

# Validate input directory
if [ ! -d "$input_dir" ]; then
    echo "Error: Input directory '$input_dir' does not exist."
    exit 1
fi

# Ensure output directory exists
mkdir -p "$output_dir"

# Post-extraction fix for freesurfer directory naming
post_extract_fix() {
    local zip_file="$1"
    local base
    base="$(basename "$zip_file")"

    # Extract sub-<id> and ses-<id> from the zip filename
    local subj
    local ses
    subj="$(echo "$base" | sed -n 's/.*\(sub-[0-9A-Za-z]\+\).*/\1/p')"
    ses="$(echo "$base" | sed -n 's/.*\(ses-[0-9A-Za-z]\+\).*/\1/p')"

    if [ -z "$subj" ] || [ -z "$ses" ]; then
        echo "Warning: Could not parse subject or session from $base. Skipping freesurfer rename."
        return
    fi

    # Locate freesurfer subject dir anywhere under output_dir (handles extra top-level folder)
    local src_dir
    src_dir="$(find "$output_dir" -type d -path "*/sourcedata/freesurfer/$subj" -print -quit)"

    if [ -z "$src_dir" ]; then
        echo "Note: No freesurfer dir found for $subj under $output_dir"
        return
    fi

    # Build destination alongside the found source dir
    local dst_dir="${src_dir%/$subj}/${subj}_${ses}"

    if [ -d "$dst_dir" ]; then
        echo "Note: Destination $dst_dir already exists. Merging contents from $src_dir."
        (
            shopt -s dotglob
            mv "$src_dir"/* "$dst_dir"/ 2>/dev/null || true
        )
        rmdir "$src_dir" 2>/dev/null || true
    else
        mv "$src_dir" "$dst_dir"
    fi
    echo "Renamed freesurfer directory to ${subj}_${ses}."
}

# Define extraction command based on whether a file pattern is provided.
if [ -n "$file_pattern" ]; then
    extract_cmd() {
        local zip_file="$1"
        echo "Extracting (with file pattern): $zip_file"
        if [ "$add_ses_freesurfer" = true ]; then
            # Preserve paths so freesurfer tree exists for rename
            7z x "$zip_file" -aoa -o"$output_dir" -r "$file_pattern"
            post_extract_fix "$zip_file"
        else
            # Default: flatten when not renaming freesurfer
            7z x "$zip_file" -aoa -o"$output_dir" -r "$file_pattern"
        fi
    }
else
    extract_cmd() {
        local zip_file="$1"
        echo "Extracting: $zip_file"
        7z x "$zip_file" -aoa -o"$output_dir"
        if [ "$add_ses_freesurfer" = true ]; then
            post_extract_fix "$zip_file"
        fi
    }
fi

# Extraction using subject list if provided
if [ -n "$sub_list" ]; then
    if [ ! -f "$sub_list" ]; then
        echo "Error: Subject list file '$sub_list' not found."
        exit 1
    fi
    while IFS= read -r subid || [ -n "$subid" ]; do
        [ -z "$subid" ] && continue
        subid=${subid%$'\r'}   # strip carriage return if present
        found=false
        for zip_file in "$input_dir"/*${subid}*.zip; do
            if [ -f "$zip_file" ]; then
                extract_cmd "$zip_file"
                found=true
            fi
        done
        if [ "$found" == false ]; then
            echo "Warning: No zip file found for subject: $subid"
        fi
    done < "$sub_list"
else
    # Extraction without a subject list (i.e. all zip files)
    found=false
    for zip_file in "$input_dir"/*.zip; do
        if [ -f "$zip_file" ]; then
            extract_cmd "$zip_file"
            found=true
        fi
    done
    if [ "$found" == false ]; then
        echo "Warning: No zip files found in $input_dir"
    fi
fi

echo "Extraction complete."
