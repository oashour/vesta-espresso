#!/bin/bash

# Authored by Omar A. Ashour, 20230717
# This script lives at https://github.com/oashour/vesta-espresso
# VESTA itself is licensed by http://jp-minerals.org/vesta/en/
#
# If the script isn't passed any .pwi or .in files, it just calls
# VESTA transparently with its arguments (if any)
# Otherwise, it creates a .cif file (myfile.pwi -> .myfile.pwi.cif) in 
# the same directory as the input, and opens it in VESTA.

if [ -z "$PMG_VENV" ]; then
	PMG_VENV="$PWD/venv"
fi
if [ -z "$VESTA" ]; then
	VESTA="/Applications/VESTA/VESTA.app/Contents/MacOS/VESTA"
fi

# Check if any of the input files are .pwi or .in (i.e., QE)
qe_found=false
for f in "$@"; do
    ext="${f##*.}"
    if [[ $ext == "pwi" || $ext == "in" ]]; then
        # Set up log file in the same directory as the the first input
        log="$(dirname "${f}")/.vesta-espresso.log"
        rm -f "${log}" # Delete left over logs
        # Activate venv
        source "${PMG_VENV}/bin/activate"
        qe_found=true
        break
    fi
done

# If not QE, just open VESTA
if [[ $qe_found == false ]]; then
    $VESTA "$@"
else
    # If QE, convert only QE files to .cif and open VESTA
    temp_files=()
    perm_files=()
    for f in "$@"; do
        dir=$(dirname "${f}")
        basename=$(basename "${f}")
        ext="${f##*.}"
        # If the file is .pwi or .in (i.e., QE), make temporary .cif
        if [[ $ext == "pwi" || $ext == "in" ]]; then
            cif=${dir}/.${basename}.cif
            temp_files+=("${cif}")
            echo -e "Input: ${f}\nOutput: ${cif}\n****" >> "${log}" 2>&1
            python -c "from pymatgen.io.espresso.inputs import PWin;\
                    struct = PWin.from_file('${f}').structure;\
                    struct.to('${cif}');" >> "${log}" 2>&1
            echo -e "\n----" >> "${log}" 2>&1
        else
            # Otherwise, just use the file as-is
            perm_files+=("${f}")
        fi
    done

    echo -e "Opening VESTA...\n" >> "${log}" 2>&1
    echo -e "Command: $VESTA ${temp_files[@]} ${perm_files[@]}\n" >> "${log}" 2>&1
    $VESTA "${temp_files[@]}" "${perm_files[@]}" && rm "${temp_files[@]}"
fi
