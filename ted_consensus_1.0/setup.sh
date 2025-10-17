#!/bin/bash

# This file is a part of TED: The Encyclopedia of Domains. If you utilize or reference any content from this file,
# please cite the following paper:
# Lau et al., 2024. Exploring structural diversity across the protein universe with The Encyclopedia of Domains.

set -eo pipefail

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

# install dir can be overridden by setting the INSTALL_DIR environment variable
TED_INSTALL_DIR=${TED_INSTALL_DIR:-$SCRIPT_DIR}

# Define the name of the virtual environment directory
VENV_DIR="${TED_INSTALL_DIR}/ted_consensus"
WEIGHTS_DIR="${TED_INSTALL_DIR}/programs/merizo/weights"
UNIDOC_DIR="${TED_INSTALL_DIR}/programs/unidoc"

# Define base URL and weights files in an array
BASE_URL="https://github.com/psipred/Merizo/raw/main/weights"
WEIGHTS_FILES=("weights_part_0.pt" "weights_part_1.pt" "weights_part_2.pt")

# UniDoc package download URL
UNIDOC_URL="https://yanglab.qd.sdu.edu.cn/UniDoc/download/UniDoc.tgz"
UNIDOC_TGZ="${TED_INSTALL_DIR}/programs/unidoc.tgz"

# Function to check Python version
check_python_version() {
    PYTHON_VERSION=$($1 -c "import sys; print('.'.join(map(str, sys.version_info[:3])))")
    if [[ $PYTHON_VERSION == 3.1[01]* ]]; then
        echo $1
        return 0
    else
        return 1
    fi
}

# Check for Python 3.11 or 3.10
PYTHON_COMMAND=""
if check_python_version python3.11; then
    PYTHON_COMMAND="python3.11"
elif check_python_version python3.10; then
    PYTHON_COMMAND="python3.10"
else
    echo "Python 3.10 or 3.11 is required but not found. Please install Python 3.10 or 3.11 and try again."
    exit 1
fi

# Check if pip is installed
if ! command -v pip3 &> /dev/null
then
    echo "pip is not installed. Please install pip and try again."
    exit 1
fi

# Create a virtual environment
if [ ! -d "$VENV_DIR" ]; then
    echo "Creating virtual environment..."
    "$PYTHON_COMMAND" -m venv "$VENV_DIR"
else
    echo "Virtual environment already exists."
fi

# Activate the virtual environment
source "$VENV_DIR/bin/activate"

# Upgrade pip to the latest version
echo "Upgrading pip..."
pip install --upgrade pip

# Install dependencies from requirements.txt
if [ -f "${SCRIPT_DIR}/requirements.txt" ]; then
    echo "Installing dependencies from requirements.txt..."
    pip install -r "${SCRIPT_DIR}/requirements.txt"
else
    echo "Requirements file '${SCRIPT_DIR}/requirements.txt' not found. Please make sure it exists in the install directory."
    exit 1
fi

if [ ! -d "${TED_INSTALL_DIR}" ]; then
    echo "Creating install directory..."
    mkdir -p "${TED_INSTALL_DIR}"
else
    echo "Install directory already exists."
fi

# Copy over dirs if they do not exist
for dirname in programs scripts; do
    if [ ! -d "${TED_INSTALL_DIR}/${dirname}" ]; then
        echo "Copying ${dirname} directory..."
        cp -r "${SCRIPT_DIR}/${dirname}" "${TED_INSTALL_DIR}/${dirname}"
    else
        echo "${dirname} directory already exists."
    fi
done

# Check if the weights directory exists
if [ -d "$WEIGHTS_DIR" ]; then
    echo "programs/merizo/weights exists. Checking for missing weights files..."
else
    echo "programs/merizo/weights directory not found. Creating directory..."
    mkdir -p "$WEIGHTS_DIR"
fi

# Download missing weights files
for WEIGHT_FILE in "${WEIGHTS_FILES[@]}"; do
    FILE_PATH="${WEIGHTS_DIR}/${WEIGHT_FILE}"
    if [ ! -f "$FILE_PATH" ]; then
        wget -O "$FILE_PATH" "${BASE_URL}/${WEIGHT_FILE}"
    fi
done

# Check if the unidoc directory exists
if [ -d "$UNIDOC_DIR" ]; then
    echo "programs/unidoc directory already exists."
else
    if test ! -f "${UNIDOC_TGZ}"; then
        echo "programs/unidoc directory not found. Downloading UniDoc ..."
        wget --no-check-certificate -O "${UNIDOC_TGZ}" "${UNIDOC_URL}"
    fi

    echo "Unpacking UniDoc package..."
    tar -xzvf "${UNIDOC_TGZ}" -C "${TED_INSTALL_DIR}/programs"
    mv "${TED_INSTALL_DIR}/programs/UniDoc" "${TED_INSTALL_DIR}/programs/unidoc"

    # macOS C/C++ compatibility: replace non-portable <malloc.h> with <stdlib.h>
    if [[ "$OSTYPE" == "darwin"* ]]; then
        echo "Patching UniDoc C/C++ sources for macOS (malloc.h -> stdlib.h)..."
        # Replace in all source files under UniDoc
        while IFS= read -r -d '' hdr; do
            sed -E -i '' 's#<malloc.h>#<stdlib.h>#g' "$hdr"
        done < <(find "${UNIDOC_DIR}" -type f \( -name "*.h" -o -name "*.hpp" -o -name "*.c" -o -name "*.cpp" \) -print0)
    fi
fi

# Copy the extra run script over to the unidoc dir
cp "${SCRIPT_DIR}/scripts/Run_UniDoc_from_scratch_structure_afdb.py" "${UNIDOC_DIR}/"

# if running on macOS install compiler tools and compile stride from source:
if [[ "$OSTYPE" == "darwin"* ]]; then
    # macOS only
    if ! command -v gcc &> /dev/null || ! command -v make &> /dev/null; then
        echo "Installing Xcode Command Line Tools..."
        xcode-select --install
        # Wait for installation to complete or prompt user to press enter after installation
        read -p "Press enter after Xcode Command Line Tools installation is complete"
    fi

    cd "${TED_INSTALL_DIR}/programs/chainsaw/stride" || exit
    # Remove all files except stride.tgz
    find . -type f ! -name 'stride.tgz' -delete
    # Extract the contents of stride.tgz
    tar -zxf stride.tgz
    # Compile stride
    echo "Compiling stride for MacOS..."
    make
    chmod +x stride
    # Change back to the original directory
    cd - || exit
fi

echo "Successfully set up ted_consensus"
