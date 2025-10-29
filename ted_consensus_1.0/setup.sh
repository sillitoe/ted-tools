#!/bin/bash

# This file is a part of TED: The Encyclopedia of Domains. If you utilize or reference any content from this file,
# please cite the following paper:
# Lau et al., 2024. Exploring structural diversity across the protein universe with The Encyclopedia of Domains.

set -e -o pipefail

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

# Define the name of the virtual environment directory
VENV_DIR="ted_consensus"
WEIGHTS_DIR="${SCRIPT_DIR}/programs/merizo/weights"

UNIDOC_LATEST_DIR="${SCRIPT_DIR}/programs/unidoc_latest"
UNIDOC_TED_V1_0_DIR="${SCRIPT_DIR}/programs/unidoc_ted_v1_0"
UNIDOC_DIR="${SCRIPT_DIR}/programs/unidoc"
TED_VERSION=${TED_VERSION:-"latest"}

if [ "$TED_VERSION" == "latest" ]; then
    echo "Using latest UniDoc version"
    UNIDOC_SRC_DIR="$UNIDOC_LATEST_DIR"
elif [ "$TED_VERSION" == "1.0" ]; then
    echo "Using TED UniDoc version 1.0"
    UNIDOC_SRC_DIR="$UNIDOC_TED_V1_0_DIR"
else
    echo "Unknown TED_VERSION: $TED_VERSION"
    exit 1
fi

# Copy UniDoc source files to unidoc directory
echo "Setting up UniDoc in programs/unidoc (TED ${TED_VERSION}) ..."
mkdir -p "$UNIDOC_DIR"
rsync -a --delete "$UNIDOC_SRC_DIR/" "$UNIDOC_DIR/"

# Define base URL and weights files in an array
BASE_URL="https://github.com/psipred/Merizo/raw/main/weights"
WEIGHTS_FILES=("weights_part_0.pt" "weights_part_1.pt" "weights_part_2.pt")

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
    $PYTHON_COMMAND -m venv $VENV_DIR
else
    echo "Virtual environment already exists."
fi

# Activate the virtual environment
source $VENV_DIR/bin/activate

# Upgrade pip to the latest version
echo "Upgrading pip..."
pip install --upgrade pip

# Install dependencies from requirements.txt
if [ -f "${SCRIPT_DIR}/requirements.txt" ]; then
    echo "Installing dependencies from requirements.txt..."
    pip install -r "${SCRIPT_DIR}/requirements.txt"
else
    echo "requirements.txt not found. Please make sure it exists in the current directory."
    exit 1
fi

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


# if running on macOS install compiler tools
if [[ "$OSTYPE" == "darwin"* ]]; then
    # macOS only
    if ! command -v gcc &> /dev/null || ! command -v make &> /dev/null; then
        echo "Installing Xcode Command Line Tools..."
        xcode-select --install
        # Wait for installation to complete or prompt user to press enter after installation
        read -p "Press enter after Xcode Command Line Tools installation is complete"
    fi
fi

# check to see if the Stride binary runs, if not compile it
STRIDE_BIN="${UNIDOC_DIR}/bin/stride"
if "${STRIDE_BIN}" > /dev/null 2>&1; then
    echo "Stride binary found."
else
    rc=$?
    echo "Stride binary failed (exit code ${rc})."

    ORIG_DIR=$(pwd)
    cd "${SCRIPT_DIR}/programs/chainsaw/stride" || exit
    # Remove all files except stride.tgz
    find . -type f ! -name 'stride.tgz' -delete
    # Extract the contents of stride.tgz
    tar -zxf stride.tgz
    # Compile stride
    echo "Compiling stride for MacOS..."
    make
    chmod +x stride

    echo "Copying compiled stride to unidoc bin directory..."
    cp stride "${STRIDE_BIN}"

    cd "${ORIG_DIR}"
fi

# check to see if the UniDoc_struct binary runs
# with a zero exit code, if not compile it
UNIDOC_STRUCT_BIN="${UNIDOC_DIR}/bin/UniDoc_struct"
if "${UNIDOC_STRUCT_BIN}" > /dev/null 2>&1; then
    echo "UniDoc_struct binary found."
else
    rc=$?
    echo "UniDoc_struct binary failed (exit code ${rc})."

    if [ $TED_VERSION == "1.0" ]; then
        echo "ERROR: cannot continue. UniDoc_struct binary is missing or not working and we cannot compile from source (since TED_VERSION is set to 1.0). If you want to run on MacOS, set TED_VERSION to 'latest'."
        exit 1
    fi

    # Compile UniDoc_struct
    echo "Compiling UniDoc_struct from src..."
    ORIG_DIR=$(pwd)
    cd "${UNIDOC_DIR}/src"
    rm -f "${UNIDOC_STRUCT_BIN}"
    # patch to remove all includes for malloc.h for macOS
    sed -i.bak '/#include <malloc.h>/d' *.h
    g++ -std=c++0x   -O3 -ffast-math -lm -o "${UNIDOC_STRUCT_BIN}" UniDoc_struct.cpp

    # Change back to the original directory
    cd "${ORIG_DIR}"
fi


echo "Successfully set up ted_consensus"
