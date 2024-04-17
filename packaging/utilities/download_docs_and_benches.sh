#!/bin/bash

# Function to prompt user for download options
prompt_download_options() {
    echo "Would you like to download additional resources?"
    echo "1. Documentation"
    echo "2. Benchmark Cases"
    echo "3. Both"
    echo "4. None"
    read -p "Enter your choice (1/2/3/4): " choice
    case $choice in
        1) download_documentation=true ;;
        2) download_benchmark=true ;;
        3) download_documentation=true
           download_benchmark=true ;;
        4) ;;
        *) echo "Invalid choice. Please enter a valid option."
           prompt_download_options ;;
    esac
}

# Function to download documentation using curl
download_documentation_with_curl() {
    echo "Downloading documentation with curl..."
    # Replace <documentation_url> with the actual URL of the documentation
    curl -sS -o documentation.zip "https://tribshms.readthedocs.io/_/downloads/en/latest/htmlzip/"
    unzip -q documentation.zip -d documentation
    rm documentation.zip
}

# Function to download benchmark cases using curl
download_benchmark_cases_with_curl() {
    echo "Downloading benchmark cases with curl..."
    curl -sS -o hj "https://zenodo.org/records/10909507/files/happy_jack.gz?download=1"
    curl -sS -o bs "https://zenodo.org/records/10951574/files/big_spring.gz?download=1"

    mkdir benchmark_cases

    tar -xzf hj -C benchmark_cases
    tar -xzf bs -C benchmark_cases

    rm hj
    rm bs
}

# Function to download documentation using wget
download_documentation_with_wget() {
    echo "Downloading documentation with wget..."
    wget -q -O documentation.zip "https://tribshms.readthedocs.io/_/downloads/en/latest/htmlzip/"
    unzip -q documentation.zip -d documentation
    rm documentation.zip
}

# Function to download benchmark cases using wget
download_benchmark_cases_with_wget() {
    echo "Downloading benchmark cases with wget..."
    wget -q -O hj "https://zenodo.org/records/10909507/files/happy_jack.gz?download=1"
    wget -q -Q bs "https://zenodo.org/records/10951574/files/big_spring.gz?download=1"

    mkdir benchmark_cases

    tar -xzf hj -C benchmark_cases
    tar -xzf bs -C benchmark_cases

    rm hj
    rm bs
}

# Main script starts here

# Check if curl is available
if command -v curl &> /dev/null; then
    download_function_documentation=download_documentation_with_curl
    download_function_benchmark=download_benchmark_cases_with_curl
else
    echo "Warning: curl is not available. Falling back to wget."
    download_function_documentation=download_documentation_with_wget
    download_function_benchmark=download_benchmark_cases_with_wget
fi

# Prompt user for download options
prompt_download_options

# Additional logic based on user's choice
if [ "$download_documentation" = true ]; then
    $download_function_documentation
fi

if [ "$download_benchmark" = true ]; then
    $download_function_benchmark
fi

echo "Packaging complete."
