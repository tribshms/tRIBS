#!/bin/zsh

# intended to download benchmarks for big_spring and happy_jack and update the config.json file if needed.
curl -sS -o hj "https://zenodo.org/records/10909507/files/happy_jack.gz?download=1"
curl -sS -o bs "https://zenodo.org/records/10951574/files/big_spring.gz?download=1"

path_to_dir="testing/black_box/benchmark_cases"

mkdir -p $path_to_dir

tar -xzf hj -C $path_to_dir
tar -xzf bs -C $path_to_dir

rm hj
rm bs
