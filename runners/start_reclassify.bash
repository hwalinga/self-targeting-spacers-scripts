#!/bin/bash

for bash_folder in ../reclassify/hielke/*p; do
    bash ./run_reclassify.bash "${bash_folder}"/raw_files/
done;
