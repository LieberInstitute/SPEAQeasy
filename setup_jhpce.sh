#!/usr/bin/env bash
 module load conda_R/4.1.x
 Rscript scripts/check_R_packages_JHPCE.R
  
 #echo "Setting up test files..."
 #Rscript scripts/make_test_manifests.R -d $(pwd)
    
 # echo "Preparing main script..."
 #   sed -i "s|ORIG_DIR=.*|ORIG_DIR=$(pwd)|" run_pipeline_jhpce.sh
    
 echo "Done."

