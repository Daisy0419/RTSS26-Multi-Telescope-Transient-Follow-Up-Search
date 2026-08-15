#!/bin/bash
#SBATCH --job-name=rtss_tests
#SBATCH --account=engr-lab-jbuhler
#SBATCH --partition=general-cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=36
#SBATCH --mem=1500G
#SBATCH --time=5-00:00:00
#SBATCH --output=rtss_tests_%j.out
#SBATCH --error=rtss_tests_%j.err

cd /home/compute/w.yanwang/Multi-Telescope-Followup-Searching-cluster/results2
export GRB_LICENSE_FILE=$HOME/.gurobi/wls/gurobi.lic

echo "Running on $(hostname)"
date

GXX_LIBSTDCPP=$(g++ -print-file-name=libstdc++.so.6)
GXX_LIBDIR=$(dirname "$GXX_LIBSTDCPP")

echo "GCC libstdc++: $GXX_LIBSTDCPP"
strings "$GXX_LIBSTDCPP" | grep GLIBCXX_3.4.29
strings "$GXX_LIBSTDCPP" | grep CXXABI_1.3.13

export LD_LIBRARY_PATH="$GXX_LIBDIR:$GUROBI_HOME/lib:$LD_LIBRARY_PATH"

ldd ../build2/ts | grep libstdc++

# python -u run_script.py
LD_PRELOAD="$GXX_LIBSTDCPP" python -u run_script.py

date