#!bin/bash
# Specify in & out directories
cd ~/Desktop/Research/PalEON_CR/env_regional/

# Script transfers & packages environmental drivers into new version
dir_in="env_paleon"
dir_out="phase2_env_drivers_v1"


# drivers=(domain_mask co2 lulcc soil biome nitrogen)
drivers=(domain_mask lulcc soil biome nitrogen)

mkdir -p ${dir_out}

cp README.txt ${dir_out}

for VAR in ${drivers[@]}
do
    echo ${VAR}
    tar -jcvf ${dir_out}/${VAR}.tar.bz2 ${dir_in}/${VAR}/
done

