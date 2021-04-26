#!/bin/sh

RUN_PREFIX="mpirun -np 1"
DFTB_PATH="../build/install/bin/"
PARA="Al-Al.acetb_fits.json"
FOLDER="Results_fits/"
DATABCC="data/BCC_MD_20K/MD-20K/"
DATAFCC="data/FCC_MD_500K/MD-500K/"

# Cutoff parameters
RCUT=12.
RENV=3.
ZENV=1

CONFIGS=10 #configurations used (BCC and FCC) - does not do anything here, this needs to be changed manually in config.json

ENV_DEG=4.
ORDER=2

for DEG in {6..7};
do
   echo "$ORDER, $DEG"
   echo "$FOLDER"
   # Including the right parameters
   >$PARA

# Get fitted model from folder with right names
cat "$FOLDER"pot_order"$ORDER"_deg"$DEG"_envdeg"$ENV_DEG"_rcut"$RCUT"_renv"$RENV"_zenv"$ZENV"_conf"$CONFIGS".json | awk '/"fit":/{gsub(/1/, "0")};{print}' > $PARA 

# Run the calculation
$RUN_PREFIX $DFTB_PATH/dftb+

$DFTB_PATH/dp_bands band.out band

cat band_tot.dat > "$FOLDER"BANDS_order"$ORDER"_deg"$DEG"_envdeg"$ENV_DEG"_rcut"$RCUT"_renv"$RENV"_zenv"$ZENV"_conf"$CONFIGS".dat

done

