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
   # Cutoff parameters
   cat >> $PARA <<EOF
   {
       "training_datasets": [
           "$DATAFCC/SK-supercell-noHS-0$CONFIGS.h5",
           "$DATABCC/SK-supercell-noHS-0$CONFIGS.h5"
       ],
       "HS_datasets": [
           "FCC_SK-unitcell-withHS-010.h5"
       ],
       "model": {
           "predict" : 1,
           "fit" : 1,
           "cutoff_params" : {
               "rcut" : $RCUT,
               "renv" : $RENV,
               "zenv" : $ZENV,
               "pcut" : 1
           },
           "fit_params" : {
               "degree" : $DEG,
               "order" : $ORDER,
               "env_deg" : $ENV_DEG
           }
       },
       "onsite-terms": {
           "Al" : [
               -5.5708730114137865e+01,
               -4.0432979354580754e+00,
               -4.9611320942331338e-01,
               -2.6243905672312331e+00,
               -2.6243905672300967e+00,
               -2.6241786637179012e+00,
               -4.6415677514380915e-01,
               -4.6415677508608638e-01,
               -4.6415677490900636e-01,
               -2.8478483849360947e-01,
               -2.8478483839241397e-01,
               -2.8478483828616641e-01,
               -1.8885840374342641e-01,
               -1.8885840355491695e-01
           ]
       }
   }
EOF

# Run the calculation
$RUN_PREFIX $DFTB_PATH/dftb+

# Move outputs to right folder and right names
cp $PARA "$FOLDER"pot_order"$ORDER"_deg"$DEG"_envdeg"$ENV_DEG"_rcut"$RCUT"_renv"$RENV"_zenv"$ZENV"_conf"$CONFIGS".json

done

