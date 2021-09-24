#!/bin/bash

echo "**** Job starts ****"
date

## Clean up older files
# rm -fr trash

# echo "Are you sure you want to delete the contents of trash/ and trash the contents of processed-data/?"
# read response
#
# if [ $response = "y" ]; then
#         echo "OK"
#         set $checkpoint = "t"
# elif [ $response = "n" ]; then
#         echo "Please set the PATH variable in your .bashrc file before running this script."
#         set $checkpoint = "f"
# else
#         echo "Please only use 'y' and 'n'."
#         set $checkpoint = "f"
# fi

# if [ $checkpoint = "t" ]; then
## Could be fancier and use the date or something
mkdir -p trash
mv logs/* trash/
rm -r trash/processed-data
mv processed-data/ trash/

  ## Create logs dir if needed
mkdir -p logs
mkdir -p processed-data/01_get_inv_quantile_norm
mkdir -p processed-data/02_prep_inputs
# fi

echo "**** Job ends ****"
date
