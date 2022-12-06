#!/bin/bash

inDIR=/mgi_storage/sk/stomics
visualSif=/mgi_storage/sk/stomics/SAW_v5.4.0.sif

export SINGULARITY_BIND=$inDIR

singularity exec ${visualSif} checkGTF \
    -i ${inDIR}/reference/Triticum_aestivum.IWGSC54/Triticum_aestivum.IWGSC.54.gff3 \
    -o ${inDIR}/reference/Triticum_aestivum.IWGSC54/genes_new.gtf

