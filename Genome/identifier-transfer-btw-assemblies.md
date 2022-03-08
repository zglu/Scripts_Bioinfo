*** gene model transfer by RATT ***
#Convert from or to EMBL format: http://ratt.sourceforge.net/transform.html (I used EMBL2GFF and it worked fine)

# Set config to eukaryotic
RATT_HOME=/nfs/pathogen003/tdo/Tools/PAGIT/64bit/PAGIT/RATT/
RATT_CONFIG=/nfs/users/nfs_t/tdo/Bin/ratt/RATT.config_euk; export RATT_CONFIG
 
# run command format: /nfs/users/nfs_t/tdo/Bin/ratt/start.ratt.sh <Directory with embl-files> <Query-fasta sequence> <Resultname> <Transfer type>
# I used the latest version from Thomas's directory. RATT was run on unmasked genome. I've tested different transfer types: Assembly/Strain/Species... and the Pacbio setting has the best accuracy.
# This is the command I used. After the transfer, take *.final.embl as the result file.
# It seems Pacbio is not a pre-set parameter any more. Try others or the Free option with c = 800, l = 30, g = 2000 
bsub.py 24 ratt_UM_Pacbio /nfs/pathogen003/tdo/Tools/PAGIT/64bit/PAGIT/RATT/start.ratt.sh /nfs/repository/working_area/SHISTO/chado_dump_08032017/artemis/EMBL/Smansoni/1 ~zl3/Smansoni/V7/Smansoni_v7.fa TransferPacbio Pacbio

# get pep sequences from gff: 
~zl3/scripts/gff2fasta.pl (this script couldn't get the sequence of the last gene) OR
~zl3/scripts/GFF2fa_argv.py 

# evaluate accuracy of transferred models by pariwise sequence alignment (or normal blast)
~zl3/scripts/bash/Blast-Array.sh (need to change the name of id pair file)

*** compare coordinates of RATT and Augustus models ***
gffcompare [-r <reference_mrna.gtf> [-R]] [-T] [-V] [-s <seq_path>]
    [-o <outprefix>] [-p <cprefix>]
    {-i <input_gtf_list> | <input1.gtf> [<input2.gtf> .. <inputN.gtf>]}

~zl3/bin/software/gffcompare/gffcompare -r ../../eval/TransferPacbio.SM_V7_1.final.gff -R -T -o rattUM_aug_C1 ../../augustus_new_masked/aug_C1.gff (-T not to produce temp files, -R not used)

# In the *annotated.gtf file gets the coordinates of Augustus, but also the corresponding gene name in RATT, this will help to find out whether this is 1 to 1, or many to 1. May need to swap the two gff files in the above command to find the 1-to-many situation (?).
 
*** change the ids in gff file ***
# e.g., change g1.t1 to Smp_xxxxx.1 and g10.t1 to Smp_xxxxxx.1
~zl3/scripts/pair_replace.pl
