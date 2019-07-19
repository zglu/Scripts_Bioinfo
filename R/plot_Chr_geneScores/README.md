## Plot gene scores (eg. DEG logFC values) on chromosomes

Example input: 

Chr startCoord endCoord geneName score strand chr_length

~~~~~~
SM_V7_1 17077443 17081024 Smp_186930 4.525 - 88881357
SM_V7_1 19844129 19845153 Smp_341340 -6.43 + 88881357
SM_V7_2 4691398 4746854 Smp_146940 5.46 - 48130368
SM_V7_4 27170077 27175632 Smp_145490 9.964 - 47279781
SM_V7_ZW 86588373 86589433 Smp_324380 -6.918 - 88385488
~~~~~~

Usage: 

    Rscript plotChr_genesScores.R [INPUTFILE]
