#pip install gefpy

## bgef 转换 bgem 
from gefpy.bgef_reader_cy import BgefR

input_gef = "./test.gef" ## input file => .gef
bin_size = 1  ## bin size
thread = 4    ## thread number
out_gem = "./test.gem"  ## output file => .gem

bgef = BgefR(filepath= input_gef, bin_size= bin_size, n_thread= thread)
bgef.to_gem(out_gem)


######################################
from gefpy.gef_to_gem_cy import gefToGem

strout = "SS200000795TL_B2.tissue.gem"
strsn = "SS200000795TL_B2"
obj = gefToGem(strout, strsn)

# generate bgem
strbgef = "SS200000795TL_B2.tissue.gef"
binsize = 1
obj.bgef2gem(strbgef, binsize)

# generate cgem by bgef and cgef
strcgef = "FP200000617TL_B6.cgef"
obj.cgef2gem(strcgef, strbgef) # here strbgef should be SN.raw.gef

# generate cgem by bgef and mask
strmask = "FP200000617TL_B6_mask.tif"
obj.bgef2cgem(strmask, strbgef)

##############################
## cgef based on bgem and mask
from gefpy.cgef_writer_cy import generate_cgef

mask_file = "FP200000617TL_B6_mask.tif"
bgef_file = "FP200000617TL_B6.raw.bgef"
cgef_file = "FP200000617TL_B6.cgef"
block_sizes = [256, 256]

# Generate cgef by bgef and mask
generate_cgef(cgef_file, bgef_file, mask_file, block_sizes)
