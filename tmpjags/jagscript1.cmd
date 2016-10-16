model in "/Users/stevefleming/Dropbox/Research/Metacognition/BN_model/Self-evaluation-paper/sampleMetaConf.txt"
data in jagsdata.R
compile, nchains(1)
parameters in jagsinit1.R
initialize
update 1000
monitor set d, thin(1)
monitor set a, thin(1)
monitor set x1, thin(1)
update 100000
coda *, stem('CODA1')
