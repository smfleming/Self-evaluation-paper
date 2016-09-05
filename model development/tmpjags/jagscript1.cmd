load dic
model in "/Users/stevefleming/Dropbox/Research/Metacognition/BN_model/development code/model development/confidence_sequential.txt"
data in jagsdata.R
compile, nchains(1)
parameters in jagsinit1.R
initialize
update 10
monitor set d, thin(1)
monitor set a, thin(1)
monitor deviance
update 10000
coda *, stem('CODA1')
