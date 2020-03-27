import pdb
from os import listdir
import subprocess


param_files = [f for f in listdir("params/")]
pdb.set_trace()

for param in param_files:
    network = param[:param.find('.param')]
    args = ("./netgen  < params/" +  param + " > " + "networks/netgen/" + network)
    popen = subprocess.Popen(args, stdout=subprocess.PIPE, shell=True)
