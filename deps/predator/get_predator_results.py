from __future__ import print_function
from subprocess import Popen, PIPE
import os

#with open('acinus.seq') as file_handle:
get_res_1 = Popen(['predator', '-h', 'acinus.seq'], stdout=PIPE)   #cmd result: ncoils-osf -f acinus.seq
print(get_res_1.communicate()[0])  #happily ignore err1,
