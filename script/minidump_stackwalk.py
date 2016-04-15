import os
import glob

os.chdir(os.path.dirname(__file__))
list= glob.glob("*.dmp")

for file in list:
    cmd = "minidump_stackwalk.exe " + file + " symbols > stack_" + file + ".txt"
    print cmd
    os.system(cmd)
raw_input()
