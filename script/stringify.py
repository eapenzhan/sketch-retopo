import sys

fname_in = sys.argv[1]
fname_out = fname_in + ".stringify.txt"

file_in  = open(fname_in, "r")
file_out = open(fname_out, "w")

for row in file_in:
    file_out.write("\"" + row.rstrip() + "\\n\",\n")
