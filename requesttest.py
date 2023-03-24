import requests
import argparse
import csv

parser = argparse.ArgumentParser()
parser.add_argument("-Z")
parser.add_argument("-I")
args = parser.parse_args()

Z = args.Z
I = args.I
x = requests.get(f'https://physics.nist.gov/cgi-bin/ASD/energy1.pl?de=0&spectrum=Z%3D{Z}+{I}&submit=Retrieve+Data&units=0&format=2&output=0&page_size=15&multiplet_ordered=1&conf_out=on&term_out=on&level_out=on&j_out=on&temp=')
print(x.text)
# content = x.text.replace("\"","")
# content = content.replace("=","")
# lines = content.splitlines()
# #print(content.splitlines())

# file = open(f"NISTtest\\NIST_data_{Z}0{I}.txt", "w")
# for l in lines[1:]:
#     if "Limit" in l:
#         break
#     file.write(l+"\n")
# file.close()
