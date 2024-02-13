import sys
import re
import os
files = ["15.0.0_O3-Os-None.txt", "15.0.6_O3-Os.txt", "16.0.0_O3-Os-None.txt", "16.0.4_O3-Os.txt", "17.0.2_O3-Os.txt", "17.0.6_O3-Os-None.txt", "15.0.6-16.0.4-17.0.2_None.txt"]
for file in files:
    print("Opening: "+file)
    f = open(file, "r")
    for line in f:
        if "timecop_fail" in line:
            data = line.split(" ", 9)
            version = data[8].split('/clang/')[1].split('/bin/clang')[0]
            flag = "None"
            if "-O3" in data[8]:
                flag = "O3"
            elif "-Os" in data[8]:
                flag = "Os"
            os.makedirs('timecop_fails/'+data[7], exist_ok=True)
            f2 = open('timecop_fails/'+data[7]+'/timecop_fail_'+version+'_'+flag+'.txt', "a")
            f2.write(data[9])