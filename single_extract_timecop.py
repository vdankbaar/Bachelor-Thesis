import sys
import re
import os
f = open(sys.argv[1], "r")
for line in f:
    if "timecop_fail" in line:
        data = line.split(" ", 9)
        # Adjust the searchers below to match your flagset
        flag = re.search('-g_-(.+?)_', data[8]).group(1)
        if flag == "lgmp":
            flag = "None"
        os.makedirs('timecop_fails/'+data[7], exist_ok=True)
        f2 = open('timecop_fails/'+data[7]+'/timecop_fail_'+flag+'.txt', "a")
        f2.write(data[9])