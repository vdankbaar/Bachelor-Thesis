import sys
import re
files = ["15.0.0_O3-Os-None.txt", "15.0.6_O3-Os.txt", "16.0.0_O3-Os-None.txt", "16.0.4_O3-Os.txt", "17.0.2_O3-Os.txt", "17.0.6_O3-Os-None.txt", "15.0.6-16.0.4-17.0.2_None.txt"]
testResults = set(())
for file in files:
    f = open(file, "r")
    for line in f:
        if "timecop_fail" in line or "timecop_pass" in line or "timecop_error" in line:
            data = line.split(" ", 9)
            flag = "None"
            if "-O3" in data[8]:
                flag = "O3"
            elif "-Os" in data[8]:
                flag = "Os"
            # Adjust the searchers below to match your flagset
            test = data[7]+'\t'+data[6]+'\t'+data[8].split('/clang/')[1].split('/bin/clang')[0]+'\t'+flag
            if not test in testResults:
                testResults.add(test)
res = list(testResults)
res.sort()
for t in res:
    print(t)