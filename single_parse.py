import sys
import re
testResults = set(())
f = open(sys.argv[1], "r")
for line in f:
    if "timecop_fail" in line or "timecop_pass" in line or "timecop_error" in line:
        data = line.split(" ", 9)
        # Adjust the searchers below to match your flagset
        flag = re.search('-g_-(.+?)_', data[8]).group(1)
        if flag == "lgmp":
            flag = "None"
        test = data[7]+'\t'+data[6]+'\t'+data[8].split('/clang/')[1].split('/bin/clang')[0]+'\t'+flag
        if not test in testResults:
            testResults.add(test)
res = list(testResults)
res.sort()
for t in res:
    print(t)