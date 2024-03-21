# Raw results
The folder 'raw_results' contains renamed 'data' files as they are provided upon completion of a SUPERCOP run.\
Which version/flags were used is listed in the file names.

# Godbolt
The folder 'godbolt' contains the files with C code that I used in godbolt to get a look at the assembly produced with and without optimization. I also included some screenshots for some.

# Parsing scripts
The python files can be used to sort and transform the raw data files.
To be of use in parsing new data files, compiler versions or flags the scripts will need to be rewritten as they were hardcoded to parse my specific raw data files and flag combinations.\
All scripts print their data to the console by default. To store the results in a file please pipe them to a file manually.\
Example: './py.exe parse.py > results.csv'

## parse.py
This script takes a list of txt files and extracts all the timecop fails. It then tab deliminates some of the interesting values and prints them out.\
single_parse.py does the same but allows the user to input a single file as a parameter and runs the parsing on specifically that file.\
results.csv was generated with this script.

## extract_timecop.py
This script takes a list of txt files and extracts all the timecop fails. It then generates a system of folders and text files that contain the full TIMECOP output of that failure.\
single_extract_timecop.py does the same but allows the user to input a single file as a parameter and runs the process on specifically that file.\
The contents of the results-timecop_fails directory was generated with this script.

# Results
## results.csv
The file 'results.csv' together with 'results_18.1.1.csv' contains the result of every benchmark including each compiler version and flag combination and whether it passed, failed or errored on the TIMECOP run.

## results-timecop_fails
The folder 'results-timecop_fails' contains a directory structure with text files for each TIMECOP failure.\
The compiler version and flags used are listed in the name of each text file.
The text files contain the full TIMECOP output of each failure. Since the debug flag '-g' was used in all my tests these files also contain the name of the source-code file and exact line number that failed.\

## regex.txt
The file 'regex.txt' contains the three regex strings I used to search through the SUPERCOP suite. There might be value in adjusting these to find more primitives with similar code to test.
