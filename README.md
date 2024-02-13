# Raw results
The folder 'raw_results' contains renamed 'data' files as they are provided upon completion of a SUPERCOP run.
Note that which version/flags used are listed in the file names.

# Parsing scripts
The python files can be used to sort and transform the raw data files.
To be of use in parsing new data files, compiler versions or flags the scripts will need to be adjusted as they were hardcoded to parse my specific raw data files and flag combinations.

# Results
## results.csv
The file 'results.csv' contains the result of every benchmark including each compiler version and flag combination and whether it passed, failed or errored on the TIMECOP run.

## results-timecop_fails
The folder 'results-timecop_fails' contains a directory structure with text files for each TIMECOP failure.
The text files contain the full TIMECOP output. If the debug flag '-g' was used then this should also contain the source file name and line number that failed.
The compiler version and flags used are listed in the file name.
