#!/usr/bin/python3

import sys
import pandas as pd
import gzip

df = pd.read_csv(sys.argv[1], compression="gzip", sep='\t', header=None)
print('average coverage: ', int(df.iloc[:,[2]].mean()))
print('median coverage: ', int(df.iloc[:,[2]].median()))
#print(int(df.iloc[:,[2]].mean()), '\t', int(df.iloc[:,[2]].median()))
