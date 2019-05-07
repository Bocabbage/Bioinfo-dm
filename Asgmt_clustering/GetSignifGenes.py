SIGNIFICANT_LEVEL = 0.01

import pandas as pd

rfile = open('E:/Programming/Dataset/GEO/GSE33532/GEO2R_Result.txt')
rfile.readline()
results = []
for line in rfile.readlines():
    temp = line.split("\"")
    if float(temp[3]) >= 0.01:
        break
    if temp[13]!='' and temp[13] not in results:
        results.append(temp[13])
rfile.close()

pd.DataFrame(results).to_csv('E:/Programming/Dataset/GEO/GSE33532/SignificantGenes.csv',sep=',')