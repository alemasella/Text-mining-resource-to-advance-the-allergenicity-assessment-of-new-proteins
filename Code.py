import csv
import subprocess
import os
from os.path import exists
import numpy as np
from collections import Counter
import pandas as pd


list_of_lists = []
counter = 0
Rows_wanted = int(input("Write the number of rows to visualize:")) + 1

with open("HPO_litsympt-morgs_0012393_add.tsv") as fd1:
    HPOs = csv.reader(fd1, delimiter="\t", quotechar='"')
    
    for row in HPOs:
        counter += 1
        
        if counter % Rows_wanted == 0:
            break
        
        HPOrow = row[0]
        MESHrow = row[2]
        
        S_HPOrow = HPOrow.split(":")
        S_MESHrow = MESHrow.split(":")
        
        D_HPOrow = S_HPOrow[0] + "_" + S_HPOrow[1]
        D_MESHrow = S_MESHrow[0] + "_" + S_MESHrow[1]
        
        lis = [D_HPOrow, D_MESHrow]
        
        try:
            with open("0012393_HPO_2024_PMIDs/" + D_HPOrow + ".pmids") as fd2:
                HPO = fd2.read().rstrip("\n\n").split("\n")
        except FileNotFoundError:
            with open("0012393_lit_sympt_u_nq_2024_PMIDs/" + HPOrow + ".pmids") as fd3:
                HPO = fd3.read().rstrip("\n\n").split("\n")
        
        with open("0012393_mesh_orgs_PMIDs/" + D_MESHrow + ".pmids") as fd4:
            Org = fd4.read().rstrip("\n\n").split("\n")
            
            for row1 in HPO:
                for row2 in Org:
                    if row1 == row2:
                        lis.append(row1)
                        break  # Assuming row1 is unique in both HPO and Org files
                else:
                    continue
        
        list_of_lists.append(lis)

#print(list_of_lists)


list_years=[]
for i in range(len(list_of_lists)):
	wanted_list=list_of_lists[i]
	years2=[wanted_list[0],wanted_list[1]]
	for k in range(2,len(wanted_list)):
		if not exists("PMIDs/"+wanted_list[k]+".pmid"):
			url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pubmed&id="+ wanted_list[k] + "&rettype=medline"
			result = subprocess.run(["wget", url, "-O", "PMIDs/"+wanted_list[k]+".pmid", "-o", "/dev/null"])

		with open("PMIDs/"+ wanted_list[k]+".pmid") as cd4:
			XML=cd4.read().lstrip("\n").split("\n")

		for j in range(min(15, len(XML))):
			XMLj=XML[j].split(" ")
			if XMLj[0]=="DP":
				Y=int(XMLj[3])
			else:
				continue
		years2.append(Y)
	list_years.append(years2)
#print(list_years)



for i in range(len(list_years)):
	list_considered=list_years[i]
	only_years=[]
	for k in range(2,len(list_considered)):
		only_years.append(list_considered[k])
	Y = np.array(only_years)
	co=Counter(Y)
	X=np.array(list(co.items()))
	total_years=list(range(min(only_years),max(only_years)+1))
	total_frec=[]
	years_cons=[]
	frecuencies_cons=[]
	for i in range(len(X)):
		rowYF=X[i]
		years_cons.append(rowYF[0])
		frecuencies_cons.append(rowYF[1])
	for year in total_years:
		if year in years_cons:
			ind=years_cons.index(year)
			total_frec.append(frecuencies_cons[ind])
		else:
			total_frec.append(0)
	x = np.array(total_years)
	y = np.array(total_frec)
	slope, intercept = np.polyfit(x, y, 1)
	if slope > 0.01:
		print("For the pair ("+ list_considered[0]+","+list_considered[1]+"), the relevancy of the relation has increased in the years")
	elif slope < -0.01:
		print("For the pair ("+ list_considered[0]+","+list_considered[1]+"), the relevancy of the relation has decreased in the years")
	else:
		print("For the pair ("+ list_considered[0]+","+list_considered[1]+"), the tendency in the relevancy of the relation is not clear")
#Note that the criterion of selection in the slope is arbitrary, but justified with the limits of naked eye vision
