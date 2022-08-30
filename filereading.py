import datetime
import csv
import pandas as pd
import numpy as np
import re
from collections import Counter

print("\n Started at --> {0} \n".format(datetime.datetime.now()))

filename = '/Users/daniel/Desktop/Master-Project/ireceptor_script-main/results_test2/TRB/data_TRA_CSARSGGNTGELFF.csv'

df = pd.read_csv(filename,usecols=['subject.diagnosis','sample','junction_aa','locus','cell_id','sequence_id'])
df = df.rename(columns={'subject.diagnosis' : 'subject_diagnosis'})
### Remove the word junction from some rows because of the merge between all csv files
df = df[~df.junction_aa.str.contains("junction_aa")]
# Remove non floating type because of the error you get when you do (if label in)
#df = df.loc[df.subject_diagnosis.apply(type) != float]
print(df.shape)
### Remove duplicates
df = df.drop_duplicates(keep = 'first')


all_junction = []

all_disease = []

all_phenotype = []

all_junction_phenotype = []

all_locus = []

all_ids = []


for i in (df['locus']):
    if not i:
        locus_stat = "None"
        all_locus.append(locus_stat)

    else:
        locus_stat = i
        all_locus.append(locus_stat)



for i in (df['cell_id']):
    if not i:
        id_stat = "None"
        all_ids.append(locus_stat)

    else:
        id_stat = i
        all_ids.append(id_stat)


for index,i in enumerate(df['subject_diagnosis']):

	try:
		if "label': '" in i:
			disease = re.search("'label': '(.+?)'", i)
			junction_amino = df["junction_aa"].iloc[index]
			all_junction.append (junction_amino)
			all_disease.append((disease.group(1)))
		else:
			disease = "None"
			junction_amino = df["junction_aa"].iloc[index]
			all_junction.append (junction_amino)
			all_disease.append(disease)
	except:
		disease = "None"
		junction_amino = df["junction_aa"].iloc[index]
		all_junction.append (junction_amino)
		all_disease.append(disease)

for index,i in enumerate(df['sample']):
	
	try:

		if "cell_phenotype':" in i:

			phenotype = re.search("'cell_phenotype': (.+?),", i)
		
			junction_amino = df["junction_aa"].iloc[index]

			all_phenotype.append((phenotype.group(1)))
			all_junction_phenotype.append(junction_amino)

		else:
			phenotype = "None"
			junction_amino = df["junction_aa"].iloc[index]

			all_phenotype.append(phenotype)
			all_junction_phenotype.append(junction_amino)

	except:
		phenotype = "None"
		junction_amino = df["junction_aa"].iloc[index]

		all_phenotype.append(phenotype)
		all_junction_phenotype.append(junction_amino)

#print(set(all_phenotype))

stat_dataframe_phenotype = pd.DataFrame({'Junction': all_junction_phenotype,'Phenotype': all_phenotype})

#print("##### Phenotype",stat_dataframe_phenotype.shape)

print(set(all_disease))

stat_dataframe = pd.DataFrame({'Junction': all_junction,'locus': all_locus,'cell_ID': all_ids,'Disease': all_disease,'Phenotype': all_phenotype})

paired_sequences = ["Not found"] * len(stat_dataframe['cell_ID'])

for index,i in enumerate(stat_dataframe['cell_ID']):

	for index1,j in enumerate(stat_dataframe['cell_ID']):
			
		if index != index1 and i == j and stat_dataframe['locus'][index] != stat_dataframe['locus'][index1]:
			paired = "{0} paired with {1} from {2}".format(stat_dataframe['Junction'][index],stat_dataframe['Junction'][index1],stat_dataframe['Disease'][index1])
			paired_sequences[index] = paired




stat_dataframe = pd.DataFrame({'Junction': all_junction,'locus': all_locus,'cell_ID': all_ids,'Paired':paired_sequences,'Disease': all_disease,'Phenotype': all_phenotype})



print("\n Dataframe dimensions are : ",stat_dataframe.shape)

#print(Counter(stat_dataframe["Junction"]))

#print(Counter(stat_dataframe["Disease"]))
print("\n")

with open('filereading00000000001.csv', 'w', newline= '') as csvfile:

	fieldnames = [ 'Disease','Junction', 'Total_hits','Hits_CD4','Hits_CD8','TRA_hits','TRB_hits','Paired']

	thewriter = csv.DictWriter(csvfile, fieldnames=fieldnames)

	thewriter.writeheader()


	for diseases in set(all_disease):

		count_sequences_disease = stat_dataframe[stat_dataframe['Disease'] == diseases].count()

	
		disease_sequence = set(stat_dataframe[stat_dataframe['Disease']==diseases]['Junction'])

		disease_phenotype = set(stat_dataframe[stat_dataframe['Disease']==diseases]['Phenotype'])
		print(diseases)
	

		#print("\n#### Disease : {0} \n {1} \n".format(diseases,count_sequences_disease))

		#print(disease_phenotype)

		#print(disease_sequence)
		thewriter.writerow({'Disease':"\t"})

		for j in disease_sequence:

		#print("Occurrences of {0}".format(j),len(stat_dataframe[(stat_dataframe['Junction']==j) & (stat_dataframe['Disease']==diseases)]))
			junction_occurrences = len(stat_dataframe[(stat_dataframe['Junction']==j) & (stat_dataframe['Disease']==diseases)])

			#pairring_tmp = len(stat_dataframe[(stat_dataframe['Junction'] == j) & (stat_dataframe['Paired'] == "None") & (stat_dataframe['Disease']==diseases)])

			pairring_tmp = stat_dataframe[(stat_dataframe['Junction'] == j) & (stat_dataframe['Paired'] != "Not found") & (stat_dataframe['Disease']==diseases)]

			pairring = pairring_tmp['Paired'].to_string(index=False)

			TRA = len(stat_dataframe[(stat_dataframe['Junction'] == j) & (stat_dataframe['locus'] == "TRA") & (stat_dataframe['Disease']==diseases)])

			TRB = len(stat_dataframe[(stat_dataframe['Junction'] == j) & (stat_dataframe['locus'] == "TRB") & (stat_dataframe['Disease']==diseases)])

			CD4 = len(stat_dataframe[(stat_dataframe['Junction'] == j) & (stat_dataframe['Phenotype'] == "'CD4+'") & (stat_dataframe['Disease']==diseases)])

			CD8 = len(stat_dataframe[(stat_dataframe['Junction'] == j) & (stat_dataframe['Phenotype'] == "'CD8+'") & (stat_dataframe['Disease']==diseases)])
			#print("{0} hits".format(junction_occurrences),"CD4: %.3f"%(CD4 / junction_occurrences),"CD8: %.3f"%(CD8 / junction_occurrences))
			thewriter.writerow({ 'Disease':diseases,'Junction':j, 'Total_hits':junction_occurrences ,'Hits_CD4':CD4, 'Hits_CD8':CD8, 'TRA_hits':TRA, 'TRB_hits':TRB,'Paired':pairring})





	

print("\n Finished at --> {0} \n ".format(datetime.datetime.now()))



#for phenotypes in set(all_phenotype):

	#count_sequences_phenotype = stat_dataframe[stat_dataframe['Phenotype'] == phenotypes].count()


	#print("\n#### Phenotype : {0} \n {1} \n".format(phenotypes,count_sequences_phenotype),set(stat_dataframe[stat_dataframe['Phenotype']==phenotypes]['Junction']))





