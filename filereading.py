import datetime
import csv
import pandas as pd
import numpy as np
import re
from collections import Counter
import pandas as pd
import glob
import os
import argparse
import sys

print("\n Started at --> {0} \n".format(datetime.datetime.now()))
def getArguments():
    # Set up the command line parser
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description="")

    parser.add_argument(
    "-v",
    "--verbose",
    action="store_true",
    help="Run the program in verbose mode.")
    
    parser.add_argument("path")
    parser.add_argument("chain")
    parser.add_argument( "--paired_sequencing", action="store_true")
    parser.add_argument( "--same_VJ", action="store_true", help="keep only seq with the same V,J genes (perfect match)")
    args = parser.parse_known_args()[0]
    if args.same_VJ:
    	parser.add_argument('--sourceSequence', dest= 'p2')
    	parser.add_argument('--Vgene', dest= 'p3')
    	parser.add_argument('--Jgene', dest= 'p4')	
    options = parser.parse_args()
    return options
# setting the path for joining multiple files


if __name__ == "__main__":

	options = getArguments()

	folder_name = options.path
	Chain = options.chain
	#print(len(os.listdir(folder_name)))

	files = os.path.join(folder_name, "*{0}*".format(Chain))

# list of merged files returned
	files = glob.glob(files)
	Master_list = pd.concat((pd.read_csv(f,low_memory = False) for f in files), ignore_index=True)
	filename = folder_name + "/Master_list_{0}.csv".format(Chain)
	if os.path.isfile(filename):
		print("ERROR! /Master_list_{0}.csv file already exists! Remove it to avoid duplicating the results!".format(Chain))
		sys.exit(1)
	#Master_list = pd.concat(map(pd.read_csv, files), ignore_index=True)
	Master_list.to_csv(filename)

try:
	sourceSequence = options.p2
except:
	print('No source sequence because no Same VJ required in the arguments')

try:
	Vgene = options.p3
except:
	print('No V gene required in the arguments')

try:
	Jgene = options.p4
except:
	print('No J gene required in the arguments')


df1 = pd.read_csv(filename,low_memory=False)
print('\n {0}'.format(df1.shape))
df1 = df1.rename(columns={'subject.diagnosis' : 'subject_diagnosis'})
df1 = df1.rename(columns={'subject.subject_id' : 'subject_subject_id'})
df1 = df1.rename(columns={'v_gene' : 'v_gene_old'})
df1 = df1.rename(columns={'j_gene' : 'j_gene_old'})

aaa =df1.assign(v_gene=df1['v_gene_old'].str.split(', or ')).explode('v_gene')

aaa =aaa.assign(j_gene=df1['j_gene_old'].str.split(', or ')).explode('j_gene')


df = aaa.drop_duplicates(subset =['subject_subject_id', 'junction_aa', 'v_gene', 'j_gene'])

patients = df
patients = patients.drop_duplicates(subset =[ "subject_subject_id"])
print('\n There are {0} different patients'.format(patients.shape[0]))

print('Dropping clonotype duplicates')
df = df.drop_duplicates(subset =[ "subject_subject_id","junction_aa","v_gene","j_gene"])

print('\n {0}'.format(df.shape))

z = []

if 'cell_id' in df.columns:

	for i in df['cell_id']:
		if not pd.isna(i):
			z.append(i)

print("{0} cell IDs were found \n".format(len(z)))

all_junction = []

all_disease = []

all_phenotype = ["None"] * len(df['sample'])

all_locus = []

all_ids = []
all_vgene = []
all_jgene = []
all_patients = []


for i in (df['locus']):
    if not i:
        locus_stat = "None"
        all_locus.append(locus_stat)

    else:
        locus_stat = i
        all_locus.append(locus_stat)


for i in df['subject_subject_id']:
	if not i:
		id_subject = "None"
		all_patients.append(id_subject)

	else:
		id_subject = i
		all_patients.append(id_subject)



for i in df['v_gene']:
	if not i:
		vgene = "None"
		all_vgene.append(vgene)

	else:
		vgene = i
		all_vgene.append(vgene)

for i in df['j_gene']:
	if not i:
		jgene = "None"
		all_jgene.append(jgene)

	else:
		jgene = i
		all_jgene.append(jgene)


if 'cell_id' in df.columns:

	for i in (df['cell_id']):
		if not i:
			id_stat = "None"
			all_ids.append(locus_stat)
		else:
			id_stat = i
			all_ids.append(id_stat)
else:
	all_ids = ["None"] * len(df['junction_aa'])


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

		if "cell_subset':" in i:
			

			phenotype_tmp = re.search("'cell_phenotype': (.+?),", i)
			
			cell_subset = phenotype_tmp.group(1)
			

			if "CD4+" in cell_subset:
				all_phenotype[index] = "CD4-positive"

			elif "CD8+" in cell_subset:
				all_phenotype[index] = "CD8-positive"

			else:
				all_phenotype[index] = "None"



		else:
			continue

	except:
		all_phenotype[index] = "None"


stat_dataframe = pd.DataFrame({'Junction': all_junction,'locus': all_locus,'cell_ID': all_ids,'Disease': all_disease,'Phenotype': all_phenotype,'V_Gene': all_vgene,'J_Gene':all_jgene,'Patients': all_patients})



paired_sequences = ["Not found"] * len(stat_dataframe['cell_ID'])

if options.paired_sequencing:
	print('\nNow searching for paired sequences')

	for index,i in enumerate(stat_dataframe['cell_ID']):
	
		for index1,j in enumerate(stat_dataframe['cell_ID']):
				
			if index != index1 and i == j and stat_dataframe['locus'][index] != stat_dataframe['locus'][index1]:
				paired = "{0} paired with {1} from {2}".format(stat_dataframe['Junction'][index],stat_dataframe['Junction'][index1],stat_dataframe['Disease'][index1])
				paired_sequences[index] = paired


stat_dataframe = pd.DataFrame({'Junction': all_junction,'locus': all_locus,'cell_ID': all_ids,'Paired':paired_sequences,'Disease': all_disease,'V_Gene': all_vgene,'J_Gene':all_jgene,'Phenotype': all_phenotype,'Patients': all_patients})

print("\n Dataframe dimensions are : ",stat_dataframe.shape)

#print(Counter(stat_dataframe["Junction"]))

#print(Counter(stat_dataframe["Disease"]))
print("\n")

if options.same_VJ:
	source_jgene = stat_dataframe[stat_dataframe['Junction']=='{0}'.format(sourceSequence)]['J_Gene']
	
	source_vgene = stat_dataframe[stat_dataframe['Junction']=='{0}'.format(sourceSequence)]['V_Gene']
	
	source_jgenelist = list(source_jgene)
	source_vgenelist = list(source_vgene)

	all_possible_vgene = []
	for i in source_vgenelist:

		try:
			if "," in i:
				i = i.replace("[","")
				if "or" in i:
					i = i.replace("or","")
				i = i.replace("]","")
				i = i.replace(" ","")
				i = i.split(',')

				for j in i:
					j = j.replace("'","")
					all_possible_vgene.append(j)
	
			else:
				all_possible_vgene.append(i)

		except:
			continue

	all_possible_vgene = set(all_possible_vgene)
	all_possible_vgene = list(all_possible_vgene)
	
	print(all_possible_vgene)

	all_possible_jgene = []
	for i in source_jgenelist:
		try:
			if "," in i:
				i = i.replace("[","")
				i = i.replace("]","")
				i = i.replace(" ","")
				i = i.split(',')
				for j in i:
					j = j.replace("'","")
					all_possible_jgene.append(j)
	
			else:
				all_possible_jgene.append(i)
		except:
			continue
	

	
	all_possible_jgene = set(all_possible_jgene)
	all_possible_jgene = list(all_possible_jgene)

	print(all_possible_jgene)

	if not all_possible_vgene:
		print("THE SOURCE SEQUENCE MIGHT NOT EXIST IN THE MASTER LIST FILE !")

	if not all_possible_vgene:
		print("THE SOURCE SEQUENCE MIGHT NOT EXIST IN THE MASTER LIST FILE !")	
	

	if options.p3:

		sameVJ_df = stat_dataframe[stat_dataframe['V_Gene'].str.contains("{0}$".format(Vgene),case=False,na=True)]
		sameVJ_df = sameVJ_df[sameVJ_df['V_Gene'].notna()]
		print("Dataframe shape after V only gene haircut {0}".format(sameVJ_df.shape))	


	if options.p4:
		if not options.p3:
			sameVJ_df = stat_dataframe
		sameVJ_df = sameVJ_df[sameVJ_df['J_Gene'].str.contains("{0}$".format(Jgene),case=False,na=True)]
		sameVJ_df = sameVJ_df[sameVJ_df['J_Gene'].notna()]
		print("Dataframe shape after J gene haircut {0}".format(sameVJ_df.shape))



	with open('{file_path}/Disease_{v_gene}_{j_gene}_{chain_type}.csv'.format(file_path = folder_name,v_gene = options.p3,j_gene = options.p4,chain_type = Chain), 'w', newline= '') as csvfile:
	
		fieldnames = [ 'Disease','Junction', 'Total_hits','Hits_CD4','Hits_CD8','TRA_hits','TRB_hits','Paired']
	
		thewriter = csv.DictWriter(csvfile, fieldnames=fieldnames)
	
		#thewriter.writeheader()
	
	
		for diseases in set(all_disease):
			
	
			count_sequences_disease = sameVJ_df[sameVJ_df['Disease'] == diseases].count()
	
		
			disease_sequence = set(sameVJ_df[sameVJ_df['Disease']==diseases]['Junction'])
	
			disease_phenotype = set(sameVJ_df[sameVJ_df['Disease']==diseases]['Phenotype'])
			
		
	
			#print("\n#### Disease : {0} \n {1} \n".format(diseases,count_sequences_disease))
	
			#print(disease_phenotype)
	
			#print(disease_sequence)
			thewriter.writeheader()
			thewriter.writerow({'Disease':"\t"})
	
			for j in disease_sequence:
	
			#print("Occurrences of {0}".format(j),len(sameVJ_df[(sameVJ_df['Junction']==j) & (sameVJ_df['Disease']==diseases)]))
				junction_occurrences = len(sameVJ_df[(sameVJ_df['Junction']==j) & (sameVJ_df['Disease']==diseases)])
	
				#pairring_tmp = len(sameVJ_df[(sameVJ_df['Junction'] == j) & (sameVJ_df['Paired'] == "None") & (sameVJ_df['Disease']==diseases)])
	
				pairring_tmp = sameVJ_df[(sameVJ_df['Junction'] == j) & (sameVJ_df['Paired'] != "Not found") & (sameVJ_df['Disease']==diseases)]
	
				pairring = pairring_tmp['Paired'].to_string(index=False)
	
				TRA = len(sameVJ_df[(sameVJ_df['Junction'] == j) & (sameVJ_df['locus'] == "TRA") & (sameVJ_df['Disease']==diseases)])
	
				TRB = len(sameVJ_df[(sameVJ_df['Junction'] == j) & (sameVJ_df['locus'] == "TRB") & (sameVJ_df['Disease']==diseases)])
	
				CD4 = len(sameVJ_df[(sameVJ_df['Junction'] == j) & (sameVJ_df['Phenotype'] == "CD4-positive") & (sameVJ_df['Disease']==diseases)])
	
				CD8 = len(sameVJ_df[(sameVJ_df['Junction'] == j) & (sameVJ_df['Phenotype'] == "CD8-positive") & (sameVJ_df['Disease']==diseases)])
				#print("{0} hits".format(junction_occurrences),"CD4: %.3f"%(CD4 / junction_occurrences),"CD8: %.3f"%(CD8 / junction_occurrences))
				if options.same_VJ:
					current_vgene = sameVJ_df[sameVJ_df['Junction']==j]['V_Gene']
					
					
					#junction_occurrences = len(sameVJ_df[(sameVJ_df['Junction']==j) & (sameVJ_df['Disease']==diseases) & (current_vgene in source_vgene)])
					#print(junction_occurrences)
	
					thewriter.writerow({ 'Disease':diseases,'Junction':j, 'Total_hits':junction_occurrences ,'Hits_CD4':CD4, 'Hits_CD8':CD8, 'TRA_hits':TRA, 'TRB_hits':TRB,'Paired':pairring})
	
				else:
					thewriter.writerow({ 'Disease':diseases,'Junction':j, 'Total_hits':junction_occurrences ,'Hits_CD4':CD4, 'Hits_CD8':CD8, 'TRA_hits':TRA, 'TRB_hits':TRB,'Paired':pairring})
	
		thewriter.writerow({'Disease':"\t"})
		
		with open('{file_path}/General_{v_gene}_{j_gene}_{chain_type}.csv'.format(file_path = folder_name,v_gene = options.p3,j_gene = options.p4,chain_type = Chain), 'w', newline= '') as csvfile:


			fieldnames_1 = [ 'Sequence', 'All_hits','CD4_count','CD8_count','Patients']
			thewriter = csv.DictWriter(csvfile, fieldnames=fieldnames_1)
			thewriter.writeheader()
			for j in set(sameVJ_df['Junction']):
				junction_occurrences_1 = len(sameVJ_df[(sameVJ_df['Junction']==j)])
				CD4 = len(sameVJ_df[(sameVJ_df['Junction'] == j) & (sameVJ_df['Phenotype'] == "CD4-positive")])
				CD8 = len(sameVJ_df[(sameVJ_df['Junction'] == j) & (sameVJ_df['Phenotype'] == "CD8-positive")])
				ww =((sameVJ_df[(sameVJ_df['Junction']==j)]))
				patients = len(set(ww.Patients))
				thewriter.writerow({ 'Sequence': j,'All_hits': junction_occurrences_1, 'CD4_count': CD4, 'CD8_count': CD8, 'Patients': patients})


else:
	with open('{file_path}/Disease_I-receptor_{chain_type}.csv'.format(file_path = folder_name,chain_type = Chain), 'w', newline= '') as csvfile:
	
		fieldnames = [ 'Disease','Junction', 'Total_hits','Hits_CD4','Hits_CD8','TRA_hits','TRB_hits','Paired']
	
		thewriter = csv.DictWriter(csvfile, fieldnames=fieldnames)
	
		#thewriter.writeheader()
	
	
		for diseases in set(all_disease):
	
			count_sequences_disease = stat_dataframe[stat_dataframe['Disease'] == diseases].count()
	
		
			disease_sequence = set(stat_dataframe[stat_dataframe['Disease']==diseases]['Junction'])
	
			disease_phenotype = set(stat_dataframe[stat_dataframe['Disease']==diseases]['Phenotype'])
			print(diseases)
		
	
			#print("\n#### Disease : {0} \n {1} \n".format(diseases,count_sequences_disease))
	
			#print(disease_phenotype)
	
			#print(disease_sequence)
			thewriter.writeheader()
			thewriter.writerow({'Disease':"\t"})
	
			for j in disease_sequence:
	
			#print("Occurrences of {0}".format(j),len(stat_dataframe[(stat_dataframe['Junction']==j) & (stat_dataframe['Disease']==diseases)]))
				junction_occurrences = len(stat_dataframe[(stat_dataframe['Junction']==j) & (stat_dataframe['Disease']==diseases)])
	
				#pairring_tmp = len(stat_dataframe[(stat_dataframe['Junction'] == j) & (stat_dataframe['Paired'] == "None") & (stat_dataframe['Disease']==diseases)])
	
				pairring_tmp = stat_dataframe[(stat_dataframe['Junction'] == j) & (stat_dataframe['Paired'] != "Not found") & (stat_dataframe['Disease']==diseases)]
	
				pairring = pairring_tmp['Paired'].to_string(index=False)
	
				TRA = len(stat_dataframe[(stat_dataframe['Junction'] == j) & (stat_dataframe['locus'] == "TRA") & (stat_dataframe['Disease']==diseases)])
	
				TRB = len(stat_dataframe[(stat_dataframe['Junction'] == j) & (stat_dataframe['locus'] == "TRB") & (stat_dataframe['Disease']==diseases)])
	
				CD4 = len(stat_dataframe[(stat_dataframe['Junction'] == j) & (stat_dataframe['Phenotype'] == "CD4-positive") & (stat_dataframe['Disease']==diseases)])
	
				CD8 = len(stat_dataframe[(stat_dataframe['Junction'] == j) & (stat_dataframe['Phenotype'] == "CD8-positive") & (stat_dataframe['Disease']==diseases)])
				#print("{0} hits".format(junction_occurrences),"CD4: %.3f"%(CD4 / junction_occurrences),"CD8: %.3f"%(CD8 / junction_occurrences))
				if options.same_VJ:
					current_vgene = stat_dataframe[stat_dataframe['Junction']==j]['V_Gene']
					
					
					#junction_occurrences = len(stat_dataframe[(stat_dataframe['Junction']==j) & (stat_dataframe['Disease']==diseases) & (current_vgene in source_vgene)])
					#print(junction_occurrences)
	
					thewriter.writerow({ 'Disease':diseases,'Junction':j, 'Total_hits':junction_occurrences ,'Hits_CD4':CD4, 'Hits_CD8':CD8, 'TRA_hits':TRA, 'TRB_hits':TRB,'Paired':pairring})
	
				else:
					thewriter.writerow({ 'Disease':diseases,'Junction':j, 'Total_hits':junction_occurrences ,'Hits_CD4':CD4, 'Hits_CD8':CD8, 'TRA_hits':TRA, 'TRB_hits':TRB,'Paired':pairring})
	
		thewriter.writerow({'Disease':"\t"})
		
		with open('{file_path}/General_I-receptor_{chain_type}.csv'.format(file_path = folder_name,chain_type = Chain), 'w', newline= '') as csvfile:



			fieldnames_1 = [ 'Sequence', 'All_hits','CD4_count','CD8_count','Patients']
			thewriter = csv.DictWriter(csvfile, fieldnames=fieldnames_1)
			thewriter.writeheader()
			sameVJ_df = stat_dataframe
			for j in set(sameVJ_df['Junction']):
				junction_occurrences_1 = len(sameVJ_df[(sameVJ_df['Junction']==j)])
				CD4 = len(sameVJ_df[(sameVJ_df['Junction'] == j) & (sameVJ_df['Phenotype'] == "CD4-positive")])
				CD8 = len(sameVJ_df[(sameVJ_df['Junction'] == j) & (sameVJ_df['Phenotype'] == "CD8-positive")])
				ww =((sameVJ_df[(sameVJ_df['Junction']==j)]))
				patients = len(set(ww.Patients))
				thewriter.writerow({ 'Sequence': j,'All_hits': junction_occurrences_1, 'CD4_count': CD4, 'CD8_count': CD8, 'Patients': patients})


		
	



	

print("\n Finished at --> {0} \n ".format(datetime.datetime.now()))



#for phenotypes in set(all_phenotype):

	#count_sequences_phenotype = stat_dataframe[stat_dataframe['Phenotype'] == phenotypes].count()


	#print("\n#### Phenotype : {0} \n {1} \n".format(phenotypes,count_sequences_phenotype),set(stat_dataframe[stat_dataframe['Phenotype']==phenotypes]['Junction']))





