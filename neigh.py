## Neigh packag v2 ##
# Joaquin Giner Lamia

from __future__ import division
from pymongo import MongoClient
from collections import Counter
import os,sys
import json

## self packages
import mongo_client
from neigh_functions import *


### connection to mongo client using mongo_client.py ###

coll_unigenes = mongo_client.mongo_connect()[0]
coll_clusters = mongo_client.mongo_connect()[1]
coll_e5 = mongo_client.mongo_connect()[2]


###### Parametters for the analysis #####

maximum_gmgc_genes = int(400) # maximun number of unigene allowed by gmgc cluster to be computed, smaller the number smaller computing time 
percentage_cut_off = int(30) # percentage of cogs (porcentaje de keggs necesario para considerar un kegg como posible candidato entre todos los neighbours keggs)
neighbours_range = 2 # number of neighbours genes to analyze around each unigene

#### input file name 
global input_unigene_file
input_unigene_file =""




## data input

gmgc_list = []

try: 
	sys.argv[1] == True
except:
	sys.exit("\nError: please use one of these options\nUsage:  neigh.py -u <unigene cluster> or \n\tneigh.py -f <file(unigenes_list).txt>")
	

if sys.argv[1] == "-f":
	try:
		argument = "-f"
		input_unigene_file= sys.argv[2]
		input_file = open(input_unigene_file,"r")
		for line in input_file:
			gmgc_list.append(line.strip())
	except:
		sys.exit("\nError: please add a unigene\nUsage: neigh.py -u <unigene> or \nneigh.py -f <file(unigenes_list).txt>")


elif sys.argv[1] == "-u":
	try:
		argument = "-u"
		input_gmgc= sys.argv[2]
		gmgc_list.append(input_gmgc)
	except:
		sys.exit("\nError: please add a unigene\nUsage: neigh.py -u <unigene> or \nneigh.py -f <file(unigenes_list).txt>")
		




##### NEIGH FUNCTIONAL ANALYSIS ####
#note gmgc means unigene cluster

def neigh_analysis(gmgc_list,neighbours_range):

	
	## RESULT LIST FOR OUTPUT WRITE 
	result_kegg_for_output_file=[]
	result_egg_for_output_file=[]
	result_unigenes_for_output_file=[]

	for gmgc in gmgc_list:

		if "GMGC" in gmgc:
			gmgc = clean_unigene(gmgc)

		## RESULT LIST
		result_list_kegg = []
		result_list_egg = []
		list_of_neigh_unigenes = []

		
		
		#Retrieve neighbours genes for every ORFs in the GMGC cluster
		#retrieve orfs sample_contig_gene that belong to unigene cluster
		gmgc_orfs_cluster = mongo_orf_find(gmgc,maximum_gmgc_genes,coll_unigenes) 
		
		#retrieve neighbours genes for every orfs cluster, by default -2-1+1+2
		neighbours_genes= retrieve_neighbours_data(gmgc_orfs_cluster,neighbours_range) 
		

		count_for_list_of_unigenes = 0 #Creates the unigenes list only one time 
		for cog in ("kegg","Egg"):
			cog_category = cog
			count_for_list_of_unigenes +=1

			#retrieve cogs assignation list and neighbourhood parametters 
			cogs_assignation_list = compute_neigh_cogs_assignation(neighbours_genes,neighbours_range,coll_unigenes,coll_clusters,cog_category)
			cogs_organization_list = cogs_assignation_list[0]
			analysed_orfs = cogs_assignation_list[1]
			number_neigh = cogs_assignation_list[2]
			number_neigh_with_cogs = cogs_assignation_list[3]
			unigenes_functions_list = cogs_assignation_list[4]
		
			# number of genes that have at least one cog in his neighobours 
			neigh_orf_with_cogs = get_neigh_orf_with_cogs(cogs_organization_list) 
		

			# unigene cluster and their neighbours unigenes for future analysis
			if count_for_list_of_unigenes < 2:
				neigh_unigenes = create_unigenes_lists_for_print(gmgc,neighbours_genes,neighbours_range,coll_unigenes)
				list_of_neigh_unigenes.append(neigh_unigenes)
			
				if list_of_neigh_unigenes[0].rstrip() == gmgc:
					list_of_neigh_unigenes= [gmgc,"singleton"]

			###
			gmgc_functional = mongo_functional_find(gmgc,coll_clusters)# retrieve functional information from GMGC query	 
			
			if cog_category == "kegg":
				query_list = gmgc_functional[0]
			else:
				query_list = gmgc_functional[1]



			# compute total number of cogs and unique cogs in the unigene cluster
			Count = Counter(unigenes_functions_list)
			unique_cogs = 0
			count_of_cogs = 0
			for k,v in Count.items():
				count_of_cogs += int(v)
				unique_cogs +=1


			# This section compute if any of the cogs are highly representated in the neighbourhood of the unigenes
			# In the case there are cogs overepresentated they are storage in kegg_description_dict
			kegg_description_dict = {}
			subject_cog_list = [] # store subject list of cog
			subject_cog_dict = {} #almacenamos todos los kegg que pasen el porcentaje para compararlos con los query_kegg_list
			for k,v in Count.items(): # v is the number of time this function has been found in the neighbours genes 
				cog = k
				KEGG_count = str(v)
				percentage = float(v/neigh_orf_with_cogs)*100 # compute percentage of number of neighbours genes with cogs
				percentage = ("{0:.2f}".format(percentage))

				if float(percentage) >= percentage_cut_off: 
					if cog_category == "kegg":   # control for cog_category output
						kegg_description_dict[cog] = get_kegg_description(cog)
					else:
						kegg_description_dict[cog] = get_egg_description(cog,coll_e5)
					
					subject_cog_dict[cog]=percentage
					subject_cog_list.append(cog)

			# Get funciontal scores and cogs description of enriched cogs
			if kegg_description_dict != {}:
				# get functional description and COG hits
				description_list = get_functional_data(subject_cog_dict, query_list, kegg_description_dict, cog_category)
				hits, query_list, description_list, hit_kegg_percentage, num_query_cog =  [description_list[i] for i in (0,1,2,3,4)]
				
				# get functional scores
				score_list = neigh_scores(hits,num_query_cog,subject_cog_list,unique_cogs,count_of_cogs,neigh_orf_with_cogs)
				positive_value, accuracy_value, func_conservation = [score_list[i] for i in (0,1,2)]

				### LIST FOR OUTPUT SECTION ###
				result =[gmgc,
						",".join(query_list),
						",".join(subject_cog_list),
						str(analysed_orfs),
						str(number_neigh),
						str(number_neigh_with_cogs),
						str(unique_cogs),
						str(count_of_cogs),
						str(func_conservation),
						",".join(hit_kegg_percentage),
						";".join(description_list)
						]

				if cog_category == "kegg":
					result_list_kegg.append(result)
				else:
					result_list_egg.append(result)
				

			else:
				result = [gmgc,"No Match"]
			
				if cog_category == "kegg":
					result_list_kegg.append(result)
				else:
					result_list_egg.append(result)

			result_kegg_for_output_file.append(result_list_kegg)
			result_egg_for_output_file.append(result_list_egg)
			result_unigenes_for_output_file.append(list_of_neigh_unigenes)

	return [result_kegg_for_output_file,result_egg_for_output_file,result_unigenes_for_output_file]





def main(gmgc_list):
	global kegg_dict
	kegg_pathways = open("data/KEGGs_pathways.txt","r")
	kegg_dict = make_kegg_dict(kegg_pathways)
	result_list = neigh_analysis(gmgc_list,neighbours_range)
	kegg_result = result_list[0]
	egg_result = result_list[1]
	list_of_neighbours_unigenes = result_list[2]

	print_results(argument, kegg_result, egg_result, list_of_neighbours_unigenes,input_unigene_file)



if __name__ == main(gmgc_list):
	main(gmgc_list)





