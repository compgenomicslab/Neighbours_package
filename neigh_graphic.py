from __future__ import division
import json
from pymongo import MongoClient
from collections import Counter
import os,sys


import mongo_client



### connection to mongo client using mongo_client.py ###
coll_unigenes = mongo_client.mongo_connect()[0]
coll_clusters = mongo_client.mongo_connect()[1]
coll_e5 = mongo_client.mongo_connect()[2]



if len(sys.argv) == 2:
	if sys.argv[1] == "-h"	or sys.argv[1] == "-help":
		sys.exit("\n"+"# unigenes cluster graphication #" +"\n" +
			 	 "\n" + "Usage: neigh_graphic.py <number_of_neighbours_genes_to_display> <unigene_cluster>"+"\n"+
			 	 "\n" + "Example: neigh_graphic.py 2 000_000_005"+"\n"
				)

elif len(sys.argv) < 2:

		sys.exit("\n"+"# unigenes cluster graphication #" +"\n" +
					 "\n" + "Error: few arguments" +
					 "\n" + "Usage: neigh_graphic.py <number_of_neighbours_genes_to_display> <unigene_cluster>"+"\n"+
					 "\n" + "Example: neigh_graphic.py 2 000_000_005"+"\n"
				)
else:
	try:
		number_of_neigh = sys.argv[1]
		unigene_cluster = sys.argv[2]

	except:
			sys.exit("\n"+"# unigenes cluster graphication #" +"\n" +
					 "\n" + "Error: many arguments" +
					 "\n" + "Usage: neigh_graphic.py <number_of_neighbours_genes_to_display> <unigene_cluster>"+"\n"+
					 "\n" + "Example: neigh_graphic.py 2 000_000_005"+"\n"
					)



### kegg pathways containing file (KEGGs_pathways.txt) 
### is parsed to generate KEEG=description hash
kegg_pathways = open("data/KEGGs_pathways.txt","r")


def make_kegg_dict(kegg_pathways):
    """generates kegg pathways hash which contains kegg descriptions"""
    global kegg_dict
    kegg_dict = {}


    for line in kegg_pathways:
       fields = line.strip("\n").split("\t")
       kegg= fields[0]
       description= " ".join(fields[1::]).rstrip(" ")
       kegg_dict[kegg]= description

    return kegg_dict



def get_kegg_description(kegg):
	"""retrieve kegg description from kegg_dict"""
	description = kegg_dict[kegg]
	return description



def mongo_orf_find(GMGC):
	""" retrieve information of an orf element (strand,locus,start,end) from a unigene cluster """

	GMGC_cluster = coll_unigenes.find({"u":GMGC})
	GMGC_dict = {}

	for orf in GMGC_cluster:
		orf = orf['o']

		for a in orf:

			gene = a['g']
			Locus = a['s']
			start = Locus[0]
			end = Locus[1]
			strand = Locus[2]

			try:
				#GMGC_dict[gene]=[contig,start,end,strand,sample]
				GMGC_dict[gene]=[start,end,strand]
				#print [contig,start,end,strand,sample]

			except:
				print("hay un problema")
	return GMGC_dict

#gmgc_cluster = mongo_orf_find("000_000_000")


def mongo_functional_find(GMGC): 
        """ retrieve functional information from mongo.cluster database """
        
        kegg_list = []
        Egg_list = []

        GMGC_function = coll_clusters.find({"u":GMGC})
     
        GMGC_function_list =[]
        for element in GMGC_function:

                kegg = element['K_P']
                kegg = kegg.split(",")

                for n in kegg: #depuro la lista de kegg para quedarme solo con los kegg pathways
                    if n in kegg_dict and n not in kegg_list:
                        kegg_list.append(n)

                Egg = element['OGs'].split(",")
                #print(Egg)
                for n in Egg:
                    if n.split("@")[1] == "1":
                        Egg.remove(n)
                #print(Egg)

                try:
                    #GMGC_function_list=[kegg_list,cog,eggnog,gene_symbol,desc]
                    GMGC_function_list=[kegg_list,Egg]
                except:
                    print("hay un problema")
        return GMGC_function_list



def retrieve_gmgc(gene):

	""" retrieve the gmgc of an unigene orf in mongo database """
	GMGC_function = coll_unigenes.find({ "o.g":gene},{"u":1}).limit(1)

	selected = []
	gmgc = ""
	for n in GMGC_function:

		try:
			gmgc = n["u"]
		except:
			gmgc = ""

	return gmgc




def retrieve_neighbours_data(gmgc_cluster,number):
	""" retrieve a list with the four genes -2-1+1+2 which 
	sourround the query unigene orf """
	
	neigh_dict = {}
	for k,v in gmgc_cluster.items():
		gene_ordered =[]

		query_gene = k
		start = v[0]
		end = v[1]
		strand = v[2]
		orf = k.split("_")
		gene = int(orf[3])
		sample_contig =orf[0:3]

		gene_list =range(-number,number+1,1)

		for gene_pos in gene_list:
			sample_cont = []
			genes =str(int(gene)+int(gene_pos))+"#"+strand
			sample_cont = sample_contig[0:3]
			sample_cont.append(str(genes))
			genes = "_".join(sample_cont)
			gene_ordered.append(genes)
			neigh_dict[query_gene]=gene_ordered


	return neigh_dict


def get_kegg_description(Egg):
	""" connect to mongo eggnog5 database 
	and retrieve description of Egnog """
	Egg = Egg.split("@")[0]
	e5_database = coll_e5.find({"e":Egg})
	description = ""

	for element in e5_database:
		description = element["d"]


	return description



def printing_genes(gmgc_list_func):
	""" this function controls ouput gene graphication """

	# draw elements for graphic output
	positive=" [{}]>"
	negative="<[{}] "
	CRED = '\033[91m'
	CYEW = '\033[33m'
	CEND = '\033[0m'
	CBLUE = '\033[44m'
	CGRE = '\033[92m'
	CGRE2 = '\033[42m'
	CGREY = '\033[100m'


	# variables
	list_for_draw=[]
	kegg_repeated=[]
	keggs_for_draw=[]
	egg_repeated=[]
	egg_for_draw=[]
	

	for n in gmgc_list_func:
		if n !="":
			

			keggs = n[0]
			strand = n[1]
			unigene = n[2]
			eggnog_list = n[3]


			final_kegg = unigene+":"
			for kegg in keggs:
				final_kegg = final_kegg+" "+kegg
				if kegg in kegg_dict:
					if kegg not in kegg_repeated:
						kegg_repeated.append(kegg)
						kegg_description=kegg_dict[kegg]
						keggs_for_draw.append((CBLUE+str(kegg)+CEND)+" : "+kegg_description)

			for egg in eggnog_list:
				if egg not in egg_repeated:

						egg_repeated.append(egg)
						egg_description = get_kegg_description(egg)
						egg_for_draw.append((CGRE2+str(egg)+CEND)+" : "+egg_description)

			egg = " ".join(eggnog_list)

			if strand == "+":
				positive =" [{}]>".format(final_kegg+" "+(CGRE+egg+CEND))
				list_for_draw.append(positive)
			
			else:
				negative="<[{}] ".format(final_kegg+" "+(CGRE+egg+CEND))
				list_for_draw.append(negative)
		else:
			list_for_draw.append(CRED+"[X]"+CEND)

	#printing the line
	line=""
	count = 0
	for element in list_for_draw:
		count +=1
		if element != CRED+"[X]"+CEND:
			if len(element.split(":")[1]) > 3:
				if count == 3:
					line+=(CGREY+" "+element+CEND)	
						
					#line+=(CRED+" "+element.replace("]","")+CEND+CRED+"]"+CEND)	
				else:
					line+=(CYEW+" "+element+CEND)	
		else:
			line+=" "+element
	print(line)

	for kegg in keggs_for_draw:
		print(kegg)
	for egg in egg_for_draw:	
		print(egg)
		
		

	print("\n")


def get_unigene_functional_data(neigh_gene_list,strand):

	#  neigh_gene_list is a depured list (contains uniq syntenies) 
	#  of ORFs and sourrounded neigh orfs

	unigene_list_description = []
	for unigene in neigh_gene_list:
		unigene_cogs = mongo_functional_find(str(unigene)) # [[keggs][Eggnogs]]	
		
		unigene_description=[]
		if unigene_cogs: 
		
			eggnog = unigene_cogs[1]
			kegg = unigene_cogs[0]
		
			# add kegg, strand, unigene and eggnogs to unigene description
			unigene_description.append(kegg)
			unigene_description.append(strand)
			unigene_description.append(unigene)
			unigene_description.append(eggnog)

		else:
			unigene_description=""

		#unigene_list_description contains the 4 neihbours genes and the orfs with their strand and cog information
		unigene_list_description.append(unigene_description)
		#print(unigene_list_description)
	return unigene_list_description





def synteny_based_on_keggs(neigh_gene_list):
	
	""" generates a list of keggs that belong to neighbours genes to filter by kegg synteny"""

	unigene_list_of_keggs = []
	for unigene in neigh_gene_list:
		unigene_cogs = mongo_functional_find(str(unigene)) # [[keggs][Eggnogs]]	
		
		unigene_description=[]
		if unigene_cogs: 
			kegg = unigene_cogs[0]
			kegg = "_".join(kegg)
			#print unigene_description
		else:
			kegg=""

		unigene_list_of_keggs.append(kegg)
	unigene_list_of_keggs.sort()	

	return unigene_list_of_keggs



def main(unigene_cluster,number_of_neigh):
	""" number of neigh is the number of genes you want to 
	show in the ouput, i.e ussing 2 you retrieve ,-2,-1,+1,+2"""
	print("\n")

	## VARIABLES AND LISTs
	unique_neigh_syntenies=[]
	keggs_synteny = []
	kegg_list = []
	gmgc_orf_dict = {}
	number_neigh = 0
	neigh_with_keggs = 0
	analysed_orfs = 0
	keggs_for_draw_list=[] #retrieve stats of genes with keggs 


	#retrieve orfs fromunigene cluster {"orf":[start,end,strand]...}
	unigene_cluster = mongo_orf_find(unigene_cluster)
	
	#retrieve neighbours for evey orfs and generate a  dictionary orf, neighbours orfs
	# {"orf":[unigene-2,unigene-1,unigen_from_query_cluster,unigene,unigen+1,unigene+2]}
	unigene_cluster_neighbours = retrieve_neighbours_data(unigene_cluster,int(number_of_neigh))


	#iteration over each orfs to analyze the synteny of Keggs
	for query_gene, neighbour_gene_list in unigene_cluster_neighbours.items():
		analysed_orfs += 1
		gmgc_list_func = []
		neigh_gene_list=[]
		neigh_gene_list_for_synteny = []
		
		#analyze each neighbour ORF
		for ORF in neighbour_gene_list:
			strand = ORF.split("#")[1] # unigenes strands 
			neigh_orf = ORF.split("#")[0] #unigene
			unigene = retrieve_gmgc(neigh_orf)
			neigh_gene_list.append(unigene)
			neigh_gene_list_for_synteny.append(unigene)

		# filtering by unigene synteny composition (first level)
		neigh_gene_list_for_synteny.sort()
		if neigh_gene_list_for_synteny not in unique_neigh_syntenies:
			unique_neigh_syntenies.append(neigh_gene_list_for_synteny)
			
			# filtering by keggs synteny composition (second level)
			list_of_keggs = synteny_based_on_keggs(neigh_gene_list)		
			if list_of_keggs not in keggs_synteny:
				keggs_synteny.append(list_of_keggs)

				#almacenamos el resultado de los cogs de todas las sintenias para graficar
				gmgc_list_func = get_unigene_functional_data(neigh_gene_list,strand)			
			
				#use printing_genes function to show funcy results in terminal
				printing_genes(gmgc_list_func)



# we create kegg_pathways dictionary
make_kegg_dict (kegg_pathways)


if __name__ == main(unigene_cluster,2):
	main(unigene_cluster)
