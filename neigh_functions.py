import os
### Parseamos el fichero con las descripciones de los KEGGs pathways
### y generamos un diccionario para almacenarlos.
kegg_pathways = open("data/KEGGs_pathways.txt","r")

def make_kegg_dict(kegg_pathways):
    """generate kegg pathway dictionary containning kegg descriptions"""
    global kegg_dict
    kegg_dict = {}


    for line in kegg_pathways:
       fields = line.strip("\n").split("\t")
       kegg= fields[0]
       description= " ".join(fields[1::]).rstrip(" ")
       kegg_dict[kegg]= description

    return kegg_dict



def get_kegg_description(kegg):
        """retrieve kegg description from kegg_dict hash """
        description = kegg_dict[kegg]
        
        return description


def get_egg_description(Egg, coll_e5):
	""" connect to mongo eggnog5 database 
	and retrieve description of Egnog """
	Egg = Egg.split("@")[0]
	e5_database = coll_e5.find({"e":Egg})
	description = ""

	for element in e5_database:
		description = element["d"]

	return description


def mongo_orf_find(GMGC,maximum_gmgc_genes,coll_unigenes):
        """ retrieve orf genomic information"""

        GMGC_cluster = coll_unigenes.find({"u":GMGC})
        GMGC_dict = {}
        gene_count = 0

        for orf in GMGC_cluster:
                orf = orf['o']     
                for a in orf:
                        
                        if  gene_count < maximum_gmgc_genes:
                                gene_count +=1
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
                                        print("something was wrong")
                else:
                        break

        return GMGC_dict










def retrieve_gmgc(gene, coll_unigenes):
        """ retrieve gmgc for every orf from unigenes db """
        GMGC_function = coll_unigenes.find({ "o.g":gene},{"u":1}).limit(1)

        selected = []
        gmgc = ""
        for n in GMGC_function:

                try:
                        gmgc = n["u"]
                except:
                        gmgc = ""

        return gmgc



def retrieve_neighbours_data(gmgc_cluster,neighbours_range):
        """ get a list containning -2-1+1+2 neoghbours genes surrounding 
        the unigenes in every contig """  
        
       
        # generate a neighbour range to analyze 2,4,6,8 or any other desired value 
        # i.e. with 2 generate this gene list [-2,-1,1,2]
        gene_list = []
        for num in range(-neighbours_range,neighbours_range+1,1):
                gene_list.append(num)
                
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

                for gene_pos in gene_list:
                        sample_cont = []
                        genes = int(gene)+int(gene_pos)
                        sample_cont = sample_contig[0:3]
                        sample_cont.append(str(genes))
                        genes = "_".join(sample_cont)
                        gene_ordered.append(genes)
                        neigh_dict[query_gene]=gene_ordered

        return neigh_dict 



def clean_unigene(gmgc):
        '''in case GMGC nomeclature was used in the input example,
         GMGC.100_000_123.UNKNOWN words after/before . are removed'''
        gmgc_clean = gmgc.split(".")[1]
        
        return gmgc_clean





def mongo_functional_find(GMGC, coll_clusters): #prueba al 350
        """ get functional information from mongo.clusters db for gmgc element"""
        kegg_list = []
        GMGC_function = coll_clusters.find({"u":GMGC})
        GMGC_function_list =[]
        for element in GMGC_function:

                kegg = element['K_P'] 
                kegg = kegg.split(",")

                for n in kegg: #depuro la lista de kegg para quedarme solo con los kegg pathways
                    if n in kegg_dict and n not in kegg_list:
                        kegg_list.append(n)

                #gene_symbol = element['g_s']
                #cog = element['COG']
                #desc = element['desc']
                Egg = element['OGs'].split(",")
                
                GMGC_function_list=[kegg_list, Egg]
            
        return GMGC_function_list





def compute_neigh_cogs_assignation(neighbours_genes,neighbours_range,coll_unigenes,coll_clusters,cog_category):
	   
        # VARIABLES
	cogs_organization_list = []  # storage cogs assgination for every unigene
	analysed_orfs = 0            # number of orfs (unigenes) that contains a GMGC cluster
	number_neigh = 0             # number of neighbours genes present in that GMGC cluster
	number_neigh_with_cogs= 0    # number of neighbours genes present in that gmgc cluster that contains cog assignation
	unigenes_functions_list = [] # storage all the cogs retrieved from neighbours genes	
	  

	for k,v in neighbours_genes.items():
                analysed_orfs += 1
                query_gene = k
                gene_list = v

                range_list = [] 
                for num in range(0,neighbours_range*2,1):
                        range_list.append(num)
                
                gmgc_list = [] # list for functional analysis developed below
                
                for n in range_list:
                        gmgc_orf = retrieve_gmgc(gene_list[n],coll_unigenes)
                        gmgc_list.append(gmgc_orf)
                
                #retrieve list of neigh orfs containing COGS instead of unigene code
                cogs_organization = []
                for unigene in gmgc_list: # unigene = cluster de GMGC al que pertenece el unigene

                        avoid_cog_repeated = [] # avoid repeated cogs in the same neigh orf
                        cog = ["NA"]
                        if unigene != '':
                                number_neigh +=1

                        try:
                                if cog_category == "kegg":
                                        cog = mongo_functional_find(unigene,coll_clusters)[0]
                                else:
                                        cog = mongo_functional_find(unigene,coll_clusters)[1]

                        except:
                                cog = ["NA"]

                        if cog == []:
                                cog=["NA"]

                        if cog != ["NA"] and cog[0] != '':
                                number_neigh_with_cogs +=1

                        if cog[0] =='':
                                cog= ["NA"]
                                cogs_organization.append(cog)
                        else:
                                cogs_organization.append(cog)
                        for n in cog: # clean the cogs to avoid future errors
                                if cog_category == "kegg":
                                        if n in kegg_dict and n not in avoid_cog_repeated and n != '':
                                                avoid_cog_repeated.append(n)
                                                unigenes_functions_list.append(n)
                                else:
                                        if n not in avoid_cog_repeated and n != '':
                                                if n != 'NA':
                                                        avoid_cog_repeated.append(n)
                                                        unigenes_functions_list.append(n)
                        
                        
                        #print cogs_organization #this print show the synteny organization of the neighbours cogs	
                        
                cogs_organization_list.append(cogs_organization)

	return 	[cogs_organization_list, 
		analysed_orfs, 
		number_neigh, 
		number_neigh_with_cogs, 
		unigenes_functions_list
		]


def neigh_scores(hits,num_query_cog,subject_cog_list,unique_cogs,count_of_cogs,neigh_orf_with_cogs):
	""" most important neigh score is functional 
        conservation developed for GMGC nature """

	#positve value (how many of preasignned keggs were retrieved keggs for that unigene)
	try:
		positive_value = ("{0:.2f}".format(int(hits)/num_query_cog)) 
	except:
		positive_value = 0

	#accuracy value (proportion of how many retrieved keggs were correct)
	try:
		accuracy_value = ("{0:.2f}".format(hits/len(subject_cog_list))) 
	except:
		accuracy_value = 0

	#operon functional conservation proxy
	try:
		func_conservation = ("{0:.3f}".format(1-(unique_cogs/count_of_cogs)/neigh_orf_with_cogs)) # neigh of with cogs its the same that the number of ORF that enter in the neigh analysis
	except:
		func_conservation = 0

	return [positive_value,
                accuracy_value,
                func_conservation
                ]

def get_neigh_orf_with_cogs(cogs_organizacion_list):
	""" calculate number of unigenes with at least 
	one neigh genes with cog assignation"""  

	cog_depured_list =[]
	for k in cogs_organizacion_list:	
		if k.count("NA") != 4:
			cog_depured_list.append(k)

	neigh_orf_with_cogs = len(cog_depured_list)
	
	return neigh_orf_with_cogs



def print_results(argument, kegg_result, egg_result, list_of_neighbours_unigenes,input_unigene_file):
	if argument == "-u":
		print("#gmgc"+"\t"+
			"query_keggs"+"\t"+
			"subject_cogs"+"\t"+
			"analysed_orfs"+"\t"+
			"number_neigh_genes"+"\t"+
			"number_neigh_with_cogs"+"\t"+
			"unique_cogs"+"\t"+
			"count_of_cogs"+"\t"+
			"cog_conservation"+"\t"+
			"hit_cog_percentage"+"\t"+
		"cog_description"
			)

		print("\t".join(kegg_result[0][0]))
		print("\t".join(egg_result[0][0]))
		print("\t".join(list_of_neighbours_unigenes[0]))

	if argument == "-f":
		
		directory = input_unigene_file+"_result"
		isdir = os.path.isdir(directory+"/") 
		if isdir == False:
                        os.mkdir(directory)

		kegg_output = open(directory+"/kegg_prediction.tsv","w")
		egg_output = open(directory+"/eggnog_prediction.tsv","w")
		unigene_output = open(directory+"/neighbours_unigenes.tsv","w")

		for line in kegg_result:
			kegg_output.write("\t".join(line[0])+"\n")

		for line in egg_result:
			egg_output.write("\t".join(line[0])+"\n")

		for line in list_of_neighbours_unigenes:
			unigene_output.write("\t".join(line)+"\n")





def get_functional_data (subject_cog_dict,query_list,kegg_description_dict,cog_category):
		
		hit_cog_percentage = [] # store percentage passed score of every COG i.e. COG@30.20
		hits = 0 # cogs present in the query unigene that were predicted by the method

		for k,v in subject_cog_dict.items():
			cog = k
			percentage = v
			k_p = k+"@"+str(percentage)
			hit_cog_percentage.append(k_p) 
			num_query_cog = len(query_list)

			if cog in query_list:
				hits +=1
						
		try:
			if query_list == [u'']:
				query_list =['NA']

			
			for cog in query_list:
				if cog_category == "kegg":	
					kegg_description_dict[cog]=get_kegg_description(cog)
				else:
					kegg_description_dict[cog]=get_egg_description(cog,coll_e5)
		except:
			query_list= ["NA"]

		try:
			description_list = []

			for cog, description in kegg_description_dict.items():
				k_d= str(cog)+"@"+str(description)
				description_list.append(k_d)

		except:
			description_list=[]
		
		return [hits, 
				query_list, 
				description_list, 
				hit_cog_percentage, 
				num_query_cog
				]



def create_unigenes_lists_for_print(gmgc,neighbours_genes,neighbours_range,coll_unigenes):
	""" create a list containing the unigenes which 
	belong every neighbours for every orf in the unigene cluster separated by @
	example: unigene_cluster \t unigene-2@unigene-1@unigene_cluster_orf@unigene+1@unigene+2"""
	
	unigenes_list_for_print = [] # list of neighbours unigenes to be save as ouput for future analysis
	analysed_orfs = 0
	for k,v in neighbours_genes.items():
			analysed_orfs += 1
			query_gene = k
			gene_list = v

			range_list = [] # [0,1,2,3,4]
			for num in range(0,(neighbours_range*2+1),1):
				range_list.append(num)

			gmgc_list = [] # list for functional analysis developed below
			

			for n in range_list:
				gmgc_orf = 	retrieve_gmgc(gene_list[n],coll_unigenes)
				gmgc_list.append(gmgc_orf)
			
			#create unigene_list_for_print, containing unigene assignation of every neigh orf
			out_list = []
			for n in gmgc_list:
				if n == '':
					n= 'NA'
				else:
					n=n
				out_list.append(n)
			
			out_line = "@".join(out_list)
			unigenes_list_for_print.append(out_line)

	line = gmgc+"\t"+",".join(unigenes_list_for_print)

	return line