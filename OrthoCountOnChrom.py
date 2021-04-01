import time
import sys,re
import argparse
from collections import defaultdict

def read_OG_dict(orthogroup_file:str) -> dict:
	"""
	:param orthogroup_file: Orthogroups.txt file
	:return: OG dict e.g: ["OG000001":[gene1, gene2],
						   "OG000002":[gene3, gene4]]
	"""
	og_dict = {}
	with open(orthogroup_file) as f:
		for i in f.readlines():
			og_dict[i.strip().split()[0]] = i.strip().split()[1:]
	return og_dict

def read_gene_list(gene_file:str) -> list:
	with open(gene_file) as target:
		target_geneID = list(target.read().splitlines())
	return target_geneID

def copyNumberOnChrom(geneID:list, ortho_dict:dict) -> dict:
	"""
	count accessory's gene orthogroups on core/accessory chromosomes

	:param geneID: Accessory chromosomes' gene id list
	:param ortho_dict: Dict from Orthogroups.txt file
	:return:
	"""
	strain_prefix = geneID[0][0:5]
	res = {}
	for gene in geneID:
		for og in ortho_dict:
			og_geneList = ortho_dict[og]
			if gene in og_geneList:
				sub_ogList = [i for i in og_geneList if i[0:5] == strain_prefix]
				geneOnCoreChrom = [i for i in sub_ogList if i not in geneID]
				geneOnAccessoryChrom = [i for i in sub_ogList if i in geneID]
				res[gene] = (len(geneOnCoreChrom), len(geneOnAccessoryChrom)-1)
				break
	return res

def main():
	orthogroup_file = "/Users/alexwang/0data/0mango/1_OrthoFinder/Orthogroups.txt"
	# accessory_gene = "/Users/alexwang/0data/0mango/accessory/gd10_accessory.geneID"
	# output = "/Users/alexwang/0data/0mango/accessory/gd10_accessoryGeneCopy.count.txt"
	accessory_gene = sys.argv[1]
	output = sys.argv[2]
	gene_list = read_gene_list(accessory_gene)
	og_list = read_OG_dict(orthogroup_file)
	res = copyNumberOnChrom(gene_list, og_list)
	with open(output,"w+") as f:
		f.write("geneID\tcore\taccessory\n")
		for gene in res:
			f.write(gene+"\t"+ str(res[gene][0]) + "\t" + str(res[gene][1]) + "\n")

if __name__ == '__main__':
	main()




