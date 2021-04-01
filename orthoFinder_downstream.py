import time
import sys,re
import argparse
from collections import defaultdict

def read_gene_list(gene_file:str) -> list:
	with open(gene_file) as target:
		target_geneID = list(target.read().splitlines())
	return target_geneID

def copyNumberOnCoreChrom(geneID:list, ortho_dict:dict) -> dict:
	"""
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
				res[gene] = len(geneOnCoreChrom)
				break
	return res


def read_OG_dict(orthogrou_file:str) -> dict:
	"""
	:param orthogrou_file: Orthogroups.txt file
	:return: OG dict e.g: ["OG000001":[gene1, gene2],
						   "OG000002":[gene3, gene4]]
	"""
	og_dict = {}
	with open("/Users/alexwang/0data/0mango/1_OrthoFinder/Orthogroups.txt") as f:
		for i in f.readlines():
			og_dict[i.strip().split()[0]] = i.strip().split()[1:]
	return og_dict


def main(args):
	if not args.go:
		N = 28
		sub_OG = []
		# with open("/Users/alexwang/0data/0mango/transcriptome/ortholog_heat/cazy.geneID") as target:
		target_geneID = read_gene_list(args.gene)
		outf1 = open(args.pre + "_count.tsv", "w+")
		outf2 = open(args.pre + "_targetGene.tsv", "w+")
	else:
		## Add GO annotation
		# with open("/Users/alexwang/0data/0mango/protein/interproscan/go_annot.txt") as f:
		with open(args.go) as g:
			go_list = g.readlines()
		outf = open("orthogroups_GO.tsv", "w+")
		og_go = defaultdict(list)

	# with open("/Users/alexwang/0data/0mango/transcriptome/ortholog_heat/Orthogroups.tsv") as f:
	with open(args.ortho) as f:
		for idx, line in enumerate(f):
			if not args.go:
				og = line.strip().split("\t")
				cc = [og[0]]
				gg = [og[0]]
				for strain_og in og[1:]:
					sub_gg = []
					target_count = 0
					gene_list = strain_og.strip().split(", ")
					for ge in gene_list:
						if ge in target_geneID:
							sub_gg.append(ge)
							target_count += 1
					cc.append(target_count)
					if sub_gg:
						gg.append(",".join(sub_gg))
				if sum(cc[1:]) != 0:
					sub_OG.append(line)
					if len(cc)!=N+1:
						cc = cc + [0]*(N+1-len(cc))
					cc = list(map(str, cc))
					outf1.writelines("\t".join(cc)+"\n")
					outf2.writelines("\t".join(gg)+"\n")
			else:
				og = re.split("[\t ,]", line.strip())
				for k in go_list:
					gene, goid, definition, ontology, term = k.strip().split("\t")
					if gene in og[1:]:
						if not term in og_go[og[0]]:
							og_go[og[0]].append(term)

	if not args.go:
		with open(args.pre + "_subOG.tsv", "w+") as outf3:
			outf3.writelines(sub_OG)
		outf2.close()
		outf1.close()
	else:
		for key in og_go:
			outf.writelines(key+"\t"+"|".join(og_go[key]) + "\n")

if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument("-og", "--orthogroup", help="input orthogroup.txt from orthoFinder2", dest="ortho")
	parser.add_argument("-target", "--target", help="input target gene id", dest="gene")
	parser.add_argument("-pre", "--prefix", help="input prefix of target genes", dest="pre")
	parser.add_argument("-go", "--GOannotation", help="add GO term to orthogroups", dest="go")
	# parser.add_argument("-s", "--sub", help="output sub orthogroup of target genes", dest="sub", action="store_true")
	args = parser.parse_args()
	main(args)
