import sys, os
import re
from functools import wraps
from collections import defaultdict
from sys import stderr
import time
import argparse

def func_timer(func):
	"""
	a decorator
	:param func: target function
	:return:
	"""
	@wraps(func)
	def timeit(*args, **kwargs):
		print("Funciton: {name} start...".format(name = func.__name__))
		t0 = time.time()
		result = func(*args, **kwargs)
		t1 = time.time()
		print('[Function: {name} finished, spent time: {time:.2f}s]'.format(name=func.__name__, time = t1 - t0))
		return result
	return timeit

class Interval(object):
	"""docstring for Interval"""

	def __init__(self, start, end, strand):
		if start > end:
			raise ValueError('Error! start is smaller than end!')
		super(Interval, self).__init__()
		self.start = start
		self.end = end
		self.strand = strand
		self.length = end - start

	def is_contain(self, start, end):
		if start > end:
			raise ValueError('Error! start is smaller than end!')
		if self.start <= start and self.end >= end:
			return True
		else:
			return False

	def is_overlap(self, start, end):
		if self.start <= start <= self.end or self.start <= end <= self.end:
			return True
		else:
			return False

	def __str__(self):
		return "interval:%i-%i" % (self.start, self.end)

	__repr__ = __str__

class ParsingError(Exception):
	pass

class Gene(object):
	""" Gene object """

	def __init__(self, field_dict):
		self.chrom = field_dict["chrom"]
		self.genomic_iv = field_dict["iv"]
		self.attr = field_dict["attr"]
		self.gene_id = self.attr["gene_id"]
		self.gene_name = self.attr.get("gene_name", self.attr.get("gene_id", None))
		self.transcripts = []
		self.principal_transcripts = []

	def __str__(self):
		return "GeneID:%s;GeneName:%s" % (self.gene_id, self.gene_name)

	__repr__ = __str__

class Transcript(object):
	""" Transcript object """

	def __init__(self, field_dict):
		self.chrom = field_dict["chrom"]
		self.genomic_iv = field_dict["iv"]
		self.attr = field_dict["attr"]
		self.gene_id = self.attr["gene_id"]
		self.gene_name = self.attr.get("gene_name", self.attr.get("gene_id", None))
		self.transcript_id = self.attr["transcript_id"]
		self.genomic_exons = []
		self.genomic_cds = []
		self.genomic_startcodon = []
		self.genomic_stopcodon = []
		self.length = 0

		self.cds = None
		self.startcodon = None
		self.stopcodon = None


	def add_feature(self, field_dict):
		feature = field_dict["feature"]
		iv = field_dict["iv"]
		if feature == "exon":
			self.genomic_exons.append(iv)
			self.length += iv.length
		elif feature == "CDS":
			self.genomic_cds.append(iv)
		elif feature == "start_codon":
			self.genomic_startcodon.append(iv)
		elif feature == "stop_codon":
			self.genomic_stopcodon.append(iv)
		else:
			raise ValueError("Error: the feature is not recognized: %s" % feature)

	def cdsLength(self):
		if self.genomic_cds:
			cds_list = self.genomic_cds
			l = sum([iv.length for iv in cds_list])
			return l
		else:
			return 0

	def __str__(self):
		return "Transcript obj: transcriptID:%s" % (self.transcript_id)

	__repr__ = __str__

def parsingGFFattr(attr_string):
	attr_dict = {}
	attr_split = re.findall(r'.*? ".*?";', attr_string)
	# print(attr_split)
	for i in attr_split:
		if i:
			k, v = i.strip().split(" ", 1)
			if v.endswith(";"):
				v = v[:-1]
			v = v.strip('"')
			attr_dict[k] = v
	return attr_dict

def parsingline(line: str)->dict:
	"""
    GTF/GFF anotations
    """
	chrom, source, feature, start, end, score, strand, frame, attr = line.strip().split("\t", 8)
	# convert to 0-based.
	iv = Interval(int(start) - 1, int(end), strand)
	field_dict = {"chrom": chrom, "source": source, "feature": feature, "iv": iv, "strand": strand, "attr": parsingGFFattr(attr)}
	return field_dict

def readGFF(filename: str, reformat = False) -> tuple:
	if not os.path.exists(filename):
		raise IOError('GFF file does not exists: %s' % filename)
	sys.stdout.write("Reading the GFF file: %s ......\n" % filename)
	gene_dict = {}
	transcript_dict = {}
	d_g = defaultdict(dict)
	with open(filename) as fn:
		for i,line in enumerate(fn):
			if line[0] == '#' or (not line.strip()):
				continue
			if line.strip().split("\t")[2] == "gene" or line.strip().split("\t")[2] == "transcript":
				line = cor_attr(line)
			field_dict = parsingline(line)
			if field_dict['feature'] == 'gene':
				gobj = Gene(field_dict)
				gene_dict[gobj.gene_id] = gobj
			elif field_dict['feature'] == 'transcript':
				tobj = Transcript(field_dict)
				transcript_dict[tobj.transcript_id] = tobj
			elif field_dict["feature"] in ["exon", "intron", "CDS", "start_codon","stop_condon"]:
				if reformat:
					g_id = field_dict['attr']['gene_id']
					t_id = field_dict['attr']['transcript_id']
					iv = field_dict['iv']
					chr = field_dict['chrom']
					source = field_dict['source']
					strand = field_dict['strand']
					if t_id not in d_g[g_id].keys():
						d_g[g_id][t_id] = [(chr, source, strand, iv, line)] ### {"g1":{"g1.t1":[()]}}
					else:
						d_g[g_id][t_id].append((chr, source, strand, iv, line))
				else:
					tid = field_dict["attr"]["transcript_id"]
					try:

						transcript_dict[tid].add_feature(field_dict)
					except KeyError:
						raise ParsingError("Error in line %i. The annotation should be three-level hierarchy of \
						                    gene => transcript => exon (or CDS)" % i)
			else:
				pass
		if reformat:
			for gene_id in d_g.keys():
				chr, source, strand = list(d_g[gene_id].values())[0][0][0:3]
				txs = list(d_g[gene_id].values())
				s = min([j[3].start for i in txs for j in i]) + 1
				e = max([j[3].end for i in txs for j in i])
				sys.stdout.write(
					"{0}\t{1}\tgene\t{2}\t{3}\t.\t{4}\t.\tgene_id \"{5}\";\n".format(chr, source, s, e, strand, gene_id))
				for tx_id in d_g[gene_id].keys():
					txs2 = d_g[gene_id][tx_id]
					ss = min([j[3].start for j in txs2]) + 1
					ee = max([j[3].end for j in txs2])
					sys.stdout.write(
						"{0}\t{1}\ttranscript\t{2}\t{3}\t.\t{4}\t.\tgene_id \"{5}\"; transcript_id \"{6}\";\n".format(chr,
						                                                                                            source,
						                                                                                            ss, ee,
						                                                                                            strand,
						                                                                                            gene_id,
						                                                                                            tx_id))
					for fe in txs2:
						sys.stdout.write(fe[4])
		else:
			return gene_dict, transcript_dict


def cor_attr(line: str) -> str:
	"""
	some GTF files' format are not standard, so correct them
	:param line: str
	:return: corrected lines
	"""
	ll = line.split("\t")
	if ll[2] == "gene":
		if not re.findall(r"gene_id \".*\";", ll[8]): ## check whether "gene" records are standard
			return "\t".join(ll[0:8]) + "\t" + "gene_id \"{0}\";\n".format(ll[8].strip())
		else:
			return line
	if ll[2] == "transcript":
		if not line.strip().endswith(";"):
			return line.strip() + ";\n"
		else:
			return line

def GroupGeneSubsets(gtf_file: str) -> tuple:
	"""
	Group gtf lines into subset that are the same gene
	Group is based on the ID attribute:
	:param gtf_data:
	:return: subsets: dict
			 gset_keys: keys of subsets
	"""
	subsets = defaultdict(list)
	gset_keys = []
	gset_keys_set = set()
	with open(gtf_file) as fin:
		for line in fin:
			if line[0] =="#" or (not line.strip()):
				continue
			if line.strip().split("\t")[2] == "gene" or line.strip().split("\t")[2] == "transcript":
				line = cor_attr(line)

			field_dict = parsingline(line)
			field_dict["line"] = line
			gid = field_dict["attr"]["gene_id"]
			if "gene_name" not in field_dict["attr"]:
				field_dict["attr"]["gene_name"] = gid
			## gid: 'Scaffold_1:SCLS00001:-'
			# chrom = field_dict["chrom"]
			# strand = field_dict["iv"].strand
			# gid0 = chrom + ":"+gid+":"+strand
			## gid: "SCLS00001"
			gid0 = gid
			subsets[gid0].append(field_dict)
			if gid0 not in gset_keys_set:
				gset_keys_set.add(gid0)
				gset_keys.append(gid0)
	#### select duplicated gene id
	# gene_id_list = [i.split(":")[1] for i in subsets.keys()]
	# gene_id_uniq = set(gene_id_list)
	# duplicated_gene_id = [i for i in gene_id_uniq if gene_id_list.count(i) >1]
	# duplicated_gene_id = set(duplicated_gene_id)
	# for i in subsets.keys():
	# 	chrom,gid,strand = i.split(":")
	# 	if gid in duplicated_gene_id:
	# 		stderr.write("Warnning, gene_id is duplicated: %s , will be renamed as: %s\n" % (gid, i))
	# 		for j in range(len(subsets[i])):
	# 			subsets[i][j]["attr"]["gene_id"] = i
	return subsets, gset_keys

def remove_record(transDict:dict, gset:dict, length = 150) -> dict:
	print("Removing transcripts whose lengths are smaller than {} bp ...".format(length))
	for tx in transDict:
		if transDict[tx].cdsLength() < length:
			gene_id, tx_id = transDict[tx].gene_id, transDict[tx].transcript_id
			# for subFeatureDict in gset[gene_id]:
			# 	attr_dict = subFeatureDict['attr']
			keep = []
			for rec in gset[gene_id]:
				try:
					if not rec['attr']['transcript_id'] == tx_id:
						keep.append(rec)
					else:
						sys.stdout.write("Removing transcript {0}, CDS length:{1}\n".format(tx_id, transDict[tx].cdsLength()))
				except KeyError:
					pass
			if len(keep) == 1 and keep[0]['feature'] == 'gene':
				del gset[gene_id]
			else:
				gset[gene_id] = keep
	return gset

def getPrimaryTranscript(transDict: dict, gset: dict) -> list:
	"""
	order transcript by CDS length
	:param transDict: Dict
	:param gset: Dict
	:return: List
	"""
	print("Extracting primary transcripts ... ")
	black_list = []
	for geneGroup in gset.values():
		transGroup = []
		for i in geneGroup:
			if i["feature"] == "transcript":
				tid = i["attr"]["transcript_id"]
				transGroup.append((tid, transDict[tid].cdsLength()))
		if len(transGroup) >= 2:
			sortByCds = sorted(transGroup, key=lambda x:x[1], reverse=True)
			for id,l in sortByCds[1:]:
				black_list.append(id)
	return black_list

def write_GTF(gset, outFileName, black_list = None):
	# with open("{0}_primaryTrans.gtf".format(gtfFile.rstrip(".gtf")), "w+") as outf:
	with open(outFileName, "w+") as outf:
		for s in gset.values():
			for f in s:
				if not black_list:
					outf.writelines(f["line"])
				else:
					if not f["feature"] == "gene":
						if not f["attr"]["transcript_id"] in black_list:
							outf.writelines(f["line"])
					else:
						outf.writelines(f["line"])

def old2new(gene_dict:dict) -> dict:
	"""
	Sort the gene dict and return an old2new mapping dict
	:param gene_dict:
	"""
	r = []
	for obj in gene_dict.keys():
		chr_num = int(re.findall(r".*_(\d+)", gene_dict[obj].chrom)[0])
		# print(chr_num)
		r.append((obj, chr_num, gene_dict[obj].genomic_iv))
	r.sort(key=lambda x: (x[1], x[2].start))
	mapping_dict = dict(zip([i[0] for i in r], range(1, len(r) + 1)))
	return mapping_dict

def renameGeneIdentifier(gset:dict, geneDict:dict, prefix:str) -> tuple:
	"""
	Used to rename gene identifier result from Augustus gene model prediction software
	:param gset:
	:param geneDict:
	:param prefix:
	:return:
	"""
	print("Rename gene identifier to {} ...".format(prefix))
	d = old2new(geneDict)
	new_dict = {}
	for geneID in gset:
		new_identifier = prefix + str(d[geneID]).zfill(5)
		for dd in gset[geneID]:
			for i in dd['attr']:  # replace old gene name with new identifier by mapping Dict
				dd['attr'][i] = dd['attr'][i].replace(geneID, new_identifier)
			dd['line'] = dd['line'].replace(geneID, new_identifier)

		new_dict[new_identifier] = gset[geneID]

	return new_dict, new_dict.keys()


def main(args):
	if args.GTFFILE is None:
		parser.error("You should at least specify a GTF file.")
	else:
		gtfFile = args.GTFFILE
		if args.REFORMAT:
			readGFF(gtfFile, reformat = True)
		else:
			geneDict, transDict = readGFF(gtfFile, reformat=False)  # return geneDict/transDict: dict of all Gene Object and Transcript Object
			gset, gset_keys = GroupGeneSubsets(gtfFile)  # return gset: dict of all  features
			if args.RENAME:
				PREFIX = args.RENAME
				gset2, gset_keys2 = renameGeneIdentifier(gset, geneDict, prefix=PREFIX)
				write_GTF(gset2, args.OUT)

			elif args.CDS_LEN_CUTOFF:
				gset2 = remove_record(transDict, gset, length=int(args.CDS_LEN_CUTOFF))
				write_GTF(gset2, args.OUT)

			elif args.KEEP_PRIMARY:
				black_list = getPrimaryTranscript(transDict, gset)
				write_GTF(gset, args.OUT, black_list=black_list)
			############### testing
			# gtfFile = "xx.gtf"
			# geneDict, transDict = readGFF(gtfFile)  # return geneDict/transDict: dict of all Gene Object and Transcript Object
			# gset, gset_keys = GroupGeneSubsets(gtfFile)  # return gset: dict of all  features
			# gset2 = remove_record(transDict,gset,150)
			# print(gset2['SCLS08366'])

if __name__ == '__main__':
	# parent_parser = argparse.ArgumentParser(add_help=False)
	# parent_parser.add_argument("-i", "--input", help="input GTF file handle", dest="GTFFILE")
	# parent_parser.add_argument("-o", "--output", help="output GTF file name", dest = "OUT")
	#
	# parser = argparse.ArgumentParser()
	# subparsers = parser.add_subparsers(help='sub-command help')
	# # add sub command "operate"
	# parser_b = subparsers.add_parser('operate', help='operate help',parents = [parent_parser])
	# parser_b.add_argument("-r", "--rename", help="rename gene identifier and order by coordinate", dest="RENAME")
	# # add sub command "filter"
	# parser_a = subparsers.add_parser('filter', help='filter help',parents = [parent_parser])
	# parser_a.add_argument("-len", "--length", help="filter transcript with CDS length < LEN", dest="CDS_LEN_CUTOFF")
	# parser_a.add_argument("-p", "--primary_transcript", help="only keep primary transcript", dest="KEEP_PRIMARY", action="store_true")


	parser = argparse.ArgumentParser()
	parser.add_argument("-i", "--input", help="input GTF file handle", dest="GTFFILE")
	parser.add_argument("-o", "--output", help="output GTF file name", dest="OUT")
	parser.add_argument("-r", "--rename", help="rename gene identifier and order by coordinate", dest="RENAME")
	parser.add_argument("-len", "--length", help="filter transcript with CDS length < LEN", dest="CDS_LEN_CUTOFF")
	parser.add_argument("-p", "--primary_transcript", help="only keep primary transcript", dest="KEEP_PRIMARY",
	                      action="store_true")
	parser.add_argument("-reformat", "--reformat", help="reformat gtf file", dest="REFORMAT",action="store_true")
	args = parser.parse_args()
	main(args)
	# main()

