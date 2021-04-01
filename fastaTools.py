#!/usr/bin/python3

from Bio import SeqIO,Seq
import sys,os,argparse
import re


def openFasta(faFile:str) -> tuple:
	if not os.path.exists(faFile):
		raise IOError('Fasta file does not exists: %s' % faFile)
	fa = SeqIO.parse(faFile,format="fasta")

	seq_dict = {}
	for seq in fa:
		seq_dict[seq.name] = seq.seq
	return faFile,seq_dict

class Fasta(object):
	def __init__(self,faInfo):
		self.fn, self.seq_dict = faInfo
		self.amount = len(self.seq_dict.keys())
		self.totalBase = sum(map(len, self.seq_dict.values()))

	def subseq(self, id:str, header=True, range:str=False):
		if not id in self.seq_dict.keys():
			raise KeyError("Error! sequence %s not in fasta file" % id)
		seq = self.seq_dict[id]
		if range:
			s,e = map(int,range.split("-"))
			seq = seq[s-1:e]
		if header:
			sys.stdout.write(">{0}\n{1}".format(id, str(seq)))
		else:
			sys.stdout.write(str(seq))
		def __str__(self):
			return "Fasta file:{fn},contains {cc} sequence.".format(fn=self.fn, cc=self.amount)
		__repr__ = __str__

def rearrange(seq_dict:dict) -> dict:
		sort_seq = sorted(seq_dict.values(),key = lambda v: len(v), reverse=True)
		sort_id = seq_dict.keys()
		sorted_seq_dict = dict(zip(sort_id, sort_seq))
		return sorted_seq_dict

def rmseq(seq_dict: dict, id: str, sort = True, keep = False) -> dict:
	# keeped_seq_dict = dict(filter(lambda item: item[0] != id, seq_dict))
	id_list = id.split(",")
	if keep:
		keeped_seq_dict = {k: v for k, v in seq_dict.items() if k in id_list}
	else:
		keeped_seq_dict = {k: v for k, v in seq_dict.items() if k not in id_list}

	if sort:
		return rearrange(keeped_seq_dict)
	else:
		return keeped_seq_dict

def filterseqBylength(seq_dict:dict, minLen=100) -> dict:
	keeped_seq_dict = {k: v for k, v in seq_dict.items() if len(v) >= minLen}
	return keeped_seq_dict


# def writeFasta(seq_dict:dict, outfile:str):
# 	with open(outfile,"w") as fh:
# 		for k,v in seq_dict.items():
# 			fh.writelines(">"+k+"\n"+str(v)+"\n")

def output(seq_dict:dict):
	for k,v in seq_dict.items():
		print(">"+k+"\n"+str(v)+"\n")

def main(args):
	if args.fa is None:
		parser.error("You should at least specify a fasta file.")
	else:
		faInfo = openFasta(args.fa)
		fa = Fasta(faInfo)
		fn, seq_dict = faInfo
		if not args.sort and not args.rm and not args.len and not args.kp:
			print("Total base:%d" % fa.totalBase)
			for id, s in seq_dict.items():
				if args.gc:
					print(id, end="\t")
					print(len(s), end="\t")
					tot = str(s).count("C") + str(s).count("c") + str(s).count("G") + str(s).count("g")
					print(tot/len(s))
				else:
					print(id, end="\t")
					print(len(s))
		else:
			if args.sort:
				dict = rearrange(seq_dict)
				output(dict)

			elif args.rm:
				dict = rmseq(seq_dict, id = args.rm, keep=False)
				output(dict)

			elif args.kp:
				dict = rmseq(seq_dict, id = args.kp, keep=True)
				output(dict)

			elif args.len:
				dict = filterseqBylength(seq_dict, minLen=args.len)
				output(dict)



if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument("-i", "--input", help="input fasta file", dest="fa")
	parser.add_argument("-gc", "--GC_content", help="output GC content", dest="gc",action="store_true")
	parser.add_argument("-s", "--sort", help="sort input by sequence length.", dest="sort", action="store_true")
	parser.add_argument("-rm", "--remove", help="remove target contigs.", dest="rm")
	parser.add_argument("-keep", "--keepID", help="only keep target contigs.", dest="kp")

	parser.add_argument("-f", "--filterBylength", help="filter sequence < 100bp.", dest="len")
	args = parser.parse_args()
	main(args)
