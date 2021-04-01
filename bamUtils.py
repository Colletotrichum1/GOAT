# dispose bam file
import pysam
import random
import multiprocessing
import shutil
import sys,os, re
import numpy as np
import argparse

def openBam(bamFile:str, returnStats:bool = False):
	try:
		bam = pysam.Samfile(bamFile, 'rb')

	except IOError:
		print("The file {} does not exist".format(bamFile))
	try:
		assert(bam.check_index() is not False)
	except:
		print("{} does not appear to have an index".format(bamFile))

	if bam.is_bam:
		mapped = bam.mapped
		unmapped = bam.unmapped
		if returnStats:
			stats = {chrom.contig: [chrom.mapped, chrom.unmapped] for chrom in bam.get_index_statistics()}
		if mapped == 0:
			sys.stderr.write("Warning:{} does not have mapped reads".format(bamFile))
	if returnStats:
		return bam, mapped, unmapped, stats
	else:
		return bam

def countReadsInInterval(args:tuple) -> tuple:
	chrom, start, end, fname, toEOF = args
	bam = openBam(fname)
	mapped = 0
	unmapped = 0
	for b in bam.fetch(chrom, start, end):
		if chrom == "*":
			unmapped += 1
			continue
		if b.pos < start:
			continue
		if not b.is_unmapped:
			mapped += 1
		else:
			unmapped += 1
	return mapped, unmapped, chrom

def getChrSize(bam: pysam.libcsamfile.Samfile):
	chrList = [(x,y) for x,y in zip(bam.references, bam.lengths)]
	return chrList

def getMapStats(bam: pysam.libcsamfile.Samfile, numThreads:int = 4)->tuple:
	header = getChrSize(bam)
	res = mapReduce([bam.filename, False], countReadsInInterval, header, numOfProcessors=numThreads)
	mapped = sum([x[0] for x in res])
	unmapped = sum([x[1] for x in res])
	stats = {x[0]: [0, 0] for x in header}
	for r in res:
		stats[r[2]][0] += r[0]
		stats[r[2]][1] += r[1]

	# We need to count the number of unmapped reads as well
	unmapped += bam.count("*")

	return mapped, unmapped, stats
# @func_timer
def mapReduce(staticArgs, func, chromSize:list, genomeChunkLegnth=None, numOfProcessors=4, self_ = None) -> list:
	"""
	:param staticArgs: tuple of the last two argument send to func
	:param func: function(chrom, start, end, staticArgs)
	:param chromSize: A list of tuple
	:param genomeChunkLegnth:
	:param numberOfProcessors:
	:return: list
	"""
	if not genomeChunkLegnth:
		genomeChunkLegnth = 1e5
	genomeChunkLegnth = int(genomeChunkLegnth)
	region_start = 0
	region_end = None
	Tasks = []
	for chrom, size in chromSize:
		start = 0 if region_start == 0 else region_start
		for startPos in range(start, size, genomeChunkLegnth):
			endPos = min(size, startPos+genomeChunkLegnth)
			regions = [[startPos, endPos]]

			for reg in regions:
				if self_ is not None:
					argsList = [self_]
				else:
					argsList = []
				argsList.extend([chrom, reg[0], reg[1]])
				argsList.extend(staticArgs)
			Tasks.append(tuple(argsList))
	if len(Tasks) > 1 and numOfProcessors > 1:
		random.shuffle(Tasks)
		pool = multiprocessing.Pool(numOfProcessors)
		res = pool.map_async(func, Tasks).get(9999999)
		pool.close()
		pool.join()    # wait for sub processing finished
	else:
		res = list(map(func, Tasks))
	return res

class ReadsCounter(object):
	def __init__(self, bamFile:str, binLength:int, numOfProcessors=4):
		self.bamFile = bamFile
		self.binLength = binLength
		self.numOfProcessors = numOfProcessors
	def getScaleFactor(self, method='rpkm') -> float:
		bam, mapped, unmapped, stats = openBam(self.bamFile, returnStats=True)
		# normalize bin coverage using RPKM
		if method == "rpkm":
			million_reads_mapped = float(mapped) / 1e6
			tile_len_in_kb = float(self.binLength) / 1000
			scale_factor = 1.0 / (million_reads_mapped * tile_len_in_kb)
			return scale_factor

	def countBinReads(self, chrom, start, end) -> np.array:
		subnum_reads_per_bin = []
		bam_handles = []
		bam_handles.append(openBam(self.bamFile))

		transcriptsToConsider = []
		transcriptsToConsider.append([(start, end, self.binLength)])
		for bam in bam_handles:
			for trans in transcriptsToConsider:
				tcov = self.getCoverage(bam, chrom, trans)
				subnum_reads_per_bin.extend(tcov)
		subnum_reads_per_bin = np.concatenate([subnum_reads_per_bin]).reshape(-1, 1, order='F')
		return subnum_reads_per_bin

	@staticmethod
	def getCoverage(bamHandle, chrom, regions:list) -> np.array:
		nbins = 0
		for reg in regions:
			nbins += (reg[1] - reg[0]) // reg[2]
			if (reg[1] - reg[0]) % reg[2] > 0:
				nbins += 1
		coverages = np.zeros(nbins, dtype='float64')

		vector_start = 0
		for idx, reg in enumerate(regions):
			if len(reg) == 3:
				tileSize = int(reg[2])
				nRegBins = (reg[1] - reg[0]) // tileSize
				if (reg[1] - reg[0]) % tileSize > 0:
					nRegBins += 1
			extension = 0
			regStart = int(max(0, reg[0] - extension))
			regEnd = reg[1] + int(extension)
			c = 0
			if chrom not in bamHandle.references:
				raise NameError("chromosome {} not found in bam file".format(chrom))

			for read in bamHandle.fetch(chrom, regStart, regEnd):
				position_blocks = read.get_blocks()
				last_eIdx = None
				for fragmentStart, fragmentEnd in position_blocks:
					if fragmentEnd is None or fragmentStart is None:
						continue
					fragmentLength = fragmentEnd - fragmentStart
					if fragmentLength == 0:
						continue
					# skip reads that are not in the region being
					# evaluated.
					if fragmentEnd <= reg[0] or fragmentStart >= reg[1]:
						continue
					if fragmentStart < reg[0]:
						fragmentStart = reg[0]
					if fragmentEnd > reg[0] + len(coverages) * tileSize:
						fragmentEnd = reg[0] + len(coverages) * tileSize

					sIdx = vector_start + max((fragmentStart - reg[0]) // tileSize, 0)
					eIdx = vector_start + min(np.ceil(float(fragmentEnd - reg[0]) / tileSize).astype('int'), nRegBins)
					if last_eIdx is not None:
						sIdx = max(last_eIdx, sIdx)
						if sIdx >= eIdx:
							continue
					sIdx = int(sIdx)
					eIdx = int(eIdx)
					coverages[sIdx:eIdx] += 1
					last_eIdx = eIdx
				c += 1
			vector_start += nRegBins
		return coverages

	@staticmethod
	def worker(self, chrom, start, end, func:str, func_args) -> tuple:
		"""
		:param chrom:
		:param start:
		:param end:
		:param func: function name to calculate bin RPKM. E.g. scaleCoverage
		:param func_args: argument for func, {'scale_factor': scale_factor}
		:return:
				A tuple of chunk information (chrom, start, end, tempfilename)
		"""
		if start > end:
			raise NameError("start position ({0}) bigger "
			                "than end position ({1})".format(start, end))
		coverage = self.countBinReads(chrom=chrom, start=start, end=end)

		line_string = "{}\t{}\t{}\t{:g}\n"
		_file = open(getTempFileName(suffix='.bg'), 'w')
		previous_value = None

		for tileIndex in range(coverage.shape[0]):
			tileCoverage = coverage[tileIndex, :]
			value = func(tileCoverage, func_args)
			if previous_value is None:
				writeStart = start + tileIndex * self.binLength
				writeEnd = min(writeStart+self.binLength, end)
				previous_value = value

			# if coverage value of two continuous bins are equal, the two bins can be merge
			elif previous_value == value:
				writeEnd = min(writeEnd+self.binLength, end)

			elif previous_value != value:
				if not np.isnan(previous_value):
					_file.write(line_string.format(chrom, writeStart, writeEnd, previous_value))
				previous_value = value
				writeStart = writeEnd
				writeEnd = min(writeStart + self.binLength, end)

		# write remaining value if not a nan
		if previous_value is not None and writeStart != end and not np.isnan(previous_value):
			_file.write(line_string.format(chrom, writeStart,
			                               end, previous_value))
		tempfilename = _file.name
		_file.close()
		return chrom, start, end, tempfilename

	def writeBedGraph(self, outFileName:str):
		scale_factor = self.getScaleFactor(method='rpkm')
		bam_handle = openBam(self.bamFile)
		chromsNameAndSize = getChrSize(bam_handle)
		genome_chunk_length = getGenomeChunkLength(bam_handle, self.binLength)
		res = mapReduce((scaleCoverage, {'scaleFactor':scale_factor}), wrapper, chromsNameAndSize, genomeChunkLegnth=genome_chunk_length, numOfProcessors=self.numOfProcessors, self_=self)

		# Determine the sorted order of the temp files
		chrom_order = dict()
		for i, _ in enumerate(chromsNameAndSize):
			chrom_order[_[0]] = i
		res = [[chrom_order[x[0]], x[1], x[2], x[3]] for x in res]
		res.sort()

		# if format == 'bedgraph':
		out_file = open(outFileName, 'wb')
		for r in res:
			if r[3]:
				_foo = open(r[3], 'rb')
				shutil.copyfileobj(_foo, out_file)
				_foo.close()
				os.remove(r[3])
		out_file.close()
		# else:
		# 	bedGraphToBigWig(chrom_names_and_size, [x[3] for x in res], out_file_name)

def wrapper(args):
	"""
	all arguments including self delivered to worker()
	"""
	return ReadsCounter.worker(*args)

def scaleCoverage(tile_coverage:list, args) -> int:
	"""
	tileCoverage should be an list with only one element
	"""
	return args['scaleFactor'] * tile_coverage[0]

def getGenomeChunkLength(bamHandle, tile_size) -> int:
	genomeLength = sum(bamHandle.lengths)
	max_reads_per_bp = float(bamHandle.mapped) / genomeLength
	# 2e6 is an empirical estimate
	genomeChunkLength = int(min(5e6, int(2e6 / max_reads_per_bp)))

	# The chunk length should be multiple fo tile_size
	genomeChunkLength -= genomeChunkLength % tile_size
	return genomeChunkLength

def getTempFileName(suffix=''):
	"""
    Return a temporary file name. The calling function is responsible for
    deleting this upon completion.
    """
	import tempfile
	_tempFile = tempfile.NamedTemporaryFile(prefix="_deeptools_",
                                            suffix=suffix,
                                            delete=False)

	memFileName = _tempFile.name
	_tempFile.close()
	return memFileName

def main(args):
	# bamfile="ph1_atac_sorted_uniq_specific.bam"
	# bam_handle = openBam(bamfile)
	# print(ReadsCounter.getCoverage(bam_handle, 'chr1', regions=[(1000, 2000, 50)]))
	# rc = ReadsCounter(bamFile=bamfile, binLength=50, numOfProcessors=4)
	# rc.writeBedGraph(outFileName="bam.bedgraph")
	if args.BAMFILE is None:
		parser.error("You should at least specify a BAM file.")
	else:
		bamfile = args.BAMFILE
		bs = int(args.bs)
		if args.REGIONS:
			chr, start, end = re.split("[-:]", args.REGIONS)
			bam_handle = openBam(bamfile)
			cov = str(ReadsCounter.getCoverage(bamHandle=bam_handle,
			                                   chrom=chr,
			                                   regions=[(int(start), int(end), bs)]))
			sys.stdout.write(cov)
		else:
			rc = ReadsCounter(bamFile=bamfile, binLength=bs, numOfProcessors=args.CPU)
			rc.writeBedGraph(outFileName=args.OUT)

if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument("-b", "--bam", help="input bam file", dest="BAMFILE")
	parser.add_argument("-r", "--regions", help="target regions", dest="REGIONS")
	parser.add_argument("-bs", "--binSize", help="bin size", dest="bs",default=50)
	parser.add_argument("-o", "--output", help="output file name", dest="OUT")
	parser.add_argument("-p", "--cpu", help="number of processors", dest="CPU", default=4)
	args = parser.parse_args()
	main(args)
