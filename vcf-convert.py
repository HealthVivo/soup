#!/usr/bin/env python
import sys, os, math, pysam, time

PHRED_THRESHOLD = -10 * math.log(0.05, 10)

# VCF class
class Vcf(object):
    def __init__(self, fasta, source):
        self.file_format = 'VCFv4.2'
        self.fasta = fasta
        self.source = source
        self.reference = fasta.filename
        self.info_list = []
        self.format_list = []
        self.alt_list = []
    
    class VcfElement(object):
    	def __init__(self, id, number, type, desc):
    		self.id = id
    		self.number = number
    		self.type = type
    		self.desc = desc
	
    def add_info(self, id, number, type, desc):
        self.info_list.append(Vcf.VcfElement(id, number, type, desc))

    def add_alt(self, id, desc):
        self.alt_list.append(Vcf.VcfElement(id, None, None, desc))

    def add_format(self, id, number, type, desc):
        self.format_list.append(Vcf.VcfElement(id, number, type, desc))
    
    # return the VCF header
    def get_header(self):
    	out = ['##fileformat=' + self.file_format, '##fileDate=' + time.strftime('%Y%m%d'), '##reference=' + self.reference]
    	out.append('##source='+self.source)
    	for x in self.info_list:
    		out.append('##INFO=<ID=' + x.id + ',Number=' + str(x.number) + ',Type=' + x.type + ',Description=\"' + x.desc + '\">')    
    	for x in self.alt_list:
    		out.append('##ALT=<ID=' + x.id + ',Description=\"' + x.desc + '\">')
    	for x in self.format_list:
    		out.append('##FORMAT=<ID=' + x.id + ',Number=' + str(x.number) + ',Type=' + x.type + ',Description=\"' + x.desc + '\">')
    	out.append('\t'.join(['#CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT']))
        return '\n'.join(out)
# end of VCF class


# SV class
class SV (object):
	def __init__(self, chrom, start, end, log2cn, phred):
		self.chrom = chrom
		self.start = start
		self.end = end
		self.log2cn = log2cn
		self.phred = phred
# end of SV class


# convert bicseq output ==> VCF
def write_vcf(vcf_fn, fasta_fn, sv_data, source):
	# fasta
	fasta = pysam.Fastafile(fasta_fn)
	
	# vcf
	myvcf = Vcf(fasta, source)
	myvcf.add_info('SVTYPE', 1, 'String', 'Type of structural variant')
	myvcf.add_info('SVLEN', '.', 'Integer', 'Difference in length between REF and ALT alleles')
	myvcf.add_info('END', 1, 'Integer', 'End position of the variant described in this record')
	myvcf.add_info('IMPRECISE', 0, 'Flag', 'Imprecise structural variation')
	myvcf.add_alt('DEL', 'Deletion')
	myvcf.add_alt('DUP', 'Duplication')
	myvcf.add_format('GT', 1, 'String', 'Genotype')
	myvcf.add_format('GQ', 1, 'String', 'Genotype quality')
	myvcf.add_format('CNQ', 1, 'Float', 'Copy number genotype quality for imprecise events')
	myvcf.add_format('CN', 1, 'Float', 'Copy number genotype for imprecise events')
	
	# create output file
	out = open(vcf_fn, 'w')
	out.write(myvcf.get_header()+'\n')
	
	# parse each line and convert to VCF format
	for x in sv_data:
		# check p-value ratio against threshold
		if x.phred < PHRED_THRESHOLD: continue
		
		# is this a duplication or a deletion?
		sv_type = ""
		sv_len = 0
		if x.log2cn > 0.0:
			sv_type = "DUP"
			sv_len = x.end - x.start
		elif x.log2cn < 0.0:
			sv_type = "DEL"
			sv_len = x.start - x.end
		
		# get FASTA ref
		ref = fasta.fetch(x.chrom, x.start-1, x.start)
		if len(ref) == 0: ref = "."
		
		# INFO
		info = "SVTYPE=%s;END=%d;SVLEN=%d;IMPRECISE" % (sv_type, x.end, sv_len)
		
		# FORMAT
		format = "GT:GQ:CN:CNQ\t./.:.:%.02f:%.02f" % (2**x.log2cn, x.phred)
				
		# output line
		# format: CHROM POS ID REF ALT QUAL FILTER INFO FORMAT
		out_line = "%s\t%d\t.\t%s\t<%s>\t%.02f\t.\t%s\t%s\n" % (x.chrom, x.start, ref, sv_type, x.phred, info, format)
		out.write(out_line)
	# end of for loop
	
	out.close()
	fasta.close()
# end of conversion


# read and parse bicseq file
def bicseq2vcf(bic_fn, fasta_fn):
	fdata = []
	try:
		with open(bic_fn, 'r') as f: fdata = f.readlines()
	except:
		print "Could not open / read the input file"
	
	# parse line
	svs = []
	for line in fdata:
		pieces = line.split()
		if len(pieces) != 7:
			continue
		chrom, start, end, log2cn, phred = "0", 0, 0, 0.0, 0.0
		try:
			chrom = pieces[0]
			start = int(pieces[1])
			end = int(pieces[2])
			log2cn = float(pieces[5])
			log10p = float(pieces[6])
			phred = -10*log10p if log10p != 0 else log10p
		except Exception, ex:
			continue
		svs.append(SV(chrom, start, end, log2cn, phred))
	# end of parse line
	
	# now create VCF output
	write_vcf(bic_fn+".vcf", fasta_fn, svs, 'bicseq')
# end of bicseq output processing


# read and parse bicseq file
def cnvnator2vcf(cnv_fn, fasta_fn):
	fdata = []
	try:
		with open(cnv_fn, 'r') as f: fdata = f.readlines()
	except:
		print "Could not open / read the input file"
	
	# parse line
	svs = []
	for line in fdata:
		pieces = line.split()
		if len(pieces) != 10:
			continue
		chrom, start, end, log2cn, phred = "0", 0, 0, 0.0, 0.0
		try:
			if float(pieces[3]) == 0: continue
			idx1 = pieces[1].find(':')
			idx2 = pieces[1].find('-')
			if idx1 <= 0 or idx2 <= 0: continue
			chrom = pieces[1][:idx1]
			if chrom.startswith('chr'): chrom = chrom[3:]
			start = int(pieces[1][idx1+1:idx2])
			end = int(pieces[1][idx2+1:])
			log2cn = math.log(float(pieces[3]), 2)
			pval = float(pieces[9])
			if pval == 0:
				phred = 1000 # some nominal big number
			else:
				phred = -10*math.log(pval, 10) if pval != 1 else pval
		except Exception, ex:
			continue
		svs.append(SV(chrom, start, end, log2cn, phred))
	# end of parse line
	
	# now create VCF output
	write_vcf(cnv_fn+".vcf", fasta_fn, svs, 'cnvnator')
# end of bicseq output processing	
	


# start of main
if __name__ == "__main__":
	def print_usage():
		print "Usage: %s {-b | -c} {input file} {path to ref genome}" % sys.argv[0]
		print "  -b: bicseq input file"
		print "  -c: cnvnator input file"
	
	if len(sys.argv) != 4:
		print_usage()
		sys.exit(1)
	
	if sys.argv[1] == '-b':
		bicseq2vcf(sys.argv[2], sys.argv[3])
	elif sys.argv[1] == '-c':
		cnvnator2vcf(sys.argv[2], sys.argv[3])
	else:
		print "Unknown option: %s" % sys.argv[1]
		print_usage()
		sys.exit(1)
# end of main