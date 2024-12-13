#GC content - output tsv with file name and GC content

import argparse
import re
import glob
import os
from Bio.SeqUtils import GC
from Bio import SeqIO
import pandas as pd

#################################

def get_argue():
	parser = argparse.ArgumentParser(description="Get GC content from fasta files in directory")
	
	fasta_dir = parser.add_argument("-f","--fasta_dir",
		help="The path to the directory containing the fasta(s) you wish to check",
		required=True)

	output_file = parser.add_argument("-o","--output_file",
		help="Output file name - '.tsv' will be added",
		required=True)

	args = parser.parse_args()

	fasta_dir = args.fasta_dir
	output_file = args.output_file
	
	return fasta_dir, output_file


#################################

def get_dir_list(path):
	exist_test = os.path.exists(path)
	if(exist_test == True):
		list_dir = os.listdir(path)
	else:
		print("### ERROR - path does not exist")

	return list_dir

#################################

def calc_gc(fasta):

	total_seq_len = 0
	total_N_count = 0
	total_GC_count = 0
	
	contig_larger_500bp = 0
	
	total_N_count_500 = 0
	total_GC_count_500 = 0
	total_seq_len_500 = 0


	for seq_record in SeqIO.parse(fasta, 'fasta'):
		seq_identifier = seq_record.id
		s = seq_record.seq
		seq_len = len(s)
		total_seq_len += seq_len
		if(seq_len > 500):
			contig_larger_500bp += 1

			total_seq_len_500 += seq_len

			for ntd in s:
				if(ntd == 'N'):
					total_N_count_500 += 1
				if(ntd == 'G' or ntd == 'C'):
					total_GC_count_500 += 1

		for ntd in s:
			if(ntd == 'N'):
				total_N_count += 1
			if(ntd == 'G' or ntd == 'C'):
				total_GC_count += 1

	total_GC_content = round((total_GC_count/total_seq_len)*100,1)

	total_GC_content_500 = round((total_GC_count_500/total_seq_len_500)*100,1)

	results_list = [fasta,total_seq_len, total_GC_count, total_GC_content, total_N_count, contig_larger_500bp, total_seq_len_500, total_GC_count_500, total_N_count_500, total_GC_content_500]

	return(results_list)


#################################

def main():
	fasta_dir, output_file = get_argue()
	print("Fasta path entered: {}\noutput_file: {}\n".format(fasta_dir,output_file))

	path_contents = get_dir_list(fasta_dir)

	lol_results = []

	for i in (path_contents):
		ass_stats = calc_gc(fasta_dir+'/'+i)
		print("{}\t{}".format(i,ass_stats))
		lol_results.append(ass_stats)


	results_df = pd.DataFrame(lol_results, columns=['isolate_name','seq_len','GC_count','GC_content','N_count','contig_larger_500bp','total_seq_len_500','total_GC_count_500','total_N_count_500','total_GC_content_500'])
	results_df.to_csv(output_file+".tsv",sep="\t", index=False)

#################################

if __name__ == '__main__':
	main()
