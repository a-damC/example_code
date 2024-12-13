#Script to locate SNP loci that contain an N from an alignment

import os
import csv
import re
import sys
import argparse
import datetime
from Bio import SeqIO


##########################################################################################################################################################################

def get_argue():
	parser = argparse.ArgumentParser(description="Script to identify isolates with an N at a loci where a SNP occurs")
	a = parser.add_argument("-a","--alignment",
		help="Path to alignment file, including the name of the file.",
		required=True)
	p =parser.add_argument("-p","--percentage", 
		help="Threshold value at which you add isolates with N at a SNP loci into your 'Offenders' bin.\tFor example:\nIf you have an alignment of 100 isolates and choose 10 as the percentage, this will only add isolates to your 'Offenders' bin if\nat a SNP loci 10 or less isolates have a N\noption should change depending on what you are looking for, in the example above if you were looking for a specific isolate you would\nwant a low -a maybe 5 or less, if you are looking for a group of isolates maybe 15 or less but again this depends on the size of your alignment",
		required=True,
		type=int)
	r =parser.add_argument("-r","--refgenome_name",
		help="Name of outlier or reference used in alignment (optional arguement)",
		default="XXXXXX")
	o =parser.add_argument("-o","--output_results",
		help="Do you wish to output the results to a file? N for No (default) and if you wish to have a results file add the name to the option, results output will be named *date*_XXXX_snpable_results.tsv",
		default="N")

	args = parser.parse_args()
	
	alignment = args.alignment
	percentage = args.percentage
	refgenome_name = args.refgenome_name
	output_name = args.output_results

	return alignment, percentage, refgenome_name, output_name

##########################################################################################################################################################################

def alignment_length(file):
	with open(file) as multi_fasta:
		for line in multi_fasta:
			if(re.search('^[ATGCN]',line)):
				DNA = line.rstrip('\n')
				alig_len = (len(DNA))
				break

	return alig_len

##########################################################################################################################################################################

def SNPable_loci(alig_len,file,refgenome_name,percentage,output_name):

	ref_genome_name = refgenome_name

	offender_dict = {}

	if(output_name != "N"):
		date = datetime.datetime.now()
		fmt_date = date.strftime("%y%m%d")
		results_file = open(fmt_date+'_'+output_name+'_snpable_results.tsv','a')
		results_file.write('Loci_position\tN_count\tATGC_breakdown\tOffenders\n')

	with open(file) as multi_fasta:
				
		seq_records = {}

		for seq_record in SeqIO.parse(multi_fasta, 'fasta'):
			i = seq_record.id
			s = seq_record.seq
			seq_records[i] = s

		for x in range(alig_len):
			#print("at range {}".format(x))
			A_count = 0
			T_count = 0
			C_count = 0
			G_count = 0
			N_count = 0
			Total_count = 0
			N_isolate = []

			for ID in seq_records.keys():
				if(str(ID) != ref_genome_name):
					if(str(seq_records[ID][x]) == 'A'):
						A_count +=1
						Total_count +=1
					if(str(seq_records[ID][x]) == 'T'):
						T_count +=1
						Total_count +=1
					if(str(seq_records[ID][x]) == 'G'):
						G_count +=1
						Total_count +=1
					if(str(seq_records[ID][x]) == 'C'):
						C_count +=1
						Total_count +=1
					if(str(seq_records[ID][x]) == 'N'):
						N_count +=1
						Total_count +=1
						N_isolate.append(ID)

			#if statement to see if there are Ns
			if(N_count > 0):
				#for a given loci, calc the percentage of isolates that have N at that loci
				perc_Ns = (N_count/Total_count)*100
				#calc top base count for the loci, if top base + N_count is less than total count it must mean there is another base and thus loci is a SNP
				top_count = max(A_count,T_count,C_count,G_count)
				top_and_N_count = top_count+N_count
				#to avoid out by 1 error
				loci = x+1
				if(perc_Ns <= percentage and top_and_N_count<Total_count):
					print("Total count: {}".format(Total_count))
					print("snpable3 position: {} has N: {}, A: {}, T: {}, G: {}, C: {}".format(loci,N_count,A_count,T_count,G_count,C_count))
					print("Offending isolate(s): {}\n ----".format(N_isolate))
					if(output_name != "N"):
						results_file.write(str(loci)+'\t'+str(N_count)+'\t'+'A: '+str(A_count)+ ' ,T: '+str(T_count)+ ' ,G: '+str(G_count)+ ' ,C: '+str(C_count)+'\t'+str(N_isolate)+'\n')

					for offender in N_isolate:
						if offender not in offender_dict.keys():
							offender_dict[str(offender)] = 1
						else:
							offender_dict[str(offender)] += 1
	if(output_name != "N"):
		results_file.close()
	return offender_dict


##########################################################################################################################################################################

def main():
	alignment, percentage, refgenome_name, output_name = get_argue()
	print("alignment path: {}\npercentage Threshold: {}\nreference genome name: {}\noutput_name: {}".format(alignment,percentage,refgenome_name,output_name))
	
	ali_length = alignment_length(alignment)
	print("\nAlignment length is: {}\n".format(ali_length))

	offenders = SNPable_loci(ali_length,alignment,refgenome_name,percentage,output_name)

	with open('snpable_count_results.tsv','a') as quick_results:
		for o in sorted(offenders, key=offenders.get, reverse=False):
			print(o, '\t' ,offenders[o])
			quick_results.write(str(o)+'\t'+str(offenders[o])+'\n')

##########################################################################################################################################################################

if __name__ == '__main__':
	main()

