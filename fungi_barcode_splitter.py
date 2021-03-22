#!/usr/bin/env python
import sys

from Bio import SeqIO, bgzf
from gzip import open as gzopen

#The old barcode splitter was far too complicated in use for the fungal dept.
#Thats why I programmed a new one just based on a hardcoded barcode list (bob10, JUL 2020)

barcode_list = ['CTCGGCCA', 'CTGCGGCA', 'GACCTTCA', 'GACGACCA', 'GATTAGCA', 'GCATCCCA', 'GCCAGCCA' , 'GCGACTCA']

def recognize_bc(records, barcode):
	filtered_records = []
	len_bc = len(barcode)
	for record in records:
		if record.seq.startswith(barcode):
			filtered_records.append(record[len_bc:])
	return filtered_records
			
#check arguments, else print help
if len(sys.argv) < 3: 
	sys.stderr.write("Usage: fungi_barcode_splitter.py *.fastq.gz\n")
	sys.stderr.write("At least two fastq.gz files should be provided.")
	sys.exit(1)
	

for i in range(1,len(sys.argv)):
	curr_file = sys.argv[i] #get filename as argument
	#R2 has the barcode!
	if 'R2' in curr_file:
		print curr_file
		if curr_file.replace('R2', 'R1') in sys.argv:
			print 'Found matching files: ' + curr_file + '; ' + curr_file.replace('R2', 'R1')
			sync_file = curr_file.replace('R2', 'R1')
			for barcode in barcode_list:
				records = SeqIO.parse(gzopen(curr_file, 'rt'), format='fastq')
				filtered_records = recognize_bc(records, barcode)
				print 'Found barcode  ' + barcode + ': ' + str(len(filtered_records))
				with bgzf.BgzfWriter(curr_file.replace('.fastq.gz', '_' +barcode + '.fastq.gz'), "wb") as outgz:
					SeqIO.write(sequences=filtered_records, handle=outgz, format="fastq")
				#Sync R1
				records = SeqIO.parse(gzopen(sync_file, 'rt'), format='fastq')
				records_sync = []
				for record in records:
					if len(filtered_records) > 0 and record.id == filtered_records[0].id:
						records_sync.append(record)
						del filtered_records[0]
				print 'Synced barcode ' + barcode + ': ' + str(len(records_sync))
				with bgzf.BgzfWriter(sync_file.replace('.fastq.gz', '_' +barcode + '.fastq.gz'), "wb") as outgz:
					SeqIO.write(sequences=records_sync, handle=outgz, format="fastq")