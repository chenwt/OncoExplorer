# BLAST cmd

###saskia###
# Unmasked version
makeblastdb -in pptd.fa -dbtype prot -parse_seqids -out pptd -title "Unmasked pptd sequences"
'''
Building a new DB, current time: 02/26/2017 11:41:30
New DB name:   /usr1/public/yifeng/Github/blast/ncbi-blast-2.6.0+/blastdb/pptd
New DB title:  Unmasked pptd sequences
Sequence type: Protein
Keep MBits: T
Maximum file size: 1000000000B
Adding sequences from FASTA; added 2532 sequences in 0.137596 seconds.
'''
# Check it
blastdbcmd -db pptd -info
'''
Database: Unmasked pptd sequences
	2,532 sequences; 1,779,511 total residues

Date: Feb 26, 2017  11:41 AM	Longest sequence: 14,507 residues

Volumes:
	/usr1/public/yifeng/Github/blast/ncbi-blast-2.6.0+/blastdb/pptd
'''

# cp pptd.fa ../workd
# cp test.fa ../workd

blastp -query test.fa -db pptd -task blastp -num_threads 30 -outfmt "0" -num_descriptions 3000 -num_alignments 0 -out out_unmask_pptd.txt

