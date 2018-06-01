'''
Project #3c :Open reading frames.
Name: Manju Yadav Akkaraboina 
Course:
Date: 09/28/2017

Description:
    Program that  will read in a specific DNA strand and creates genes using the different frames and prints the start ans stop counts for a particular protien found and largest stop block and the corresponding frame number. and print the smallest and largest genes in the genebank file and print the common genes from genebank and the fasta file

'''
import sys
import re
import collections

# Declarations of the global variables 
start_codon="ATG"
stop_codon=["TGA","TAG","TAA"]
dna=""
c_dna=""
fasta = {}
c_fasta={}
gb={}
c_gb={}


def complimentDNA(cdna):
	'''
		This method is used to create the complementary DNA strand and reverse it
		Input: CDNA - DNA string
		Return: CDNA1 - Complimetary DNA string
	'''
	cdna1=""
	for i in range(0,len(cdna)):
			# checing to replace with the complimentary string
			if cdna[i]=='A':
				cdna1=cdna1+"T"
				
			elif cdna[i]=='C':
				cdna1=cdna1+"G"
				
			elif cdna[i]=='G':
				cdna1=cdna1+"C"
				
			else:
				cdna1=cdna1+"A"
				
	return cdna1[::-1]
	
def loadDNA(file):
	'''
		This method is used to load the DNA from the fasta file
		Input: file - file pointer
		Return: None
	'''
	global dna,c_dna
	# skipping the header line
	next(file)
	# loading the dna string
	for line in file:
		dna=dna+line.rstrip()
	# loading the complimentary dna strand
	c_dna=complimentDNA(dna)

	

def loadfasta(dna1,f):
	'''
		This method is used to load the fasta dictionary with the start and stops
		Input:  dna1- DNA,
				f-file pointer
		Return: None
	'''
	frame=0
	while ( frame<=2):
		st_cnt=0
		stp_cnt=0
		prev_j=0
		l_st_blk=0
		for i in range(frame,len(dna1),3):
			start=[]
			# Check for start codon
			if dna1[i:i+3]==start_codon:
				start.append(i+1+frame)
				st_cnt+=1
				for j in range(i+3,len(dna1),3):
					# Check for intermediate starts
					if dna1[j:j+3]==start_codon:
						start.append(j+1)
						st_cnt+=1
					# check for stop codon
					if dna1[j:j+3] in stop_codon:
						stp_cnt+=1
						n_len=j+3-i
						# check to update the largest protien size
						if j-prev_j>600:
							l_st_blk+=1
						if j not in f.keys():
							f[j+3]=i 	
							f[j+3]=(start)
							print(start)
						i=j
						prev_j=j
						break
		# printing the output
		print ("Frame:",frame," startct=",st_cnt,"stopct=",stp_cnt,"large stop blocks=",l_st_blk)
		frame+=1
		
		
		
def loadgb(file):
	'''
	This method is used to load the GB file
	Input: file-file pointer for GB file
	Output: None
	'''
	for line in file:
		r1=re.search('\s+gene+\s+([0-9]+)+\.\.([0-9]+)+\s',str(line))
		if r1 is not None:
			start=r1.group(1)
			stop=r1.group(2)
			gb[stop]=start
	return

		
def loadgbCompliment(file):
	'''
	This method is used to load the complimetary GB file
	Input: file-file pointer for GB file
	Output: None
	'''
	for line in file:
		r1=re.search('\s+gene+\s+complement+\(+([0-9]+)+\.\.([0-9]+)+\)+\s',str(line))
		if r1 is not None:
			
			start=r1.group(1)
			stop=r1.group(2)
			c_gb[stop]=start
	return

def getSmallest(gene_dict):
	'''
	This method is used to get the smallest gene information
	Input: gene_dict- gb dictionary
	Output: k- key value from gb dictionary
	'''
	l=len(dna)
	for key in gene_dict:
		if l>=int(key)-int(gene_dict[key]):
			k=key
			l=int(key)-int(gene_dict[key])		
	return k

def getLargest(gene_dict):
	'''
	This method is used to get the largest gene information
	Input: gene_dict- gb dictionary
	Output: k- key value from gb dictionary
	'''
	l=0
	for key in gene_dict:
		if l<int(key)-int(gene_dict[key]):
			k=key
			l=int(key)-int(gene_dict[key])
	return k
		
if __name__ == "__main__":
	'''
		The main method begins here
	'''
	fasta_file=open("dna.fasta",'r') #open the fasta file and load the DNA
	loadDNA(fasta_file)
	fasta_file.close()
	
	
	gb_file=open("gene.gb",'r')
	loadgb(gb_file)
	gb_file.close()
	c_gb_file=open("gene.gb",'r')
	loadgbCompliment(c_gb_file)
	c_gb_file.close()
	
	# Print the output
	
	print ("Part A"+"\n")
	print("For DNA")
	loadfasta(dna,fasta)
	print("For complimentary DNA")
	loadfasta(c_dna,c_fasta)
	
	print ("\n"+"Part B"+"\n")
	s_len=getSmallest(gb)
	l_len=getLargest(gb)
	num_gene=len(gb)
	print("For forward strand:\n", "the smallest gene",gb[s_len],"..",s_len," the largest gene",gb[l_len],"..",l_len,"the number of genes", num_gene)
	s_len=getSmallest(c_gb)
	l_len=getLargest(c_gb)
	num_gene=len(c_gb)
	print("For backward strand:\n", "the smallest gene",c_gb[s_len],"..",s_len," the largest gene",c_gb[l_len],"..",l_len,"the number of genes", num_gene)
	
	print ("\n"+"Part C"+"\n")
	
	d1_keys = set(fasta.keys())
	d2_keys = set(gb.keys())
	intersect_keys = d1_keys.intersection(d2_keys)
	print(intersect_keys)
	
	d1_keys = set(c_fasta.keys())
	d2_keys = set(c_gb.keys())
	intersect_keys = d1_keys.intersection(d2_keys)
	output=open("output.txt","w")
	for keys in gb:
		output.write(keys+"\t"+ gb[keys]+"\n")
	output.close()
	output=open("output1.txt","w")
	#print(fasta)
	for k,v in fasta.items(): 
		output.write(str(k)+"\t")
		#print(len(v))
		for i in  v:
			output.write(str(i)+" ")
		output.write("\n")
	output.close()
	print(intersect_keys)


	