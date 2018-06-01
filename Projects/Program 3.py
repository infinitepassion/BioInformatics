'''
Project #3 :Open reading frames.
Name: Manju Yadav Akkaraboina 
Course:
Date: 09/28/2017

Description:
    Program that  will read in a specific DNA strand and creates genes using the different frames and prints the start ans stop counts for a particular protien found and largest stop block and the corresponding frame number. and print the smallest and largest genes in the genebank file and print the common genes from genebank and the fasta file

'''
import sys
import re
import collections
import heapq


# Declarations of the global variables 
start_codon="ATG"
stop_codon=["TGA","TAG","TAA"]
dna=""
c_dna=""
fasta = {}
c_fasta={}
gb={}
c_gb={}
fasta_genes={}
c_fasta_genes={}


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

	

def countStartStops(dna1,f):
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
		j=len(dna)
		for i in range(frame,len(dna1),3):
			if dna1[i:i+3]==start_codon:
				st_cnt+=1
			if dna1[i:i+3] in stop_codon:
				stp_cnt+=1
				if i-j>600:
					l_st_blk+=1
				j=i
				
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
			start=int(r1.group(1))
			stop=int(r1.group(2))
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
			
			start=int(r1.group(1))
			stop=int(r1.group(2))
			c_gb[start]=stop
	return

def getSmallest(gene_dict):
	'''
	This method is used to get the smallest gene information
	Input: gene_dict- gb dictionary
	Output: k- key value from gb dictionary
	'''
	l=len(dna)
	for key in gene_dict:
		if l>=key-gene_dict[key]:
			k=key
			l=key-gene_dict[key]		
	return k

def getLargest(gene_dict):
	'''
	This method is used to get the largest gene information
	Input: gene_dict- gb dictionary
	Output: k- key value from gb dictionary
	'''
	l=0
	for key in gene_dict:
		
		if l<abs(key-gene_dict[key]):
			k=key
			l=abs(key-gene_dict[key])
	return k

def loadfasta(d,f):
	'''
		This method is used to load the fasta dictionary with the start and stops
		Input:  dna1- DNA,
				f-file pointer
		Return: None
	'''
	frame=0
	while frame<=2:
		i=frame
		flag=True
		while(flag):
			start=[]
			# Check for start codon
			if d[i:i+3] ==start_codon:
				start.append(i+1)
				for j in range(i+3,len(d)-i,3):
					# Check for intermediate starts
					if d[j:j+3]== start_codon:
						start.append(j+1)
					# check for stop codon
					if d[j:j+3] in stop_codon:
						f[j+3]=start
						i=j+3
						break
			i+=3
			if(i+frame>=len(dna)):
				flag=False	
		frame+=1
def loadcfasta(d,f):
	'''
		This method is used to load the fasta dictionary with the start and stops
		Input:  dna1- DNA,
				f-file pointer
		Return: None
	'''
	frame=0
	while frame<=2:
		i=frame
		flag=True
		while(flag):
			start=[]
			# Check for start codon
			if d[i:i+3] ==start_codon:
				start.append(len(dna)-i)
				for j in range(i+3,len(d)-i,3):
					# Check for intermediate starts
					if d[j:j+3]== start_codon:
						start.append(len(dna)-j)
					# check for stop codon
					if d[j:j+3] in stop_codon:
						f[len(dna)-j-2]=start
						i=j+3
						break
			i+=3
			if(i+frame>=len(dna)):
				flag=False	
		frame+=1
						
		
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
	countStartStops(dna,fasta)
	print("For complimentary DNA")
	countStartStops(c_dna,c_fasta)
	
	print ("\n"+"Part B"+"\n")
	s_len=getSmallest(gb)
	l_len=getLargest(gb)
	num_gene=len(gb)
	print("For forward strand:\n","the smallest gene",gb[s_len],"..",s_len,"the largest gene",gb[l_len],"..",l_len,"the number of genes", num_gene)
	s_len=getSmallest(c_gb)
	l_len=getLargest(c_gb)
	num_gene=len(c_gb)
	print("For backward strand:\n","the smallest gene",c_gb[s_len],"..",s_len,"the largest gene",c_gb[l_len],"..",l_len,"the number of genes", num_gene)
	
	
	print ("\n"+"Part C"+"\n")
	
	loadfasta(dna,fasta)
	loadcfasta(c_dna,c_fasta)
	d1_keys=set(gb.keys())
	d2_keys = set(fasta.keys())
	intersect_keys = d2_keys.intersection(d1_keys)
	for i in intersect_keys:
		for v in fasta[i]:
			if gb[i]==v:
				fasta_genes[i-v]=str(v)+".."+str(i)
	k=sorted(fasta_genes.items(), reverse=True)[:5]
	print("The largest 5 genes in the forward strand are:")
	for i in k:
		s=str(i)
		print(s.split(",")[1].rstrip(")")[2:-1])
	
	
	d1_keys=set(c_gb.keys())
	d2_keys = set(c_fasta.keys())
	intersect_keys = d2_keys.intersection(d1_keys)
	for i in intersect_keys:
		for v in c_fasta[i]:
			if c_gb[i]==v:
				c_fasta_genes[v-i]=str(i)+".."+str(v)
	k=sorted(c_fasta_genes.items(), reverse=True)[:5]
	print("\nThe largest 5 genes in the complimentary strand are:")
	for i in k:
		s=str(i)
		print(s.split(",")[1].rstrip(")")[2:-1])
	
	
	
'''        
I have written the entire program as turned in and have not copied this code, or parts of this code from the internet or another student.      
Signature____________________
'''	
	
	
	
	