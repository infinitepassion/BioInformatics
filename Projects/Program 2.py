'''
Project #1 :DNA Cutting-Restriction Enzymes.
Name: Manju Yadav Akkaraboina 
Course:
Date: 09/05/2017

Description:
    Program that  will read in a specific DNA strand and cut it into multiple pieces (substrings) 
    determined by the HaeIII  restriction enzyme

'''

import sys,os


def HaeIII(dna):
    
    '''
        This method is used to cut the input strand dna into pieces at every GGCC.
        Input : dna
        Return : dna substrands
    '''
    dna=dna.split("GGCC")
    i=0
    for p in range(len(dna)):
        if i==0:
            dna[p]=dna[p]+"GG"
            i+=1
        elif i==len(dna)-1:
            dna[p]="CC"+dna[p]
        else:
            
            dna[p]="CC"+dna[p]+"GG"
            i+=1
    
    return dna

#the main method begins here
if __name__ == "__main__":
    #read the dna
    dna_strand=input("Enter DNA strand ")
    #call to the HaeIII 
    dna_strand=HaeIII(dna_strand)
        
    #iterate to number of strands to print the substrand with its complementary strands
    i=0
    for dna in dna_strand:
        print()        
        i=0
        print(dna)
        #loop to iterate  every character in the substrand and print its complementary
        while i < len(dna):
            
            if dna[i]=='A':
                sys.stdout.write("T")
                
            elif dna[i]=='C':
                sys.stdout.write("G")
                
            elif dna[i]=='G':
                sys.stdout.write("C")
                
            else:
                sys.stdout.write("A")
                
            i+=1
        print()  
        

'''
The complexity of the program is linear, O(n) and n depends upon the number of substrands
that result from cutting the input dna strand  into pieces
'''

    
'''        
I have written the entire program as turned in and have not copied this code, or parts of this code from the internet or another student.      
Signature____________________
'''
    