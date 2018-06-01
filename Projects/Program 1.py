'''
Project #1 Fibonacci Processing: Python warm up.
Name: Manju Yadav Akkaraboina 
Course:
Date: 08/29/2017

Description:
    Program that creates a list of Fibonacci numbers and then we read a locatio input from a file and then print the corresponding number of that location




'''
import sys,os

#the main method begins here

if __name__ == "__main__":
    Fibonacci =[]
    #adding  the first two numbers in the fibanocci series
    Fibonacci.append(0)
    Fibonacci.append(1)
    i=0
    #loop to generate the fibonaco for 5000 numbers after appending the first two numbers
    for x in range(0, 4999):
        Fibonacci.append(Fibonacci[i]+Fibonacci[i+1])
        i+=1
    
    
    path=os.getcwd()+"\\inputfile.txt"
    with open(path, "r") as lines:
        
        for line in lines:
            loc=int(line)
            
            print ("The Fibonacci number for",line,"is",Fibonacci[loc])
            
'''        
I have written the entire program as turned in and have not copied this code, or parts of this code from the internet or another student.      
Signature____________________
'''