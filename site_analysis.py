import copy
import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


score_matrix=pd.DataFrame(\
                          {"G":[5,1,1,0,-1,-10],"A":[1,2,1,1,1,-10],"S":[1,1,2,1,1,-10],"T":[0,1,1,3,0,-10],"P":[-1,1,1,0,6,-10],"-":[-10,-10,-10,-10,-10,-10]},\
                              index=["G","A","S","T","P","-"])

    
aminoacid_records=list(SeqIO.parse("aminoacid_seq.fasta", "fasta"))


Threshold=1
Group=[0,[1,[2,3]]]



def get_unique_list(seq):
    seen = []
    return [x for x in seq if x not in seen and not seen.append(x)]


def f(G):
    Group_4=copy.deepcopy(G)
    Group_re=[]
    for x in G:
        if isinstance(x, list):
            Group_re.append(x)
        else:
            True
    for x in Group_re:       
        if len(x)>1:
            for k in x:
                Group_4.append(k)
        else:
            True

    Group_3=get_unique_list(Group_4)
    return(Group_3)

def F(H):
    h=[]
    for n in range(100):
        H=f(H)
    for z in H:
        if isinstance(z,int):
            h.append(z)
    return(h)


GG=copy.deepcopy(Group)
for n in range(100):
    GG=f(GG)
    
GM=[]
for z in GG:
    if isinstance(z,list):
        GM.append(F(z))
GM.append(F(Group))

#単系統群の決定
s=[]
for x in GM:
    s.append(set(x))



alignment=[]
for j in range(len(aminoacid_records)):
    alignment.append(aminoacid_records[j].seq)



def check_conservation(s1,x):
    alignment_2=[]
    checker=[]
    for i in s1:
        alignment_2.append(alignment[i])
    for j in range(len(s1)-1):
        for k in range(j+1,len(s1)):
            checker.append(score_matrix[alignment_2[j][x]][alignment_2[k][x]])
    return(all([q >=Threshold for q in checker]))


def upper_set(s2):
    S2=[]
    for ss in s:
        if ss >s2:
            S2.append(ss)
        else:
            continue
    return(S2)



def check_nonconservation(s3,y):
    return(all([check_conservation(r,y)==False for r in upper_set(s3)]))


A={}
for SS in s:
    A[str(SS)]=0
    
with open("site_analysis_output.txt","w") as g:
    for n in range(len(alignment[0])):
        for S in s:
            if check_conservation(S,n)==True and check_nonconservation(S,n)==True:
                g.write("site:"+str(n+1)+"\n")
                g.write(str(S))
                g.write("\n"+"--------------------------------------"+"\n")
                A[str(S)]+=1
            else:
                continue
    g.write("============================================"+"\n")
    for SSS in s:
        g.write(str(SSS)+"\n")
        g.write(str(A[str(SSS)]))
        g.write("\n"+"--------------------------------------"+"\n")

        
        
            
                
            

        


        
            


