import copy
import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


PAM250=pd.DataFrame(\
                          {"G":[5,1,1,0,-1,-4,-3,-3,-1,1,0,0,-1,-5,-5,-7,-2,-3,-2,-3,-10],\
                           "A":[1,2,1,1,1,-2,-1,-1,0,0,0,0,0,-4,-3,-6,-1,-2,-1,-2,-10],\
                               "S":[1,1,2,1,1,-3,-1,-2,-1,0,1,0,-1,-3,-3,-2,0,0,-1,0,-10],\
                                   "T":[0,1,1,3,0,-2,0,-1,0,0,0,0,-1,-3,-3,-5,0,-1,-1,-2,-10],\
                                       "P":[-1,1,1,0,6,-3,-2,-2,-1,-1,-1,-1,0,-5,-5,-6,-1,0,0,-3,-10],\
                                           "L":[-4,-2,-3,-2,-3,6,2,4,2,-4,-3,-3,-2,2,-1,-2,-3,-3,-2,-6,-10],\
                                               "I":[-3,-1,-1,0,-2,2,5,2,4,-2,-2,-2,-2,1,-1,-5,-2,-2,-2,-2,-10],\
                                                   "M":[-3,-1,-2,-1,-2,4,2,6,2,-3,-2,-2,-1,0,-2,-4,0,0,-2,-5,-10],\
                                                       "V":[-1,0,-1,0,-1,2,4,2,4,-2,-2,-2,-2,-1,-2,-6,-2,-2,-2,-2,-10],\
                                                           "D":[1,0,0,0,-1,-4,-2,-3,-2,4,2,3,2,-6,-4,-7,0,-1,1,-5,-10],\
                                                               "N":[0,0,1,0,-1,-3,-2,-2,-2,2,2,1,1,-4,-2,-4,1,0,2,-4,-10],\
                                                                   "E":[0,0,0,0,-1,-3,-2,-2,-2,3,1,4,2,-5,-4,-7,0,-1,1,-5,-10],\
                                                                       "Q":[-1,0,-1,-1,0,-2,-2,-1,-2,2,1,2,4,-5,-4,-5,1,1,3,-5,-10],\
                                                                           "F":[-5,-4,-3,-3,-5,2,1,0,-1,-6,-4,-5,-5,9,7,0,-5,-4,-2,-4,-10],\
                                                                               "Y":[-5,-3,-3,-3,-5,-1,-1,-2,-2,-4,-2,-4,-4,7,10,0,-4,-4,0,0,-10],\
                                                                                   "W":[-7,-6,-2,-5,-6,-2,-5,-4,-6,-7,-4,-7,-5,0,0,17,-3,2,-3,-8,-10],\
                                                                                       "K":[-2,-1,0,0,-1,-3,-2,0,-2,0,1,0,1,-5,-4,-3,5,3,0,-5,-10],\
                                                                                           "R":[-3,-2,0,-1,0,-3,-2,0,-2,-1,0,-1,1,-4,-4,2,3,6,2,-4,-10],\
                                                                                               "H":[-2,-1,-1,-1,0,-2,-2,-2,-2,1,2,1,3,-2,0,-3,0,2,6,-3,-10],\
                                                                                                   "C":[-3,-2,0,-2,-3,-6,-2,-5,-2,-5,-4,-5,-5,-4,0,-8,-5,-4,-3,12,-10],\
                                                                                                       "-":[-10,-10,-10,-10,-10,-10,-10,-10,-10,-10,-10,-10,-10,-10,-10,-10,-10,-10,-10,-10,-10]},\
                              index=["G","A","S","T","P","L","I","M","V","D","N","E","Q","F","Y","W","K","R","H","C","-"])

 
aminoacid_records=list(SeqIO.parse("aminoacid_seq.fasta", "fasta"))


#アミノ酸の置換行列、保存度の閾値、系統関係を設定

score_matrix=PAM250
Threshold=6
Group=[[3,4],[2,[0,1]]]



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



def value_function(s4,z):
    alignment3=[]
    total_point=0
    mass=0
    for i in s4:
        alignment3.append(alignment[i])
    for j in range(len(s4)-1):
        for k in range(j+1,len(s4)):
            total_point+=score_matrix[alignment3[j][z]][alignment3[k][z]]
            mass+=1
    return(total_point/mass)



A={}
for SS in s:
    A[str(SS)]=0
    
with open("site_analysis_output.txt","w") as g:
    g.write("<<site_analysis>>"+"\n")
    for i in range(3):
        g.write("\n")
    for n in range(len(alignment[0])):
        for S in s:
            if check_conservation(S,n)==True and check_nonconservation(S,n)==True:
                g.write("site:"+str(n+1)+"\n")
                g.write("conservation_group:"+str(S)+"\n")
                g.write("conservation_value:"+str(value_function(S,n)))
                g.write("\n"+"--------------------------------------"+"\n")
                A[str(S)]+=1
            else:
                continue
    for i in range(6):
        g.write("\n")
    g.write("============================================"+"\n")
    g.write("<<group_analysis>>"+"\n")
    for i in range(3):
        g.write("\n")
    for SSS in s:
        g.write(str(SSS)+"\n")
        g.write("conservation_count:"+str(A[str(SSS)]))
        g.write("\n"+"--------------------------------------"+"\n")

        
        
            
                

        


        
            



