import copy
import numpy as np
import pandas as pd
from pandas import Series
from Bio import SeqIO
import matplotlib.pyplot as plt
import os



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


#アミノ酸の置換行列,系統関係を設定

score_matrix=PAM250
Group=[[9,10],[[[5,6],[7,8]],[4,[[2,3],[0,1]]]]]
site_range=True




if site_range==True:
    test_table=pd.read_csv("select_site.csv")
    test_table=test_table.dropna(how='all')

    site_list=[]

    for i in test_table.index.values:
        x=test_table.at[i,"start"]
        y=test_table.at[i,"stop"]
        for j in range(x-1,y):
            site_list.append(j)



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

Threshold=6


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

def downer_set(s7):
    S10=[]
    for ssss in s:
        if ssss <s7:
            S10.append(ssss)
        else:
            continue
    return(S10)
        

def check_nonconservation(s3,y):
    return(all([check_conservation(r,y)==False for r in upper_set(s3)]))



def value_function(s4,z):
    alignment3=[]
    total_point=[]
    for i in s4:
        alignment3.append(alignment[i])
    for j in range(len(s4)-1):
        for k in range(j+1,len(s4)):
            total_point.append(score_matrix[alignment3[j][z]][alignment3[k][z]])
    total_point_2=np.array(total_point)
    return(total_point_2)


if site_range==False:
    site_list=range(len(alignment[0]))

def CS(s5,site_1):
    alignment4=[]
    checker_2=[]
    for j in s5:
        alignment4.append(alignment[j])
    for k in range(len(s5)-1):
        for l in range(k+1,len(s5)):
            checker_2.append(score_matrix[alignment4[k][site_1]][alignment4[l][site_1]])
    return(min(checker_2))

def check_super_conservation(s6,site_2):
    return(all([CS(stutter,site_2)<CS(s6,site_2) for stutter in upper_set(s6)]))

def check_ultra_conservation(s13,site_3):
    return(all([CS(splite,site_3)==CS(s13,site_3) for splite in downer_set(s13)]))

def super_upper_set(s1):
    s4=[]
    if upper_set(s1)==[]:
        s3=[]
    else:
        s3=upper_set(s1)[0]
        for k in upper_set(s1):
            s3=s3.intersection(k)
        s3=s3-s1
        
        for l in upper_set(s1):
            s4.append(l-s3)
        s4.remove(s1)
    return(s4)
    
os.mkdir("text_file")
os.mkdir("Histogram")
os.mkdir("csv_file")


A={}
for SS in s:
    A[str(SS)]=[]
   
with open("text_file/conservative_site_analysis.txt","w") as g:
    g.write("<<site_analysis>>"+"\n")
    for i in range(3):
        g.write("\n")
    for n in site_list:
        for S in s:
            if check_super_conservation(S,n)==True:
                A[str(S)].append(n+1)
                g.write("site:"+str(n+1)+"\n")
                g.write("conservative_group:"+str(S)+"\n")
                g.write("Threshold_value:"+str(CS(S,n))+"\n")
                g.write("std:"+str(np.std(value_function(S,n))))
                g.write("\n"+"--------------------------------------"+"\n")
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
        g.write("conservative_site:"+str(A[str(SSS)])+"\n")
        g.write("conservation_count:"+str(len(A[str(SSS)]))+"\n")
        C=[]
        for spell in A[str(SSS)]:
            C.append(CS(SSS,spell-1))
        g.write("conservation_value:"+str(sum(C)/len(C))+"\n")
        g.write("\n"+"--------------------------------------"+"\n")


for focus_group in s:
    B=[]
    for x in A[str(focus_group)]:
        B.append(CS(focus_group,x-1))
    C=np.array(B)
    D=[]
    D_2=[]
    D_3=[]
    d=[]
    for spell_stutter in A[str(focus_group)]:
        if check_ultra_conservation(focus_group,spell_stutter-1)==True and all([CS(MT,spell_stutter-1)<CS(focus_group,spell_stutter-1) for MT in super_upper_set(focus_group)])==True:
            D.append(CS(focus_group,spell_stutter-1))
            D_2.append(spell_stutter)
        elif check_ultra_conservation(focus_group,spell_stutter-1)==True and all([CS(MT,spell_stutter-1)<CS(focus_group,spell_stutter-1) for MT in super_upper_set(focus_group)])==False:
            D_3.append(spell_stutter)
            d.append(CS(focus_group,spell_stutter-1))
        else:
            continue
    E=np.array(D)
    e=np.array(d)
    fig=plt.figure()
    plt.hist([E,e],label=["ultra","super"],ec="k")
    plt.legend()
    plt.title("conservation_value")
    plt.grid()
    fig.savefig("Histogram/histogram_"+str(focus_group)+".png")
    
    with open("text_file/group_analysis_"+str(focus_group)+".txt","w") as I:
        I.write("ultra_conservative_site:"+"\n")
        I.write(str(D_2)+"\n"+"\n"+"\n")
        I.write("super_conservaive_site:"+"\n")
        I.write(str(D_3)+"\n"+"\n"+"\n")
        I.write("ultra_conservation_count:"+str(len(D_2))+"\n"+"\n"+"\n")
        I.write("super_conservation_count:"+str(len(D_3))+"\n")
    F=Series(E,name="ultra")
    f=Series(e,name="super")
    g=f.describe()
    G=F.describe()
    H=pd.concat([G,g],axis=1)
    H.to_csv("csv_file/conservation_site_analysis_"+str(focus_group)+".csv")
    
    




