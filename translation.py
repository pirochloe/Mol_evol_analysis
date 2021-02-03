codontable={"TTT":"F","TTC":"F","TTA":"L","TTG":"L","TCT":"S","TCC":"S","TCA":"S",\
            "TCG":"S","TAT":"Y","TAC":"Y","TAA":".","TAG":".","TGT":"C","TGC":"C",\
            "TGA":".","TGG":"W","CTT":"L","CTC":"L","CTA":"L","CTG":"L","CCT":"P",\
            "CCC":"P","CCA":"P","CCG":"P","CAT":"H","CAC":"H","CAA":"Q","CAG":"Q",\
            "CGT":"R","CGC":"R","CGA":"R","CGG":"R","ATT":"I","ATC":"I","ATA":"I",\
            "ATG":"M","ACT":"T","ACC":"T","ACA":"T","ACG":"T","AAT":"N","AAC":"N",\
            "AAA":"K","AAG":"K","AGT":"S","AGC":"S","AGA":"R","AGG":"R","GTT":"V",\
            "GTC":"V","GTA":"V","GTG":"V","GCT":"A","GCC":"A","GCA":"A","GCG":"A",\
            "GAT":"D","GAC":"D","GAA":"E","GAG":"E","GGT":"G","GGC":"G","GGA":"G",\
            "GGG":"G"}
with open("sequence.txt","r") as f:
    result=f.read()
Z=""
for i in range(len(result)):
    if result[i] in ["A","T","G","C"]:
        Z+=result[i]
B=[]
for j in range(0,len(Z),3):
    B+=[Z[j:j+3]]
X=""
for k in range(len(B)):
    X+=codontable[B[k]]
with open("trans4.txt","w") as g:
    g.write(X)
