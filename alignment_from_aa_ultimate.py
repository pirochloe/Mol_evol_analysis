from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

nucleotide_records = list(SeqIO.parse("nucleotide_seq.fasta", "fasta"))
aminoacid_records=list(SeqIO.parse("aminoacid_seq.fasta", "fasta"))

#CDSの開始は-1することに注意
a00=392
a01=3185

a10=189
a11=2997

a20=231
a21=3030

a30=502
a31=3283

a40=513
a41=3294

a50=0
a51=100

a60=0
a61=100

a70=0
a71=100

A=[a00,a01,a10,a11,a20,a21,a30,a31,a40,a41,a50,a51,a60,a61,a70,a71]


CDS=[]
alignment=[]
nuc_alignment=[]
out=[]

for i in range(len(nucleotide_records)):
    x=A[2*i]
    y=A[2*i+1]
    CDS.append(nucleotide_records[i].seq[x:y])
    
for j in range(len(aminoacid_records)):
    alignment.append(aminoacid_records[j].seq)
    

for i in range(len(CDS)):
    B=[]
    Y=alignment[i]
    for j in range(0,len(CDS[i]),3):
        B+=([CDS[i][j:j+3]])
    for k in range(len(Y)):
        if Y[k]==Seq("-"):
            B.insert(k,Seq("---"))
        else:
            True
    s=Seq("")
    for l in B:
        s+=l
    out.append(s)
  

records=[]  

for i in range(len(nucleotide_records)):
    records.append(SeqRecord (out[i],id=nucleotide_records[i].id,description =nucleotide_records[i].description))
    
SeqIO.write (records, 'alignment_CDS.fasta', 'fasta')

    

    




    