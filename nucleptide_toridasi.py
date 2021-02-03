with open("nucleotide.txt","r") as f:
    result=f.read()
Y=""
for i in range(len(result)):
    if result[i] in ["A","T","G","C"]:
        Y+=result[i]
x=1480
y=3390

#x番目からy番目までの配列を取り出し

B=Y[x-1:y]
with open("outputfile.txt","w") as g:
    g.write(B)
