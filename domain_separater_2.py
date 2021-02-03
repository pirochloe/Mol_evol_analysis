with open("sequence.txt","r") as f:
    pre_pre_sequence=f.read()
pre_sequence=pre_pre_sequence+"---X"

#空白、改行を除く
sequence="".join(pre_sequence.split())

marker=[-3]
for i in range(1,len(sequence)-3):
    if sequence[i]=="-" and sequence[i+1]=="-" and sequence[i+2]=="-" and sequence[i+3]!="=":
        marker.append(i)
    else:
        pass
    
dom_seq=[]
for i in range(15):
    dom_seq.append(sequence[marker[i]+3:marker[i+1]])


#dom_seqの各配列をfasta形式でテキストファイルに出力する
X="Propithecus"
#Xに生物種名を入力
with open("domain.txt","a") as g:
    for i in range(15):
        g.write("\n>"+X+"_dom"+str(i+1)+"\n"+dom_seq[i])
        


