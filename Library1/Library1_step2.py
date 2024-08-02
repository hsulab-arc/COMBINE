from Bio import SeqIO
import hamdist
UMI_per_effector=[]
UMI_per_effector.append(0)
UMIs=[]
UMI_ID=[]

for seq_record in SeqIO.parse("kmercollection_R1.fasta", "fasta"):
    temp=len(seq_record)
    UMI_num=int(temp/15)
    UMI_per_effector.append(UMI_num)
    for n in range(UMI_num):
        UMIs.append(seq_record.seq[n*15:(n+1)*15])
    UMI_ID.append(seq_record.id)

dmin=len(UMIs[0])
count=0

UMI_sum_effector=[]
UMI_sum_effector.append(0)

for k in range(len(UMI_per_effector)):
    UMI_sum_effector.append(UMI_sum_effector[k]+UMI_per_effector[k])
del UMI_sum_effector[0]

UMIs_to_be_deleted=[]
for i in range(len(UMI_sum_effector)-1):
        for i2 in range(UMI_per_effector[i+1]):
            for j in range(i+1,len(UMI_sum_effector)-1):
                for j2 in range(UMI_per_effector[j+1]):       
                    dist=hamdist.hamdist(UMIs[UMI_sum_effector[i]+i2],UMIs[UMI_sum_effector[j]+j2])
                    if dist<2:
                        if UMI_per_effector[i+1]>UMI_per_effector[j+1]:
                            UMIs_to_be_deleted.append(UMI_sum_effector[i]+i2)
                        else:
                            UMIs_to_be_deleted.append(UMI_sum_effector[j]+j2)
unique_UMIs_to_be_deleted=list(set(UMIs_to_be_deleted))

for index in sorted(unique_UMIs_to_be_deleted, reverse=True):
    for l in reversed(range((len(UMI_sum_effector)-1))): # (len(UMI_sum_effector)-1)=156
        if UMI_sum_effector[l]-1<index:
            if index<=(UMI_sum_effector[l+1]-1):
                UMI_per_effector[l+1]=UMI_per_effector[l+1]-1
                del UMIs[index] 
                break

Improved_UMI_sum_effector=[]
Improved_UMI_sum_effector.append(0)

for k in range(len(UMI_per_effector)-1):
    Improved_UMI_sum_effector.append(Improved_UMI_sum_effector[k]+UMI_per_effector[k+1])

ofile = open("R1_15mer_HD2_UMIs.fasta", "w")

for i in range(len(UMI_ID)):
    ofile.write(">" +UMI_ID[i]+ "\n")
    for j in range(Improved_UMI_sum_effector[i],Improved_UMI_sum_effector[i+1],1):
        ofile.write(str(UMIs[j]))
    ofile.write("\n")
ofile.close()

from Bio import SeqIO
import hamdist
UMI_per_effector=[]
UMI_per_effector.append(0)
UMIs=[]
UMI_ID=[]

for seq_record in SeqIO.parse("kmercollection_R2.fasta", "fasta"):
    temp=len(seq_record)
    UMI_num=int(temp/15)
    UMI_per_effector.append(UMI_num)
    for n in range(UMI_num):
        UMIs.append(seq_record.seq[n*15:(n+1)*15])
    UMI_ID.append(seq_record.id)
dmin=len(UMIs[0])
count=0

UMI_sum_effector=[]
UMI_sum_effector.append(0)

for k in range(len(UMI_per_effector)):
    UMI_sum_effector.append(UMI_sum_effector[k]+UMI_per_effector[k])
del UMI_sum_effector[0]

UMIs_to_be_deleted=[]
for i in range(len(UMI_sum_effector)-1):
        for i2 in range(UMI_per_effector[i+1]):
            for j in range(i+1,len(UMI_sum_effector)-1):
                for j2 in range(UMI_per_effector[j+1]):       
                    dist=hamdist.hamdist(UMIs[UMI_sum_effector[i]+i2],UMIs[UMI_sum_effector[j]+j2])
                    if dist<2:
                        if UMI_per_effector[i+1]>UMI_per_effector[j+1]:
                            UMIs_to_be_deleted.append(UMI_sum_effector[i]+i2)
                        else:
                            UMIs_to_be_deleted.append(UMI_sum_effector[j]+j2)

unique_UMIs_to_be_deleted=list(set(UMIs_to_be_deleted))

for index in sorted(unique_UMIs_to_be_deleted, reverse=True):
    for l in reversed(range((len(UMI_sum_effector)-1))): 
        if UMI_sum_effector[l]-1<index:
            if index<=(UMI_sum_effector[l+1]-1):
                UMI_per_effector[l+1]=UMI_per_effector[l+1]-1
                del UMIs[index] 
                break

Improved_UMI_sum_effector=[]
Improved_UMI_sum_effector.append(0)

for k in range(len(UMI_per_effector)-1):
    Improved_UMI_sum_effector.append(Improved_UMI_sum_effector[k]+UMI_per_effector[k+1])

ofile = open("R2_15mer_HD2_UMIs.fasta", "w")

for i in range(len(UMI_ID)):
    ofile.write(">" +UMI_ID[i]+ "\n")
    for j in range(Improved_UMI_sum_effector[i],Improved_UMI_sum_effector[i+1],1):
        ofile.write(str(UMIs[j]))
    ofile.write("\n")
ofile.close()
