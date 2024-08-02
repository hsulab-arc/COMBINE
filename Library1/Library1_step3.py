from Bio import pairwise2
from Bio.pairwise2 import format_alignment
from Bio import SeqIO

UMI_per_effector=[]
UMI_per_effector.append(0)
UMIs=[]
UMI_ID=[]

for seq_record in SeqIO.parse("R1_15mer_HD2_UMIs_sorted.fasta", "fasta"):
    temp=len(seq_record)
    UMI_num=int(temp/15)
    UMI_per_effector.append(UMI_num)
    for n in range(UMI_num):
        UMIs.append(seq_record.seq[n*15:(n+1)*15])
    UMI_ID.append(seq_record.id)

Effector_seq=[]
Effector_ID=[]
for seq_record2 in SeqIO.parse("R1_Forward_P5_301cycles.fasta", "fasta"):
    temp=len(seq_record2)
    Effector_seq.append(seq_record2.seq)
    Effector_ID.append(seq_record2.id)

UMI_sum_effector=[]
UMI_sum_effector.append(0)


for k in range(len(UMI_per_effector)):
    UMI_sum_effector.append(UMI_sum_effector[k]+UMI_per_effector[k])
del UMI_sum_effector[0]


UMIs_to_be_deleted=[]

for i in range(len(Effector_ID)-1):
    for j in range(UMI_per_effector[i+1]):
        for k in [x for x in range(len(Effector_ID)-1) if x!=i]:
            if pairwise2.align.localms(str(UMIs[UMI_sum_effector[i]+j]),str(Effector_seq[k]),1,0,-1,-1,score_only=True)>=14:
                UMIs_to_be_deleted.append(UMI_sum_effector[i]+j)

                
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

ofile = open("R1_15mer_HD2_UMIs_score13.fasta", "w")

for i in range(len(UMI_ID)):
    ofile.write(">" +UMI_ID[i]+ "\n")
    for j in range(Improved_UMI_sum_effector[i],Improved_UMI_sum_effector[i+1],1):
        ofile.write(str(UMIs[j]))
    ofile.write("\n")
ofile.close()

UMI_per_effector=[]
UMI_per_effector.append(0)
UMIs=[]
UMI_ID=[]

for seq_record in SeqIO.parse("R1_15mer_HD2_UMIs_sorted.fasta", "fasta"):
    temp=len(seq_record)
    UMI_num=int(temp/15)
    UMI_per_effector.append(UMI_num)
    for n in range(UMI_num):
        UMIs.append(seq_record.seq[n*15:(n+1)*15])
    UMI_ID.append(seq_record.id)

UMIs_to_be_deleted=[]

for i in range(len(Effector_ID)-1):
    for j in range(UMI_per_effector[i+1]):
        for k in [x for x in range(len(Effector_ID)-1) if x!=i]:
            if pairwise2.align.localms(str(UMIs[UMI_sum_effector[i]+j]),str(Effector_seq[k]),1,0,-1,-1,score_only=True)>=13:
                UMIs_to_be_deleted.append(UMI_sum_effector[i]+j)

                
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

ofile = open("R1_15mer_HD2_UMIs_score12.fasta", "w")

for i in range(len(UMI_ID)):
    ofile.write(">" +UMI_ID[i]+ "\n")
    for j in range(Improved_UMI_sum_effector[i],Improved_UMI_sum_effector[i+1],1):
        ofile.write(str(UMIs[j]))
    ofile.write("\n")
ofile.close()


############### R 2 ################

UMI_per_effector=[]
UMI_per_effector.append(0)
UMIs=[]
UMI_ID=[]

for seq_record in SeqIO.parse("R2_15mer_HD2_UMIs_sorted.fasta", "fasta"):
    temp=len(seq_record)
    UMI_num=int(temp/15)
    UMI_per_effector.append(UMI_num)
    for n in range(UMI_num):
        UMIs.append(seq_record.seq[n*15:(n+1)*15])
    UMI_ID.append(seq_record.id)

Effector_seq=[]
Effector_ID=[]

for seq_record2 in SeqIO.parse("R2_Reverse_P7_301cycles.fasta", "fasta"):

    temp=len(seq_record2)
    Effector_seq.append(seq_record2.seq)
    Effector_ID.append(seq_record2.id)

UMI_sum_effector=[]
UMI_sum_effector.append(0)

for k in range(len(UMI_per_effector)):
    UMI_sum_effector.append(UMI_sum_effector[k]+UMI_per_effector[k])
del UMI_sum_effector[0]

UMIs_to_be_deleted=[]

for i in range(len(Effector_ID)-1):
    for j in range(UMI_per_effector[i+1]):
        for k in [x for x in range(len(Effector_ID)-1) if x!=i]:
            if pairwise2.align.localms(str(UMIs[UMI_sum_effector[i]+j]),str(Effector_seq[k]),1,0,-1,-1,score_only=True)>=14:
                UMIs_to_be_deleted.append(UMI_sum_effector[i]+j)

                
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

ofile = open("R2_15mer_HD2_UMIs_score13.fasta", "w")

for i in range(len(UMI_ID)):
    ofile.write(">" +UMI_ID[i]+ "\n")
    for j in range(Improved_UMI_sum_effector[i],Improved_UMI_sum_effector[i+1],1):
        ofile.write(str(UMIs[j]))
    ofile.write("\n")
ofile.close()

UMI_per_effector=[]
UMI_per_effector.append(0)
UMIs=[]
UMI_ID=[]

for seq_record in SeqIO.parse("R2_15mer_HD2_UMIs_sorted.fasta", "fasta"):
    temp=len(seq_record)
    UMI_num=int(temp/15)
    UMI_per_effector.append(UMI_num)
    for n in range(UMI_num):
        UMIs.append(seq_record.seq[n*15:(n+1)*15])
    UMI_ID.append(seq_record.id)

UMIs_to_be_deleted=[]

for i in range(len(Effector_ID)-1):
    for j in range(UMI_per_effector[i+1]):
        for k in [x for x in range(len(Effector_ID)-1) if x!=i]:
            if pairwise2.align.localms(str(UMIs[UMI_sum_effector[i]+j]),str(Effector_seq[k]),1,0,-1,-1,score_only=True)>=13:
                UMIs_to_be_deleted.append(UMI_sum_effector[i]+j)

                
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

ofile = open("R2_15mer_HD2_UMIs_score12.fasta", "w")

for i in range(len(UMI_ID)):
    ofile.write(">" +UMI_ID[i]+ "\n")
    for j in range(Improved_UMI_sum_effector[i],Improved_UMI_sum_effector[i+1],1):
        ofile.write(str(UMIs[j]))
    ofile.write("\n")
ofile.close()
