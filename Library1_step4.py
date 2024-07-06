from Bio import SeqIO
import math
import os

UMI_per_effector=[]
UMI_per_effector.append(0)
UMIs=[]
UMI_ID=[]

for seq_record in SeqIO.parse("R1_15mer_HD2_UMIs_score13.fasta", "fasta"):
    temp=len(seq_record)
    UMI_num=int(temp/15)
    UMI_per_effector.append(UMI_num)
    for n in range(UMI_num):
        UMIs.append(seq_record.seq[n*15:(n+1)*15])
    UMI_ID.append(seq_record.id)

UMI_sum_effector=[]
UMI_sum_effector.append(0)

for k in range(len(UMI_per_effector)):
    UMI_sum_effector.append(UMI_sum_effector[k]+UMI_per_effector[k])
del UMI_sum_effector[0]

UMIs_10=[]
UMIs_10_per_effector=[]
UMIs_10_per_effector.append(0)

for i in range(1,len(UMI_per_effector),1):
   
    if UMI_per_effector[i]>9:
        for j in range(10):
            UMIs_10.append(UMIs[UMI_sum_effector[i-1]+math.floor(0.1*j*UMI_per_effector[i])])
        UMIs_10_per_effector.append(10)
    else:
        for j in range(UMI_per_effector[i]):
            UMIs_10.append(UMIs[UMI_sum_effector[i-1]+j])
        UMIs_10_per_effector.append(UMI_per_effector[i])

UMIs_10_sum_effector=[]
UMIs_10_sum_effector.append(0)

for k in range(len(UMIs_10_per_effector)):
    UMIs_10_sum_effector.append(UMIs_10_sum_effector[k]+UMIs_10_per_effector[k])
del UMIs_10_sum_effector[0]        

ofile = open("R1_15mer_HD2_UMIs_score1415_deleted_10per_effector_temp.txt", "w")
for i in range(len(UMI_ID)):
    ofile.write(">" +UMI_ID[i]+ "\n")
    for j in range(UMIs_10_per_effector[i+1]):
        ofile.write(str(UMIs_10[UMIs_10_sum_effector[i]+j]))
    ofile.write("\n")
#do not forget to close it
ofile.close()

thisFile = "R1_15mer_HD2_UMIs_score1415_deleted_10per_effector_temp.txt"
base = os.path.splitext(thisFile)[0]
os.rename(thisFile, base + ".fasta")

UMI_per_effector=[]
UMI_per_effector.append(0)
UMIs=[]
UMI_ID=[]

for seq_record in SeqIO.parse("R1_15mer_HD2_UMIs_score12.fasta", "fasta"):
    temp=len(seq_record)
    UMI_num=int(temp/15)
    UMI_per_effector.append(UMI_num)
    for n in range(UMI_num):
        UMIs.append(seq_record.seq[n*15:(n+1)*15])
    UMI_ID.append(seq_record.id)

UMI_sum_effector=[]
UMI_sum_effector.append(0)

for k in range(len(UMI_per_effector)):
    UMI_sum_effector.append(UMI_sum_effector[k]+UMI_per_effector[k])
del UMI_sum_effector[0]

UMIs_10=[]
UMIs_10_per_effector=[]
UMIs_10_per_effector.append(0)

for i in range(1,len(UMI_per_effector),1):
   
    if UMI_per_effector[i]>9:
        for j in range(10):
            UMIs_10.append(UMIs[UMI_sum_effector[i-1]+math.floor(0.1*j*UMI_per_effector[i])])
        UMIs_10_per_effector.append(10)
    else:
        for j in range(UMI_per_effector[i]):
            UMIs_10.append(UMIs[UMI_sum_effector[i-1]+j])
        UMIs_10_per_effector.append(UMI_per_effector[i])

UMIs_10_sum_effector=[]
UMIs_10_sum_effector.append(0)

for k in range(len(UMIs_10_per_effector)):
    UMIs_10_sum_effector.append(UMIs_10_sum_effector[k]+UMIs_10_per_effector[k])
del UMIs_10_sum_effector[0]        

UMI_per_effector_fasta1415=[]
UMI_per_effector_fasta1415.append(0)
UMIs_fasta1415=[]

for seq_record in SeqIO.parse("R1_15mer_HD2_UMIs_score1415_deleted_10per_effector_temp.fasta", "fasta"):
    temp=len(seq_record)
    UMI_num_fasta1415=int(temp/15)
    UMI_per_effector_fasta1415.append(UMI_num_fasta1415)
    for n in range(UMI_num_fasta1415):
        UMIs_fasta1415.append(seq_record.seq[n*15:(n+1)*15])
os.remove("R1_15mer_HD2_UMIs_score1415_deleted_10per_effector_temp.fasta")

UMI_sum_effector_fasta1415=[]
UMI_sum_effector_fasta1415.append(0)

for k in range(len(UMI_per_effector_fasta1415)):
    UMI_sum_effector_fasta1415.append(UMI_sum_effector_fasta1415[k]+UMI_per_effector_fasta1415[k])
del UMI_sum_effector_fasta1415[0]

ofile = open("R1_15mer_UMIs_10per_effector.fasta", "w")
for i in range(len(UMI_ID)):
    ofile.write(">" +UMI_ID[i]+ "\n")
    if UMIs_10_per_effector[i+1]==10:
        for j in range(UMIs_10_per_effector[i+1]):
            ofile.write(str(UMIs_10[UMIs_10_sum_effector[i]+j]))
        ofile.write("\n")
    else:
        if UMIs_10_per_effector[i+1]!=0:
            count3=0
            switch2=0
            for j in range(UMIs_10_per_effector[i+1]):
                ofile.write(str(UMIs_10[UMIs_10_sum_effector[i]+j]))

            for k in range(UMI_per_effector_fasta1415[i+1]):
                for j in range(UMIs_10_per_effector[i+1]):
                    if UMIs_fasta1415[UMI_sum_effector_fasta1415[i]+k]!=UMIs_10[UMIs_10_sum_effector[i]+j]:
                        if j==(UMIs_10_per_effector[i+1]-1):
                            ofile.write(str(UMIs_fasta1415[UMI_sum_effector_fasta1415[i]+k]))
                            count3=count3+1
                        if count3==10-UMIs_10_per_effector[i+1]:
                            switch2=1
                    else:
                        break
                if switch2==1:
                    break
            ofile.write("\n")
        else:
            for k in range(UMI_per_effector_fasta1415[i+1]):
                ofile.write(str(UMIs_fasta1415[UMI_sum_effector_fasta1415[i]+k]))
            ofile.write("\n")
ofile.close()


########### R 2 ############


UMI_per_effector=[]
UMI_per_effector.append(0)
UMIs=[]
UMI_ID=[]

for seq_record in SeqIO.parse("R2_15mer_HD2_UMIs_score13.fasta", "fasta"):
    temp=len(seq_record)
    UMI_num=int(temp/15)
    UMI_per_effector.append(UMI_num)
    for n in range(UMI_num):
        UMIs.append(seq_record.seq[n*15:(n+1)*15])
    UMI_ID.append(seq_record.id)

UMI_sum_effector=[]
UMI_sum_effector.append(0)

for k in range(len(UMI_per_effector)):
    UMI_sum_effector.append(UMI_sum_effector[k]+UMI_per_effector[k])
del UMI_sum_effector[0]

UMIs_10=[]
UMIs_10_per_effector=[]
UMIs_10_per_effector.append(0)

for i in range(1,len(UMI_per_effector),1):
   
    if UMI_per_effector[i]>9:
        for j in range(10):
            UMIs_10.append(UMIs[UMI_sum_effector[i-1]+math.floor(0.1*j*UMI_per_effector[i])])
        UMIs_10_per_effector.append(10)
    else:
        for j in range(UMI_per_effector[i]):
            UMIs_10.append(UMIs[UMI_sum_effector[i-1]+j])
        UMIs_10_per_effector.append(UMI_per_effector[i])

UMIs_10_sum_effector=[]
UMIs_10_sum_effector.append(0)

for k in range(len(UMIs_10_per_effector)):
    UMIs_10_sum_effector.append(UMIs_10_sum_effector[k]+UMIs_10_per_effector[k])
del UMIs_10_sum_effector[0]        

ofile = open("R2_15mer_HD2_UMIs_score1415_deleted_10per_effector_temp.txt", "w")
for i in range(len(UMI_ID)):
    ofile.write(">" +UMI_ID[i]+ "\n")
    for j in range(UMIs_10_per_effector[i+1]):
        ofile.write(str(UMIs_10[UMIs_10_sum_effector[i]+j]))
    ofile.write("\n")
#do not forget to close it
ofile.close()

thisFile = "R2_15mer_HD2_UMIs_score1415_deleted_10per_effector_temp.txt"
base = os.path.splitext(thisFile)[0]
os.rename(thisFile, base + ".fasta")

UMI_per_effector=[]
UMI_per_effector.append(0)
UMIs=[]
UMI_ID=[]

for seq_record in SeqIO.parse("R2_15mer_HD2_UMIs_score12.fasta", "fasta"):
    temp=len(seq_record)
    UMI_num=int(temp/15)
    UMI_per_effector.append(UMI_num)
    for n in range(UMI_num):
        UMIs.append(seq_record.seq[n*15:(n+1)*15])
    UMI_ID.append(seq_record.id)

UMI_sum_effector=[]
UMI_sum_effector.append(0)

for k in range(len(UMI_per_effector)):
    UMI_sum_effector.append(UMI_sum_effector[k]+UMI_per_effector[k])
del UMI_sum_effector[0]

UMIs_10=[]
UMIs_10_per_effector=[]
UMIs_10_per_effector.append(0)

for i in range(1,len(UMI_per_effector),1):
   
    if UMI_per_effector[i]>9:
        for j in range(10):
            UMIs_10.append(UMIs[UMI_sum_effector[i-1]+math.floor(0.1*j*UMI_per_effector[i])])
        UMIs_10_per_effector.append(10)
    else:
        for j in range(UMI_per_effector[i]):
            UMIs_10.append(UMIs[UMI_sum_effector[i-1]+j])
        UMIs_10_per_effector.append(UMI_per_effector[i])

UMIs_10_sum_effector=[]
UMIs_10_sum_effector.append(0)

for k in range(len(UMIs_10_per_effector)):
    UMIs_10_sum_effector.append(UMIs_10_sum_effector[k]+UMIs_10_per_effector[k])
del UMIs_10_sum_effector[0]        

UMI_per_effector_fasta1415=[]
UMI_per_effector_fasta1415.append(0)
UMIs_fasta1415=[]

for seq_record in SeqIO.parse("R2_15mer_HD2_UMIs_score1415_deleted_10per_effector_temp.fasta", "fasta"):
    temp=len(seq_record)
    UMI_num_fasta1415=int(temp/15)
    UMI_per_effector_fasta1415.append(UMI_num_fasta1415)
    for n in range(UMI_num_fasta1415):
        UMIs_fasta1415.append(seq_record.seq[n*15:(n+1)*15])
os.remove("R2_15mer_HD2_UMIs_score1415_deleted_10per_effector_temp.fasta")

UMI_sum_effector_fasta1415=[]
UMI_sum_effector_fasta1415.append(0)

for k in range(len(UMI_per_effector_fasta1415)):
    UMI_sum_effector_fasta1415.append(UMI_sum_effector_fasta1415[k]+UMI_per_effector_fasta1415[k])
del UMI_sum_effector_fasta1415[0]

ofile = open("R2_15mer_UMIs_10per_effector.fasta", "w")
for i in range(len(UMI_ID)):
    ofile.write(">" +UMI_ID[i]+ "\n")
    if UMIs_10_per_effector[i+1]==10:
        for j in range(UMIs_10_per_effector[i+1]):
            ofile.write(str(UMIs_10[UMIs_10_sum_effector[i]+j]))
        ofile.write("\n")
    else:
        if UMIs_10_per_effector[i+1]!=0:
            count3=0
            switch2=0
            for j in range(UMIs_10_per_effector[i+1]):
                ofile.write(str(UMIs_10[UMIs_10_sum_effector[i]+j]))

            for k in range(UMI_per_effector_fasta1415[i+1]):
                for j in range(UMIs_10_per_effector[i+1]):
                    if UMIs_fasta1415[UMI_sum_effector_fasta1415[i]+k]!=UMIs_10[UMIs_10_sum_effector[i]+j]:
                        if j==(UMIs_10_per_effector[i+1]-1):
                            ofile.write(str(UMIs_fasta1415[UMI_sum_effector_fasta1415[i]+k]))
                            count3=count3+1
                        if count3==10-UMIs_10_per_effector[i+1]:
                            switch2=1
                    else:
                        break
                if switch2==1:
                    break
            ofile.write("\n")
        else:
            for k in range(UMI_per_effector_fasta1415[i+1]):
                ofile.write(str(UMIs_fasta1415[UMI_sum_effector_fasta1415[i]+k]))
            ofile.write("\n")
ofile.close()
