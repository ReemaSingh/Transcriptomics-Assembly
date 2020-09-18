import re
with open("P_Pal.gff", "r") as f:
    for line in f:
       if "gene" in line:
            line1 = line.split("\t")      
            print line1[0]+"\t"+line1[1]+"\t"+line1[2]+"\t"+line1[3]+"\t"+line1[4]+"\t"+line1[5]+"\t"+line1[6]+"\t"+line1[7]+"\t"+"ID="+line1[8].rstrip()
       if "mRNA" in line:
           line2 = line.split("\t")
           print line2[0]+"\t"+line2[1]+"\t"+line2[2]+"\t"+line2[3]+"\t"+line2[4]+"\t"+line2[5]+"\t"+line2[6]+"\t"+line2[7]+"\t"+"ID="+line2[8].rstrip()
       if "CDS" in line:
           print line.rstrip()
       if "CDS" in line:
            exon = re.sub("CDS","exon",line)
            #exon = re.sub("ID=","ID=exon_",exon)
            print exon.rstrip()

