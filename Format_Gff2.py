import re
with open("P_Pal_1.gff", "r") as f:
    for line in f:
       if "gene" in line:
          for part in line.split():
              if "ID=" in part:
                 Name = part
                 Name = re.sub("ID=","Name=",Name)
                 #print part+";"+Name
          line1 = line.split("\t")
          print line1[0]+"\t"+line1[1]+"\t"+line1[2]+"\t",line1[3]+"\t"+line1[4]+"\t"+line1[5]+"\t"+line1[6]+"\t"+line1[7]+"\t"+part+";"+Name

       if "mRNA" in line:
          for part1 in line.split():
             if "ID=" in part1:
               Name1 = part
               Name1 = re.sub("ID=","Name=mRNA_",Name1)
               Parent = part
               Parent = re.sub("ID=","Parent=",Parent)
              # print part1+";"+Name1+";"+Parent
          line2 = line.split("\t")
          print line2[0]+"\t"+line2[1]+"\t"+line2[2]+"\t"+line2[3]+"\t"+line2[4]+"\t"+line2[5]+"\t"+line2[6]+"\t"+line2[7]+"\t"+part1+";"+Name1+";"+Parent

       if "CDS" in line:
           for part2 in line.split():
              if "ID=" in part2:
                part2 = re.sub(";Name=.*","",part2)
                Name2 = part
                Name2 = re.sub("ID=","Name=CDS_",Name1)
                Parent1 = part1
                Parent1 = re.sub("ID=","Parent=",Parent1)
           line3 = line.split("\t")
           print line3[0]+"\t"+line3[1]+"\t"+line3[2]+"\t"+line3[3]+"\t"+line3[4]+"\t"+line3[5]+"\t"+line3[6]+"\t"+line3[7]+"\t"+part2+";"+Name2+";"+Parent1
          
       if "exon" in line:
           for part3 in line.split():
              if "ID=" in part3:
                part3= re.sub(";Name=.*","",part3)
                Name3= part1
                Name3= re.sub("ID=","Name=exon_",Name1)
                Parent2 = part1
                Parent2 = re.sub("ID=","Parent=",Parent2)
           line4 = line.split("\t")
           print line4[0]+"\t"+line4[1]+"\t"+line4[2]+"\t"+line4[3]+"\t"+line4[4]+"\t"+line4[5]+"\t"+line4[6]+"\t"+line4[7]+"\t"+part3+";"+Name3+";"+Parent2
