
from optparse import OptionParser
import os

import re

import datetime
import glob
import time
# import configparser
import gzip

import Drugseqini2 as Drugseqini


version="1.0"

datename=time.strftime('%Y%m_%d%H%M',time.localtime(time.time()))

parser = OptionParser(usage="usage: %prog -f readFold [-p <int>] [-e <int>] ", version="%prog 1")
parser.add_option("-f",
        action="store", 
        dest="readFolder",
        help="Folder with RNA-seq data (gzipped fastq files)")

parser.add_option("-b",
        action="store", 
        dest="barcode",
        help="barcode file")

parser.add_option("-o",
        action="store", 
        dest="output",
        
        help="output fold")

parser.add_option("-m",
        action="store", 
        dest="module",
        default="a",
        help="module:such as m/a(all)")


# (options, args) = parser.parse_args()

# def transformer(FRdata,FRdata2,WRdata1,WRdata2,same_seq="AAGCAGTGGTATCAACGCAGAGT"):
def transformer(FRdata1,FRdata2,WRdata1,same_seq="CAGTGGTATCAACGCAGA",WRdata1D="",BClist=[],Outdir=""):
    index=0
    datanew=open(Outdir+"/"+"datanew_R1.fastq","w+")
    datanew2=open(Outdir+"/"+"datanew_R2.fastq","w+")

    Obread1=open(Outdir+"/"+"Obread1.fastq","w+")
    Obread2=open(Outdir+"/"+"Obread2.fastq","w+")


    while True:
        line1=FRdata1.readline().rstrip().decode('utf-8')
        if not line1:
        # if index==10:
            print(index)
            break
        line2=FRdata1.readline().rstrip().decode('utf-8')

        lineUMI_before=line2.split(same_seq)
        line3=FRdata1.readline().rstrip().decode('utf-8')
        line4=FRdata1.readline().rstrip().decode('utf-8')
        line2_1=FRdata2.readline().rstrip().decode('utf-8') 
        line2_2=FRdata2.readline().rstrip().decode('utf-8')
        line2_3=FRdata2.readline().rstrip().decode('utf-8')
        line2_4=FRdata2.readline().rstrip().decode('utf-8')  
        if len(lineUMI_before)>1:
            lineUMI=lineUMI_before[1]
            lineUMI2=str(lineUMI)[0:12]
            if lineUMI2 in BClist:
                # print(lineUMI_before)
                # print(lineUMI2)
                if len(lineUMI)>23:
                    # datanew.write(line1+"\n"+lineUMI_before[0][0:20]+"\n"+line3+"\n"+line4+"\n")
                    # datanew2.write(line2_1+"\n"+line2_2+"\n"+line2_3+"\n"+line2_4+"\n")
                    line1BC=line1.split(" ")[0]+"_"+lineUMI2+"_"+lineUMI[12:22]+" "+line1.split(" ")[1]
                    line2_1BC=line2_1.split(" ")[0]+"_"+lineUMI2+"_"+lineUMI[12:22]+" "+line2_1.split(" ")[1]
                    line2toBC=line2.split(same_seq)[1][22:]
                    # [0:20]
                    # print(line2toBC)
                    line4toBC=line4[len(line4)-len(line2toBC):len(line4)]
                    # datanew.write(line1+"\n"+line2toBC+"\n"+line3+"\n"+line4toBC+"\n")
                    datanew.write(line1BC+"\n"+line2toBC+"\n"+line3+"\n"+line4toBC+"\n")
                    datanew2.write(line2_1BC+"\n"+line2_2+"\n"+line2_3+"\n"+line2_4+"\n")

                index+=1

            else:
                Obread1.write(line1+"\n"+line2+"\n"+line3+"\n"+line4+"\n")
                Obread2.write(line2_1+"\n"+line2_2+"\n"+line2_3+"\n"+line2_4+"\n")
                WRdata1.write(lineUMI2+"\n")  
        else :
            # WRdata1.write("NA"+"\n")
            WRdata1D.write(line2+"\n")


def countUMI(gtfname="",genomedir="",fastqname="datanew_R2.fastq ",featurecountthread="4",starthread="8",outdirname=".",outname="",featureCount_g='gene_name'):
    infile=os.path.join(outdirname,fastqname)
    outfile=os.path.join(outdirname,outname)
    outbam=os.path.join(outdirname,outname+"_Aligned.sortedByCoord.out.bam")
    os.system("STAR --runThreadN "+starthread+" --genomeDir "+genomedir+"\
                        --readFilesIn   "+infile+"  \
                        --outFileNamePrefix "+outdirname+"/"+outname+"_"+" \
                        --outFilterMultimapNmax 1  \
                        --outSAMtype BAM SortedByCoordinate")
    print(outbam)
    # @date 2024/11/15 增加了"-g gene_name,用于使得定量的矩阵默认行名是基因名称
    os.system("featureCounts -a "+gtfname+" -o gene_assigned -R BAM "+outbam+"  -T "+featurecountthread + " -g " + featureCount_g)
    os.system("mv gene_assigned* "+outdirname)
    os.system("mv *featureCounts* "+outdirname)



    os.system("samtools sort "+outbam+".featureCounts.bam -o "+outfile+"_assigned_sorted.bam")
    os.system("samtools index "+outfile+"_assigned_sorted.bam")
    # os.system("umi_tools count --per-gene --gene-tag=XT --assigned-status-tag=XS --per-cell -I "+outfile \
        # +"assigned_sorted.bam  "+outdirname+"/"+outname+"_percell.csv")"assigned_sorted.bam"+" -S "+outfile+"_counts.tsv.gz")
    os.system("umi_tools count --per-gene --gene-tag=XT --assigned-status-tag=XS --per-cell --wide-format-cell-counts -I  "+outfile \
        +"_assigned_sorted.bam -S "+outfile+"_counts.tsv.gz")

        # index+=1

def BCfunc(bcfile):
    namedir=[]
    BCname=open(bcfile,"r+")
    for i in BCname:
        # print(i.rstrip("\n"))
        if  i.rstrip("\n") !='':
            namedir.append(i.rstrip("\n"))
    return namedir

if __name__ == '__main__':
    (options, args) = parser.parse_args()
    FQ=options.readFolder
    BCfile=options.barcode
    outdir=options.output
    gtfname=str(Drugseqini.gtfname2)
    genomedir=str(Drugseqini.genomedir2)
    featurecountthread=str(Drugseqini.featurecountthread2)
    starthread=str(Drugseqini.starthread2)
    print("Note: Pepline start")
    print("Note:the gtf file is", gtfname)
    print("Note:the genomedir file is",genomedir)
    print("Note:the featurecountthread  is",featurecountthread)
    print("Note:the starthread  is",starthread)
    # os.system("rm -rf "+outdir)

    try:
        if not options.readFolder :
            parser.error('Soryy,Folder with  data (gzipped fastq files) is not given.')
        else:
            # print(glob.glob(FQ+"/*_1*.gz"))
            if options.module=="a":
                print(glob.glob(FQ+"/*R1.f*.gz"))
                FQ1=glob.glob(FQ+"/*R1.f*.gz")[0]
                FQ2=glob.glob(FQ+"/*R2.f*.gz")[0]
                print(FQ1)
                print(FQ2)
                outname2=FQ1.split("/")[-1].split("_")[0]
                BCfilelist=BCfunc(BCfile)
                os.system("rm -rf "+outdir)
                os.mkdir(outdir)
                transformer(FRdata1=gzip.open(FQ1,"r"),FRdata2=gzip.open(FQ2,"r"),WRdata1=open(outdir+"/"+FQ+".fq","w"),WRdata1D=open(outdir+"/"+FQ+".fqD","w"),BClist=BCfilelist,Outdir=outdir)
                print("-----------------Transform done-----------------")
                print("-----------------Barcode-UMI-Count analysis start-------------")

                countUMI(outname=outname2,outdirname=outdir,gtfname=gtfname,genomedir=genomedir,featurecountthread=featurecountthread,starthread=starthread)
            elif options.module=="m":
                Indata=FQ
                
                # outdir="Datatest"+FQ.split("_")[0]
                print(Indata)
                print(outdir)
                os.system("rm -rf "+outdir)
                os.mkdir(outdir)
                temin=os.path.join(os.getcwd(),Indata)
                os.system("ln -s "+temin+" "+outdir)
                outname2="Data_"+FQ.split("_")[0]
                
                countUMI(fastqname=Indata,outname=outname2,outdirname=outdir,gtfname=gtfname,genomedir=genomedir,featurecountthread=featurecountthread,starthread=starthread)


            
    except KeyboardInterrupt:
        print("\n")
        print("Note:Interrupt..")
        print("\n")
      
