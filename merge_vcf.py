#!/usr/bin/env python
# -*- coding: ASCII -*-

###This program merges individual emit all vcf files created by the aDNA_GenoCaller and GenoCaller_indent scripts for biallelic SNPs
###pysam needs to be installed. Two input files are needed.
##The first input file needs to be a tab delimited bed file containing the list of SNPs to be included in the output vcf, one per line, of the form <chrom><start><end><all1><all2><snp_name>
##e.g. 2\t136608645\t136608646\tG\tA\t2:136608646\n
##One of two alleles should be the reference allele, though the order does not matter
##The second input file contains a list of individual emit_all vcf files to be merged, followed by the assigned sample name, tab seperated.
##e.g.
##../2072_LCT.RG.LCT.bed.aDNA.emit_all.vcf	samp_2072
##../84001_LCT.RG.LCT.bed.aDNA.emit_all.vcf	samp_84001
##../84005_LCT.RG.LCT.bed.aDNA.emit_all.vcf	samp_84005
###to run type:
##merge_vcf.py <SNP_bedfile> list_of_emit_alls> <reference genome>
##Output will be a merged vcf.gz file.

import gzip
import string
import pysam
import numpy as np
import math
import gzip
from sys import argv


snp_bedfile=argv[1]#'1000G_EURASIAN_wYRI_1.2m.bedfile'#argv[1]
samp_list_file=argv[2]#'AllMed_1000G'#argv[2]
ref_file=argv[3]#'/vault/public/GRCh37/GRCh37.fa'

ref=pysam.FastaFile(ref_file)

def phred2prob(x):
    return 10.0**(-x/10.0)

def prob2phred(x):
    return -10*math.log10(x)

PL_snp_dic={}

#AA,AC,AG,AT,CC,CG,CT,GG,GT,TT

PL_snp_dic['AC']=np.array([0,1,4])
PL_snp_dic['AG']=np.array([0,2,7])
PL_snp_dic['AT']=np.array([0,3,9])
PL_snp_dic['CA']=np.array([4,1,0])
PL_snp_dic['CG']=np.array([4,5,7])
PL_snp_dic['CT']=np.array([4,6,9])
PL_snp_dic['GA']=np.array([7,2,0])
PL_snp_dic['GC']=np.array([7,5,4])
PL_snp_dic['GT']=np.array([7,8,9])
PL_snp_dic['TA']=np.array([9,3,0])
PL_snp_dic['TC']=np.array([9,6,4])
PL_snp_dic['TG']=np.array([9,8,7])

geno_list=['0/0','0/1','1/1']

file=open(samp_list_file,'r')
data=file.read()
data=string.split(data,'\n')
if data[-1]=='':
    del(data[-1])

files_use=[]
samps=[]

for g in range(len(data)):
    file_use,samp=string.split(data[g],'\t')
    files_use.append(file_use)
    samps.append(samp)

file=open(snp_bedfile,'r')

data=file.read()
data=string.split(data,'\n')
if data[-1]=='':
    del(data[-1])

snps={}
for g in range(len(data)):
    k=string.split(data[g],'\t')
    key=k[0]+'_'+k[2]
    all1=k[3]
    all2=k[4]
    if (all1<>'0') and (all2<>'0'):
        ref_all=ref.fetch(k[0],int(k[1]),int(k[2]))
        if all1==ref_all:
            snps[key]=[k[0],k[2],[all1,all2],k[5]]
        elif all2==ref_all:
            snps[key]=[k[0],k[2],[all2,all1],k[5]]

fileout=gzip.open(samp_list_file+'.vcf.gz','w')


###set up output headers
header='##fileformat=VCFv4.1\n'
for i in range(len(ref.lengths)):
    header=header+'##contig=<ID='+ref.references[i]+',length='+str(ref.lengths[i])+'>\n'
header=header+'##reference=file:'+ref_file+'\n'

header=header+'#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT'

for g in range(len(samps)):
    header=header+'\t'+samps[g]

header=header+'\n'

fileout.write(header)


file_read={}
for g in range(len(files_use)):
    file_read[g]=open(files_use[g],'r')


X=[]
Y=[]


for g in range(len(file_read)):
    x=file_read[g].readline()
    while x[0]=='#':
        x=file_read[g].readline()
    X.append(x)
    Y.append(string.split(x[:-1],'\t'))    
    x='#'

count=0


wrong_ver=len(files_use)*[0]
multi_ver=len(files_use)*[0]

while X[0]<>'':
    if X[0][0]<>'#':
        count+=1
        Y_t=zip(*Y)
        if Y_t[1].count(Y_t[1][0])==len(samps):

            chrom=Y_t[0][0].replace('chr','')
            pos=Y_t[1][0]
            key=chrom+'_'+Y[0][1]
            
            refs=list(set(Y_t[3]))
            alts_list=[]
            for g in range(len(Y_t[4])):
                alts_list.extend(string.split(Y_t[4][g],','))
            alts=list(set(alts_list))
            alts.sort()
            try:
                refs.remove('.')
            except:
                clean=1
            try:
                alts.remove('.')
            except:
                clean=1
            
            formats=Y_t[8]
            vals=Y_t[9]
            GTs=[]
            DPs=[]
            ADs=[]
            GQs=[]
            PL10s=[]
            VQs=[]
            VQs=list(Y_t[5])
            for g in range(len(formats)):
                if formats[g]=='.':
                    GTs.append('./.')
                    DPs.append('.')
                    ADs.append('.')
                    GQs.append('.')
                    PL10s.append('.')
                    
                else:
                    format_ls=string.split(formats[g],':')
                    vals_ls=string.split(vals[g],':')
                    if 'GT' in format_ls:
                        GT=vals_ls[format_ls.index('GT')]
                        if GT=='0/0':
                            GTs.append(GT)
                        else:
                            if len(Y_t[4][g])==1:
                                GT=GT.replace('1',str(alts.index(Y_t[4][g])+1))                          
                            else:
                                alt_alls=string.split(Y_t[4][g],',')
                                GT=GT.replace('1',str(alts.index(alt_alls[0])+1))
                                GT=GT.replace('2',str(alts.index(alt_alls[1])+1)) 
                            GTs.append(GT)
                    else:
                        GTs.append('./.')

                    if 'DP' in format_ls:
                        DP=vals_ls[format_ls.index('DP')]
                        DPs.append(DP)
                    else:
                        DPs.append('.')
                    
                    alt_DP=['0']*(len(alts)+1)
                    if 'AD' in format_ls:
                        AD=vals_ls[format_ls.index('AD')]
                        AD_sep=string.split(AD,',')
                        alt_DP[0]=AD_sep[0]
                        if len(Y_t[4][g])==1:
                            alt_DP[alts.index(Y_t[4][g])+1]=AD_sep[1]
                        else:
                            alt_alls=string.split(Y_t[4][g],',')
                            
                            alt_DP[alts.index(alt_alls[0])+1]=AD_sep[1]
                            alt_DP[alts.index(alt_alls[1])+1]=AD_sep[2]                            
                        ADs.append(string.join(alt_DP,','))

                    else:
                        if 'DP' in format_ls:
                            alt_DP[0]=DP
                            ADs.append(string.join(alt_DP,','))
                        else:
                            ADs.append('.')

                    if 'GQ' in format_ls:
                        GQ=vals_ls[format_ls.index('GQ')]
                        GQs.append(GQ)
                    else:
                        GQs.append('.')

                    if 'PL' in format_ls:
                        PL=vals_ls[format_ls.index('PL')]
                        PL10s.append(PL)
                    else:
                        PL10s.append('.')


            DP2=np.zeros((len(DPs)),dtype='int32')
            for g in range(len(DPs)):
                if DPs[g]=='.':
                    DP2[g]=0
                else:
                    DP2[g]=int(DPs[g])

            depth_string='DP_mean='+str(int(np.mean(DP2)))+';DP_min='+str(np.min(DP2))+';DP_max='+str(np.max(DP2))

            if len(refs)>0:                       
                fields=[GTs,DPs,ADs,GQs,PL10s,VQs]
                fields_T=zip(*fields)
                if len(alts)==0:
                    out=chrom+'\t'+pos+'\t'+key+'\t'+refs[0]+'\t.\t.\t.\t'+depth_string+'\tGT:DP:AD:GQ:PL:VQ'
                    for g in range(len(fields_T)):
                        out=out+'\t'+string.join(fields_T[g],':')
                if len(alts)>0:
                    out=chrom+'\t'+pos+'\t'+key+'\t'+refs[0]+'\t'+string.join(alts,',')+'\t.\t.\t'+depth_string+'\tGT:DP:AD:GQ:PL:VQ'
                    for g in range(len(fields_T)):
                        out=out+'\t'+string.join(fields_T[g],':')              
            
                if snps.has_key(key):
                    snp_code=snps[key][2][0]+snps[key][2][1]
                    PLs=[]
                    ADs2=[]
                    DPs3=[]
                    for g in range(len(PL10s)):
                        AD_targ=string.split(ADs[g],',')
                        if PL10s[g]=='.':
                            PLs.append('.')
                            DPs3.append('.')
                            ADs2.append('.')
                        else:
                            AD_cat=AD_targ[0]
                            DP_tot=int(AD_targ[0])
                            try:
                                DP_tot+=int(AD_targ[alts.index(snps[key][2][1])+1])
                                AD_cat=AD_cat+','+str(int(AD_targ[alts.index(snps[key][2][1])+1]))
                            except:
                                AD_cat=AD_cat+',0'

                                
                            if DP_tot==0:
                                PLs.append('.')
                                DPs3.append('.')
                                ADs2.append('.')
                                wrong_ver[g]+=1
                                #print count,wrong_ver
                            else:
                                PL=list(np.array(string.split(PL10s[g],','))[PL_snp_dic[snp_code]])
                                if '0' not in PL:
                                    prob=np.zeros(3)
                                    for gg in range(len(PL)):
                                        prob[gg]=phred2prob(float(PL[gg]))
                                    PL_new=-10*np.log10(prob/np.max(prob))
                                    PL_new2=np.nan_to_num(PL_new)
                                    PL=[str(int(PL_new2[0])),str(int(PL_new2[1])),str(int(PL_new2[2]))]
                                if PL.count('0')>1:
                                    multi_ver[g]+=1
                                    PLs.append('.')
                                    DPs3.append('.')
                                    ADs2.append('.')                                    
                                else:
                                    PLs.append(PL)
                                    DPs3.append(str(DP_tot))
                                    ADs2.append(str(AD_cat))

                    GT3s=[]
                    GQ2s=[]
                    for g in range(len(PLs)):
                        if PLs[g]=='.':
                            GT3s.append('./.')
                            GQ2s.append('.')
                        else:

                            GT3s.append(geno_list[PLs[g].index('0')])
                            GQ=[]
                            for gg in range(len(PLs[g])):
                                GQ.append(int(PLs[g][gg]))
                            GQ.sort()
                            GQ2s.append(str(GQ[1]))
                            PLs[g]=string.join(PLs[g],',')



                    fields2=[GT3s,DPs3,ADs2,GQ2s,PLs]
                    fields2_T=zip(*fields2)
                    
                    out2=chrom+'\t'+pos+'\t'+snps[key][3]+'\t'+snps[key][2][0]+'\t'+snps[key][2][1]+'\t.\t.\t'+depth_string+'\tGT:DP:AD:GQ:PL'
                    for g in range(len(fields2_T)):
                        out2=out2+'\t'+string.join(fields2_T[g],':')
                    fileout.write(out2+'\n')
                
        if (count%10000)==0:
            print count,chrom,pos
##            print wrong_ver
##            print multi_ver

    X=[]
    Y=[]
    for g in range(len(file_read)):
        x=file_read[g].readline()
        X.append(x)
        Y.append(string.split(x[:-1],'\t'))

    
fileout.close()


                      
