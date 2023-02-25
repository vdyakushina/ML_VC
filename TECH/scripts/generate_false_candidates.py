#!/usr/bin/env python3

import random
import pandas as pd
from subprocess import run

samples=[]
with open('/home/gkhvorykh/samples/meta.csv', 'r') as in_s:
  for line in in_s:
    samples.append(line.split(',')[0])

N=160

VAF={'0-0.01':1, '0.01-0.05':2, '0.05-0.15':1, '0.15-1':1}
SB={'0-0.3':1, '0.3-0.8':1, '0.8-1.1':2}

v=sum(VAF.values())*sum(SB.values())

## SNV stratification

for freq in VAF.keys():
  for strd in SB.keys():
    var_n=(N/v)*(SB[strd]*VAF[freq])
    vars=[]
    samples_s=[]
    samples_s+=samples
    freq_bottom=float(freq.split('-')[0])
    freq_top=float(freq.split('-')[1])
    strd_bottom=float(strd.split('-')[0])
    strd_top=float(strd.split('-')[1])
    i=0
    while i < var_n and len(samples_s)>0:  
      sample=random.choices(samples_s, k=1)[0]
      print(sample)
      samples_s.remove(sample)
      tsv=pd.read_csv('/home/gkhvorykh/samples/%s/analysis/variant.raw.SNV.tsv' % sample, sep='\t')
      tsv=pd.concat([tsv[tsv['Mutect2_TO_FILTER']< 40.0], tsv[tsv['Mutect2_TO_FILTER'].isnull()]], 0)
      tsv=tsv[(tsv['SGA_VAF']>freq_bottom) & (tsv['SGA_VAF']<=freq_top) & (tsv['SGA_SB']>=strd_bottom) & (tsv['SGA_SB']< strd_top)]
      if len(tsv) > 0:
        select=random.randint(0,len(tsv)-1)
        variant=tsv.iloc[select]['variant_name']
        print(variant)
        if variant not in vars:
          nsamples=run("grep \'%s\' /home/gkhvorykh/samples/*/analysis/variant.raw.SNV.tsv " "|" "wc -l" % variant, shell=True, capture_output=True, text=True).stdout.strip()
          with open('/home/gkhvorykh/train_false/false_SNV.tsv', 'a') as result:
            result.write('\t'.join([freq, strd, variant, sample, nsamples, '\t'.join(list(tsv.iloc[select][1:].astype(str)))])+'\n')
          i+=1
          vars.append(variant)
        else:
          continue
      else:
        fail='\t'.join([freq, strd, sample])
        with open('/home/gkhvorykh/train_false/failed_SNV.tsv', 'a') as result:
          result.write(fail+'\n')

run('sed -i 1i\'frequncy\tSB\tvariant\tsample\tNumber_positive_samples\tMutect2_TO_FILTER\tMutect2_TO_STRQ\tMutect2_TO_TLOD\tMutect2_TNP_FILTER\tStrelka_TO_FILTER\tStrelka_TO_QUAL\tStrelka_TNP_FILTER\tSC_FILTER\tSC_PENALTY\tSINVICT\tVD_FILTER\tVD_Q\tSGA_VAF\tSGA_SB\tSGA_RepeatUnit\tSGA_RepeatRefCount\tSGA_DP\' /home/gkhvorykh/train_false/false_SNV.tsv', shell=True)

## SNVs with VAF=0.0

freq=0.0
var_n=20
vars=[]
samples_s=[]
samples_s+=samples
i=0
while i < var_n and len(samples_s)>0:  
  sample=random.choices(samples_s, k=1)[0]
  print(sample)
  samples_s.remove(sample)
  tsv=pd.read_csv('/home/gkhvorykh/samples/%s/analysis/variant.raw.SNV.tsv' % sample, sep='\t')
  tsv=pd.concat([tsv[tsv['Mutect2_TO_FILTER']< 40.0], tsv[tsv['Mutect2_TO_FILTER'].isnull()]], 0)
  tsv=tsv[tsv['SGA_VAF']==0.0]
  if len(tsv) > 0:
    select=random.randint(0,len(tsv)-1)
    variant=tsv.iloc[select]['variant_name']
    print(variant)
    if variant not in vars:
      nsamples=run("grep \'%s\' /home/gkhvorykh/samples/*/analysis/variant.raw.SNV.tsv " "|" "wc -l" % variant, shell=True, capture_output=True, text=True).stdout.strip()
      with open('/home/gkhvorykh/train_false/false_SNV_VAF0.tsv', 'a') as result:
        result.write('\t'.join([str(freq), '-', variant, sample, nsamples, '\t'.join(list(tsv.iloc[select][1:].astype(str)))])+'\n')
      i+=1
      vars.append(variant)
    else:
      continue
  else:
    with open('/home/gkhvorykh/train_false/failed_SNV_VAF0.tsv', 'a') as result:
      result.write(sample+'\n')

run('sed -i 1i\'frequncy\tSB\tvariant\tsample\tNumber_positive_samples\tMutect2_TO_FILTER\tMutect2_TO_STRQ\tMutect2_TO_TLOD\tMutect2_TNP_FILTER\tStrelka_TO_FILTER\tStrelka_TO_QUAL\tStrelka_TNP_FILTER\tSC_FILTER\tSC_PENALTY\tSINVICT\tVD_FILTER\tVD_Q\tSGA_VAF\tSGA_SB\tSGA_RepeatUnit\tSGA_RepeatRefCount\tSGA_DP\' /home/gkhvorykh/train_false/false_SNV_VAF0.tsv', shell=True)

## SNVs with MutectFilter=40.0-100.0

MutF_bottom=40
MutF_top=100
var_n=20
vars=[]
samples_s=[]
samples_s+=samples
i=0
while i < var_n and len(samples_s)>0:  
  sample=random.choices(samples_s, k=1)[0]
  print(sample)
  samples_s.remove(sample)
  tsv=pd.read_csv('/home/gkhvorykh/samples/%s/analysis/variant.raw.SNV.tsv' % sample, sep='\t')
  tsv=tsv[(tsv['Mutect2_TO_FILTER']>=MutF_bottom) & (tsv['Mutect2_TO_FILTER']<=MutF_top)]
  if len(tsv) > 0:
    select=random.randint(0,len(tsv)-1)
    variant=tsv.iloc[select]['variant_name']
    print(variant)
    if variant not in vars:
      nsamples=run("grep \'%s\' /home/gkhvorykh/samples/*/analysis/variant.raw.SNV.tsv " "|" "wc -l" % variant, shell=True, capture_output=True, text=True).stdout.strip()
      with open('/home/gkhvorykh/train_false/false_SNV_MutF_40_100.tsv', 'a') as result:
        result.write('\t'.join(['-', '-', variant, sample, nsamples, '\t'.join(list(tsv.iloc[select][1:].astype(str)))])+'\n')
      i+=1
      vars.append(variant)
    else:
      continue
  else:
    with open('/home/gkhvorykh/train_false/failed_SNV_MutF_40_100.tsv', 'a') as result:
      result.write(sample+'\n')

run('sed -i 1i\'frequncy\tSB\tvariant\tsample\tNumber_positive_samples\tMutect2_TO_FILTER\tMutect2_TO_STRQ\tMutect2_TO_TLOD\tMutect2_TNP_FILTER\tStrelka_TO_FILTER\tStrelka_TO_QUAL\tStrelka_TNP_FILTER\tSC_FILTER\tSC_PENALTY\tSINVICT\tVD_FILTER\tVD_Q\tSGA_VAF\tSGA_SB\tSGA_RepeatUnit\tSGA_RepeatRefCount\tSGA_DP\' /home/gkhvorykh/train_false/false_SNV_MutF_40_100.tsv', shell=True)

## INDEL stratification

for freq in VAF.keys():
  for strd in SB.keys():
    var_n=(N/v)*(SB[strd]*VAF[freq])
    vars=[]
    samples_s=[]
    samples_s+=samples
    freq_bottom=float(freq.split('-')[0])
    freq_top=float(freq.split('-')[1])
    strd_bottom=float(strd.split('-')[0])
    strd_top=float(strd.split('-')[1])
    i=0
    while i < var_n and len(samples_s)>0:  
      sample=random.choices(samples_s, k=1)[0]
      print(sample)
      samples_s.remove(sample)
      tsv=pd.read_csv('/home/gkhvorykh/samples/%s/analysis/variant.raw.INDEL.tsv' % sample, sep='\t')
      tsv=pd.concat([tsv[tsv['Mutect2_TO_FILTER']< 40.0], tsv[tsv['Mutect2_TO_FILTER'].isnull()]], 0)     
      tsv=tsv[(tsv['SGA_VAF']>freq_bottom) & (tsv['SGA_VAF']<=freq_top) & (tsv['SGA_SB']>=strd_bottom) & (tsv['SGA_SB']< strd_top)]
      if len(tsv) > 0:
        select=random.randint(0,len(tsv)-1)
        variant=tsv.iloc[select]['variant_name']
        print(variant)
        if variant not in vars:
          nsamples=run("grep \'%s\' /home/gkhvorykh/samples/*/analysis/variant.raw.INDEL.tsv " "|" "wc -l" % variant, shell=True, capture_output=True, text=True).stdout.strip()
          with open('/home/gkhvorykh/train_false/false_INDEL.tsv', 'a') as result:
            result.write('\t'.join([freq, strd, variant, sample, nsamples, '\t'.join(list(tsv.iloc[select][1:].astype(str)))])+'\n')
          i+=1
          vars.append(variant)
        else:
          continue
      else:
        fail='\t'.join([freq, strd, sample])
        with open('/home/gkhvorykh/train_false/failed_INDEL.tsv', 'a') as result:
          result.write(fail+'\n')

run('sed -i 1i\'frequncy\tSB\tvariant\tsample\tNumber_positive_samples\tMutect2_TO_FILTER\tMutect2_TO_STRQ\tMutect2_TO_TLOD\tMutect2_TNP_FILTER\tStrelka_TO_FILTER\tStrelka_TO_QUAL\tStrelka_TNP_FILTER\tSC_FILTER\tSC_PENALTY\tSINVICT\tVD_FILTER\tVD_Q\tSGA_VAF\tSGA_SB\tSGA_RepeatUnit\tSGA_RepeatRefCount\tSGA_DP\' /home/gkhvorykh/train_false/false_INDEL.tsv', shell=True)

## INDEL with VAF=0.0

freq=0.0
var_n=20
vars=[]
samples_s=[]
samples_s+=samples
i=0
while i < var_n and len(samples_s)>0:  
  sample=random.choices(samples_s, k=1)[0]
  print(sample)
  samples_s.remove(sample)
  tsv=pd.read_csv('/home/gkhvorykh/samples/%s/analysis/variant.raw.INDEL.tsv' % sample, sep='\t')
  tsv=pd.concat([tsv[tsv['Mutect2_TO_FILTER']< 40.0], tsv[tsv['Mutect2_TO_FILTER'].isnull()]], 0)
  tsv=tsv[tsv['SGA_VAF']==0.0]
  if len(tsv) > 0:
    select=random.randint(0,len(tsv)-1)
    variant=tsv.iloc[select]['variant_name']
    print(variant)
    if variant not in vars:
      nsamples=run("grep \'%s\' /home/gkhvorykh/samples/*/analysis/variant.raw.INDEL.tsv " "|" "wc -l" % variant, shell=True, capture_output=True, text=True).stdout.strip()
      with open('/home/gkhvorykh/train_false/false_INDEL_VAF0.tsv', 'a') as result:
        result.write('\t'.join([str(freq), '-', variant, sample, nsamples, '\t'.join(list(tsv.iloc[select][1:].astype(str)))])+'\n')
      i+=1
      vars.append(variant)
    else:
      continue
  else:
    with open('/home/gkhvorykh/train_false/failed_INDEL_VAF0.tsv', 'a') as result:
      result.write(sample+'\n')

run('sed -i 1i\'frequncy\tSB\tvariant\tsample\tNumber_positive_samples\tMutect2_TO_FILTER\tMutect2_TO_STRQ\tMutect2_TO_TLOD\tMutect2_TNP_FILTER\tStrelka_TO_FILTER\tStrelka_TO_QUAL\tStrelka_TNP_FILTER\tSC_FILTER\tSC_PENALTY\tSINVICT\tVD_FILTER\tVD_Q\tSGA_VAF\tSGA_SB\tSGA_RepeatUnit\tSGA_RepeatRefCount\tSGA_DP\' /home/gkhvorykh/train_false/false_INDEL_VAF0.tsv', shell=True)

## INDEL with MutectFilter=40.0-100.0

MutF_bottom=40
MutF_top=100
var_n=20
vars=[]
samples_s=[]
samples_s+=samples
i=0
while i < var_n and len(samples_s)>0:  
  sample=random.choices(samples_s, k=1)[0]
  print(sample)
  samples_s.remove(sample)
  tsv=pd.read_csv('/home/gkhvorykh/samples/%s/analysis/variant.raw.INDEL.tsv' % sample, sep='\t')
  tsv=tsv[(tsv['Mutect2_TO_FILTER']>=MutF_bottom) & (tsv['Mutect2_TO_FILTER']<=MutF_top)]
  if len(tsv) > 0:
    select=random.randint(0,len(tsv)-1)
    variant=tsv.iloc[select]['variant_name']
    print(variant)
    if variant not in vars:
      nsamples=run("grep \'%s\' /home/gkhvorykh/samples/*/analysis/variant.raw.INDEL.tsv " "|" "wc -l" % variant, shell=True, capture_output=True, text=True).stdout.strip()
      with open('/home/gkhvorykh/train_false/false_INDEL_MutF_40_100.tsv', 'a') as result:
        result.write('\t'.join(['-', '-', variant, sample, nsamples, '\t'.join(list(tsv.iloc[select][1:].astype(str)))])+'\n')
      i+=1
      vars.append(variant)
    else:
      continue
  else:
    with open('/home/gkhvorykh/train_false/failed_INDEL_MutF_40_100.tsv', 'a') as result:
      result.write(sample+'\n')

run('sed -i 1i\'frequncy\tSB\tvariant\tsample\tNumber_positive_samples\tMutect2_TO_FILTER\tMutect2_TO_STRQ\tMutect2_TO_TLOD\tMutect2_TNP_FILTER\tStrelka_TO_FILTER\tStrelka_TO_QUAL\tStrelka_TNP_FILTER\tSC_FILTER\tSC_PENALTY\tSINVICT\tVD_FILTER\tVD_Q\tSGA_VAF\tSGA_SB\tSGA_RepeatUnit\tSGA_RepeatRefCount\tSGA_DP\' /home/gkhvorykh/train_false/false_INDEL_MutF_40_100.tsv', shell=True)
