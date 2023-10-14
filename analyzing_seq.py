"""In this project I am going to analyze RNA sequences fro mthe covid veriants
I am going to use the biopython package to do so"""

import pandas as pd
import numpy as np
sequences=pd.read_csv('C:\myPythonProjects\covid RNA sequences analyze\\ncbi_datasets.csv')
print(sequences.info())
print(sequences.shape)
#lets remover US State and Host Name columns since they are empty

sequences.drop(['US State','Host Name'],axis=1,inplace=True)
print(sequences.info())
#now we dont have any null columns

#lets convert the dates into date format 

sequences['Collection Date']=pd.to_datetime(sequences['Collection Date'],format='mixed')
print(sequences.info())
#now it is a date time object
#lets answers some questions
# 1. When did the first case in each continent occur?
sequences['Geo Location']=sequences['Geo Location'].astype(str)
sequences['Geo Location']=sequences['Geo Location'].str.split()

def first_word(sen):
    if sen[0]=='North'or sen[0]=='South':
        return sen[0]+ " " +sen[1]
    return sen[0]
sequences['continent']=sequences['Geo Location'].apply(first_word)
continents=pd.pivot_table(data=sequences,values=['Nucleotide Accession','Collection Date'],index='continent',aggfunc=np.min)
print(continents)
"""we have a list of the first sequence collected in each continent. The first sample was collected in China
so this checks out. Also the last continent the the virus hit was Africa"""

#2.how many sequences collected in each continent?
print(sequences['continent'].value_counts())

"""The continent that the most amount of sequences were collected from It is North America with over half a million sequences
And the continent that the least amount of seqences were collected is South America with just 755"""
#3. what are the shortest and longest sequences that were collected?
#first lets look at some bins so make sure there are'nt any outliers

print(sequences['Nucleotide Length'].value_counts(bins=15))
#all of the values lay somewhere between 28200 and 30000. there is one outlier that has length lestg than 4700, lets delete it

sequences.drop(sequences[sequences['Nucleotide Length']<5000].index,inplace=True)
print(sequences['Nucleotide Length'].value_counts(bins=10))#we got rid of the outlier


print(sequences['Nucleotide Length'].max())#30018
print(sequences['Nucleotide Length'].min())#27366
#4. how many samples werw taked each month? lets move to tableau to make this much easier

"""The trends are clear, when the Delta has spread in the late 2020 and early 2021 there was a huge surge
in sample collection. the same trend is true for the Omicron, when there was the outbreak in late 2021,
there was a hugh surge in sample collctions"""
#now lets pull some accual sequences- the reference sequence, the base sequence and one of delta's and omicron's sequences

print(sequences[sequences['Sequence Type']=='RefSeq'])#reference sequence
print(sequences[sequences['Isolate Name'].str.contains('Delta').fillna(False)])#delta
print(sequences[sequences['Isolate Name'].str.contains('Omicron').fillna(False)])#Omicron
ids=['NC_045512.2','OM061695.1','OM095411.1','MN985325.1']#now we have the list of the mentioned sequences ids
names=['reference','delta','omicron','base']
#lets get the metadate for those ids and then download them
selected_ids=sequences[sequences['Nucleotide Accession'].isin(ids)]
print(selected_ids)#now we have a list of the 4 base sequences
#lets download the sequences using the Entrez module from bio python
from Bio import Entrez
Entrez.email='amitradin888@gmail.com'

def download_seq(id1):
    handle=Entrez.esearch(db='nucleotide',term=id1,retmax='1')#searching the data base
    record=Entrez.read(handle)
    handle2=Entrez.efetch(db='nucleotide',id=record['IdList'][0],rettype='fasta',retmode='text')#second search, it is required with this module
    return handle2.read()#the result
#lets get all of the sequences
data={}
for id in ids:
    data[id]={'fasta':download_seq(id)}
print(data)
#lets extract the actual sequences to make it workable
from Bio import SeqIO
import io
actual={}
for key,value in data.items():
    fake_file=io.StringIO(value['fasta'])#creating a fake file to fool biopython
    actual[key]=list(SeqIO.parse(fake_file,'fasta'))[0]

print(actual)#we have the actual sequences now! we can start analyzing them
"""Now we want to see how to virus mutated over time.
We have the reference code, and we can compare it to the veriants"""
from Bio import Align
aligner=Align.PairwiseAligner()
aligner.algorithm #we have out algorithm, now lets compare some pairs

#first lets compare the refrence and the omicron
score=aligner.score(actual['NC_045512.2'].seq,actual['OM095411.1'].seq)

print(score/len(actual['NC_045512.2'].seq))#0.994 


#now the reference and the delta
score=aligner.score(actual['NC_045512.2'].seq,actual['OM061695.1'].seq)
print(score/len(actual['NC_045512.2'].seq))#0.997
#finally delta vs omicron
score=aligner.score(actual['OM095411.1'].seq,actual['OM061695.1'].seq)
print(score/len(actual['OM095411.1'].seq))#0.996
"""We can conclude that the veriants are actually very similer in terms of the RNA 
sequence itself, and yet, the verients behave totally differently from each other which is pretty crazy"""
#now lets create a dataframe that compares each veriant with the other
matrix=np.empty([4,4])
print(matrix)
for i,first in enumerate(ids):
    for j,second in enumerate(ids):
        score=aligner.score(actual[first].seq,actual[second].seq)
        matrix[i][j]=score/len(actual[first].seq)
print(matrix)#now we have a comparison matrix, lets convert it into a DataFrame

scores=pd.DataFrame(matrix,columns=names,index=names)
print(scores)
"""Now we have a complete view of the mutation of the virus
we can see that the reference and the base are almost identical, because the 
first case in China and the US is very close in terms of dates
As the dates go by the virus changes more and more. and the difference 
between the reference and the omicron veriant is the largets with a 6 percent 
change"""
"""Now I want so see the exact mutation points
We can do this yet again with the align module"""

seq1=actual['NC_045512.2'].seq #reference
seq2=actual['OM061695.1'].seq #delta
delta_alignments=aligner.align(seq1,seq2)

delta_alignment = delta_alignments[0]
print(delta_alignment.shape)
print(delta_alignment.aligned)#those are the alignments, but lets make it more readable
seq_1_end=None
seq_2_end=None
for alignments in zip(delta_alignment.aligned[0],delta_alignment.aligned[1]):
    if seq_1_end and seq_2_end:
        seq1_not_aligned=seq1[seq_1_end:alignments[0][0]]
        seq2_not_aligned=seq1[seq_2_end:alignments[1][0]]
        #lets print the mismatches
        print('1: {}'.format(seq1_not_aligned))
        print('2: {}'.format(seq2_not_aligned))
    seq_1_end=alignments[0][1]
    seq_2_end=alignments[1][1]

#thats cool but it does not tell us the type of mismatch and the location of it
#lets add those things
from IPython.display import HTML

def color_print(s, color='black'):
    return "<span style='color:{}'>{}</span>".format(color, s)
seq1_end = None
seq2_end = None
display_seq = []
for alignments in zip(delta_alignment.aligned[0], delta_alignment.aligned[1]):
    
    if seq1_end and seq2_end:
        seq1_mismatch = seq1[seq1_end:alignments[0][0]]
        seq2_mismatch = seq2[seq2_end:alignments[1][0]]
        if len(seq2_mismatch)==0:
            display_seq.append(color_print(seq1[seq1_end:alignments[0][0]], "red"))#deletion
        elif len(seq1_mismatch)==0:
            display_seq.append(color_print(seq2[seq2_end:alignments[1][0]], "green"))#insertion
        else:
            display_seq.append(color_print(seq2[seq2_end:alignments[1][0]], "blue"))#substitution
    
    display_seq.append(seq1[alignments[0][0]:alignments[0][1]])#adding the next part of the sequence without color
    
    seq1_end = alignments[0][1]
    seq2_end = alignments[1][1]
#now we need to join the display sequences
display_seq = [str(i) for i in display_seq]

print(HTML('<br>'.join(display_seq)))
"""Great! now we can see exactly which kind of mismatches occured
 and where they occured. This iis really cool to see
 We can clearly see that there aren't many mismatches 
 and yet the veriants are completly different"""