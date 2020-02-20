#!/usr/bin/env python
import sys
import re
import math

'''
HMM - CpG Islands

Given a long sequence of emits (fasta sequence), try to decode the states for each position
of the sequence.

Also implemented if we want to parse a generated data (output model_to_data).

'''
class Seq:
    def __init__(self):
        self.id = 0
        self.seq=''
        self.features=''

    def count_freq_seq(self):
        first_order = dict()

        total = 0
        for aa in self.seq:
            first_order[aa] = first_order.get(aa, 0) + 1
            total += 1

        for a in first_order:
            first_order[a] /= total

        return first_order


def read_fasta(fasta_file):
    '''
    This function returns a list containing a Seq() object for
    a fasta file.
    '''
    record=[]
    nrec=-1
    inseq=0
    with open(fasta_file) as f:
        for line in f:
            if re.match ( r'^>', line):
                nrec+=1
                record.append(Seq())
                mobj=re.match ( r'^>(\S*)\s*(.*)', line)

                if (mobj):
                    record[nrec].id=mobj.group(1)
                    record[nrec].features=mobj.group(2)
                inseq=0
            else :
                if inseq==0 :
                    inseq=1
                    record[nrec].seq=line
                else:
                    cstring=record[nrec].seq+line
                    record[nrec].seq=cstring


    for x in range (0,nrec+1):
        record[x].seq=re.sub (r'[ \n\t\r]',"",record[x].seq)
        sys.stdout.write (">%s\n%s\n"%(record[x].id, record[x].seq))
    return record

LOG_ZERO=-999999999
SMALL=0.00000001

def log_multiply (v1, v2):
    if v1==LOG_ZERO or v2 ==LOG_ZERO:
        return LOG_ZERO
    return v1+v2

def viterbi (Data, Model):
  max_k=LOG_ZERO
  ptr_k=""
  L=len(Data)
  States=Model.keys()
  print(States)
  V   ={}#This will contain your Viterbi Matrix - the best scores
  PTR ={}#This will contain the traceback, THE MEMORY WHERE WE COME FROM WILL BE OUR PREDICTION
  for i in range (0, L+1):
      V[i]={} #enters V[0], v[1], ... and as value is and empty dictionary
      PTR[i]={}
      '''
      >>> print(PTR) if L==6
      {0: {}, 1: {}, 2: {}, 3: {}, 4: {}, 5: {}, 6: {}}
      '''
  for k in (States):
      V[0][k]=0

  for i in range (1,L+1):
      symbol=Data[i-1]
      for l in (States):
        max_k=LOG_ZERO
        ptr_k="" #keep tracking previous pointer
        for k in (States):
            v = log_multiply(V[i-1][k], Model[k][l])
            #V[i-1][k]
            if v>max_k or max_k==LOG_ZERO: #who is the parent of v
                #use log multiply rather than multiplying probabilities
                max_k = v #now v is is the max score
                ptr_k = k #we update the States

        V[i][l] = log_multiply(Model[l][symbol], max_k)
        PTR[i][l] = ptr_k

  max_k=LOG_ZERO
  ptr_k=""

  for k in (States):
    vv=V[L][k]
    if vv>max_k or max_k==LOG_ZERO:
        max_k=vv
        ptr_k=k
  vd={}#this will contain the list of decoded states
  for i in range (L, 0,-1):
      vd[i-1]=ptr_k
      ptr_k=PTR[i][ptr_k]

  return vd

def readmodel (name):
    Model={}
    with open(name) as f:
        for line in f:
            mobj=re.match ( r'(\S*)::(\S*)::(\S*)', line)
            if mobj:
                s1=mobj.group(1)
                s2=mobj.group(2)
                p=float (mobj.group(3))
                if (p<SMALL): p=LOG_ZERO
                else: p=math.log(p)

                if s1 not in Model:
                    Model[s1]={}

                Model[s1][s2]=p

        return Model

def readdata (name):
    Data={}
    Ref={}
    i=0
    with open(name) as f:
        for line in f:
            Data[i]=line[:1]
            Ref [i]=line[1:-1]

            i+=1
        return (Data,Ref)

if __name__ == '__main__':
    if sys.argv[1].endswith('.txt'):
        (Data,Ref) = readdata (sys.argv[1])
        Model = readmodel (sys.argv[2])
        ViterbiD = viterbi (Data, Model)

        L=len(Data)
        for i in range (0,L):
          #we will print the original state as provided in data along with the prediction
          sys.stdout.write("{}{}{}\n".format(Data[i], Ref[i], ViterbiD[i]))
    if sys.argv[1].endswith('.fasta'):
        Data = read_fasta (sys.argv[1])
        Model = readmodel (sys.argv[2])
        ViterbiD = viterbi (Data[0].seq, Model)

        L=len(Data[0].seq)
        for i in range (0,L):
          #we will print the original state as provided in data along with the prediction
          sys.stdout.write("{}{}\n".format(Data[0].seq[i], ViterbiD[i]))
