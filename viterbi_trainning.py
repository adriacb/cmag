'''
HMM - CpG Islands

Viterbi trainning Model

'''
import sys
import re
import time
import math
import random

#sfasta_file = sys.argv[1] #'/Users/Adria/Desktop/cpg/brca1.fasta'
data = sys.argv[1]
model = sys.argv[2] #import model


LOG_ZERO=-999999999
SMALL=0.00000001

def log_multiply (v1, v2):
    if v1==LOG_ZERO or v2 ==LOG_ZERO:
        return LOG_ZERO
    return v1+v2


def readmodel(name):
    Model={}
    with open(name) as f:
        for line in f:
            mobj=re.match ( r'(\S*)::(\S*)::(\S*)', line)
            if mobj:
                s1=mobj.group(1)
                s2=mobj.group(2)
                p=float (mobj.group(3))
                if s1 not in Model:
                    Model[s1]={}

                Model[s1][s2]=p

        States=Model.keys()


        Emits=[]
        EmitsH={}

        for k in (States):
            k1=Model[k].keys()
            for l in (k1):
                if l not in Model and l not in EmitsH:
                    Emits.append(l)
                    EmitsH[l]=1

    return (Model, States, Emits)

def readdata (data):
    Data={}
    Ref={}
    i=0
    with open(data) as f:
        for line in f:
            Data[i]=line[:1]
            Ref[i]=line[1:-1]

            i+=1
    return (Data, Ref)

def viterbi (Data, Model):
    maxk=LOG_ZERO
    ptr_k=""
    L=len(Data)
    States=Model.keys()
    Ns=len (States)
    V   ={}
    PTR ={}
    for i in range (0, L+1):
        V[i]={}
        PTR[i]={}
    for k in (States):
        V[0][k]=0
        #print(V[0])
    for i in range (1,L+1):
        symbol=Data[i-1] 
    #print(V)
        for l in (States):
            max_k=LOG_ZERO
            ptr_k=""
            for k in (States):
                v=log_multiply(V[i-1][k],Model[k][l]);
                #print('-->',v)
                if v>max_k or max_k==LOG_ZERO:
                    max_k=v
                    ptr_k=k
            V[i][l]=log_multiply(Model[l][symbol],max_k)
            PTR[i][l]=ptr_k

    max_k=LOG_ZERO
    ptr_k=""

    for k in (States):
        vv=V[L][k]
        if vv>max_k or max_k==LOG_ZERO:
            max_k=vv
            ptr_k=k
    vd={}
    for i in range (L, 0,-1):
        vd[i-1]=ptr_k
        ptr_k=PTR[i][ptr_k]

    return (vd, max_k)

def viterbi_trainning (Data, Model, States, Emits): #restimate the parameters and update the model

    (v,score)=viterbi(Data,Model)
    L=len(Data)

    A={}#will be used to re-estimate the transitions
    E={}#will be used to re-estimate the emissions
    #print('--- 4 viterbi training ---')
    for k in (States):
        A[k]={}
        E[k]={}

        for l in (States):
            A[k][l]=0
        for l in (Emits):
            E[k][l]=0
        #print('Transitions dict {} \n Emissions dict {}'.format(A, E))

    for i in range (1,L):
        A[v[i-1]][v[i]] += 1 #This is the transition of the current state from the previous
        #print('v[i-1]',v[i-1],'v[i]', v[i])
        E[v[i]][Data[i]] += 1 #V[i] is the STATE and Data[i] is the emit
        #print('v[i]',v[i],'Data[i]', Data[i])
    for k in (States):    #just count
        for l in (States):A[k][l]+=1
        for l in (Emits) :E[k][l]+=1

    for k in (States):
        num=float(0)
        for l in (States): num+=A[k][l]
        for l in (States): Model[k][l]=math.log(A[k][l]/num)

    for k in (States):
        num=float (0)
        for l in (Emits): num+=E[k][l]
        for l in (Emits): Model[k][l]=math.log(E[k][l]/num)

    (v,score)=viterbi(Data,Model)
    return (Model,score)

def model2flat (Model, States, Emits):
    #print('--- 3 model2flat ---')
    print()
    for k in States:

        tr = 0
        for l in States:
            #print('State: ', l)
            r = random.random()
            tr += r
            Model[k][l] = r
            #print('random', Model[k][l])
        for l in States:
            Model[k][l] /= tr
            #print('/=tr',Model[k][l])

        tr = 0
        for l in Emits:
            #print('Emit: ', l)
            r = random.random()
            tr += r
            Model[k][l] = r
            #print('random', Model[k][l])
        for l in Emits:
            Model[k][l] /= tr
            #print('/tr', Model[k][l])
    return Model



if __name__ == '__main__':
    Model1, S, E = readmodel(model)
    D, R = readdata(data)
    #fasta = read_fasta(fasta_file)
    M=model2flat(Model1,S,E)
    pscore=0
    for i in range (0, 100):
        (Model,score)=viterbi_trainning (D, M,S,E)
        sys.stderr.write("iteration %d Score=%d\n"%(i,score))
        if (score==pscore):
            break
        pscore=score

    print('Writting Emits and States')
    for l in (S):

        for k in (S):
            sys.stderr.write("{}::{}::{:.5f}\n".format(l,k,math.exp(Model[l][k])))
        for k in (E):
            sys.stderr.write("{}::{}::{:.5f}\n".format(l,k,math.exp(Model[l][k])))
