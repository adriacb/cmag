'''
HMM - CpG Islands

This generates a Data given a transitions and emission probabilities.

Usage: $ model_to_data.py model.txt

'''
import sys
import re
import time
import math
import random

file_model = sys.argv[1] #'/model.txt'


def readmodel (name):
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

def state2state (Model, cstate):
    f=random.random()
    t=0
    States=['A+','C+','G+','T+', 'A-','C-', 'G-', 'T-']
    #print(States)
    for i in (States):
        t += Model[cstate][i]
        if (t>=f):
            return i

def state2emission (Model, cstate):
    Emits=['A','C','G','T']
    f = random.random()
    t = 0
    for i in Emits:
       #print(Model[cstate])
        t += Model[cstate][i]
        if t >= f:
            return i


if __name__ == '__main__':
    M, S, E= readmodel(file_model)
    cstate='A+'
    f = open('data_test.txt', "w+")
    for i in range (0, 10000):
        v=state2emission(M, cstate)
        print("{}{}".format(v, cstate))
        f.write("{}{}\n".format(v, cstate))
        cstate=state2state(M, cstate)
    f.close()
