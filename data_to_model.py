'''
HMM - CpG Islands

'''
import sys
import re
import time
import math
import random


data = sys.argv[1]



pstate = 0
Model={}
stateL={}
emitL={}
Ttrans={}
Temit={}

with open(sys.argv[1]) as f:
    for line in f:
        state=line[1:-1]
        emit=line[:1]
        stateL[state] = 0
        emitL[emit] = 0

for s1 in stateL.keys():
    #Initialise the Model and the Totalcount arrays Ttrans and Temit
	Ttrans[s1] = 0
	Temit[s1] = 0
	Model[s1] = dict()
	for s2 in stateL.keys():
		Model[s1][s2] = 0
	for s3 in emitL.keys():
		Model[s1][s3] = 0
#print(Ttrans)
#print(Temit)
#print(Model)


with open(sys.argv[1]) as f:
	for line in f:
		state = line[1:-1]
		emitL[line[:1]] += 1
		Temit[state] +=1
		Model[state][line[:1]] += 1
		if pstate:
			Model[pstate][state] += 1
			Ttrans[state] += 1
		pstate = state

print(Temit)
print(emitL)
print ("#Transitions")
for s1 in stateL.keys():
    for s2 in stateL.keys():
            sys.stdout.write ("%s::%s::%.3f\n"%(s1,s2, Model[s1][s2]/float(Ttrans[s1])))

print ("#Emission")
for s1 in (stateL.keys()):
	print('{} Dice: '.format(s1))
	for s2 in (emitL.keys()):
		sys.stdout.write ("%s::%s::%.3f\n"%(s1,s2,Model[s1][s2]/float(Temit[s1])))
