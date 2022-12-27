#! /usr/bin/python
import numpy as np

def revcomp(seq):
    sdict={'A':'T','T':'A','G':'C','C':'G'}
    newseq = [sdict[s] for s in seq]
    return ''.join(newseq[::-1])

def cbccorrect(cbcs,qcbcs,cbcfreq_dict):
    newcbcs=[]
    if len(cbcs)>0:
        qscores = ['!','"','#','$','%','&','\'','(',')','*','+',',','-','.','/','0','1','2','3','4','5','6','7','8','9',':',';','<','=','>','?','@','A','B','C','D','E','F','G','H','I','J','K','L','M','N']
        qdict = {q:float(i) for i,q in enumerate(qscores)}
        nts=['A','G','T','C']
        cbclen=len(cbcs[0])
        corrected=0
        for cbc,qcbc in zip(cbcs,qcbcs):
            if qcbc=='0':
                newcbcs.append(cbc)
            else:
                T=0
                k=0
                tcbcs=[]
                Llist=np.zeros(4*cbclen)
                for i in range(cbclen):
                    pedit=max([0.0005,pow(10.0,-qdict[qcbc[i]]/10.0)])
                    for nt in nts:
                        tcbc=''.join([n if j!=i else nt for j,n in enumerate(cbc)])
                        tcbcs.append(tcbc)
                        if tcbc in cbcfreq_dict:
                            L = cbcfreq_dict[tcbc]*pedit
                            Llist[k]=L
                            T+=L
                        k+=1
                if T>0:
                    Llist=np.array(Llist)/T
                    if np.max(Llist)>0.975:
                        newcbcs.append(tcbcs[np.argmax(Llist)])
                        corrected+=1
                    else:
                        newcbcs.append(cbc)
                else:
                    newcbcs.append(cbc)
    return newcbcs
