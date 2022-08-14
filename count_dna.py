#! /usr/bin/python
import gzip
import io
import numpy as np
from os.path import exists

def revcomp(seq):
    sdict={'A':'T','T':'A','G':'C','C':'G'}
    newseq = [sdict[s] for s in seq]
    return ''.join(newseq[::-1])


def cbcfreq(address,bcset):
    cbcfreq_dict = {}
    cbcs = []
    qcbcs = []
    tot=0
    with io.BufferedReader(gzip.open(address,'rb')) as f:
        for line in f:
            dlist = line.decode().split()
            cbc=dlist[1]
            if cbc in bcset:
                tot+=1
                if cbc in cbcfreq_dict:
                    cbcfreq_dict[cbc]+=1
                else:
                    cbcfreq_dict[cbc]=1
                qcbcs.append('0')
            else:
                qcbcs.append(dlist[2])
            cbcs.append(cbc)
    tot=float(tot)
    for cbc in cbcfreq_dict:
        cbcfreq_dict[cbc]=float(cbcfreq_dict[cbc])/tot
    return cbcs,qcbcs,cbcfreq_dict

def cbccorrect(cbcs,qcbcs,cbcfreq_dict):
    qscores = ['!','"','#','$','%','&','\'','(',')','*','+',',','-','.','/','0','1','2','3','4','5','6','7','8','9',':',';','<','=','>','?','@','A','B','C','D','E','F','G','H','I']
    qdict = {q:float(i) for i,q in enumerate(qscores)}
    nts=['A','G','T','C']
    newcbcs = []
    cbclen=len(cbcs[0])
    print(qcbcs.count('0'),len(qcbcs))
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
    print('corrected',corrected)
    return newcbcs

def fragments(sample,ref,address,fragfile,whitelist,rc,chrs):
    if exists(whitelist):
        if rc==1:
            bcset = set([revcomp(line.split()[0]) for line in open(whitelist)])
        else:
            bcset = set([line.split()[0] for line in open(whitelist)])
    else:
        print("Error: can't find barcode whitelist file.")
        exit()    
    cbcs,qcbcs,cbcfreq_dict = cbcfreq(address,bcset)
    newcbcs = cbccorrect(cbcs,qcbcs,cbcfreq_dict)
    frags={ch:{} for ch in chrs}
    with io.BufferedReader(gzip.open(address,'rb')) as f:
        for i,line in enumerate(f):
            dlist = line.decode().split()
            ch=dlist[3]
            address='_'.join([dlist[4],dlist[5],newcbcs[i]])
            if address not in frags[ch]:
                frags[ch][address]=1
            else:
                frags[ch][address]+=1
    with open(fragfile,'w') as g:
        g.write('# id=%(sample)s\n' % vars())
        g.write('# description=\n')
        g.write('#\n')
        g.write('# pipeline_name=dna10x\n')
        g.write('# pipeline_version=dna10x-1.0.0\n')
        g.write('#\n')
        g.write('# reference_path=%(ref)s\n' % vars())
        g.write('# reference_fasta_hash=\n')
        g.write('# reference_gtf_hash=\n')
        g.write('# reference_version=\n')
        g.write('# mkref_version=\n')
        g.write('#\n')
        for ch in chrs:
            g.write('# primary_contig=%(ch)s\n' % vars())
        for ch in chrs:
            fdict=frags[ch]
            p1s=np.array([int(k.split('_')[0]) for k in fdict])
            ind=np.argsort(p1s)
            srtkey=np.array(list(fdict.keys()))[ind]
            for sk in srtkey:
                address=sk.split('_')
                st=ch+'\t'+address[0]+'\t'+address[1]+'\t'+address[2]+'\t'+str(fdict[sk])+'\n'
                g.write(st)
