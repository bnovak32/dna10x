#! /usr/bin/python
import gzip
import io

def addressct(address,addressct,stats):
    r1addresses = set()
    r2addresses = set()
    addresses = {}
    readcts={}
    with io.BufferedReader(gzip.open(address,'rb')) as f:
        for line in f:
            dlist = line.decode().split()
            cbc=dlist[1]
            ch=dlist[2]
            p1add = cbc,ch,dlist[3],dlist[7]
            p2add = cbc,ch,dlist[4],dlist[7]
            if cbc in readcts:
                readcts[cbc]+=1
            else:
                readcts[cbc]=1
            if p1add not in r1addresses and p2add not in r2addresses:
                r1addresses.add(p1add)
                r2addresses.add(p2add)
                if cbc in addresses:
                    addresses[cbc].append((ch,dlist[3],dlist[4],dlist[7]))
                else:
                    addresses[cbc] = [(ch,dlist[3],dlist[4],dlist[7])]  
    with open(addressct,'w') as g1, open(stats,'w') as g2:
        for cbc in addresses:
            address=addresses[cbc]
            for tup in address:
                st=cbc+'\t'+'\t'.join([p for p in tup])+'\n'
                g1.write(st)
            reads = readcts[cbc]
            frags = len(address)
            st = cbc+'\t'+str(reads)+'\t'+str(frags)+'\n'
            g2.write(st)
    return 0
    
            

