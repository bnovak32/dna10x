#! /usr/bin/python
import gzip
import io

def addressct(address,addressct,stats):
    addresses = {}
    with io.BufferedReader(gzip.open(address,'rb')) as f:
        for line in f:
            dlist = line.decode().split()
            cbc = dlist[1]
            tup = dlist[2],dlist[3],dlist[4],dlist[7]
            if cbc in addresses:
                addresses[cbc].append(tup)
            else:
                addresses[cbc] = [tup]
    readcts = {cbc:len(addresses[cbc]) for cbc in addresses}
    addresses = {cbc:set(addresses[cbc]) for cbc in addresses}
    with open(addressct,'w') as g1, open(stats,'w') as g2:
        for cbc in addresses:
            for tup in addresses[cbc]:
                st=cbc+'\t'+'\t'.join([p for p in tup])+'\n'
                g1.write(st)
            reads = readcts[cbc]
            frags = len(addresses[cbc])
            st = cbc+'\t'+str(reads)+'\t'+str(frags)+'\n'
            g2.write(st)
    return 0
    
            

