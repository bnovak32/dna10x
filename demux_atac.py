#! /usr/bin/python
import gzip
import io
import sys


def demux(r1,r2,r3,p,sra=False):
    i=0
    end=p+16
    with io.BufferedReader(gzip.open(r1,'rb')) as f1, io.BufferedReader(gzip.open(r2,'rb')) as f2, io.BufferedReader(gzip.open(r3,'rb')) as f3:
        for line1,line2,line3 in zip(f1,f2,f3):
            if i == 0:
                store11=line1.decode()
                store13=line3.decode()
                i+=1
            elif i==1:
                store21 = line1.decode()
                store23 = line3.decode()
                bc = line2.decode()[p:end]
                i+=1
            elif i==2:
                store31=line1.decode()
                store33=line3.decode()
                i+=1
            elif i==3:
                bcq = line2.decode()[p:end]
                if not sra:
                    newlines=store11.split()[0]+' BC:Z:'+bc+bcq+'\n'+store21+store31+line1.decode()+store13.split()[0]+' BC:Z:'+bc+bcq+'\n'+store23+store33+line3.decode()
                else:
                    bc_tag = 'CR:Z:'+bc
                    bcq_tag = 'CY:Z:'+bcq
                    read_name = "\t".join((store11.split(" ")[1],bc_tag, bcq_tag))
                    newlines="@"+read_name+'\n'+store21+"+"+'\n'+line1.decode()+"@"+read_name+'\n'+store23+"+"+'\n'+line3.decode()
                sys.stdout.write(newlines)
                i=0
    return 0
