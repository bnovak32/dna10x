#! /usr/bin/python
from pysam import AlignmentFile
import gzip
import io

def address(address,bamout,bclen,bcset,alignment_score,insert_size):
    cbcs = []
    qcbcs = []
    cbcfreq_dict = {}   
    tot=0 
    with gzip.open(address, 'wb') as g:
        with AlignmentFile(bamout,'rb') as f:
            for read in f:
                score = int(read.get_tag('AS'))
                isize = int(read.tlen)
                if read.is_read1:
                    if score>alignment_score and abs(isize)<insert_size and isize!=0:
                        readid = ':'.join(read.qname.split(':')[3:7])
                        cbc = read.get_tag('BC')[0:16]
                        ch = f.getrname(read.reference_id)
                        cbcs.append(cbc)
                        if cbc in bcset:
                            if cbc in cbcfreq_dict:
                                cbcfreq_dict[cbc]+=1
                            else:
                                cbcfreq_dict[cbc]=1
                            qcbcs.append('0')
                            tot+=1
                        else:
                            qcbcs.append(read.get_tag('BC')[16::])
                        if isize>0:
                            p1 = read.reference_start+4
                            p2 = p1+isize-5
                        else:
                            p1 = read.next_reference_start+4
                            p2 = p1-isize-5
                        newline=readid+'\t'+cbc+'\t'+ch+'\t'+str(p1)+'\t'+str(p2)+'\t'+str(score)+'\t'+str(isize)+'\n'
                        g.write(newline.encode())   
    for cbc in cbcfreq_dict:
        cbcfreq_dict[cbc]=cbcfreq_dict[cbc]/tot
    return cbcs,qcbcs,cbcfreq_dict

def contig_address(address,chrs,bamout,bclen,bcset,alignment_score,insert_size):
    outfiles = [address+'.'+ch+'.address.txt.gz' for ch in chrs]
    outputs = {ch:gzip.open(outfile,'wb') for ch,outfile in zip(chrs,outfiles)}
    ch_cbcdict={ch:[] for ch in chrs}
    ch_qcbcdict={ch:[] for ch in chrs}
    cbcfreq_dict = {}
    tot=0
    with AlignmentFile(bamout,'rb') as f:
        for read in f:
            if read.is_read1:
                score = int(read.get_tag('AS'))
                isize = int(read.tlen)
                if score>alignment_score and abs(isize)<insert_size and isize!=0:
                    readid = ':'.join(read.qname.split(':')[3:7])
                    cbc = read.get_tag('BC')[0:bclen]
                    ch = f.getrname(read.reference_id)
                    if isize>0:
                        p1 = read.reference_start+4
                        p2 = p1+isize-5
                    else:
                        p1 = read.next_reference_start+4
                        p2 = p1-isize-5
                    newline=readid+'\t'+cbc+'\t'+ch+'\t'+str(p1)+'\t'+str(p2)+'\t'+str(score)+'\t'+str(isize)+'\n'
                    outputs[ch].write(newline.encode())
                    ch_cbcdict[ch].append(cbc)
                    if cbc in bcset:
                        if cbc in cbcfreq_dict:
                            cbcfreq_dict[cbc]+=1
                        else:
                            cbcfreq_dict[cbc]=1
                        ch_qcbcdict[ch].append('0')
                        tot+=1
                    else:
                        ch_qcbcdict[ch].append(read.get_tag('BC')[bclen::])
    for cbc in cbcfreq_dict:
        cbcfreq_dict[cbc]=cbcfreq_dict[cbc]/tot
    return ch_cbcdict,ch_qcbcdict,cbcfreq_dict

