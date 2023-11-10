#! /usr/bin/python
from pysam import AlignmentFile
import gzip


def get_raw_cbcfreq_dict(inbam):
    cbcfreq_dict = {}
    tot = 0
    with AlignmentFile(inbam, 'rb') as f:
        for read in f:
            if read.is_read1 and not read.is_secondary and not read.is_supplementary:
                cbc = read.get_tag('CR')
                qcbc = read.get_tag('CY')
                if cbc in cbcfreq_dict:
                    cbcfreq_dict[cbc]+=1
                else:
                    cbcfreq_dict[cbc]=1 
                tot+=1
    for cbc in cbcfreq_dict:
        cbcfreq_dict[cbc]=cbcfreq_dict[cbc]/tot  
    return cbcfreq_dict


def contig_address(address,chrs,bamout,bcset,insert_size,minscore):
    outfiles = [address+'.'+ch+'.address.txt.gz' for ch in chrs]
    print(bamout)
    print(len(outfiles))
    outputs = {ch:gzip.open(outfile,'wb') for ch,outfile in zip(chrs,outfiles)}
    ch_cbcdict={ch:[] for ch in chrs}
    ch_qcbcdict={ch:[] for ch in chrs}
    cbcfreq_dict = {}
    rdict={}
    tot=0
    with AlignmentFile(bamout,'rb') as f:
        print(f.filename)
        for read in f:
            score = read.get_tag('AS')
            isize = read.tlen
            rlen= read.rlen
            if read.is_read1 and score>0 and isize!=0 and score/rlen>minscore and abs(isize)<insert_size and abs(isize) > 6:
                readid = ':'.join(read.qname.split(':')[3:7])
                # cbc_qcbc = read.get_tag('CR').split(' CY:Z:')
                # cbc = cbc_qcbc[0]
                # qcbc = cbc_qcbc[1]
                cbc = read.get_tag('CR')
                qcbc = read.get_tag('CY')
                ch = f.getrname(read.reference_id)
                if isize>0:
                    p1 = read.reference_start+4
                    p2 = p1+isize-5
                else:
                    p1 = read.next_reference_start+4
                    p2 = p1-isize-5
                newline=readid+'\t'+cbc+'\t'+qcbc+'\t'+ch+'\t'+str(p1)+'\t'+str(p2)+'\t'+str(rlen)+'\t'+str(score)+'\t'+str(isize)
                rdict[readid]=newline
            else:
                readid = ':'.join(read.qname.split(':')[3:7])
                if readid in rdict:
                    score = read.get_tag('AS')
                    if score>0 and score/rlen>minscore:
                        rlen=read.rlen
                        nl=rdict[readid]
                        cbc=nl.split()[1]
                        ch=nl.split()[3]
                        newline = nl+'\t'+str(rlen)+'\t'+str(score)+'\n'
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
                            ch_qcbcdict[ch].append(nl.split()[2])
                        del rdict[readid]
    for cbc in cbcfreq_dict:
        cbcfreq_dict[cbc]=cbcfreq_dict[cbc]/tot
    return ch_cbcdict,ch_qcbcdict,cbcfreq_dict
