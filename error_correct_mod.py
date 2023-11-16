#! /usr/bin/python
import numpy as np

qscores = ['!', '"', '#', '$', '%', '&', '\'', '(', ')', '*', '+', ',', '-', '.', '/', '0', '1', '2', '3', '4', '5',
           '6', '7', '8', '9', ':', ';', '<', '=', '>', '?', '@', 'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K',
           'L', 'M', 'N']
qdict = {q: float(i) for i, q in enumerate(qscores)}
nts = ['A', 'G', 'T', 'C']


def revcomp(seq):
    sdict = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
    newseq = [sdict[s] for s in seq]
    return ''.join(newseq[::-1])


def correct_cbc(cbc, qcbc, cbcfreq_dict):
    cbclen = len(cbc)

    T = 0
    k = 0
    tcbcs = []
    Llist = np.zeros(4 * cbclen)
    corrected = False
    corrected_cbc = cbc

    for i in range(cbclen):
        pedit = max([0.0005, pow(10.0, -qdict[qcbc[i]] / 10.0)])
        for nt in nts:
            tcbc = ''.join([n if j != i else nt for j, n in enumerate(cbc)])
            tcbcs.append(tcbc)
            if tcbc in cbcfreq_dict:
                L = cbcfreq_dict[tcbc] * pedit
                Llist[k] = L
                T += L
            k += 1
    if T > 0:
        Llist = np.array(Llist) / T
        if np.max(Llist) > 0.975:
            corrected = True
            corrected_cbc = tcbcs[np.argmax(Llist)]

    return corrected, corrected_cbc


def cbccorrect(cbcs, qcbcs, cbcfreq_dict):
    newcbcs = []
    if len(cbcs) > 0:
        # cbclen=len(cbcs[0])
        corrected = 0
        for cbc, qcbc in zip(cbcs, qcbcs):
            if qcbc == '0':
                newcbcs.append(cbc)
            else:
                is_corrected, corrected_cbc = correct_cbc(cbc, qcbc, cbcfreq_dict)
                if is_corrected:
                    corrected += 1
                newcbcs.append(corrected_cbc)
    print("Corrected CBC: {corrected}")
    return newcbcs
