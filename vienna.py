import RNA

def getBPPM(sequence, structure = "", bppm_cutoff = 0.00001):
    """
        Requires ViennaRNAtools Python module
        Returns the base pair probability matrix using Vienna pf_fold, get_pr and free_pf_arrays functions.
        returns upper triangular matrix, whose entries exceed a threshold
    """
    bppm = {}
    
    
     #'--noPS', '-d 2', t, P
     
     
     
    if structure != "":
        RNA.cvar.fold_constrained = 1
    else:
        RNA.cvar.fold_constrained = 0
    #print "Before", structure
    RNA.pf_fold(sequence, structure)
    #print "After", structure
    seq_len = len(sequence)+1
    for i in xrange(1, seq_len):
        for j in xrange(1, seq_len):
            if i<j:
                bpp = RNA.get_pr(i,j)
                if bpp > bppm_cutoff:
                    bppm[str(i) + "_" + str(j)] = bpp
                else:
                    bppm[str(i) + "_" + str(j)] = 0
    RNA.free_pf_arrays()
    #print bppm
    #exit(1)
    return bppm
    
    
def getAccuracy(struct_dict, bppm):
    """
        Calculate average structuredness of given structure(stack) within bppm
        struct_dict is a dictionary of form  struct_dict[i] = j indicating the base pair (i,j)|i<j wihtin the bppm 

    """
    acc = 0
    for sq_i in struct_dict.keys():
        v = str(sq_i) + "_" + str(struct_dict[sq_i])
        if v in bppm:
            acc += bppm[v]
            #acc += math.pow(bppm[v], 2) / len(struct_stack)
    return acc

