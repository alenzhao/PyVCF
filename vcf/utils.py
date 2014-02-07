"""
Utilities for VCF files.
"""

def walk_together(*readers, **kwargs):
    """
    Simultaneously iteratate over two or more VCF readers. For each 
    genomic position with a variant, return a list of size equal to the number 
    of VCF readers. This list contains the VCF record from readers that have
    this variant, and None for readers that don't have it. 

    Args:
        vcf_record_sort_key: function that takes a VCF record and returns a 
            tuple that can be used as a key for comparing and sorting VCF 
            records across all readers. This tuple defines what it means for two 
            variants to be equal (eg. whether it's only their position or also 
            their allele values), and implicitly determines the chromosome 
            ordering since the tuple's 1st element is typically the chromosome 
            name (or calculated from it).
        ignore_chr: (default True) 'chrX' and 'X' will be considered equal. Only
            used if vcf_record_sort_key is not set.
        check_contig_order: (default True) raise an exception if readers have 
            incompatible contig orderings (eg. 1, 10, 2 .. vs.  1, 2, 10 ..)
    """
    check_contig_order = kwargs.get('check_contig_order', True)
    is_user_defined_key = 'vcf_record_sort_key' in kwargs
    if is_user_defined_key:
        get_key = kwargs['vcf_record_sort_key']
    else:
        if kwargs.get('ignore_chr', True):
            get_key = lambda r: (r.CHROM.replace('chr', ''), r.POS)
        else:
            get_key = lambda r: (r.CHROM, r.POS)
    
    nexts = []
    for reader in readers:
        try:
            nexts.append(reader.next())
        except StopIteration:
            nexts.append(None)

    min_k = (None,)   # keep track of the previous min key's contig
    finished_contigs = []   # used for checking contig order
    while any([r is not None for r in nexts]):
        next_idx_to_k = dict(
            (i, get_key(r)) for i, r in enumerate(nexts) if r is not None)
        keys_with_prev_contig = [
            k for k in next_idx_to_k.values() if k[0] == min_k[0]]

        if any(keys_with_prev_contig):
            min_k = min(keys_with_prev_contig)   # finish previous contig
        else:
            min_k = min(next_idx_to_k.values())   # move on to next contig

            if check_contig_order and len(finished_contigs) > 1:
                prev_contig = finished_contigs[-1]
                err = None
                if is_user_defined_key and prev_contig > min_k[0]:
                    err = ("Order of '%s', '%s' records in %s (or other VCFs) "
                        "contradicts the order defined by vcf_record_sort_key")
                elif min_k[0] in finished_contigs:
                    # Don't enforce the default contig ordering, since it may 
                    # not match what's in the VCF files. Just make sure the 
                    # ordering of emitted contigs is monotonic / without repeats
                    err = ("The order of contigs '%s' and '%s' is ambiguous or "
                        "contradictory in %s (or other VCFs). To define the "
                        "contig order, set the 'vcf_record_sort_key' arg.")
                if err:
                    bad_i = [i for i, k in next_idx_to_k.items() if k == min_k]
                    bad_file = readers[bad_i[0]].filename
                    raise Exception(err % (prev_contig, min_k[0], bad_file))
            finished_contigs.append(min_k[0])

        min_k_idxs = set([i for i, k in next_idx_to_k.items() if k == min_k])
        yield [nexts[i] if i in min_k_idxs else None for i in range(len(nexts))]

        for i in min_k_idxs:
            try:
                nexts[i] = readers[i].next()
            except StopIteration:
                nexts[i] = None


def trim_common_suffix(*sequences):
    """
    Trim a list of sequences by removing the longest common suffix while
    leaving all of them at least one character in length.

    Standard convention with VCF is to place an indel at the left-most
    position, but some tools add additional context to the right of the
    sequences (e.g. samtools). These common suffixes are undesirable when
    comparing variants, for example in variant databases.

        >>> trim_common_suffix('TATATATA', 'TATATA')
        ['TAT', 'T']

        >>> trim_common_suffix('ACCCCC', 'ACCCCCCCC', 'ACCCCCCC', 'ACCCCCCCCC')
        ['A', 'ACCC', 'ACC', 'ACCCC']

    """
    if not sequences:
        return []
    reverses = [seq[::-1] for seq in sequences]
    rev_min = min(reverses)
    rev_max = max(reverses)
    if len(rev_min) < 2:
        return sequences
    for i, c in enumerate(rev_min[:-1]):
        if c != rev_max[i]:
            if i == 0:
                return sequences
            return [seq[:-i] for seq in sequences]
    return [seq[:-(i + 1)] for seq in sequences]
