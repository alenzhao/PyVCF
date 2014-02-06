"""
Utilities for VCF files.
"""

def walk_together(*readers, **kwargs):	
	""" Simultaneously iteratate two or more VCF readers and return
		lists of concurrent records from each
		reader, with None if no record present.  Caller must check the
		inputs are sorted in the same way and use the same reference
		otherwise behaviour is undefined.
		
		Args:
			vcf_record_sort_key: function that takes a VCF record and returns a tuple that can be used as the key for comparing and sorting VCF records across all given VCFReaders. The tuple's 1st element should be the contig name.
	"""
	if 'vcf_record_sort_key' in kwargs:
		get_key = kwargs['vcf_record_sort_key']
	else:
		get_key = lambda r: (r.CHROM, r.POS)
	
	nexts = []
	for reader in readers:
		try:
			nexts.append(reader.next())
		except StopIteration:
			nexts.append(None)

	min_k = (None,)   # keep track of the previous min key's contig
	while True:
		kdict = {i: get_key(x) for i,x in enumerate(nexts) if x is not None}
		keys_with_prev_contig = [k for k in kdict.values() if k[0] == min_k[0]]
		if any(keys_with_prev_contig):
			# finish all records from previous contig	
			min_k = min(keys_with_prev_contig) 
		else:			
			# move on to the next contig
			min_k = min(kdict.values())  
		
		min_k_idxs = set([i for i, k in kdict.items() if k == min_k])
		yield [nexts[i] if i in min_k_idxs else None for i in range(len(nexts))]

		for i in min_k_idxs:
			try:
				nexts[i] = readers[i].next()
			except StopIteration:
				nexts[i] = None
				
		if all([x is None for x in nexts]):
			break


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
