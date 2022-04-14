# -*- coding: utf-8 -*-
'''

'''


from gimmemotifs.scanner import as_fasta
from tqdm.auto import tqdm

# !!! This function is still under development. May be removed in the future version.
def custom_scan(self, seqs, nreport=100, scan_rc=True, zscore=False, gc=False, batch_size=50000, verbose=True):
    """
    This is the modified version of gimmemotifs's scanner.scan function.

    1) This function uses custom batch size.
    2) Verbose can be silenced.
    """
    seqs = as_fasta(seqs, genome=self.genome)
    if zscore:
        if gc:
            if len(self.meanstd) <= 1:
                self.set_meanstd(gc=gc)
        else:
            if len(self.meanstd) != 1:
                self.set_meanstd(gc=gc)

    # progress bar
    pbar = tqdm(
        desc="scanning",
        unit=" sequences",
        total=len(seqs),
        disable=(verbose==False),  # can be silenced
    )

    #logger.debug("Scanning")
    for batch_idx in range(0, len(seqs), batch_size):
        it = self._scan_sequences(
            seqs.seqs[batch_idx : batch_idx + batch_size],
            nreport,
            scan_rc,
            #zscore=zscore,
        )
        for result in it:
            yield result
            pbar.update(1)
    pbar.close()
