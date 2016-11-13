#ifndef __RUN_HMM_H
#define __RUN_HMM_H

void writeDNA();
void writeMeta();
void writeAminoAcids(FILE* aa_outfile_fp, thread_data* td, unsigned int buffer);

#endif
