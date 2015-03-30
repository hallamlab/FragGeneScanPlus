#!/bin/bash
make clean
make
./FragGeneScan -s sample_inputs/SRS063040_mod.scaffolds-10000.fa -m 1024 -o out -w 0 -t illumina_1
diff out sample_outputs/10000_out
diff out.faa sample_outputs/10000_out.faa
diff out.ffn sample_outputs/10000_out.ffn
rm out
rm out.faa
rm out.ffn
