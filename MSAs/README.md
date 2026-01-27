To generate multiple sequence alignments (MSAs) for Potts analysis, readers should follow the protocol described in the “Methods and Protocols” section of the main text.  

For a detailed discussion of MSA generation and processing, see [Sternke et al., Methods in Enzymology 643:149-179 (2020)](https://doi.org/10.1016/bs.mie.2020.06.001).  Scripts for generating and processing an MSA using the outlined workflow, along with detailed instructions, can be found at https://github.com/msternke/protein-consensus-sequence.  

Briefly, sequences can be obtained a sequence database such as the PFAM or Interpro from a sequence data bases.  If sequences in the database are not aligned they can be aligned with a program such as [MAFFT](https://mafft.cbrc.jp/alignment/server/index.html).  Prior to alignment, it can be helpful to remove sequences that deviate from the median sequence length by greater than a threshold value (e.g. +/- 30% away from median value).

To reduce the effect of gaps in the MSA on the downstream analysis, positions with gap characters in >50% area sequences can be removed from the MSA. Furthermore, sequences in the MSA that are composed of > a threshold percentage gap characters (values ranging from 10%-50% were used in this study) can be removed.  
