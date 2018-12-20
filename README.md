# CONAN - Co-Variation Network Analyzer

This software analyzes a residue co-variation network with the aim to emphasize local evolutionary constraints and consequently detect functional and specificity determinant sites. Its only mandatory input is a multiple sequence alignment in one of the following formats: Pfam, Stockholm or FASTA; and the output consists of several networks and communities files from several cuts of the network.

This software was still not pubblished, but further information about this methodology can be accessed at:

[da Fonseca, N. J., Afonso, M. Q. L., de Oliveira, L. C., & Bleicher, L. (2018). A new method bridging graph 
theory and residue co-evolutionary networks for specificity determinant positions detection. Bioinformatics](https://academic.oup.com/bioinformatics/advance-article-abstract/doi/10.1093/bioinformatics/bty846/5123433?redirectedFrom=fulltext).

**Required external libraries**
- [Python-Levenshtein](https://pypi.org/project/python-Levenshtein/)
- [Numpy](https://pypi.org/project/numpy/)
- [Scipy](https://www.scipy.org/install.html)
- [Pandas](https://pandas.pydata.org/pandas-docs/stable/install.html)
- [NetworkX](https://networkx.github.io/documentation/latest/install.html)

**The mandatory inputs consits of:**
+ -i <filename> - A multiple sequence alignment file
+ -o <directory> - An output directory path

**Optional parameters:**
* -O <value> - Minimum Occupancy value. It is used to remove fragments of sequences. (0-1)
* -I <value> - Maximum Identity value. It is used to remove high identity sequences. (0-1)
* -f <value> - Minimum node frequency. It removes nodes (residues) that rarely occurs in the alignment. (0-1)
* -F <value> - Maximum node frequency. It removes higly conserved nodes (residues). (0-1)
* -m <0 or 1> - Method to Statistically validate the nodes.
  * 0 - Tumminello (Network based validation)
  * 1 - DRCN (Frequency based validation)
* -e <0 or 1> - Include marginally conservation properties.
  * 0 - Consider only co-variation between amino acids.
  * 1 - Also include stereochemical and structural amino acids properties.

Use -h for access the Help text.


You can also find more information at http://www.biocomp.icb.ufmg.br/conan

