# Discovery2
A pipeline that identifies common and uncommon variants of microorganisms and viruses. 
This is possible by translating nucleotide sequences into 6 frame protein sequences and matching the sequences against different databases. 
It also produces a classification of each sequence and an interactive table of all results. 

We also have other results that I have not worked with or developed further. 
In addition, a list of sequences that do not match anything in the databases used in the pipeline is produced for each sample.

This pipeline is divided into two parallel runs. One is for reads analysis and others are after assembly data is produced and analyzed.

![alt text](/UML_diagram/discovery2.png)
