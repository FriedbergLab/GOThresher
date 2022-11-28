## An expanded explanation of GOThresher

Here we provide a basic glossary and an expanded explanation for some of the command line options in GOThresher.

#### Glossary:
+ _Aspect_: The Gene Ontology comprises three apects: Molecular Function (MFO, F), Biological Process (BPO, P), and Cellular Component (CCO, C). these are also known as _Categories_ in the Gene Ontology. Each protein may have one or more annotations from each category. For more information see: [[1]](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3037419/)
    + _Molecular function_ is defined as the biochemical activity (including specific binding to ligands or structures) of a gene product. 
    + _Biological process_ refers to a biological objective to which the gene or gene product contributes. A process is accomplished via one or more ordered assemblies of molecular functions.
    + _Cellular component_ refers to the place in the cell where a gene product is active.

+ _Evidence Code_ : A GO annotation is a statement about the function of a particular gene. Each annotation includes an _evidence code_ to indicate how the annotation to a particular term is supported. There are several types of evidence codes, and GOThresher allows for filtering based on them.  Evidence codes can be grouped by their type. For example, the Experimental Evidence Codes  indicate that there is evidence from an experiment directly supporting the annotation of the gene. The Computatioanl Analysis evidence codes indicate that the annotation is based on a computational analysis of the gene sequence and/or other data as described in the cited reference.  For more information see the relevant [GO pages](http://geneontology.org/docs/guide-go-evidence-codes/).

+ _Information content_: the information content (IC) of a GO term is a numerical representation of how specific that term is. The term ["Catalytic Activity"](https://www.ebi.ac.uk/QuickGO/term/GO:0003824) has a lower information content than ["Hydrolase activity"](https://www.ebi.ac.uk/QuickGO/term/GO:0016787) which, in turn, has a lower information content than ["Alpha Amylase Activity"](https://www.ebi.ac.uk/QuickGO/term/GO:0004556). 
    + The _Phillip_lord _ information content is calculated as `IC = -log(P(i))` where `P(i)` is the frequency of the GO term `i` in the corpus. For further reading see: [[2]](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1000443).
    + The _Wyatt_Clark_ information content  is calculated by treating the Gene Ontolgoy as a Bayesian network, with each  as defined in [[3]](https://academic.oup.com/bioinformatics/article/29/13/i53/195366).
+ _References_: a "Reference" is the source of the annotation: that can be a paper with a PubmedID used to identify it, or a knowledgebase such as [Reactome](https://reactome.org/).


#### Command line options expanded 

+ `--cutoff_prot=n`
If a reference annotates more that `n` proteins, do not include its annotations in the analysis. This is because we have observed that annotations coming from references annotating many proteins usually describe high throughput experiments and contain annotations with a low information content [[4]](https://journals.plos.org/ploscompbiol/article/authors?id=10.1371/journal.pcbi.1003063)

+ `--cutoff_attn=n`
Similar to the above, only removing annotations from a reference that provide more than a certain number of annotations (as opposed to a threshold based on a number of proteins).

+ `--output`
Output directory, other than the default.

+ `--evidence`
A GO annotation is a statement about the function of a particular gene. Each annotation includes an [_evidence code_](http://geneontology.org/docs/guide-go-evidence-codes/) to indicate how the annotation to a particular term is supported. if `--evidence` is used, then only proteins annotated with the support of evidence codes listed after this argument will be included in the output. There are also special arguments that are sets of several evidence codes:

    + EXPEC: all experimental evidence codes
    + COMPEC: all computational evidence codes
    + AUTHREC: all evidence codes derived from author statements
    + CUREC: all evidence codes derived by curator

+ `--evidence_inverse` excludes the evidece codes given in this argument.
+ `--aspect` filter by any combination of Biological Process, Cellular Component, and Molecular Function by using the letters P C and F respectively.
+ `--assigned_by` filter in which database assigns the annotation. Default: all 
+ `--assigned_by_inverse` filter out annotations assigned by theses databases
+ `--recalculate`: calculating the information accretion (aka Wyatt Clark) for a GO DAG takes time. This only needs to be used if running on a new corpus of data.
+ `--info_threshold`: annotations with an information content (Wyatt Clark or Phillip Lord)  that is below a certain percentile or absolute value will be filtered out. Example: -WCTHRESH 20 heeps the top 20% of annotation

 



