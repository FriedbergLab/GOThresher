## Command line options

#### Glossary:

_References_: a "Reference" is the source of the annotation: that can be a paper with a PubmedID used to identify it, or a knowledgebase such as Reactome

#### Here we explain the 

+ `--cutoff_prot=n`
If a reference annottates nore that `n` proteins, do not include it's annotations in the analysis. This is because we have observed that annotations coming from references annotating many proteins tend to have a low information content [1](https://journals.plos.org/ploscompbiol/article/authors?id=10.1371/journal.pcbi.1003063)

+ `--cutoff_attn=n`
Similar to the above, only removing annotations from a reference that provide more than a certain number of annotations.

+ `--output`
Output directory. (default?)

+ `--evidence`
A GO annotation is a statement about the function of a particular gene. Each annotation includes an [_evidence code_](http://geneontology.org/docs/guide-go-evidence-codes/) to indicate how the annotation to a particular term is supported. if `--evidence` is used, then only proteins annotated with the support of evidence codes listed after this argument will be included in the output. Special arguments include
++ EXPEC: all experimental evidence codes
++ COMPEC: all computational evidence codes

