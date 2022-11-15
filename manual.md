## An expanded explanation of GOThresher

#### Glossary:
+ _Information content_: the information content (IC) of a GO term is a numerical representation of how specific that term is. The term ["Catalytic Activity"](https://www.ebi.ac.uk/QuickGO/term/GO:0003824) has a lower information content than ["Hydrolase activity"](https://www.ebi.ac.uk/QuickGO/term/GO:0016787) which, in turn, has a lower information content than ["Alpha Amylase Activity"](https://www.ebi.ac.uk/QuickGO/term/GO:0004556). The inforamtion content is calculates as `IC = -log(P(i))` where `P(i)` is the frequency of the GO term `i` in the corpus. For further reading see: [[1]](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1000443)
+ _References_: a "Reference" is the source of the annotation: that can be a paper with a PubmedID used to identify it, or a knowledgebase such as Reactome
+ _Wyatt_Clark_: 

#### Command line options expanded 

+ `--cutoff_prot=n`
If a reference annottates nore that `n` proteins, do not include it's annotations in the analysis. This is because we have observed that annotations coming from references annotating many proteins tend to have a low information content [[2]](https://journals.plos.org/ploscompbiol/article/authors?id=10.1371/journal.pcbi.1003063)

+ `--cutoff_attn=n`
Similar to the above, only removing annotations from a reference that provide more than a certain number of annotations.

+ `--output`
Output directory. (default?)

+ `--evidence`
A GO annotation is a statement about the function of a particular gene. Each annotation includes an [_evidence code_](http://geneontology.org/docs/guide-go-evidence-codes/) to indicate how the annotation to a particular term is supported. if `--evidence` is used, then only proteins annotated with the support of evidence codes listed after this argument will be included in the output. There are also special arguments that are sets of several evidence codes:

    + EXPEC: all experimental evidence codes
    + COMPEC: all computational evidence codes
    + AUTHREC: all evidence codes derived from author statements
    + CUREC: all evidence codes derived by curator

+ `--evidence_inverse`: excludes the evidece codes in this argument



