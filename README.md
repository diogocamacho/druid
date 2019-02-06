# DRUID: DRUg Indication Discoverer

DRUID, or DRUg Indication Discoverer, is an algorithm that identifies drug profiles that revert or mimic a condition of interest.  For example, DRUID can be used to repurpose compounds for novel indications given a gene expression profile, it can be used to prioritize compounds for a compendium of disease states, and it can be incorporated into computational chemistry pipelines for the identification and characterization of drug properties for novel design. 

## Term frequency and inverse document frequency
At the core, DRUID builds on the notions of natural language processing and text mining, and the concepts of term frequency and inverse document frequency.  Briefly, given a corpus of documents, the term frequency (tf) and the inverse document frequency (idf) of any given word is given by:

\[
tf_{w, d} = \frac{word w in document d}{total number of words in document d}
\]

and 

\[
idf_{w, D} = log(\frac{N}{|\{w \in D : w \in d\}|})
\]

where N is the total number of documents in the corpus, \[|\{w \in D : w \in d\}|\] is the number of documents where the term w appears, and d is a given document.  From this comes that 

\[
tfidf_{w, d, D} = tf_{w, d} * idf_{w, D}
\]

We have taken this concept to explore the identification of drugs (ie, drug profiles) to match a given phenotypic profile of interest.  

## TF-IDF for gene expression data
Consider a gene expression profile for a disease of interest.  Also, let's consider gene expression profiles of cell lines, individuals, that correspond to effects of drugs treatments (like those found the the [CMAP](https://www.broadinstitute.org/connectivity-map-cmap) or [LINCS](https://clue.io) repositories.)  We can conceptualize that each phenotype (drug or disease) is a collection of words that are represented by the effects on the effect of a given gene.  As an example, the connectivity map data for methotrexate on HL-60 cells at a concentration of 8.8 uM for 6h leads to 200 genes showing differential expression, including repression of CDC20, CDC37, or CDKN3, and activation of GATA2, MTRR, or NEU1.  Therefore, for the entire collection of drugs in the CMAP data (our compendium or corpus), we have a total of 6,100 documents (the drugs tested) and 19,672 words (the genes and how these are affected.)  We can then use the tf-idf calculation to define the relevance of each gene and each document in the entire corpus.  

Because drugs can have promiscuous effects on genes, we calculate 2 tf-idf matrices: one on the DxG space, and one in the GxD space.  These corrections account for the promiscuity of both genes (DxG matrix) and drugs (GxD space) and a corrected tf-idf matrix can be computed:

\[
ctfidf = tfidf_{DxG} x tfidf_{GxD}
\]

## DRUID
Now that the tf-idf matrix is defined, a similarity between a phenotype of interest and our compendium/corpus can be computed using a simple vector similarity.  Implemented in DRUID is a cosine similarity between the query vector and any given corpus vector (drug profile).  In order to assess the significance of a similarity score, we also define a set of random vectors of the same size as our query vector, to quantify the probability of the query vector similarity being better than a random similarity.  Finally, a DRUID score for any given drug is calculated taking into account all of these terms:

\[
DS_{d} = 1 + \frac{cs_{d}}{max(cs)} - log10(random_p)
\]

As can easily be seen, the DRUID score (DS) has a lower bound of 1 and an upper bound of 2 + log10(random_sets). 

## Considerations
DRUID is built agnostically, which implies that any given omics data type can be used, though there needs to be consistency between the query vector and the compendium.
