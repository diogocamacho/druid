# DRUID: DRUg Indication Discoverer

DRUID, or DRUg Indication Discoverer, is an algorithm that identifies drug profiles that revert or mimic a condition of interest.  For example, DRUID can be used to repurpose compounds for novel indications given a gene expression profile, it can be used to prioritize compounds for a compendium of disease states, and it can be incorporated into computational chemistry pipelines for the identification and characterization of drug properties for novel design. 

At the core, DRUID builds on the notions of natural language processing and text mining, and the concepts of term frequency and inverse document frequency.  Briefly, given a corpus of documents, the term frequency (tf) and the inverse document frequency (idf) of any given word is given by:

$tf_{w, d} = \frac{word w in document d}{total number of words in document d}$

and 

$idf_{w, D} = log(\frac{N}{|\{w \in D : w \in d\}|})$

where N is the total number of documents in the corpus, $|\{w \in D : w \in d\}|$ is the number of documents where the term w appears, and d is a given document.  From this comes that 

$tfidf_{w, d, D} = tf_{w, d} * idf_{w, D}$

We have taken this concept to explore the identification of drugs (ie, drug profiles) to match a given phenotypic profile of interest.  

DRUID is built agnostically, which implies that any given omics data type can be used, though there needs to be consistency between the query vector and the compendium.
