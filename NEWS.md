# NEWS

Version 1.1.0 (09/08/2020)
+ filtering out drug matches with less than `min_matches`
+ added `AnnotationDbi` dependency
+ added `org.Hs.eg.db` dependency
+ added gene symbols to genes that are matched from the query to a given drug signature

Version 1.0.0 (09/08/2020)
+ re-versioning number for easier following of updates;
+ updated DESCRIPTION file to include the `biocViews:` call-out for Bioconductor installs;
+ fixed a bug on score function where signatures with 0 matches where getting good DRUID scores. 