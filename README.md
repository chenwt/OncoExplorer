# Onco Knowledge Explorer
Explore cancer data interactively.

## Metrics
* # patients: 4453/4468

|can| blca | brca | coad | esca | gbm | hnsc | kirc | kirp | lihc | luad | lusc | ov | prad | read | stad | ucec |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
|#pat|  200|  841| 182 |  149 |  201 |  458 |  424 |  168 |  147 |  383 |  136 |  319 |  398 |  77|  176|  193|




## Connected pathway graph
* SGA sets: 559
* DEG sets: 8,841
* Nodes(SGA | DEG): 9,158 
* Directed edges(TDI pairs) = 1,496,128(?????)

## Noisy-or > 0.9
* train edges = 31,304
* test edges = 31,293 (overlap: 29,871)
* base graph = 32,726 (overlap: 31,172)

## Raw dataset
* Directed edges(TDI pairs, same path in different patients count as different ones): 11,769,129(filtered,...3,279,967).

Note: we do not remove the duplicated records, they account for 200/3,279,967.
We use the aggregate of SGA and DEG in advance.



* Check the overlap of SGA & DEG w/ KEGG pathway.

*TDI_Results_new.csv: 1 + 11,769,129 records (No duplication, containing unit SGA).
lowercase the record: tr A-Z a-z < TDI_Results_new.csv > out (Also No duplication, containing unit SGA)

min posterior = 0.011143. (sort -t$',' -k 4g TDI_Results_new.csv).

*TDI_Results_filter_no_unit.csv (1+3,281,866) [remove duplicated entities]->  (1+3,281,018) [overlap with TDI_Results_new.csv]-> ensemble.txt (1+3,281,018)

tr A-Z a-z < TDI_Results_filter_no_unit.csv > out

min posterior = 0.1

* SGAs.txt (559 filtered SGAs); DEGs.txt (8,844 filtered DEGs)

__Overlap:242; SGAs & hsa05200(cancer) = 25/179; DEGs & hsa05200 = 100/179;__

cut -f5 ensemble.txt | tail -n +2 | sort | uniq > SGAs.txt

cut -f6 ensemble.txt | tail -n +2 | sort | uniq > DEGs.txt


X (The sort here will change MARCH9 to 9-mar)cut -f2 -d',' TDI_Results_new_lc.csv | tail -n +2 | sort | uniq > SGAl.txt

SGAl.txt (19,850 unfiltered SGAs); DEGl.txt (19,414 unfiltered DEGs);

__Overlap:14,984; SGAs & hsa05200(cancer) = 168/179; DEGs & hsa05200 = 171/179;__

cut -f2 -d',' TDI_Results_new_lc.csv | tail -n +2 | sort -g | uniq > SGAl.txt

cut -f3 -d',' TDI_Results_new_lc.csv | tail -n +2 | sort -g| uniq > DEGl.txt


* Test with GO database for community detection.

* try different scoring methods.

* try different methods of clustering.

* grant.

EOF.
