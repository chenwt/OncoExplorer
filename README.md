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



* Test with organelles structure, other than biological process.
* Consider the directions of the perturbations? (Future direction)
* (Try with the unfiltered datasets of TDI pairs)


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


* __Test with GO database for community detection.__
* examine the GO terms:
* SGA
cdkn2b  go:0000079(B)-
cdkn2b  go:0000086(B)
cdkn2b  go:0004861(M)
cdkn2b  go:0005515(M)
cdkn2b  go:0005634(C)
cdkn2b  go:0005654(C)
cdkn2b  go:0005737(C)
cdkn2b  go:0005829(C)
cdkn2b  go:0007050(B)-+
cdkn2b  go:0007093(B)
cdkn2b  go:0008285(B)-
cdkn2b  go:0019901(M)
cdkn2b  go:0030219(B)
cdkn2b  go:0030511(B)-+
cdkn2b  go:0031668(M)
cdkn2b  go:0031670(B)
cdkn2b  go:0042326(B)-
cdkn2b  go:0045944(B)
cdkn2b  go:0048536(B)-+
cdkn2b  go:0050680(B)
cdkn2b  go:0071901(B)-
cdkn2b  go:2000134(B)

 It seems that the nodes do not contain each other in the relationship of 'is_a'..

rbm17   go:0000166
rbm17   go:0000380
rbm17   go:0003723
rbm17   go:0005515
rbm17   go:0005681
rbm17   go:0043234

* __We could use the 'is_a' relationship to cluster the DEGs. Further, we can propogate the SGA->DEG->GOTerm->subsets goterm of 8150.__

* Examine the class of SGA through GOTerm, or through DEG and then GOTerm, compare the differences.

![alt tag](https://github.com/yifengtao/OncoExplorer/blob/master/figure/fig8_rbm17.jpg)


* try different methods of clustering, e.g. using spetral clustering to check the differences between cancers.

* Test with BRCA, GBM, OV.

python analysis02.py --inputData /usr1/public/yifeng/Github/outputData --outputData /usr1/public/yifeng/Github/outputData --cancer gbm


BRCA: 480,577 pat,sga,deg records -> aggregate 22,499, 21,141 SGA edges, 529 nodes.

GBM: 151,030 -> 7,150, 3,043 SGA edges, 313 nodes.

OV: 110,125 -> 9,397, 6,251 SGA edges, 418 nodes.

![alt tag](https://github.com/yifengtao/OncoExplorer/blob/master/figure/fig1_hist_weight.jpg)

It seems that merely aggregate the total number of overlap is not enough to distinguish the clusters...

![alt tag](https://github.com/yifengtao/OncoExplorer/blob/master/figure/fig2_ov_cluster1.jpg)

* try different scoring methods, e.g., Jaccard.

python analysis03.py --inputData /usr1/public/yifeng/Github/outputData --outputData /usr1/public/yifeng/Github/outputData --cancer gbm

The parameter of edges and nodes is the same with the raw methods.

It seems to be better.
![alt tag](https://github.com/yifengtao/OncoExplorer/blob/master/figure/fig3_ov_cluster2.jpg)

This can also be seen from the distribution of the (SGA,SGA) weight, which shows some edges have much larger weight, other than the noise in simple overlap version.
![alt tag](https://github.com/yifengtao/OncoExplorer/blob/master/figure/fig6_hist_SGAweight_jaccard.jpg)
Simple overlap version:
![alt tag](https://github.com/yifengtao/OncoExplorer/blob/master/figure/fig7_hist_SGAweight_simple.jpg)
But may become better if we can set a threshold?
Well... Actually I am really not sure if this is better after the experiment... if merely from the figure. Might need some other source?

Remember to check duplicated edges among patients.


* Our methods is able to analyze within different types of cancers with different number of patients.

* grant.

* Test within one subset:

molecular function
molecular activities of gene products
__cellular component__
where gene products are active
__biological process__
pathways and larger processes made up of the activities of multiple gene products.

* t-SNE.
It looks crazy, even if I tried to carefully tune the perplexity, the graph in 2D is hard to cluster into different classes. The example of BRCA is shown below:
![alt tag](https://github.com/yifengtao/OncoExplorer/blob/master/figure/fig9_tsne_brca.jpg)

There might be some other parameters to tune, e.g., the calculation of weight and distance.







* Also we can consider the role within through patient-wise rule.


Note that we may want to filter the GoTerm which is not related to biological process in the goAnn.cfacts file: Already filltered.



* Continue to test with pathway, although they overlap little.

* It would be pretty interesting to check the results from three parts: simplified (SGA,SGA) weighted.
* Check with SGA->GoTerm classification.
* Check with SGA->DEG->GoTerm classification.

EOF.
