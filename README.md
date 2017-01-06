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


* __Test with GO database for community detection.__
* examine the GO terms:
* SGA
cdkn2b  go:0000079
cdkn2b  go:0000086
cdkn2b  go:0004861
cdkn2b  go:0005515
cdkn2b  go:0005634
cdkn2b  go:0005654
cdkn2b  go:0005737
cdkn2b  go:0005829
cdkn2b  go:0007050
cdkn2b  go:0007093
cdkn2b  go:0008285
cdkn2b  go:0019901
cdkn2b  go:0030219
cdkn2b  go:0030511
cdkn2b  go:0031668
cdkn2b  go:0031670
cdkn2b  go:0042326
cdkn2b  go:0045944
cdkn2b  go:0048536
cdkn2b  go:0050680
cdkn2b  go:0071901
cdkn2b  go:2000134

grin2b  go:0000165
grin2b  go:0004972
grin2b  go:0005088
grin2b  go:0005234
grin2b  go:0005515
grin2b  go:0005622
grin2b  go:0005886
grin2b  go:0005887
grin2b  go:0006810
grin2b  go:0007215
grin2b  go:0007268
grin2b  go:0007611
grin2b  go:0008270
grin2b  go:0009986
grin2b  go:0016594
grin2b  go:0017146
grin2b  go:0030054
grin2b  go:0034220
grin2b  go:0035235
grin2b  go:0043005
grin2b  go:0043547
grin2b  go:0045211
grin2b  go:0045471
grin2b  go:0048013

mrpl27  go:0003735
mrpl27  go:0005515
mrpl27  go:0005743
mrpl27  go:0005762
mrpl27  go:0006412
mrpl27  go:0044822
mrpl27  go:0070125
mrpl27  go:0070126

rbm17   go:0000166
rbm17   go:0000380
rbm17   go:0003723
rbm17   go:0005515
rbm17   go:0005681
rbm17   go:0043234


tvp23c  go:0009306
id: GO:0009306
name: protein secretion
namespace: biological_process
alt_id: GO:0045166
alt_id: GO:0045731
def: "The controlled release of proteins from a cell." [GOC:ai]
subset: gosubset_prok
synonym: "glycoprotein secretion" NARROW []
synonym: "protein secretion during cell fate commitment" NARROW []
synonym: "protein secretion resulting in cell fate commitment" NARROW []
is_a: GO:0015031 ! protein transport
is_a: GO:0032940 ! secretion by cell

tvp23c  go:0016192
id: GO:0016192
name: vesicle-mediated transport
namespace: biological_process
alt_id: GO:0006899
def: "A cellular transport process in which transported substances are moved in membrane-bounded vesicles; transported substances are enclosed in the vesicle lumen or located in the vesicle membrane. The process begins with a step that directs a substance to the forming vesicle, and includes vesicle budding and coating. Vesicles are then targeted to, and fuse with, an acceptor membrane." [GOC:ai, GOC:mah, ISBN:08789310662000]
subset: goslim_aspergillus
subset: goslim_candida
subset: goslim_chembl
subset: goslim_generic
subset: goslim_pir
subset: goslim_pombe
subset: gosubset_prok
synonym: "nonselective vesicle transport" NARROW []
synonym: "protein sorting along secretory pathway" RELATED []
synonym: "vesicle trafficking" RELATED []
synonym: "vesicle transport" EXACT []
synonym: "vesicular transport" EXACT [GOC:mah]
is_a: GO:0006810 ! transport

tvp23c  go:0030173
id: GO:0030173
name: integral component of Golgi membrane
namespace: cellular_component
def: "The component of the Golgi membrane consisting of the gene products and protein complexes having at least some part of their peptide sequence embedded in the hydrophobic region of the membrane." [GOC:go_curators]
synonym: "Golgi integral membrane protein" RELATED []
synonym: "integral to Golgi membrane" NARROW []
is_a: GO:0031228 ! intrinsic component of Golgi membrane
is_a: GO:0031301 ! integral component of organelle membrane
relationship: part_of GO:0000139 ! Golgi membrane

* DEG

adamts5 go:0004222
adamts5 go:0005178
adamts5 go:0005515
adamts5 go:0005576
adamts5 go:0005578
adamts5 go:0005615
adamts5 go:0005788
adamts5 go:0006508
adamts5 go:0008201
adamts5 go:0008237
adamts5 go:0008270
adamts5 go:0022617
adamts5 go:0036066
adamts5 go:0042742
adamts5 go:0044691
adamts5 go:0050840

cyp4f2  go:0000038
cyp4f2  go:0001676
cyp4f2  go:0003091
cyp4f2  go:0003095
cyp4f2  go:0004497
cyp4f2  go:0005506
cyp4f2  go:0005515
cyp4f2  go:0005737
cyp4f2  go:0005789
cyp4f2  go:0006690
cyp4f2  go:0006691
cyp4f2  go:0007596
cyp4f2  go:0008217
cyp4f2  go:0008392
cyp4f2  go:0016324
cyp4f2  go:0016709
cyp4f2  go:0017144
cyp4f2  go:0018685
cyp4f2  go:0019369
cyp4f2  go:0019373
cyp4f2  go:0020037
cyp4f2  go:0031090
cyp4f2  go:0032304
cyp4f2  go:0032305
cyp4f2  go:0036101
cyp4f2  go:0042360
cyp4f2  go:0042361
cyp4f2  go:0042376
cyp4f2  go:0042377
cyp4f2  go:0043231
cyp4f2  go:0050051
cyp4f2  go:0052869
cyp4f2  go:0052871
cyp4f2  go:0052872
cyp4f2  go:0055078
cyp4f2  go:0055114
cyp4f2  go:0097258
cyp4f2  go:0097259
cyp4f2  go:0097267

kcnj13  go:0005242
kcnj13  go:0005887
kcnj13  go:0006813
kcnj13  go:0010107
kcnj13  go:0034765

pmepa1  go:0000139
pmepa1  go:0005515
pmepa1  go:0005654
pmepa1  go:0005886
pmepa1  go:0010008
pmepa1  go:0010991
pmepa1  go:0016021
pmepa1  go:0030512
pmepa1  go:0030521
pmepa1  go:0031901
pmepa1  go:0043231
pmepa1  go:0050699
pmepa1  go:0060394
pmepa1  go:0070412

tnni3k  go:0002027
tnni3k  go:0004672
tnni3k  go:0004674
tnni3k  go:0004871
tnni3k  go:0005515
tnni3k  go:0005524
tnni3k  go:0005634
tnni3k  go:0005737
tnni3k  go:0006468
tnni3k  go:0008022
tnni3k  go:0031013
tnni3k  go:0035556
tnni3k  go:0046872
tnni3k  go:0055117
tnni3k  go:0086069
tnni3k  go:1903779


* try different scoring methods, e.g., Jaccard.

* try different methods of clustering, e.g. using spetral clustering to check the differences between cancers.

* grant.

* Test within one subset:

molecular function
molecular activities of gene products
__cellular component__
where gene products are active
__biological process__
pathways and larger processes made up of the activities of multiple gene products.

* (Try with the unfiltered datasets of TDI pairs)

EOF.
