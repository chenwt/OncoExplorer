# Onco Knowledge Explorer
Explore cancer data interactively.

## Metrics
* # patients: 4453/4468
* blca 200, brca 1041, cancer 1042, coad 1224, esca 1373, gbm 1574, hnsc 2032, kirc 2456, kirp 2624, lihc 2771, luad 3154, lusc 3290, ov 3609, prad 4007, read 4084, stad 4260, ucec 4453

| blca 200| brca 1041| cancer 1042| coad 1224| esca 1373| gbm 1574| hnsc 2032| kirc 2456| kirp 2624| lihc 2771| luad 3154| lusc 3290| ov 3609| prad 4007| read 4084| stad 4260| ucec 4453|
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| blca 200| brca 841| cancer 1| coad 182 | esca 149 | gbm 201 | hnsc 458 | kirc 424 | kirp 168 | lihc 147 | luad 383 | lusc 136 | ov 319 | prad 398 | read 77| stad 176| ucec 193|


| cancer    | blca |
|--------|----------------------|
|# patients | 200  |

| Tables        | Are           | Cool  |
| ------------- |:-------------:| -----:|
| col 3 is      | right-aligned | $1600 |
| col 2 is      | centered      |   $12 |
| zebra stripes | are neat      |    $1 |



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
