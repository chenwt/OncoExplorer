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
