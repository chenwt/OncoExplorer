# deepwalkcmd
# simple version
# deepwalk --input data/graph_brca.adjlist --output brca_64.embeddings

wl=2
output_path=result

echo "----deepwalk----$filename ------- number walks = $wl ----------"

for filename in brca gbm do
  echo "----deepwalk----$filename ------- number walks = $nw ----------"
done

for filename in brca gbm
do
  for nw in 10 20 50 100 200
  do
    echo "----deepwalk----$filename ------- number walks = $nw ----------"
    time deepwalk --input data/graph_${filename}.mat --output $output_path/${filename}_km40_wl_${wl}_nw${nw}.embeddings --representation-size 40 --number-walks $nw --walk-length $wl --workers 30
  done

done





# Actually using version:
# deepwalk --input data/graph_brca.adjlist --output brca_km40_wl2_nw10.embeddings --representation-size 40 --number-walks 10 --walk-length 2 --workers 30

# deepwalk --input data/graph_brca.adjlist --output brca_km40_wl2_nw0.embeddings --representation-size 40 --number-walks 200 --walk-length 2 --workers 30