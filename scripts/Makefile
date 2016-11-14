export PROPPR_JAVA_ARGS=-Xmx450g
proppr compile pathway.ppr
proppr set --programFiles pathway.wam:isDEG.cfacts:sga2deg.graph

proppr answer test.examples --apr eps=1e-4 --threads 20 > answer_test.out
echo 'calculating roc...'
python ~/Github/OncoExplorer/src/predict_roc.py --filename 'normal_eps=1e-4'
