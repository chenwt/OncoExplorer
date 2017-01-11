export PROPPR_JAVA_ARGS=-Xmx150g
proppr compile classify.ppr
proppr set --programFiles classify.wam:goAnn.cfacts:isA.cfacts

# proppr answer test_linked.examples --apr eps=1e-3 --threads 20 > answer_test.out
proppr answer test.examples --apr eps=1e-7 --threads 60 > answer_test.out
#proppr answer test_linked.examples --apr eps=1e-8 --threads 60 --unnormalized > answer_test.out

echo 'Done!'
# python ~/Github/OncoExplorer/src/predict_roc.py --filename 'normal_eps=1e-4'
