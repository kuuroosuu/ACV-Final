function computedClass=myClassifier(DS, qcParam)
% myClassifier: My classifier for SU/V detection
% In this example, we use the single-gaussian classifier. You can change it to other classifiers of your choice.

computedClass=qcEval(DS, qcParam);