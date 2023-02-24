#!/usr/bin/env python3
# Wrapper script to process the whole pipeline

import sys
from pathlib import Path
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.datasets import make_classification
from sklearn.model_selection import train_test_split
from sklearn.metrics import roc_auc_score, accuracy_score, roc_curve, plot_confusion_matrix, confusion_matrix, ConfusionMatrixDisplay
import matplotlib.pyplot as plt
import joblib
from ML_func import prepare_f, accuracy_f

### Initiate
(set, out, vtp) = sys.argv[1:4]
Path(out).mkdir(exist_ok = True)

ds_start=pd.read_csv(set, sep = '\t')
ds=ds_start
ds=prepare_f(ds, vtp)
vrbls=['Mutect2_TO_FILTER', 'Mutect2_TO_TLOD','Mutect2_TNP_FILTER','Strelka_TO_FILTER','Strelka_TO_QUAL','Strelka_TNP_FILTER','SINVICT','VD_FILTER','VD_Q','SGA_VAF','SGA_SB','SGA_RepeatRefCount','SGA_DP']

vals=ds[vrbls]  # Features
lbls=ds['cls']  # Labels

# Split dataset into training set and test set
vals_train, vals_test, lbls_train, lbls_test = train_test_split(vals, lbls, test_size=0.3) # 70% training and 30% test

#Define N trees and features
ntrees_acc={}
for ntrs in range(100, 550, 50):
  model=RandomForestClassifier(n_estimators=ntrs,
                             bootstrap = True,
                             max_features = 'sqrt')
  model.fit(vals_train,lbls_train)
  lbls_pred=model.predict(vals_test)
  ntrees_acc[ntrs]=accuracy_score(lbls_test, lbls_pred)

plt.plot(
    list(ntrees_acc.keys()),
    [round(100 * i, 3) for i in ntrees_acc.values()],
    color="black",
    lw=1.5)
plt.ylabel("Accuracy")
plt.xlabel("Number of trees")
plt.title("%s. Numer of features = 'sqrt'" % vtp)
plt.margins(y=0.5)
plt.savefig('%s/Number_of_trees_%s.png' %(out, vtp))
plt.close()

ntrs=500 if vtp=='SNV' else 300
mtry_acc={}
for mtrs in range(2, 7):
  model=RandomForestClassifier(n_estimators=ntrs,
                             bootstrap = True,
                             max_features = mtrs)
  model.fit(vals_train,lbls_train)
  lbls_pred=model.predict(vals_test)
  mtry_acc[mtrs]=accuracy_score(lbls_test, lbls_pred)

plt.plot(
    list(mtry_acc.keys()),
    [round(100 * i, 3) for i in mtry_acc.values()],
    color="black",
    lw=1.5)
plt.ylabel("Accuracy")
plt.xlabel("Number of features")
plt.title("%s. Numer of trees = %s" % (vtp, ntrs))
plt.margins(y=0.5)
plt.savefig('%s/Number_of_features_%s.png' %(out, vtp)) 
plt.close()

#Train a Model, Save
(ntrs, mtrs)= [500, 3] if vtp=='SNV' else [300,4]
model=RandomForestClassifier(n_estimators=ntrs,
                             bootstrap = True,
                             max_features = mtrs)
model.fit(vals_train,lbls_train)

joblib.dump(model, '/pipeline/scripts/RFmodel_%s.joblib' %vtp) #Save model

##Add Predictions to source dataset
lbls_pred=model.predict(vals_train)
set_prediction=ds_start.iloc[list(vals_train.index)]
set_prediction['Predicted']=lbls_pred
set_prediction.to_csv('%s/set_train_prediction_%s.tsv' %(out, vtp), sep='\t', index=False)

##Accuracy assesment
lbls_pred=model.predict(vals_test) #Classification on test sample
lbls_probs = model.predict_proba(vals_test)[:, 1] # Probabilies for each clas
accuracy_f(lbls_test, lbls_pred, lbls_probs, out, vtp, model)

fi = pd.DataFrame({'feature': vrbls,
                   'importance': model.feature_importances_}).\
                    sort_values('importance', ascending = False)
fi.to_csv('%s/feature_importance_%s.tsv' %(out, vtp), sep='\t', index=False)

set_prediction=ds_start.iloc[list(vals_test.index)]
set_prediction['Predicted']=lbls_pred
set_prediction.to_csv('%s/set_test_prediction_%s.tsv' %(out, vtp), sep='\t', index=False)