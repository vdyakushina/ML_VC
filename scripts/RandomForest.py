#!/usr/bin/env python3
# Wrapper script to process the whole pipeline

import sys
import os
from pathlib import Path
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import roc_auc_score, accuracy_score, roc_curve, plot_confusion_matrix, confusion_matrix, ConfusionMatrixDisplay
import matplotlib.pyplot as plt
import joblib
from ML_func import prepare_f, accuracy_f
from func import read_resources

### Initiate
(set, vtp) = sys.argv[1:3]

mode='apply'
if len(sys.argv) >3:
	(mode, out) = sys.argv[3:5]
	Path(out).mkdir(exist_ok = True)

res = os.getcwd()+"/resources.csv"
path = read_resources(res)
model_path=path['model']
model = joblib.load('%sRFmodel_%s.joblib' %(model_path,vtp))

ds_start=pd.read_csv(set, sep = '\t')
ds=ds_start
ds=prepare_f(ds, vtp)
vrbls=['Mutect2_TO_FILTER', 'Mutect2_TO_TLOD','Mutect2_TNP_FILTER','Strelka_TO_FILTER','Strelka_TO_QUAL','Strelka_TNP_FILTER','SINVICT','VD_FILTER','VD_Q','SGA_VAF','SGA_SB','SGA_RepeatRefCount','SGA_DP']

vals=ds[vrbls]  # Features
if mode == 'validate':
	lbls=ds['cls']  # Labels

lbls_pred=model.predict(vals)
lbls_probs = model.predict_proba(vals)[:, 1] # Probabilies for each clas

set_prediction=pd.read_csv(set, sep = '\t')
set_prediction['variant_TRUE_ML']=lbls_pred
set_prediction.to_csv(set, sep='\t', index=False)

##Accuracy assesment. Only in case of mode == 'validate'
if mode == 'validate':
	accuracy_f(lbls, lbls_pred, lbls_probs, out, vtp, model)