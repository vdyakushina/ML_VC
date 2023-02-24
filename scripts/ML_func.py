# The module contains usefull functions for the pipeline

from sklearn.ensemble import RandomForestClassifier
from sklearn.datasets import make_classification
from sklearn.model_selection import train_test_split
from sklearn.metrics import roc_auc_score, accuracy_score, roc_curve, plot_confusion_matrix, confusion_matrix, ConfusionMatrixDisplay
import pandas as pd
import matplotlib.pyplot as plt
import joblib

def prepare_f(tbl, var_tp):
  vrbls=['Mutect2_TO_FILTER', 'Mutect2_TO_TLOD','Mutect2_TNP_FILTER','Strelka_TO_FILTER','Strelka_TO_QUAL','Strelka_TNP_FILTER','SINVICT','VD_FILTER','VD_Q','SGA_VAF','SGA_SB','SGA_RepeatRefCount','SGA_DP']
  varbls_dict={'Mutect2_TO_TLOD':0, 'Strelka_TO_QUAL':0,'SINVICT':0, 'VD_Q':0, 'SGA_VAF':0, 'SGA_SB':0, 'SGA_DP':0, 
               'Mutect2_TO_FILTER':274, 'Mutect2_TNP_FILTER': 1140,'Strelka_TO_FILTER':84, 'Strelka_TNP_FILTER':640,
               'VD_FILTER': 36, 'SGA_RepeatRefCount':200}

  if var_tp=="INDEL":
    vrbls+=['Mutect2_TO_STRQ','Scalpel_TO_FILTER','Scalpel_TNP_FILTER']
    varbls_dict.update({'Mutect2_TO_STRQ':0, 'Scalpel_TO_FILTER':84,'Scalpel_TNP_FILTER':200})
  for (key,value) in varbls_dict.items():
    tbl[key] = tbl[key].fillna(value)

  return(tbl)


def accuracy_f(lables_test, lables_pred, lables_probs, out_folder, var_tp, md):
  accuracy=round(accuracy_score(lables_test, lables_pred), 4)
  roc_value = round(roc_auc_score(lables_test, lables_probs), 4)
  TP=len(lables_pred[(lables_pred==lables_test) & (lables_pred=='tr')])
  TN=len(lables_pred[(lables_pred==lables_test) & (lables_pred=='fls')])
  FN=len(lables_pred[(lables_pred!=lables_test) & (lables_pred=='tr')])
  FP=len(lables_pred[(lables_pred!=lables_test) & (lables_pred=='fls')])
  TPR=round(TP/(TP + FN), 4)*100
  TNR=round(TN/(TN + FP), 4)*100
  PPV=round(TP/(TP + FP), 4)*100
  NPV=round(TN/(TN + FN), 4)*100
  acc=pd.DataFrame({'accuracy':[accuracy], 'AUC':[roc_value], 'TPR':[TPR], 'TNR':[TNR], 'PPV':[PPV], 'NPV':[NPV]})
  acc.to_csv('%s/accuracy_%s.tsv' % (out_folder, var_tp), sep='\t')
  acc.to_csv('%s/accuracy_%s.tsv' % (out_folder, var_tp), sep='\t', index=False)
  cm = confusion_matrix(lables_test, lables_pred, labels=md.classes_) # Confusion_matrix
  disp = ConfusionMatrixDisplay(confusion_matrix=cm, display_labels=md.classes_)
  disp.plot()
  plt.savefig('%s/confusion_matrix_%s.png' %(out_folder, var_tp)) 
  plt.close()
  #Plot ROC curve
  fpr, tpr, thresholds = roc_curve(lables_test, lables_probs, pos_label='tr')
  plt.figure()
  plt.plot(fpr, tpr, color="darkred", lw=1.5, label="ROC curve (area = %0.3f)" % roc_value)
  plt.plot([0, 1], [0, 1], color="darkgreen", lw=1.5, linestyle="--")
  plt.xlim([-0.05, 1.05])
  plt.ylim([-0.05, 1.05])
  plt.xlabel("False Positive Rate")
  plt.ylabel("True Positive Rate")
  plt.title("Receiver operating characteristic example")
  plt.legend(loc="lower right")
  plt.savefig('%s/ROC_%s.png' %(out_folder, var_tp)) 
  plt.close()