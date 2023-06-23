import DeepPurpose.DTI as models
from DeepPurpose.utils import *
from DeepPurpose.dataset import *
import pandas as pd
import numpy as np

X_pred = pd.read_csv("D:\\Py\\Cmap\\CMap_Rproj_1\\bindingDB_traindat/bindingDB_train_dat_2.csv")
# X_pred_1 = X_pred[X_pred['test_label'] == True]
# X_pred = X_pred[X_pred['test_label'] == False]
X_drug = X_pred['smiles_rdkit']
X_target = X_pred['sequence']
y = X_pred['label']

drug_encoding = 'CNN'
target_encoding = 'CNN'
train, val, test = data_process(list(X_drug), list(X_target), list(y),
                                drug_encoding, target_encoding,
                                split_method='HTS',frac=[0.7,0.1,0.2])

test.to_csv("result/test_2.csv")

config = generate_config(drug_encoding = drug_encoding,
                         target_encoding = target_encoding,
                         cls_hidden_dims = [1024,1024,512],
                         train_epoch = 100,
                         LR = 0.0001,
                         batch_size = 256,
                         cnn_drug_filters = [32,64,96],
                         cnn_target_filters = [32,64,96],
                         cnn_drug_kernels = [4,6,8],
                         cnn_target_kernels = [4,8,12],
                        )
model = models.model_initialize(**config)

model.train(train, val, test)
model.save_model('bindingDB_train_sample_2')







