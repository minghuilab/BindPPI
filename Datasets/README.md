# Datasets

## About

<font size=4>

S802.txt: The training dataset for parameterizing BindPPI models and it contains 802 complexes.

S192.txt: The independent testing dataset, which contains 192 protein-protein heterodimer complexes.

S108.txt: The independent testing dataset, which contains 108 multimers.

S365.txt: The independent testing dataset, which contains 365 protein-peptide complexes.

</font> 

## Terms in the text files

<font size=4>

PDB ID: The PDB entry of the protein-protein complex.

Partner1: The protein chains of Partner1 in the protein-protein complex.

Partner2: The protein chains of Partner2 in the protein-protein complex.

dG: Experimental binding affinity of the protein-protein complex(in kcal/mol).

RF_{13}: The predicted values from model RF_{13}.RF_{13} is the Random forest model using 13 selected handcrafted features.

RF_{174}: The predicted values from model RF_{174}.RF_{174} is the Random forest model using 174 handcrafted features.

RF_{1428}: The predicted values from model RF_{1428}.RF_{1428} is the Random forest model using 1428 handcrafted features.

XGBoost_{13}: The predicted values from model XGBoost_{13}.XGBoost_{13} is the eXtreme Gradient Boosting model using 1428 handcrafted features.

MLP_{13}: The predicted values from model MLP_{13}.MLP_{13} is the Multilayer perceptron model using 13 selected handcrafted features.

MLP_{174}: The predicted values from model MLP_{174}.MLP_{174} is the Multilayer perceptron model using 174 handcrafted features.	

MLP_{1428}: The predicted values from model MLP_{1428}.MLP_{1428} is the Multilayer perceptron model using 1428 handcrafted features.	

RF_{5120}: The predicted values from model RF_{5120}.RF_{5120} is the Random forest model using 5120-dimensional embedding feature.

XGBoost_{5120}: The predicted values from model XGBoost_{5120}.XGBoost_{5120} is the eXtreme Gradient Boosting model using 5120-dimensional embedding feature.

MLP_{5120}: The predicted values from model MLP_{5120}.MLP_{5120} is the Multilayer perceptron model using 5120-dimensional embedding feature.

SC_{5120}: The predicted values from model SC_{5120}.SC_{5120} is the Multilayer perceptron model with skip connections using a 5120-dimensional embedding feature.

MHA_{5120}: The predicted values from model MHA_{5120}.MHA_{5120} is the Multilayer perceptron model with multi-head attention and skip connections using a 5120-dimensional embedding feature.

AvgEns: The predicted values from model AvgEns.AvgEns is the Ensemble model composed of RF_{13} and MLP_{5120} using the averaging method.

WtdAvgEns: The predicted values from model WtdAvgEns.WtdAvgEns is the Ensemble model composed of RF_{13} and MLP_{5120} using the weighted averaging method.

LREns: The predicted values from model LREns.LREns is the Ensemble model composed of RF_{13} and MLP_{5120} using linear regression.

RFEns: The predicted values from model RFEns.RFEns is the Ensemble model composed of RF_{13} and MLP_{5120} using random forest.

RF_{5120, 13}: The predicted values from model RF_{5120, 13}.RF_{5120, 13} is the Random forest model using a 5120-dimensional embedding feature and 13 selected handcrafted features.

RF_{5120, 174}:The predicted values from model RF_{5120, 174}.RF_{5120, 174} is the Random forest model using a 5120-dimensional embedding feature and 174 handcrafted features.

RF_{5120, 1428}: The predicted values from model RF_{5120, 1428}.RF_{5120, 1428} is the Random forest model using a 5120-dimensional embedding feature and 1428 handcrafted features.

MLP_{5120, 13}: The predicted values from model MLP_{5120, 13}.MLP_{5120, 13} is the Multilayer perceptron model using a 5120-dimensional embedding feature and 13 selected handcrafted features.

MLP_{5120, 174}: The predicted values from model MLP_{5120, 174}.MLP_{5120, 174} is the Multilayer perceptron model using a 5120-dimensional embedding feature and 174 handcrafted features.

MLP_{5120, 1428}: The predicted values from model MLP_{5120, 1428}.MLP_{5120, 1428} is the Multilayer perceptron model using a 5120-dimensional embedding feature and 1428 handcrafted features.	

MLP_{5120////256, 13/208}: The predicted values from model MLP_{5120////256, 13/208}.MLP_{5120////256, 13/208} is the Multilayer perceptron model involves increasing the dimensionality of each handcrafted feature from 1 to 16 using one MLP layers, and reducing the dimensionality of the 5120-dimensional complex embedding to 256 using four MLP layers. Subsequently, the up-sampled and down-sampled features are concatenated as inputs for the MLP architecture. 

MLP_{5120///512, 13//416}: The predicted values from model MLP_{5120///512, 13//416}.MLP_{5120///512, 13//416} is the Multilayer perceptron model involves increasing the dimensionality of each handcrafted feature from 1 to 32 using two MLP layers, and reducing the dimensionality of the 5120-dimensional complex embedding to 512 using three MLP layers. Subsequently, the up-sampled and down-sampled features are concatenated as inputs for the MLP architecture.

PPI-Affinity: The predicted values from PPI-Affinity.

PPA_Pred2: The predicted values from PPA_Pred2.

Minpredictor: The predicted values from Minpredictor.	

PRODIGY: The predicted values from PRODIGY.

ISLAND: The predicted values from ISLAND.

FoldX: The predicted values from FoldX.	

MMPBSA: The predicted values from MMPBSA.	

Rosetta: The predicted values from Rosetta.	

ZRANK: The predicted values from ZRANK.	

ZRANK2: The predicted values from ZRANK2.	

RosettaDock: The predicted values from RosettaDock.	

pyDock_Tot: The predicted values from pyDock_Tot.	

Sipper: The predicted values from Sipper.	

AP_Pisa: The predicted values from AP_Pisa.

FireDock: The predicted values from FireDock.	

FireDock_AB: The predicted values from FireDock_AB.	

FireDock_EI: The predicted values from FireDock_EI.

CP_PIE: The predicted values from CP_PIE.

</font> 