#!/usr/bin/env python
# coding: utf-8
import random
import pandas as pd
import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F
from torch import Tensor
from torch.utils.data import DataLoader 
from models import *
import sys
esm_path = './esm/'
sys.path.append(esm_path)
import esm
import json
import argparse
import os
import warnings
warnings.filterwarnings('ignore')


parser = argparse.ArgumentParser(description='Process input_sequence JSON file path.')
parser.add_argument('-f', '--file', type=str, help='Path to input_sequence JSON file')
parser.add_argument('-d', '--dir', type=str, help='Path to save the embedding JSON file')
parser.add_argument('-o', '--output', type=str, help='Path to save the predition file')

args = parser.parse_args()

json_file_path = args.file
embedding_file_path = args.dir
output_file = args.output
pathinput = os.path.join(os.getcwd(), 'inputfiles/')

if json_file_path is None:
    print('Please provide the JSON file path.')
    exit()
    
if embedding_file_path is None:
    embedding_file_path = os.path.join(os.getcwd(), 'embedding/')

os.makedirs(embedding_file_path, exist_ok=True)

if output_file is None:
    output_file = os.path.join(os.getcwd(), 'sequence_based/MLP5120_prediction.txt')

with open(json_file_path, 'r') as fw:
    seq_dict = json.load(fw)

model_name = "esm2_t36_3B_UR50D"
model, alphabet = esm.pretrained.esm2_t36_3B_UR50D()
batch_converter = alphabet.get_batch_converter()
model.eval()
def get_embed(datatmp):    
    batch_labels, batch_strs, batch_tokens = batch_converter(datatmp)
    with torch.no_grad():
        results = model(batch_tokens, repr_layers = [36], return_contacts=True)
    token_representations = results["representations"][36]
    sequence_representations = []
    for i, (_, seq) in enumerate(datatmp):
        sequence_representations.append(token_representations[i, 1 : len(seq) + 1].mean(0))
    final1 = {}
    for pdb_chain,representations in zip(datatmp,sequence_representations):
        final1[pdb_chain[0]] = representations
    return final1

complex_data = []
for comp in seq_dict:
    for partner in seq_dict[comp]:
        for seq,seq_id in zip(seq_dict[comp][partner],range(len(seq_dict[comp][partner]))):
            complex_data.append((comp+'_'+partner+'_'+str(seq_id),seq[0]))

for i in complex_data:
    embedding_dict = get_embed([i])
    embedding_dict[i[0]] = embedding_dict[i[0]].numpy().tolist()
    with open(f'{embedding_file_path}{i[0]}.json','w') as f:
        json.dump(embedding_dict,f)
all_torch = {}
for comp in seq_dict:
    comp_tensors = []
    for partner in seq_dict[comp]:
        partner_tensors = []
        for seq_id in range(len(seq_dict[comp][partner])):
            id_use = comp+'_'+partner+'_'+str(seq_id)
            with open(f'{embedding_file_path}{id_use}.json','r') as fw:
                partner_tensor = json.load(fw)
                partner_tensors.append(torch.tensor(list(partner_tensor.values())).squeeze(0))
        comp_tensors.append(torch.mean(torch.stack(partner_tensors), dim=0))
    all_torch[comp] = torch.cat(comp_tensors,axis = 0) 
test_keys = [key for key in all_torch]
test_dataset = Dataset(test_keys,all_torch)
test_dataloader = DataLoader(test_dataset, batch_size = 32, shuffle=True) 
label_list = {}
node_dims = 5120
device = torch.device('cuda:0' if torch.cuda.is_available() else 'cpu')
loss_fn = nn.MSELoss()
df_list = []
for num in range(50):
    seed = 0
    torch.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)
    np.random.seed(seed)
    random.seed(seed) 
    torch.backends.cudnn.deterministic = True
    model_test = torch.load(pathinput + '/mlp5120_models/mlp5120_{}.pkl'.format(num),map_location=device)
    model_test.eval()
    pred_data = []
    lable_data = []
    key_data = []
    loss_data = []
    for i_batch, sample_test in enumerate(test_dataloader):
        graph = sample_test['data'].to(device)
        graph = torch.squeeze(graph , 1)
        key = sample_test['key']
        logits = model_test.forward(graph)
        pred = torch.squeeze(logits , 1)
        pred_value = torch.squeeze(logits, 1)
        for i in list(pred_value.cpu().detach().numpy()):
            pred_data.append(i)
        for i in key:
            key_data.append(i)
    tmp_df = pd.DataFrame({'key':key_data,'pred':pred_data})
    df = tmp_df.sort_values(['key']).reset_index(drop=True)
    df_list.append(df)
df_mean = pd.concat(df_list,axis = 1)
df_mean['MLP_{5120}'] = np.average(df_mean['pred'], axis=1)
df_mean['Complex'] = df_mean['key'].iloc[:,:1]
df_mean = df_mean[['Complex','MLP_{5120}']]
df_mean.to_csv(f'{output_file}',sep='\t',index = 0)



