import torch
import torch.nn as nn
import torch.nn.functional as F
from torch_geometric.nn import GCNConv, GATConv
from torch_geometric.nn.conv import MessagePassing
from torch_geometric.nn.inits import glorot, uniform
from torch_geometric.utils import softmax as Softmax
from torchmetrics.functional import pairwise_cosine_similarity
import math

# General convolution class for different types of GNN layers
class GeneralConv(nn.Module):
    def __init__(self, conv_name, in_hid, out_hid, num_types, num_relations, n_heads, dropout, use_norm=True):
        super(GeneralConv, self).__init__()
        self.conv_name = conv_name
        self.res_att = None
        # Choose graph convolution layers based on conv_name
        if self.conv_name == 'hgt':
            self.base_conv = HGTConv(in_hid, out_hid, num_types, num_relations, n_heads, dropout, use_norm)
        elif self.conv_name == 'gcn':
            self.base_conv = GCNConv(in_hid, out_hid)
        elif self.conv_name == 'gat':
            self.base_conv = GATConv(in_hid, out_hid // n_heads, heads=n_heads)

    def forward(self, meta_xs, node_type, edge_index, edge_type):
        if self.conv_name == 'hgt':
            a = self.base_conv(meta_xs, node_type, edge_index, edge_type)
            self.edge_attention_weights = self.base_conv.edge_attention_weights
            return a
        elif self.conv_name == 'gcn':
            return self.base_conv(meta_xs, edge_index)
        elif self.conv_name == 'gat':
            return self.base_conv(meta_xs, edge_index)
        elif self.conv_name == 'dense_hgt':
            return self.base_conv(meta_xs, node_type, edge_index, edge_type)

# HGTConv: Heterogeneous Graph Transformer Convolution Layer
class HGTConv(MessagePassing):
    def __init__(self, in_dim, out_dim, num_types, num_relations, n_heads, dropout=0.2, use_norm=True, distribution='uniform', **kwargs):
        super(HGTConv, self).__init__(node_dim=0, aggr='add', **kwargs)

        self.in_dim = in_dim
        self.out_dim = out_dim
        self.num_types = num_types
        self.num_relations = num_relations
        self.total_rel = num_types * num_relations * num_types
        self.n_heads = n_heads
        self.d_k = out_dim // n_heads
        self.sqrt_dk = math.sqrt(self.d_k)
        self.use_norm = use_norm
        self.distribution = distribution
        self.att = None
        self.res_att = None
        self.res = None

        # Initialize layers for the graph transformation
        self.k_linears = nn.ModuleList()
        self.q_linears = nn.ModuleList()
        self.v_linears = nn.ModuleList()
        self.a_linears = nn.ModuleList()
        self.norms = nn.ModuleList()

        for t in range(num_types):
            self.k_linears.append(nn.Linear(in_dim, out_dim))
            self.q_linears.append(nn.Linear(in_dim, out_dim))
            self.v_linears.append(nn.Linear(in_dim, out_dim))
            self.a_linears.append(nn.Linear(out_dim, out_dim))
            if use_norm:
                self.norms.append(nn.LayerNorm(out_dim))

        # Parameters for attention weights and message passing
        self.relation_pri = nn.Parameter(torch.ones(num_relations, self.n_heads))
        self.relation_att = nn.Parameter(torch.Tensor(num_relations, n_heads, self.d_k, self.d_k))
        self.relation_msg = nn.Parameter(torch.Tensor(num_relations, n_heads, self.d_k, self.d_k))
        self.skip = nn.Parameter(torch.ones(num_types))
        self.drop = nn.Dropout(dropout)

        # Initialize weight matrices
        glorot(self.relation_att)
        glorot(self.relation_msg)

    def forward(self, node_inp, node_type, edge_index, edge_type):
        return self.propagate(edge_index, node_inp=node_inp, node_type=node_type, edge_type=edge_type)

    def message(self, edge_index_i, node_inp_i, node_inp_j, node_type_i, node_type_j, edge_type, edge_index):
        self.edge_attention_weights = {}
        data_size = edge_index_i.size(0)
        self.res_att = torch.zeros(data_size, self.n_heads).to(node_inp_i.device)
        res_msg = torch.zeros(data_size, self.n_heads, self.d_k).to(node_inp_i.device)

        # Iterate over node types and relations to compute attention and messages
        for source_type in range(self.num_types):
            sb = (node_type_j == int(source_type))
            k_linear = self.k_linears[source_type]
            v_linear = self.v_linears[source_type]
            for target_type in range(self.num_types):
                tb = (node_type_i == int(target_type)) & sb
                q_linear = self.q_linears[target_type]
                for relation_type in range(self.num_relations):
                    idx = (edge_type == int(relation_type)) & tb
                    if idx.sum() == 0:
                        continue

                    # Compute attention weights and messages
                    target_node_vec = node_inp_i[idx]
                    source_node_vec = node_inp_j[idx]

                    q_mat = q_linear(target_node_vec).view(-1, self.n_heads, self.d_k)
                    k_mat = k_linear(source_node_vec).view(-1, self.n_heads, self.d_k)
                    k_mat = torch.bmm(k_mat.transpose(1, 0), self.relation_att[relation_type]).transpose(1, 0)
                    self.res_att[idx] = (q_mat * k_mat).sum(dim=-1) * self.relation_pri[relation_type] / self.sqrt_dk
                    
                    v_mat = v_linear(source_node_vec).view(-1, self.n_heads, self.d_k)
                    res_msg[idx] = torch.bmm(v_mat.transpose(1, 0), self.relation_msg[relation_type]).transpose(1, 0)

                    attention_weights = (q_mat * k_mat).sum(dim=-1) * self.relation_pri[relation_type] / self.sqrt_dk
                    
                    # Save attention weights for each edge
                    for index, weight in zip(idx.nonzero().squeeze(1), attention_weights):
                        source_node_index = edge_index[0][index]  # Source node
                        target_node_index = edge_index[1][index]  # Target node
                        cell_gene_pair = (int(source_node_index), int(target_node_index))
                        self.edge_attention_weights[cell_gene_pair] = weight.cpu().detach().numpy()

        res = res_msg * Softmax(self.res_att, edge_index_i).view(-1, self.n_heads, 1)
        return res.view(-1, self.out_dim)

    def update(self, aggr_out, node_inp, node_type):
        aggr_out = F.gelu(aggr_out)
        res = torch.zeros(aggr_out.size(0), self.out_dim).to(node_inp.device)
        for target_type in range(self.num_types):
            idx = (node_type == int(target_type))
            if idx.sum() == 0:
                continue
            trans_out = self.drop(self.a_linears[target_type](aggr_out[idx]))
            alpha = torch.sigmoid(self.skip[target_type])
            if self.use_norm:
                res[idx] = self.norms[target_type](trans_out * alpha + node_inp[idx] * (1 - alpha))
            else:
                res[idx] = trans_out * alpha + node_inp[idx] * (1 - alpha)
        self.res = res
        return res

    def __repr__(self):
        return '{}(in_dim={}, out_dim={}, num_types={}, num_types={})'.format(
            self.__class__.__name__, self.in_dim, self.out_dim,
            self.num_types, self.num_relations)


# Define similarity function and contrastive loss function
def sim(z_i, z_j):
    z_i_norm = F.normalize(z_i, p=2, dim=1)
    z_j_norm = F.normalize(z_j, p=2, dim=1)
    return torch.mm(z_i_norm, z_j_norm.t())

def contrastive_loss(z_i, z_j, z_k, temperature):
    pos_sim = sim(z_i, z_j) / temperature
    neg_sim = sim(z_i, z_k) / temperature
    exp_pos_sim = torch.exp(pos_sim)
    exp_neg_sim = torch.exp(neg_sim).sum(dim=1, keepdim=True)
    loss = -torch.log(exp_pos_sim / (exp_pos_sim + exp_neg_sim))
    return loss.mean()

# Label smoothing loss for classification
class LabelSmoothing(nn.Module):
    def __init__(self, smoothing=0.0):
        super(LabelSmoothing, self).__init__()
        self.confidence = 1.0 - smoothing
        self.smoothing = smoothing

    def forward(self, x, target):
        logprobs = torch.nn.functional.log_softmax(x, dim=-1)
        nll_loss = -logprobs.gather(dim=-1, index=target.unsqueeze(1))
        nll_loss = nll_loss.squeeze(1)
        smooth_loss = -logprobs.mean(dim=-1)
        loss = self.confidence * nll_loss + self.smoothing * smooth_loss
        return loss.mean()

# Compute contrastive losses for different temperatures
def compute_losses(z_reference, high_score_features, low_score_features, z_other_sampled, temperature_d1=0.1, temperature_d2=0.5, temperature_d3=1.0):
    loss_d1 = contrastive_loss(z_reference, z_reference, high_score_features, temperature_d1)
    loss_d2 = contrastive_loss(high_score_features, z_reference, low_score_features, temperature_d2)
    loss_d3 = contrastive_loss(low_score_features, high_score_features, z_other_sampled, temperature_d3)
    return loss_d1, loss_d2, loss_d3


import os
import torch
import torch.nn as nn
import torch.nn.functional as F
import torch.optim as optim
import numpy as np
from tqdm import tqdm
from torch.utils.data import SubsetRandomSampler
from scipy.sparse import csc_matrix, csr_matrix
from scipy.io import mmwrite, mmread
from collections import defaultdict

# Define the PathwayUCellLoss class for computing pathway-based loss
class PathwayUCellLoss(nn.Module):
    def __init__(self, pathway_genes, gene_name_to_index, pathway_directions, device):
        super(PathwayUCellLoss, self).__init__()
        self.pathway_genes = pathway_genes
        self.gene_name_to_index = gene_name_to_index
        self.pathway_directions = pathway_directions  # +1 for positive pathways, -1 for negative pathways
        self.device = device

    def forward(self, RNA_matrix, labels, sample_type, target_sample_type):
        # Map pathways to gene indices
        pathway_indices = {
            pathway: [self.gene_name_to_index[gene] for gene in genes if gene in self.gene_name_to_index] 
            for pathway, genes in self.pathway_genes.items()
        }

        # Compute pathway scores
        pathway_scores = defaultdict(list)
        for i, (pathway, indices) in enumerate(pathway_indices.items()):
            if indices:
                pathway_expression = RNA_matrix[:, indices]
                pathway_expression_tensor = torch.tensor(
                    pathway_expression.todense(), dtype=torch.float32
                ).to(self.device)
                pathway_scores[pathway] = pathway_expression_tensor.mean(dim=1)

        # Identify underrepresented labels
        unique_labels, counts = np.unique(labels.cpu(), return_counts=True)
        underrepresented_labels = unique_labels[counts < 0.05 * len(labels)]

        # Initialize ucell_loss and indices
        ucell_loss = 0
        high_score_indices = {}
        low_score_indices = {}

        for label in underrepresented_labels:
            label_mask = np.array((labels == label).cpu())
            target_mask = (np.array(sample_type) == target_sample_type) & label_mask
            target_mask_label = target_mask[label_mask]

            # Compute total pathway scores for the label with directions applied
            total_pathway_scores = torch.sum(
                torch.stack([pathway_scores[pathway][label_mask] * self.pathway_directions[i]
                            for i, pathway in enumerate(pathway_scores)]), dim=0
            )

            # Z-score normalization
            z_scores = zscore(total_pathway_scores.cpu().numpy())

            # Classify indices based on z-scores
            high_score_indices[label] = np.where((z_scores > 1.64) & target_mask_label)[0]
            low_score_indices[label] = np.where((z_scores <= 1.64) & target_mask_label)[0]

            # Add to the loss
            ucell_loss += torch.mean(total_pathway_scores)

        # Final ucell loss tensor
        ucell_loss_tensor = -torch.tensor(ucell_loss, dtype=torch.float32) / len(underrepresented_labels) if underrepresented_labels.size > 0 else torch.tensor(0.0, dtype=torch.float32)

        return ucell_loss_tensor, high_score_indices, low_score_indices


# Graph Neural Network (GNN) class to process raw data
class GNN_from_raw(nn.Module):
    def __init__(self, in_dim, n_hid, num_types, num_relations, n_heads, n_layers, dropout=0.2, conv_name='hgt',
                 prev_norm=True, last_norm=True):
        super(GNN_from_raw, self).__init__()
        self.gcs = nn.ModuleList()
        self.num_types = num_types
        self.in_dim = in_dim
        self.n_hid = n_hid
        self.adapt_ws = nn.ModuleList()
        self.drop = nn.Dropout(dropout)
        self.embedding1 = nn.ModuleList()

        # Initialize MLP weight matrices
        for ti in range(num_types):
            self.embedding1.append(nn.Linear(in_dim[ti], 256))

        for t in range(num_types):
            self.adapt_ws.append(nn.Linear(256, n_hid))

        # Initialize graph convolution layers
        for l in range(n_layers - 1):
            self.gcs.append(
                GeneralConv(conv_name, n_hid, n_hid, num_types, num_relations, n_heads, dropout, use_norm=prev_norm))
        self.gcs.append(
            GeneralConv(conv_name, n_hid, n_hid, num_types, num_relations, n_heads, dropout, use_norm=last_norm))

    def encode(self, x, t_id):
        h1 = F.relu(self.embedding1[t_id](x))
        return h1

    def forward(self, node_feature, node_type, edge_index, edge_type):
        node_embedding = []
        for t_id in range(self.num_types):
            node_embedding += list(self.encode(node_feature[t_id], t_id))

        node_embedding = torch.stack(node_embedding)
        # Initialize result matrix
        res = torch.zeros(node_embedding.size(0), self.n_hid).to(node_feature[0].device)

        # Process each node type
        for t_id in range(self.num_types):
            idx = (node_type == int(t_id))
            if idx.sum() == 0:
                continue
            # Update result matrix
            res[idx] = torch.tanh(self.adapt_ws[t_id](node_embedding[idx]))

        # Apply dropout to the result matrix
        meta_xs = self.drop(res)
        del res

        # Initialize a container for the final layer's attention weights
        self.final_layer_attention_weights = {}
        # Iterate through graph convolution layers
        for i, gc in enumerate(self.gcs):
            meta_xs = gc(meta_xs, node_type, edge_index, edge_type)
            # Only store the attention weights from the final layer
            if i == len(self.gcs) - 1:
                self.final_layer_attention_weights = gc.edge_attention_weights

        return meta_xs


# Define a simple fully connected network (used later in the model)
class Net(nn.Module):
    def __init__(self, dim_in, dim_out):
        super(Net, self).__init__()
        self.fc1 = nn.Linear(dim_in, dim_out)

    def forward(self, x):
        x = F.relu(self.fc1(x))
        return x


# Node dimension reduction class with training loop
class NodeDimensionReduction(nn.Module):
    def __init__(self, RNA_matrix, indices, ini_p1, n_hid, n_heads, n_layers, labsm, lr, 
                 wd, device, num_types=2, num_relations=1, epochs=1):
        super(NodeDimensionReduction, self).__init__()
        self.RNA_matrix = RNA_matrix
        self.indices = indices
        self.ini_p1 = ini_p1
        self.in_dim = [RNA_matrix.shape[0], RNA_matrix.shape[1]]
        self.n_hid = n_hid
        self.num_types = num_types
        self.num_relations = num_relations
        self.n_heads = n_heads
        self.n_layers = n_layers
        self.labsm = labsm
        self.lr = lr
        self.wd = wd
        self.device = device
        self.epochs = epochs
        self.LabSm = LabelSmoothing(self.labsm)

        self.gnn = GNN_from_raw(in_dim=self.in_dim,
                                n_hid=self.n_hid,
                                num_types=self.num_types,
                                num_relations=self.num_relations,
                                n_heads=self.n_heads,
                                n_layers=self.n_layers,
                                dropout=0.3).to(self.device)

        self.optimizer = torch.optim.AdamW(self.gnn.parameters(), lr=self.lr, weight_decay=self.wd)
        self.scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(self.optimizer, 'min', factor=0.5, patience=5,
                                                                    verbose=True)

    def train_model(self, n_batch):
        print('The training process for the NodeDimensionReduction model has started. Please wait.')
        self.subgraph_attention_weights = {}

        for epoch in tqdm(range(self.epochs)):
            total_loss_kl1 = 0
            total_loss_cluster = 0
            all_cell_emb = []

            for batch_id in np.arange(n_batch):
                gene_index = self.indices[batch_id]['gene_index']
                cell_index = self.indices[batch_id]['cell_index']
                gene_feature = self.RNA_matrix[list(gene_index),]
                cell_feature = self.RNA_matrix[:, list(cell_index)].T
                gene_feature = torch.tensor(np.array(gene_feature.todense()), dtype=torch.float32).to(self.device)
                cell_feature = torch.tensor(np.array(cell_feature.todense()), dtype=torch.float32).to(self.device)

                node_feature = [cell_feature, gene_feature]
                gene_cell_sub = self.RNA_matrix[list(gene_index),][:, list(cell_index)]
                gene_cell_edge_index1 = list(np.nonzero(gene_cell_sub)[0] + gene_cell_sub.shape[1]) + list(
                    np.nonzero(gene_cell_sub)[1])
                gene_cell_edge_index2 = list(np.nonzero(gene_cell_sub)[1]) + list(
                    np.nonzero(gene_cell_sub)[0] + gene_cell_sub.shape[1])
                gene_cell_edge_index = torch.LongTensor([gene_cell_edge_index1, gene_cell_edge_index2]).to(self.device)

                edge_index = gene_cell_edge_index
                node_type = torch.LongTensor(np.array(
                    list(np.zeros(len(cell_index))) + list(np.ones(len(gene_index))))).to(self.device)
                edge_type = torch.LongTensor(np.array(list(np.zeros(np.nonzero(gene_cell_sub)[0].shape[0])) + list(
                    np.ones(np.nonzero(gene_cell_sub)[1].shape[0])))).to(self.device)
                l = torch.LongTensor(np.array(self.ini_p1)[[cell_index]]).to(self.device)

                node_rep = self.gnn.forward(node_feature, node_type,
                                            edge_index,
                                            edge_type).to(self.device)
                
                self.subgraph_attention_weights[batch_id] = self.gnn.final_layer_attention_weights
                
                cell_emb = node_rep[node_type == 0]
                gene_emb = node_rep[node_type == 1]
                
                all_cell_emb.append(cell_emb)
                
                decoder1 = torch.mm(gene_emb, cell_emb.t())
                gene_cell_sub = torch.tensor(np.array(gene_cell_sub.todense()), dtype=torch.float32).to(self.device)

                logp_x1 = F.log_softmax(decoder1, dim=-1)
                p_y1 = F.softmax(gene_cell_sub, dim=-1)

                loss_kl1 = F.kl_div(logp_x1, p_y1, reduction='mean')
                loss_cluster = self.LabSm(cell_emb, l)

                batch_loss = loss_kl1 + loss_cluster

                self.optimizer.zero_grad()
                batch_loss.backward()
                self.optimizer.step()

                total_loss_kl1 += loss_kl1.item()
                total_loss_cluster += loss_cluster.item()

            total_loss_kl1 /= n_batch
            total_loss_cluster /= n_batch

            all_cell_emb = torch.cat(all_cell_emb, dim=0)
            all_cell_pre = all_cell_emb.argmax(dim=1)

            print(f'Epoch {epoch + 1}:')
            print(f'  KL Loss: {total_loss_kl1}')
            print(f'  Cluster Loss: {total_loss_cluster}')

        print('The training for the NodeDimensionReduction model has been completed.')
        return self.gnn
        
        
import torch
import torch.nn as nn
import torch.optim as optim
import torch.nn.functional as F
import numpy as np
from tqdm import tqdm
import torch.cuda as cuda
import os
from scipy.sparse import csr_matrix, csc_matrix
from bioservices import KEGG
import h5py
import anndata as ad
import scanpy as sc
import pandas as pd
from sklearn.metrics import adjusted_rand_score
from sklearn.metrics.cluster import normalized_mutual_info_score
import random
import math
import matplotlib.pyplot as plt
from scipy.io import mmwrite


# Define the EmitGCL Class for Graph Neural Network training with pathway-based loss
class EmitGCL(nn.Module):
    def __init__(self, gnn, labsm, n_hid, n_batch, device, lr, wd, pathway_genes, pathway_directions, gene_names, sample_type, sample_list, num_epochs=1):
        super(EmitGCL, self).__init__()
        self.lr = lr
        self.wd = wd
        self.gnn = gnn
        self.n_hid = n_hid
        self.n_batch = n_batch
        self.device = device
        self.num_epochs = num_epochs
        self.net = Net(2 * self.n_hid, self.n_hid).to(self.device)
        self.gnn_optimizer = torch.optim.AdamW(self.gnn.parameters(), lr=self.lr, weight_decay=self.wd)
        self.net_optimizer = torch.optim.AdamW(self.net.parameters(), lr=1e-2)
        self.labsm = labsm
        self.LabSm = LabelSmoothing(self.labsm)
        
        # Initialize PathwayUCellLoss
        gene_name_to_index = {gene: i for i, gene in enumerate(gene_names)}
        self.pathway_loss_fn = PathwayUCellLoss(pathway_genes, gene_name_to_index, pathway_directions, device)  
        
        # Store the sample list for sequential sampling
        self.sample_list = sample_list
        
    def forward(self, indices, RNA_matrix, ini_p1, sample_type, nodes_id):
        cluster_l = list()
        cluster_kl_l = list()
        sim_l = list()
        nmi_l = list()
        ini_p1 = np.array(ini_p1)
        
        self.subgraph_attention_weights = {}
        
        for epoch in range(self.num_epochs):
            total_loss_kl1 = 0
            total_loss_cluster = 0
            total_ucell_loss = 0
            total_contrastive_loss = torch.tensor(0.0, device=self.device)  # Ensure the initial value is a tensor
            valid_contrastive_loss_count = 0
            all_cell_emb = []
            
            # Process each batch
            for batch_id in tqdm(np.arange(self.n_batch)):
                gene_index = indices[batch_id]['gene_index']
                cell_index = indices[batch_id]['cell_index']
                gene_feature = RNA_matrix[list(gene_index),]
                cell_feature = RNA_matrix[:, list(cell_index)].T
                gene_feature = torch.tensor(np.array(gene_feature.todense()), dtype=torch.float32).to(self.device)
                cell_feature = torch.tensor(np.array(cell_feature.todense()), dtype=torch.float32).to(self.device)
                node_feature = [cell_feature, gene_feature]
                gene_cell_sub = RNA_matrix[list(gene_index),][:, list(cell_index)]
                gene_cell_edge_index1 = list(np.nonzero(gene_cell_sub)[0] + gene_cell_sub.shape[1]) + list(np.nonzero(gene_cell_sub)[1])
                gene_cell_edge_index2 = list(np.nonzero(gene_cell_sub)[1]) + list(np.nonzero(gene_cell_sub)[0] + gene_cell_sub.shape[1])
                gene_cell_edge_index = torch.LongTensor([gene_cell_edge_index1, gene_cell_edge_index2]).to(self.device)
                node_type = torch.LongTensor(np.array(list(np.zeros(len(cell_index))) + list(np.ones(len(gene_index))))).to(self.device)
                edge_type = torch.LongTensor(np.array(list(np.zeros(np.nonzero(gene_cell_sub)[0].shape[0])) + list(np.ones(np.nonzero(gene_cell_sub)[1].shape[0])))).to(self.device)
                l = torch.LongTensor(np.array(ini_p1)[[cell_index]]).to(self.device)
                node_rep = self.gnn.forward(node_feature, node_type, gene_cell_edge_index, edge_type).to(self.device)
                
                self.subgraph_attention_weights[batch_id] = self.gnn.final_layer_attention_weights
                
                cell_emb = node_rep[node_type == 0]
                gene_emb = node_rep[node_type == 1]

                all_cell_emb.append(cell_emb)

                decoder1 = torch.mm(gene_emb, cell_emb.t())
                gene_cell_sub = torch.tensor(np.array(gene_cell_sub.todense()), dtype=torch.float32).to(self.device)
                logp_x1 = F.log_softmax(decoder1, dim=-1)
                p_y1 = F.softmax(gene_cell_sub, dim=-1)

                loss_kl1 = F.kl_div(logp_x1, p_y1, reduction='mean')
                total_loss_kl1 += loss_kl1
                loss_cluster = self.LabSm(cell_emb, l)
                total_loss_cluster += loss_cluster

            all_cell_emb = torch.cat(all_cell_emb, dim=0)
            all_cell_pre = all_cell_emb.argmax(dim=1)
            
            # Sequentially calculate ucell_loss and contrastive loss for sample pairs
            for i in range(len(self.sample_list) - 1):
                sample_1 = self.sample_list[i]
                sample_2 = self.sample_list[i + 1]
                
                # Create masks for the two samples in the current pair
                mask_1 = (np.array(sample_type) == sample_1)
                mask_2 = (np.array(sample_type) == sample_2)
                
                # Combine the two masks
                sample_pair_mask = mask_1 | mask_2
                
                print("RNA_matrix shape:", RNA_matrix.shape)
                
                RNA_matrix_pair = RNA_matrix[:,nodes_id][:, sample_pair_mask]
                labels_pair = all_cell_pre[sample_pair_mask]
                sample_type_pair = np.array(sample_type)[sample_pair_mask]
                sample_pair_emb = all_cell_emb[sample_pair_mask]
                
                # Compute ucell_loss and indices for high and low scores
                ucell_loss_tensor, high_score_indices, low_score_indices = self.pathway_loss_fn(
                    RNA_matrix_pair.transpose(), labels_pair, sample_type_pair, sample_2
                )
                
                total_ucell_loss += ucell_loss_tensor
                
                # Compute contrastive loss for each label in the current sample pair
                for label in high_score_indices:
                    label_mask = (labels_pair == label)
                    reference_mask = (sample_type_pair == sample_1) & label_mask.cpu().numpy()

                    z_reference = sample_pair_emb[reference_mask]
                    
                    high_score_features = sample_pair_emb[label_mask][high_score_indices[label]]
                    low_score_features = sample_pair_emb[label_mask][low_score_indices[label]]

                    # Define the mask to select cells that are not of the current label and belong to sample_2
                    other_cells_mask = (labels_pair != label) & (sample_type_pair == sample_2)
                    z_other = sample_pair_emb[other_cells_mask]

                    reference_indices = np.where(reference_mask)[0]
                    other_cells_indices = np.where(other_cells_mask)[0]
                    overlapping_indices = np.intersect1d(reference_indices, other_cells_indices)
                    print(f"Overlapping indices: {overlapping_indices}")
                    
                    # Sample other cells
                    sample_size = int(len(z_other) * 1.0)
                    sampler = torch.utils.data.SubsetRandomSampler(torch.randperm(len(z_other))[:sample_size])
                    z_other_sampled = z_other[list(sampler)]
                    
                    print(f"Label: {label}")
                    print("Values in sample_type_pair[label_mask][high_score_indices[label]]:", sample_type_pair[label_mask][high_score_indices[label]])
    
                    # Compute contrastive losses
                    loss_d1, loss_d2, loss_d3 = compute_losses(z_reference, high_score_features, low_score_features, z_other_sampled)
                    
                    print(f"\tLosses: loss_d1={loss_d1}, loss_d2={loss_d2}, loss_d3={loss_d3}")
                    
                    print("\tz_reference shape:", z_reference.shape)
                    print("\thigh_score_features shape:", high_score_features.shape)
                    print("\tlow_score_features shape:", low_score_features.shape)
                    print("\tz_other_sampled shape:", z_other_sampled.shape)

                    # Add valid losses to total
                    if not torch.isnan(loss_d1) and not torch.isnan(loss_d2) and not torch.isnan(loss_d3):
                        total_contrastive_loss += loss_d1 + loss_d2 + loss_d3
                        valid_contrastive_loss_count += 1
            
            # Average total contrastive loss if there are valid counts
            if valid_contrastive_loss_count > 0:
                total_contrastive_loss /= valid_contrastive_loss_count

            # Compute weighted losses
            weighted_contrastive_loss = total_contrastive_loss * 0.3
            weighted_ucell_loss = total_ucell_loss * 0.3
            weighted_kl_loss = (total_loss_kl1 / self.n_batch) * 10
            weighted_cluster_loss = (total_loss_cluster / self.n_batch) * 1.0
            
            # Compute total loss
            total_loss = weighted_kl_loss + weighted_cluster_loss + weighted_ucell_loss + weighted_contrastive_loss

            print(f'Epoch {epoch + 1}:')
            print(f'  KL Loss: {weighted_kl_loss}')
            print(f'  Cluster Loss: {weighted_cluster_loss}')
            print(f'  UCell Loss: {weighted_ucell_loss}')
            print(f'  Total Contrastive Loss: {weighted_contrastive_loss}')
            print(f'  Total Loss: {total_loss}')
            
            self.gnn_optimizer.zero_grad()
            self.net_optimizer.zero_grad()
            total_loss.backward()
            self.gnn_optimizer.step()
            self.net_optimizer.step()

        return self.gnn

    def train_model(self, indices, RNA_matrix, ini_p1, sample_type, nodes_id):
        self.train()
        print('The training process for the EmitGCL model has started. Please wait.')
        EmitGCL_gnn = self.forward(indices, RNA_matrix, ini_p1, sample_type, nodes_id)
        print('The training for the EmitGCL model has been completed.')
        return EmitGCL_gnn

# Prediction function for EmitGCL model
def EmitGCL_pred(RNA_matrix, EmitGCL_gnn, indices, nheads, nodes_id, cell_size, device, gene_names, node_dim_reduction_model):
    n_batch = len(indices)
    embedding = []
    l_pre = []
    EmitGCL_result = {}
        
    # Assume global counts of cells and genes are known
    num_cells = RNA_matrix.shape[1]  # Total global cell count
    num_genes = RNA_matrix.shape[0]  # Total global gene count

    # Initialize global attention matrices, creating one matrix for each head
    global_attention_matrices = {f'head-{i+1}': torch.zeros(num_genes, num_cells, device=device) for i in range(nheads)}

    with torch.no_grad():
        for batch_id in tqdm(range(n_batch)):
            gene_index = indices[batch_id]['gene_index']
            cell_index = indices[batch_id]['cell_index']
            gene_feature = torch.tensor(np.array(RNA_matrix[list(gene_index),].todense()), dtype=torch.float32).to(device)
            cell_feature = torch.tensor(np.array(RNA_matrix[:, list(cell_index)].T.todense()), dtype=torch.float32).to(device)
            node_feature = [cell_feature, gene_feature]
            
            gene_cell_sub = RNA_matrix[list(gene_index),][:, list(cell_index)]
            gene_cell_edge_index = torch.LongTensor([np.nonzero(gene_cell_sub)[0] + gene_cell_sub.shape[1], np.nonzero(gene_cell_sub)[1]]).to(device)
            
            edge_index = gene_cell_edge_index
            node_type = torch.LongTensor(np.array(list(np.zeros(len(cell_index))) + list(np.ones(len(gene_index))))).to(device)
            edge_type = torch.LongTensor(np.array(list(np.zeros(gene_cell_edge_index.shape[1])))).to(device)
            
            node_rep = EmitGCL_gnn.forward(node_feature, node_type, edge_index, edge_type).to(device)
            
            # Retrieve current subgraph's attention weights
            attention_weights = node_dim_reduction_model.subgraph_attention_weights[batch_id]
            
            # Get global cell and gene indices for the current subgraph
            global_gene_indices = indices[batch_id]['gene_index']
            global_cell_indices = indices[batch_id]['cell_index']
            
            # Iterate through each edge in the current subgraph and update global attention matrices
            for (local_gene_idx, local_cell_idx), weights_vector in attention_weights.items():
                global_gene_idx = global_gene_indices[local_gene_idx-len(indices[batch_id]['cell_index'])] 
                global_cell_idx = global_cell_indices[local_cell_idx]

                # Update global attention matrices for each head
                for head_index in range(nheads):
                    head_key = f'head-{head_index+1}'
                    weight = weights_vector[head_index]
                    global_attention_matrices[head_key][global_gene_idx, global_cell_idx] = torch.tensor(weight, device=device)


            cell_emb = node_rep[node_type == 0]
            gene_emb = node_rep[node_type == 1]

            if device.type == "cuda":
                cell_emb = cell_emb.cpu()
            embedding.append(cell_emb.detach().numpy())

            cell_pre = list(cell_emb.argmax(dim=1).detach().numpy())
            l_pre.extend(cell_pre)

    cell_embedding = np.vstack(embedding)
    cell_clu = np.array(l_pre)
    
    # Create a mapping from the original indices to the shuffled indices
    index_mapping = {original_index: node_id for original_index, node_id in enumerate(list(nodes_id))}

    # Reorder cell_clu
    sorted_cell_clu = [None] * len(cell_clu)
    for original_index, node_id in index_mapping.items():
        sorted_cell_clu[node_id] = cell_clu[original_index]

    # Reorder cell_embedding
    sorted_cell_embedding = np.zeros_like(cell_embedding)
    for original_index, node_id in index_mapping.items():
        sorted_cell_embedding[node_id, :] = cell_embedding[original_index, :]


    EmitGCL_result = {'pred_label': sorted_cell_clu, 'cell_embedding': sorted_cell_embedding, 'attention_weights': global_attention_matrices}
    return EmitGCL_result
    
    
import os
import torch
import numpy as np
import scipy
from scipy.sparse import csr_matrix
from scipy.io import mmwrite
from bioservices import KEGG
from torch.utils.data import Dataset, DataLoader
import torch.nn as nn
import torch.nn.functional as F
import torch.cuda as cuda
import argparse
import random
import math
import anndata as ad
import scanpy as sc
import matplotlib.pyplot as plt
from collections import Counter
from tqdm import tqdm
import multiprocessing as mp

# Function to save a single attention matrix in .mtx format
def save_attention_matrix_to_mtx(tensor_matrix, save_path, file_name):
    """
    Saves a PyTorch tensor as a sparse matrix in Matrix Market (.mtx) format.
    
    Parameters:
    - tensor_matrix: PyTorch tensor to be saved.
    - save_path: Directory where the matrix will be saved.
    - file_name: Name of the saved file (without extension).
    """
    # Convert PyTorch tensor to numpy array
    if isinstance(tensor_matrix, torch.Tensor):
        tensor_matrix = tensor_matrix.cpu().numpy()
    # Convert numpy array to CSR format sparse matrix
    sparse_matrix = csr_matrix(tensor_matrix)
    # Save as .mtx format using mmwrite
    mmwrite(os.path.join(save_path, f"{file_name}.mtx"), sparse_matrix)


# Function to save EmitGCL results, including attention matrices, cell embeddings, and predictions
def save_result(EmitGCL_result, output_file, attention_file, save_attention=True):
    """
    Save the results of EmitGCL including attention matrices and other data.
    
    Parameters:
    - EmitGCL_result: dict containing the results of the EmitGCL model, including 'pred_label', 'cell_embedding', and 'attention_weights'.
    - output_file: str, the file path to save the results.
    - attention_file: str, the file path to save the attention matrices.
    - save_attention: bool, whether to save the attention matrices.
    """
    
    # Save attention matrices if save_attention is True
    if save_attention:
        attention_matrices = EmitGCL_result['attention_weights']
        # Iterate through the dictionary and save each attention matrix
        for head, matrix in attention_matrices.items():
            save_attention_matrix_to_mtx(matrix, attention_file, head)
    
    # Save additional results to .npy files
    np.save(output_file + "pred.npy", EmitGCL_result['pred_label'])
    np.save(output_file + "cell_embedding.npy", EmitGCL_result['cell_embedding'])


# Set random seed for reproducibility
def set_random_seed(seed=0):
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    torch.cuda.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)
    os.environ['PYTHONHASHSEED'] = str(seed)
    torch.backends.cudnn.deterministic = True
    torch.backends.cudnn.benchmark = False


# Argument parser for training configurations
def parse_arguments(labsm=0.3, wd=0.1, lr=0.0005,
                   n_hid=104, nheads=8, nlayers=3, 
                   cell_size=30, neighbor=20, egrn=True, 
                   output_file='/users/PCON0022/duanmaoteng/Metastasis/GSE197177'):
    parser = argparse.ArgumentParser(description='Training GNN on gene cell graph')
    parser.add_argument('--labsm', type=float, default=labsm)  # Label Smoothing rate
    parser.add_argument('--wd', type=float, default=wd)  # Weight decay
    parser.add_argument('--lr', type=float, default=lr)  # Learning rate
    parser.add_argument('--n_hid', type=int, default=n_hid)  # Hidden layer size
    parser.add_argument('--nheads', type=int, default=nheads)  # Number of attention heads
    parser.add_argument('--nlayers', type=int, default=nlayers)  # Number of graph convolution layers
    parser.add_argument('--cell_size', type=int, default=cell_size)  # Cell feature dimension
    parser.add_argument('--neighbor', type=int, default=neighbor)  # Number of neighboring nodes considered
    parser.add_argument('--egrn', type=bool, default=egrn) 
    parser.add_argument('--output_file', type=str, default=output_file)
    
    # Return parsed arguments
    args = parser.parse_args([])
    return (args.output_file, args.labsm, args.lr, args.wd, 
            args.n_hid, args.nheads, args.nlayers, args.cell_size, 
            args.neighbor, args.egrn)


# Ensure the directory exists, create it if it doesn't
def ensure_dir(file_path):
    if not os.path.exists(file_path):
        os.makedirs(file_path)


# Function to fetch cancer metastasis genes from KEGG pathways
def get_cancer_metastasis_genes():
    kegg = KEGG()

    # Define cancer metastasis pathways
    cancer_metastasis_pathways = {
        'VEGF signaling pathway': 'hsa04370',
        'Focal adhesion': 'hsa04510',
        'Pathways in cancer': 'hsa05200',
        'PI3K-Akt signaling pathway': 'hsa04151',
        'MAPK signaling pathway': 'hsa04010',
        "Hippo signaling pathway": "hsa04390",
        "ECM-receptor interaction": "hsa04512"
    }

    pathway_genes = {}

    # Iterate over each pathway and fetch gene symbols
    for pathway_name, pathway_id in cancer_metastasis_pathways.items():
        # Retrieve pathway information
        pathway_info = kegg.get(pathway_id)
        parsed_pathway = kegg.parse(pathway_info)

        # Get and store the gene list
        genes = parsed_pathway['GENE']
        gene_symbols = []
        for gene_id, gene_info in genes.items():
            gene_symbol = gene_info.split(' ')[0].split(';')[0]
            gene_symbols.append(gene_symbol)

        pathway_genes[pathway_name] = gene_symbols

    return pathway_genes


# Function to perform initial clustering using Scanpy
def initial_clustering(RNA_matrix, custom_n_neighbors=None, n_pcs=40, custom_resolution=None, use_rep=None, random_seed=0):
    """
    Perform initial clustering using the Leiden algorithm with Scanpy.

    Parameters:
    - RNA_matrix: Gene expression matrix (cells x genes).
    - custom_n_neighbors: Custom number of neighbors for graph construction.
    - n_pcs: Number of principal components to use for clustering.
    - custom_resolution: Custom resolution for the Leiden algorithm.
    - use_rep: Use precomputed embeddings if available.
    - random_seed: Seed for reproducibility.

    Returns:
    - initial clusters from the Leiden algorithm.
    """
    print(
        '\tWhen the number of cells is less than or equal to 500, it is recommended to set the resolution value to 0.2.')
    print('\tWhen the number of cells is within the range of 500 to 5000, the resolution value should be set to 0.5.')
    print('\tWhen the number of cells is greater than 5000, the resolution value should be set to 0.8.')

    def segment_function(x):
        if x <= 500:
            return 0.2, 5
        elif x <= 5000:
            return 0.5, 10
        else:
            return 0.8, 15

    adata = ad.AnnData(RNA_matrix.transpose(), dtype='int32')

    # If the user did not provide a custom resolution or n_neighbors value, use the values calculated by segment_function
    if custom_resolution is None or custom_n_neighbors is None:
        resolution, n_neighbors = segment_function(adata.shape[0])
    else:
        resolution = custom_resolution
        n_neighbors = custom_n_neighbors

    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    
    # Use the user-provided embedding if available, otherwise use n_pcs
    if use_rep is not None:
        adata.obsm['use_rep'] = use_rep
        sc.pp.neighbors(adata, use_rep='use_rep', n_neighbors=n_neighbors)
    else:
        sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=n_pcs)

    sc.tl.leiden(adata, resolution, random_state=random_seed)
    return adata.obs['leiden']  # Return the clusters


# Function for subgraph sampling, used in constructing the graph for the GNN
def subgraph_1(graph, seed, n_neighbors, node_sele_prob):
    """
    Sample a subgraph for the GNN from a starting node (seed), expanding it layer by layer.

    Parameters:
    - graph: The graph in sparse format (cells x genes).
    - seed: The starting node for subgraph sampling.
    - n_neighbors: The number of neighbors to consider in each layer.
    - node_sele_prob: Probability distribution for selecting neighbors.

    Returns:
    - List of nodes selected for the subgraph.
    """
    total_matrix_size = 1 + np.cumprod(n_neighbors).sum()  # Number of nodes in the subgraph
    picked_nodes = {seed}  # One node in the batch
    last_layer_nodes = {seed}

    # Number of nodes selected in each layer. Initially, only the seed node is selected.
    to_pick = 1
    for n_neighbors_current in n_neighbors:  # Current layer neighbors
        to_pick = to_pick * n_neighbors_current
        neighbors = graph[list(last_layer_nodes), :].nonzero()[1]  # Find neighbors of last_layer_nodes

        neighbors_prob = node_sele_prob[list(neighbors)]
        neighbors = list(set(neighbors))  # Make all nodes from the last layer part of the neighbors set
        n_neigbors_real = min(to_pick, len(neighbors))  # Handle the case where the required number of neighbors is less than the actual number of neighbors
        if len(neighbors_prob) == 0:
            continue
        last_layer_nodes = set(np.random.choice(neighbors, n_neigbors_real, replace=False, p=softmax(neighbors_prob)))  # Select non-repeated nodes from neighbors
        picked_nodes |= last_layer_nodes  # Update picked_nodes as last_layer_nodes ∪ picked_nodes
    indices = list(sorted(picked_nodes - {seed}))
    return indices


# Softmax function for node selection
def softmax(x):
    return np.exp(x) / np.sum(np.exp(x))


# Batch processing functions to create subgraphs
def batch_process_1(args):
    """
    Processes one batch of data to create subgraphs for the GNN.
    
    Parameters:
    - args: Tuple containing index, node ids, RNA matrix, neighbors, and cell size for the batch.
    
    Returns:
    - Batch information including selected gene and cell indices.
    """
    i, node_ids, RNA_matrix, neighbor, cell_size = args
    gene_indices_all = []
    dic = {}

    start_index = i * cell_size
    end_index = min(start_index + cell_size, len(node_ids))

    for node_index in node_ids[start_index:end_index]:
        rna_ = RNA_matrix[:, node_index].todense()
        rna_[rna_ < 5] = 0
        gene_indices = subgraph_1(RNA_matrix.transpose(), node_index, neighbor, np.squeeze(np.array(np.log(rna_ + 1))))
        dic[node_index] = {'g': gene_indices}
        gene_indices_all.extend(gene_indices)

    gene_indices_all = list(set(gene_indices_all))
    node_indices_all = node_ids[start_index:end_index]
    h = {'gene_index': gene_indices_all, 'cell_index': node_indices_all}
    return (h, dic)

def batch_process_2(args):
    """
    Processes one batch of data to create subgraphs for the GNN (alternative method).
    
    Parameters:
    - args: Tuple containing index, node ids, RNA matrix, neighbors, and cell size for the batch.
    
    Returns:
    - Batch information including selected gene and cell indices.
    """
    i, node_ids, RNA_matrix, neighbor, cell_size = args
    gene_indices_all = []
    dic = {}

    start_index = i * cell_size
    end_index = min(start_index + cell_size, len(node_ids))

    for node_index in node_ids[start_index:end_index]:
        gene_indices = subgraph_2(RNA_matrix.transpose(), node_index, neighbor)
        dic[node_index] = {'g': gene_indices}
        gene_indices_all.extend(gene_indices)

    gene_indices_all = list(set(gene_indices_all))
    node_indices_all = node_ids[start_index:end_index]
    h = {'gene_index': gene_indices_all, 'cell_index': node_indices_all}
    return (h, dic)


# Main function to partition data and select subgraphs based on cell type
def batch_select_whole(RNA_matrix, label, neighbor=[20], cell_size=30):
    """
    Partitions the RNA matrix into batches based on sample type (e.g., 'P' for primary and 'M' for metastasis).
    
    Parameters:
    - RNA_matrix: Gene expression matrix (cells x genes).
    - label: List of labels corresponding to each cell.
    - neighbor: Number of neighboring cells to consider for subgraph creation.
    - cell_size: The size of each batch.
    
    Returns:
    - Indices of selected cells and genes for each batch.
    """
    print('Partitioning data into batches based on sample type.')
    dic = {}

    # Randomly shuffle cell indices
    shuffled_indices = np.random.choice(len(label), size=len(label), replace=False)

    # Create a dictionary mapping shuffled indices to original indices
    index_mapping = {new_idx: original_idx for original_idx, new_idx in enumerate(shuffled_indices)}

    # Get the corresponding labels based on the shuffled indices
    shuffled_labels = [label[index_mapping[i]] for i in range(len(label))]

    # Initialize dictionary to hold IDs for each sample type
    sample_type_ids = {sample_type: [] for sample_type in ['P', 'M']}
    
    # Populate sample_type_ids based on shuffled labels
    for i, l in enumerate(shuffled_labels):
        if l in sample_type_ids:
            sample_type_ids[l].append(index_mapping[i])

    # Calculate batch numbers for each sample type
    n_batches = {sample_type: math.ceil(len(ids) / cell_size) for sample_type, ids in sample_type_ids.items()}

    with mp.Pool(processes=48) as pool:
        results = []

        # Process each sample type with batch_process_1 and a progress bar
        for sample_type, ids in sample_type_ids.items():
            tasks = [(i, ids, RNA_matrix, neighbor, cell_size) for i in range(n_batches[sample_type])]
            results.extend(
                list(tqdm(pool.imap_unordered(batch_process_1, tasks), total=n_batches[sample_type], desc=f"Processing {sample_type} samples"))
            )

    # Merge results
    indices_ss = [res[0] for res in results]
    for res in results:
        dic.update(res[1])
    
    all_cell_indices = [index for batch in indices_ss for index in batch['cell_index']]
    # The returned node_ids are original indices
    return indices_ss, all_cell_indices, dic


from bioservices import KEGG
import anndata as ad
import torch
import numpy as np
import pandas as pd
import math
from scipy.sparse import csr_matrix
from scipy.io import mmwrite
import scanpy as sc
import matplotlib.pyplot as plt
from collections import Counter
from tqdm import tqdm
import multiprocessing as mp


# Function to fetch cancer metastasis-related genes from KEGG pathways
def get_cancer_metastasis_genes():
    """
    Fetches genes associated with cancer metastasis pathways from KEGG database.

    Returns:
    - pathway_genes: A dictionary containing the pathways and their associated genes.
    """
    kegg = KEGG()

    # Define cancer metastasis pathways with their KEGG IDs
    cancer_metastasis_pathways = {
        'Protein processing in endoplasmic reticulum': 'hsa04141',
        'mTOR signaling pathway': 'hsa04150',
        'NF-kappa B signaling pathway': 'hsa04064',
        'Autophagy': 'hsa04140',
        'p53 signaling pathway': 'hsa04115',
        'Apoptosis': 'hsa04210'
    }

    pathway_genes = {}

    # Retrieve genes associated with each pathway
    for pathway_name, pathway_id in cancer_metastasis_pathways.items():
        pathway_info = kegg.get(pathway_id)
        parsed_pathway = kegg.parse(pathway_info)

        # Extract gene symbols and store them in the dictionary
        genes = parsed_pathway['GENE']
        gene_symbols = []
        for gene_id, gene_info in genes.items():
            gene_symbol = gene_info.split(' ')[0].split(';')[0]
            gene_symbols.append(gene_symbol)

        pathway_genes[pathway_name] = gene_symbols

    return pathway_genes


# Loading and Preprocessing Data
directory = 'pn0325'
file_path = f"/fs/ess/PAS1475/Maoteng/Pancreas_meta/{directory}/scRNA_harmony.h5ad"

# Reading the .h5ad file using Scanpy
try:
    scRNA_harmony = ad.read_h5ad(file_path)
    print(f"Successfully read file: {file_path}")

    # Display AnnData object information
    display(scRNA_harmony)

    # Check cell metadata (obs)
    print("\nCell metadata (obs):")
    display(scRNA_harmony.obs.head())

    # Check gene information (var)
    print("\nGene information (var):")
    display(scRNA_harmony.var.head())

    # Check dimensions of expression matrix
    print("\nDimensions of expression matrix (genes, cells):")
    print(scRNA_harmony.shape)
except FileNotFoundError:
    print(f"File not found, please check the path: {file_path}")

# Remove 'M' labeled records and replace 'NA' and NaN with 'M'
scRNA_harmony = scRNA_harmony[scRNA_harmony.obs['Pid'] != 'M', :]
scRNA_harmony.obs['Pid'] = scRNA_harmony.obs['Pid'].replace(['NA', np.nan], 'M')


# Main Execution
if __name__ == "__main__":
    ignore_warnings()
    set_random_seed()

    # Define pathway directions
    pathway_directions = [1, 1, 1, 1, 1, 1]

    # Parse command-line arguments
    (output_file, labsm, lr, wd, n_hid, nheads, nlayers, cell_size, neighbor, egrn) = parse_arguments(
        egrn=True, output_file=f'/fs/ess/PAS1475/Maoteng/Pancreas_meta/{directory}/')

    # Get the genes associated with the pathways
    pathway_genes = get_cancer_metastasis_genes()

    # Set output paths
    output_file = f'/fs/ess/PAS1475/Maoteng/Pancreas_meta/{directory}/Metastasis_Result_new/'
    attention_file = f'/fs/ess/PAS1475/Maoteng/Pancreas_meta/{directory}/Metastasis_Result_new/Attention/'

    # Ensure the output directories exist
    ensure_dir(output_file)
    ensure_dir(attention_file)

    # Extract cell labels and gene/cell names
    labels = list(scRNA_harmony.obs['Pid'])
    gene_names = list(scRNA_harmony.var['genes'])
    cell_names = list(scRNA_harmony.obs_names)

    # Create sparse matrix for RNA expression data
    cell_counts_matrix = scRNA_harmony.X
    cell_counts_matrix_csr = csr_matrix(cell_counts_matrix)

    # Create RNA matrix (genes as rows, cells as columns)
    gene_cell = cell_counts_matrix_csr.transpose()

    # Assign gene and cell names to the matrix
    gene_cell.obs_names = gene_names
    gene_cell.var_names = cell_names
    RNA_matrix = gene_cell

    # Output matrix dimensions
    cell_num = RNA_matrix.shape[1]
    gene_num = RNA_matrix.shape[0]

    print(f"Cell number: {cell_num}")
    print(f"Gene number: {gene_num}")

    # Set device to CUDA if available
    device = torch.device("cuda" if cuda.is_available() else "cpu")
    print('You will use: ', device)

    # Perform clustering with Scanpy
    initial_pre = initial_clustering(RNA_matrix.copy())
    cluster_ini_num = len(set(initial_pre))
    ini_p1 = [int(i) for i in initial_pre]

    # Call the batch_select_whole function to get batches for subgraph construction
    indices, Node_Ids, dic = batch_select_whole(RNA_matrix, labels, neighbor=[neighbor], cell_size=cell_size)
    sample_type = list(np.array(labels)[Node_Ids])
    n_batch = len(indices)


# In[29]:

from sklearn.metrics import adjusted_rand_score

# Assuming the initial clusters and Seurat clusters are provided
seurat_clusters = list(scRNA_harmony.obs['seurat_clusters'])
ini_p1 = seurat_clusters  # Initial cluster assignments
ini_p1 = [int(i) for i in ini_p1]


# In[ ]:

n_batch = len(indices)

# Reduce the dimensionality of features for cell, gene, and peak data
node_model = NodeDimensionReduction(
    RNA_matrix,  # Use the RNA matrix subset
    indices,  # Use the subset indices
    ini_p1,  # Use initial clustering results
    n_hid=n_hid,  # Hidden layer size
    n_heads=nheads,  # Number of attention heads
    n_layers=nlayers,  # Number of layers in the model
    labsm=labsm,  # Label smoothing parameter
    lr=lr,  # Learning rate
    wd=wd,  # Weight decay
    device=device,  # Device (CPU or GPU)
    num_types=2,  # Number of node types (e.g., cells and genes)
    num_relations=1,  # Number of relations (e.g., between cells and genes)
    epochs=10  # Number of training epochs
)

# Train the node model
gnn = node_model.train_model(n_batch=n_batch)


# In[ ]:

# Define sample types for training
sample_list = ['M', 'P']

# Train the EmitGCL model using the node model
EmitGCL_model = EmitGCL(
    gnn=gnn,  # Use the trained graph neural network model
    labsm=labsm,  # Label smoothing coefficient
    n_hid=n_hid,  # Hidden layer size
    n_batch=n_batch,  # Number of batches
    device=device,  # Computing device
    lr=lr,  # Learning rate
    wd=wd,  # Weight decay coefficient
    pathway_genes=pathway_genes,  # Pathway-related genes
    pathway_directions=pathway_directions,  # Directions of pathways
    gene_names=gene_names,  # Gene names for the dataset
    sample_type=labels,  # Cell type labels
    sample_list=sample_list,  # List of sample types (e.g., 'M' for metastasis, 'P' for primary)
    num_epochs=10  # Number of training epochs
)

# Train the EmitGCL model
EmitGCL_gnn = EmitGCL_model.train_model(indices=indices, RNA_matrix=RNA_matrix, ini_p1=ini_p1, sample_type=labels, nodes_id=Node_Ids)


# In[ ]:

# Generate predictions using the trained EmitGCL model
EmitGCL_result = EmitGCL_pred(
    RNA_matrix, EmitGCL_gnn=EmitGCL_gnn, indices=indices, nheads=nheads,
    nodes_id=Node_Ids, cell_size=cell_size, device=device,
    gene_names=gene_names, node_dim_reduction_model=node_model)

# Assuming EmitGCL_result['attention_weights'] contains your attention matrices
attention_matrices = EmitGCL_result['attention_weights']


# In[ ]:

# Iterate over the attention matrices and save each one
for head, matrix in attention_matrices.items():
    save_attention_matrix_to_mtx(matrix, attention_file, head)

# Save the result files for Node IDs, predictions, and cell embeddings
np.save(output_file + "Node_Ids.npy", Node_Ids)
np.save(output_file + "pred.npy", EmitGCL_result['pred_label'])
np.save(output_file + "cell_embedding.npy", EmitGCL_result['cell_embedding'])
