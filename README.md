# EXP2SL: a Machine Learning Framework for Cell-Line Specific Synthetic Lethality Prediction

# Requirements
Pytorch 1.0.1

# How to run
1. Download L1000 expression profiles (i.e., input features) from . Put it in folder L1000.
2. To reproduce our results, run the following command:
   <code>python train.py [cell_line] [test scenario] [semi-supervised loss weight] [dnn layer] [hidden dimension] [l2 weight]</code>
   
   Here, the range of cell line is [A549, HT29, A375]. The range of test scenario is [edge, gene].
# Contacts
If you have any questions or comments, please feel free to email Fangping Wan (wanfangping92[at]gmail[dot]com) and/or Jianyang Zeng (zengjy321[at]tsinghua[dot]edu[dot]cn).
