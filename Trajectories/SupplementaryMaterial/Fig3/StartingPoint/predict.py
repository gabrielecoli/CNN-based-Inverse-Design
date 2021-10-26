# import libraries
import os
import numpy as np
import sys
import linecache

import torch
from torch import nn
import torch.nn.functional as F
from torch.autograd import Variable


from utils import Snapshot
from utils import read_snapshot
from utils import remove_center_pick 
from utils import process_snap
 
 
def get_var_from_np(np_array, cuda=False, requires_grad=False):
     temp = Variable(torch.from_numpy(np_array), requires_grad=requires_grad).type(torch.FloatTensor)
     if cuda: temp = temp.cuda()
     return temp

class Net(nn.Module):
     def __init__(self):
         super(Net, self).__init__()
         self.conv1 = nn.Sequential(
             nn.Conv2d(1, 9, kernel_size=4, stride=1, padding=1),
             nn.ReLU(),
             nn.MaxPool2d(kernel_size=2, stride=2))
         self.conv2 = nn.Sequential(
             nn.Conv2d(9, 4, kernel_size=3, stride=1, padding=1),
             nn.ReLU(),
             nn.MaxPool2d(kernel_size=2, stride=2))
         #self.drop_out = nn.Dropout()
         self.fc1 = nn.Linear(8 * 8 * 4, 20)
         self.fc2 = nn.Linear(20, 6)

     def forward(self, x):
         out = self.conv1(x)
         out = self.conv2(out)
         out = out.reshape(out.size(0), -1)
         #out = self.drop_out(out)
         out = F.relu(self.fc1(out))
         out = self.fc2(out)
         return out

def predict_from_file(file_name, model):
     # read snapshot from file
     file_in = open(file_name, 'r')
     snap = read_snapshot(file_in)
     file_in.close()
     # get model output
     diff = process_snap(snap, flatten=False)
     out = net(get_var_from_np(np.array([[diff]])))
     # get class label and vector of probabilities
     _, label = torch.max(out, 1)
     label = label.item()
     prob = F.softmax(out, dim=1).detach().numpy().squeeze()
     return label, prob
     
     
     
# classes dictionary
class_list = ['fluid', 'square', 'hexagonal', 'qc12', 'qc10', 'qc18']
class_dic = {'fluid':0, 'square':1, 'hexagonal':2, 'qc12':3, 'qc10':4, 'qc18':5}

# load model
model_path = './convNet/net.pyt'
net = Net()
net.load_state_dict(torch.load('./convNet/net.pyt'))

nJump = int(linecache.getline('in.dat',5).split()[-1])
nTotal = int(float(int(linecache.getline('in.dat',4).split()[-1]))/nJump)

mean_prob_dict = {'fluid':0, 'square':0, 'hexagonal':0, 'qc12':0, 'qc10':0, 'qc18':0}
for snapIndex in range(nTotal):
        # predict example
        file_name = 'conf' + '/conf_' + str((snapIndex+1)*nJump).zfill(6) + '.dat'
        label, prob = predict_from_file(file_name, net)
        for phase in mean_prob_dict.keys():
                mean_prob_dict[phase] += (prob[class_dic[phase]]/nTotal)
# print prediction
fileOut = open('probs_cnn.dat','w')

phase = 'fluid'
fileOut.write(str(mean_prob_dict[phase]) + ' ')
phase = 'square'
fileOut.write(str(mean_prob_dict[phase]) + ' ')
phase = 'hexagonal'
fileOut.write(str(mean_prob_dict[phase]) + ' ')
phase = 'qc12'
fileOut.write(str(mean_prob_dict[phase]) + ' ')
phase = 'qc10'
fileOut.write(str(mean_prob_dict[phase]) + ' ')
phase = 'qc18'
fileOut.write(str(mean_prob_dict[phase]) + ' ')

fileOut.write("\n")
fileOut.close()

