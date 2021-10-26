# import libraries
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
import matplotlib.path as mpath
from matplotlib.patches import Circle
import torch
from torch import nn
from torch.autograd import Variable

class Snapshot:
    def __init__(self, particles, box):
        self.particles = particles
        self.box = box
    def get_npart(self):
        return len(self.particles)
    def get_volume(self):
        return self.box[0]*self.box[1]*self.box[2]
    def get_density(self):
        return float(self.get_npart())/self.get_volume()



def read_snapshot(file_in):
    # read number of particles
    Npart = int(file_in.readline())
    # read box
    temp = [float(n) for n in file_in.readline().split()]
    box = [temp[0], temp[1], temp[2]]
    # read coordinates
    particles = []
    for i in range(Npart):
        temp = [n for n in file_in.readline().split()]
        xyz = [float(temp[1]), float(temp[2]), float(temp[3])]
        particles.append(xyz)
    # convert to np arrays
    box = np.array(box)
    particles = np.array(particles)
    # create snapshot
    snap = Snapshot(particles, box)
    return snap

def rotate_snap(snap, angle=45):
    theta = angle/180.0*np.pi
    rot_matrix = np.array([[np.cos(theta), -np.sin(theta), 0.],[np.sin(theta), np.cos(theta), 0.], [0., 0., 1.]])
    box = snap.box
    half_box = box*0.5
    par = snap.particles
    # copy original snap
    v1 = np.array([box[0], 0, 0])
    v2 = np.array([box[0], box[1], 0])
    v3 = np.array([0, box[1], 0])
    v4 = np.array([2*box[0], 0, 0])
    v5 = np.array([2*box[0], box[1], 0])
    v6 = np.array([box[0], 2*box[1], 0])
    v7 = np.array([2*box[0], 2*box[1], 0])
    v8 = np.array([0, 2*box[1], 0])
    par_larger = np.concatenate((par, par+v1, par+v2, par+v3, par+v4, par+v5, par+v6, par+v7, par+v8), axis=0)
    par_larger = par_larger - 1.5*np.array([box[0], box[1], 0])
    # rotate snap
    par_larger = np.dot(par_larger, rot_matrix)
    par_rot = np.array([p for p in par_larger if (p[0]<half_box[0] and p[0]>-half_box[0] and p[1]<half_box[1] and p[1]>-half_box[1])])
    half_box[2] = 0
    par_rot = par_rot + half_box
    rot_snap = Snapshot(par_rot, box)
    return rot_snap

def diffraction_pattern3d(snap, size=64):
    par = snap.particles
    box = snap.box
    half_box = box*0.5
    par = par - half_box
    maxb = np.max(half_box[:2])
    maxq = int(size/2)
    nx = np.pi*np.arange(-maxq, maxq+1, 1)/maxb
    ny = np.pi*np.arange(-maxq, maxq+1, 1)/maxb
    nz = np.pi*np.arange(-maxq, maxq+1, 1)/box[2]
    qs = np.array([[x,y,z] for x in nx for y in ny for z in nz])
    dot = np.dot(qs,par.T)
    realpart = np.sum(np.cos(dot), axis=1)
    impart = np.sum(np.sin(dot), axis=1)
    diffraction = np.sqrt(realpart*realpart + impart*impart)/len(par)
    diffraction = diffraction.reshape((size+1,size+1,size+1))
    return diffraction


def remove_center_pick(diff, stride=0):
    # removes the center pick in the diffraction pattern
    i = int(len(diff)/2)
    j = int(len(diff[0])/2)
    k = int(len(diff[0][0])/2)
    diff_new = np.copy(diff)
    diff_new[i-stride:i+stride+1, j-stride:j+stride+1, k-stride:k+stride+1] = 0.0
    return diff_new

def split_snap(snap):
    box = snap.box
    par = snap.particles
    par1 = np.array([p[:2] for p in par if p[2]<box[2]/3.])
    par2 = np.array([p[:2] for p in par if p[2]>box[2]/3. and p[2]<2.*box[2]/3])
    par3 = np.array([p[:2] for p in par if p[2]>2.*box[2]/3.])
    boxlayer = box[:2]
    layer1 = Snapshot(par1, boxlayer)
    layer2 = Snapshot(par2, boxlayer)
    layer3 = Snapshot(par3, boxlayer)
    return layer1, layer2, layer3

def diffraction_pattern2d(snap, size=64):
    par = snap.particles
    box = snap.box
    half_box = box*0.5
    par = par - half_box
    maxb = np.max(half_box)
    maxq = int(size/2)
    nx = np.arange(-maxq, maxq+1, 1)
    ny = np.arange(-maxq, maxq+1, 1)
    ns = np.array([[x,y] for x in nx for y in ny])
    qs = np.pi*ns/maxb
    dot = np.dot(qs,par.T)
    realpart = np.sum(np.cos(dot), axis=1)
    impart = np.sum(np.sin(dot), axis=1)
    diffraction = np.sqrt(realpart*realpart + impart*impart)/len(par)
    diffraction = diffraction.reshape((size+1,size+1))
    return diffraction

def process_snap(snap, flatten=False):
    """ Compute and process the diffraction pattern of a snapshot
    Args:
        snap (Snapshot): a snapshot to analyze.
        flatten (bool): if True the returned diffraction pattern is a 1d array.
                        Default is False.
    Returns:
        diffraction (np.array): the 33x33 processed diffraction pattern.
    """
    diffraction = diffraction_pattern3d(snap, size=64)
    pool = nn.MaxPool3d(kernel_size=4, stride=4, padding=2)
    diffraction = pool(get_var_from_np(np.array([diffraction])))
    diffraction = np.array(diffraction.squeeze())
    if flatten:
        diffraction = diffraction.flatten()
    return diffraction



def get_var_from_np(np_array, cuda=False, requires_grad=False):
    temp = Variable(torch.from_numpy(np_array), requires_grad=requires_grad).type(torch.FloatTensor)
    if cuda: temp = temp.cuda()
    return temp
