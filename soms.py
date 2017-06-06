# -*- coding: utf-8 -*-
"""
Created on Tue Mar 21 16:46:53 2017

@author: Alexander
"""
from __future__ import division
import numpy as np

def distance(v1,v2):
    # return the distance between the two input vectors/arrays
    # scaled by size of array
    d = np.sqrt(np.sum((v1-v2)**2))/np.sqrt(v1.size)
    return d

class Som:
    def __init__(self, data, mapshape, l_rate=0.1, map_rad=1/2, weight_inds=[], do_shuffle=True):
        # data should be an n-dimensional array-like, with the first dim
        # the one to be trained along
        self.data = data
        self.norm_data = self.normalize(self.data)
        self.mapshape = mapshape    # tuple
        self.mapsize = np.prod(mapshape)
        self.weightshape = data.shape[1:]
        
        self.t = 0    #iteration number
        self.maxit = data.shape[0]
        self.map_rad = map_rad
        self.nbh_shrink_rate = abs(self.maxit/np.log(self.map_rad)) * 0.1
        
        self.init_learning_rate = l_rate
        self.do_shuffle = do_shuffle
        
        if weight_inds:
            self.weights = self.norm_data[weight_inds]
        else:
            self.weights = np.random.uniform(0, 1, np.append([self.mapsize],[self.weightshape]))
        # define the position of box centres across the map grid
        x = np.linspace(1/(2*self.mapshape[0]), 1-1/(2*self.mapshape[0]), self.mapshape[0])
        y = np.linspace(1/(2*self.mapshape[1]), 1-1/(2*self.mapshape[1]), self.mapshape[1])
        #x = np.linspace(0,1,self.mapshape[0])
        #y = np.linspace(0,1,self.mapshape[1])
        self.x_pos, self.y_pos = np.meshgrid(x,y)
        self.x_pos = self.x_pos.reshape(self.mapsize)
        self.y_pos = self.y_pos.reshape(self.mapsize)
        
    def train(self):        
        shuffled_data = self.norm_data.copy()
        if self.do_shuffle:
            np.random.shuffle(shuffled_data)
        for input_vector in shuffled_data:
            bmu = self.calc_bmu(input_vector)            
            nearby = self.find_nearby_nodes(bmu)
            
            for near in nearby:
                learning_rate = self.init_learning_rate*np.exp(-self.t/self.maxit)
                theta = np.exp(-self.distances[near]**2/(2*self.nbh_rad**2))

                self.weights[near] += theta*learning_rate*(input_vector-self.weights[near])

            self.t += 1
        self.update_bmus()
        
     
    def find_nearby_nodes(self, bmu):
        # return a list of node indices within a certain neighborhood
        self.nbh_rad = self.map_rad * np.exp(-self.t/self.nbh_shrink_rate)
        nearby = [i if ((self.x_pos[i] - self.x_pos[bmu])**2 + (self.y_pos[i] - self.y_pos[bmu])**2 <= self.nbh_rad) else [] for i in range(self.mapsize)]
        nearby = [x for x in nearby if x != []]
        return nearby
    

    def update_bmus(self):
        self.bmus = self.calc_bmus(self.norm_data)
        
        
    def calc_bmus(self, input_vector):
        #update a list of bmus for the SOM
        input_vector = self.normalize(input_vector)
        bmus = np.zeros(input_vector.shape[0])
        for i,vector in enumerate(input_vector):
            bmus[i] += [self.calc_bmu(vector)]
        return bmus
            
    def calc_bmu(self, input_vector):
        # find node with weights closest to the input vector
        distances = np.zeros(self.mapsize)
        for i in range(self.mapsize):
            distances[i] = distance(input_vector, self.weights[i])
        bmu = np.argmin(distances)
        self.distances = distances
        return bmu
    
    def is_normalized(self, input_vector):
        is_normalized = False
        if (np.amax(input_vector) < 1) and (np.amin(input_vector) > 0):
            is_normalized = True
        return is_normalized
    
    def normalize(self, input_vector):
        if self.is_normalized(input_vector):
            normed_vector = input_vector
        else:
            dmax = np.amax(input_vector)
            dmin = np.amin(input_vector)
            normed_vector = (input_vector - dmin)/(dmax - dmin)
        return normed_vector

    def err(self):
        # update and return quantization error and topographic error
        self.calc_qe()
        self.calc_te()
    
    
        return (self.qe, self.te)
        
    def calc_qe(self):
        # calculate the quantitative error: the average distance between each 
        # data vector and its bmu
        self.qes = []
        for (i,input_vector) in enumerate(self.norm_data):
            self.qes.append(distance(input_vector,self.weights[self.bmus[i]]))
        self.qe = np.mean(self.qes)
        
    def calc_te(self):
        #calculate the topographic error: the proportion of all data vectors
        # for which first and second bmus are not adjacent
        self.te = 0
        return
        