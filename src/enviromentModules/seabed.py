"""
-------------------------------------------\n
-         University of Stavanger          \n
-         Hui Cheng (PhD student)          \n
-          Lin Li (Medveileder)            \n
-      Prof. Muk Chen Ong (Supervisor)     \n
-------------------------------------------\n
Any questions about this code,
please email: hui.cheng@uis.no \n
wave direction is x+
"""
import numpy as np

water_depth=-1.0 # m
buffer_layer=0.1

def reaction_element(node_position,list_element,dt):
    force=np.zeros_like(node_position)
    for each_element in list_element:
        count=0
        for i in each_element.node_index:
            if node_position[i][2]<water_depth:
                count+=1
        upward_force=count/len(each_element.node_index)*each_element.weight
        for i in each_element.node_index:
            force[i]=[0,0,upward_force]
    
        
        
