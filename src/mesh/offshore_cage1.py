'''
/--------------------------------\
|    University of Stavanger     |
|           Hui Cheng            |
\--------------------------------/
Any questions about this code, please email: hui.cheng@uis.no
The center of the floating collar is (0,0,0)
Fish cage is along the Z- direction
Z=0 is the free surface
Z<0 is the water zone
Z>0 is the air zone
The sinkers are attached to the floating collar
'''

from killSalomeWithPort import killMyPort
from salome.smesh import smeshBuilder
import SMESH
import salome

import os
import sys
import json
import numpy as np
from numpy import pi



# #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Parameter template start
# #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# default value
parameters={
'Environment':
{
    'current':[0.5,0,0],   # current velcity [m/s]
    'fluidDensity':1000,  # density of fresh water/sea water [kg/m3]
    'waterDepth':60,
    'draft':0.36,
    'model_ratio':1,
},
'NetShape':
  {
    'origin':[[0,0,0]],
    'cageDiameter':0.98, # [m] diameter of the fish cage
    'cageHeight':0.28,  # [m] height of the fish cage
    'cageConeHeight':0.02,
    'elementOverCir':16,
    'elementOverHeight':2,
    'elementOverCone':2,
  },
'Frame':
    {
    'Diameter':1,  #[m]
    'height':0.3,
    },    
'Pontoon':
    {
        'diameter':0.1,
        'height':0.05,
        'length':0.06,
    },
'Net':
  {
    'nettingType':'square',
    'Sn': 0.25,   # solidity ratio
    'twineDiameter':0.6e-3, # [m]the twine diameter of the physical net
    'meshLength': 8e-3, # [m]the half mesh length
    'netYoungmodule':2e9, # [Pa]
    'netRho':1140.0, #[kg/m3] density of the net material
  },
}


# #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Parameter template finish
# #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

path_to_setting=os.path.join(os.getcwd(),'setting.json')
if os.path.isfile(path_to_setting):
    print('\nYes, the setting file exists. The default parameters will be overwritted. \n')
    with open(path_to_setting) as json_file:
        parameters = json.load(json_file)
    json_file.close()
else:
    print('\n The default parameters will be used. \n')


floater_center = parameters['NetShape']['origin']
cr = float(parameters['NetShape']['cageDiameter']/2.0)
cage_height = float(parameters['NetShape']['cageHeight'])
cage_cone_height = float(parameters['NetShape']['cageConeHeight'])

NT = int(parameters['NetShape']['elementOverCir'])     # Number of the nodes in circumference
NN = int(parameters['NetShape']['elementOverHeight'])  # number of section in the height, thus, the nodes should be NN+1
BN = int(parameters['NetShape']['elementOverCone'])    # number of section along the cone, thus, the nodes should be NN+1


element={}
point={}

# Step 1. Net 
# skip, net on frame
point['initial']=[]

element['netting']=[]

#  step 2 frame
###  point
cr=float(parameters['Frame']['Diameter']/2.0)

for i in range(16):
    for j in range(2): 
        point['initial'].append(
            [cr * np.cos(i * 2 * pi / float(NT)),
             cr * np.sin(i * 2 * pi / float(NT)),
              - j * cage_height / float(NN)])

    for j in range(2): 
        point['initial'].append(
            [cr * ((BN - j) / BN) * np.cos(i * 2 * pi / float(NT)),
             cr * ((BN - j) / BN) * np.sin(i * 2 * pi / float(NT)),
             - cage_height - j * (cage_cone_height) / float(BN)])
        
point['initial'].append([0, 0, -cage_cone_height-cage_height])  # the last point should be at the cone tip
point['initial'].append([0, 0, 0])  # the last point should be at the cone tip

# triangular element
for i in range(NT-1):
    element['netting'].append([len(point['initial'])-2,i*4+3,     (i+1)*4+3]) # bottom tri net
element['netting'].append(    [len(point['initial'])-2,(i+1)*4+3, 3])    # bottom tri net

# quad element
for j in range(NT-1):
    for i in range(3):
        element['netting'].append([4*j+i,4*j+i+1,4*j+i+4,4*j+i+5])   
for i in range(3):
    element['netting'].append(    [4*j+i+4,4*j+i+5, i,   i+1     ]) 
      
print(element['netting'])    


### element
element['thickColumn']=[]
element['thinColumn']=[]
element['inclined']=[]  

element['arcpipe_up']=[]     # thick/upper and lower
element['arcpipe_mid']=[]    # thin/middle
element['arcpipe_bottom']=[]  # thick/upper and lower

element['brace']=[] 
element['brace_up']=[]  


element['thickColumn'].append([64,65])
for i in range(8):
    element['thickColumn'].append([2+i*8,66+i*34])
    for j in range(2):
        element['thickColumn'].append([i*8+j,1+i*8+j])
        element['thinColumn'].append([4+i*8+j,4+1+i*8+j])
        
for i in range(16-1):
    element['arcpipe_up'].append(    [i*4,4+i*4])
    element['arcpipe_mid'].append(   [1+i*4,5+i*4])
    element['arcpipe_bottom'].append([2+i*4,6+i*4])

element['arcpipe_up'].append(    [0,4+i*4])    
element['arcpipe_mid'].append(   [1,5+i*4])
element['arcpipe_bottom'].append([2,6+i*4])
  
  
  
for i in range(16):
    element['brace'].append([64,3+i*4])         
    element['brace'].append([3+i*4,2+i*4])

element['brace_up'].append([65,0])
element['brace_up'].append([65,8])
element['brace_up'].append([65,24])
element['brace_up'].append([65,32])
element['brace_up'].append([65,40])
element['brace_up'].append([65,56])


for i in range(8-1):
    element['inclined'].append([2+8*i,4+8*i])
    element['inclined'].append([4+8*i,10+8*i])
element['inclined'].append([58,60])
element['inclined'].append([60,2])



# step 3 pontoon
#### point
Column_height=float(parameters['Frame']['height']) #0.3m
depth1=-float(parameters['Pontoon']['length']) -Column_height
depth2=-float(parameters['Pontoon']['length']) -Column_height- float(parameters['Pontoon']['height'])
pr=float(parameters['Pontoon']['diameter'])/2.0
  
for i in range(8):# side pontoon
    point['initial'].append(
                [cr * np.cos(i * 2 * pi / float(8)),
                 cr * np.sin(i * 2 * pi / float(8)),
                  -Column_height])
    point['initial'].append(
                [cr * np.cos(i * 2 * pi / float(8)),
                 cr * np.sin(i * 2 * pi / float(8)),
                  depth2])
    for j in range(16):
        point['initial'].append(
                [cr * np.cos(i * 2 * pi / float(8))+pr* np.cos(j * 2 * pi / float(16)),
                 cr * np.sin(i * 2 * pi / float(8))+pr* np.sin(j * 2 * pi / float(16)),
                  depth1])

        point['initial'].append(
                [cr * np.cos(i * 2 * pi / float(8))+pr* np.cos(j * 2 * pi / float(16)),
                 cr * np.sin(i * 2 * pi / float(8))+pr* np.sin(j * 2 * pi / float(16)),
                  depth2 ])
# central pontoon
for j in range(16):        
        point['initial'].append(
                [pr* np.cos(j * 2 * pi / float(16)),
                 pr* np.sin(j * 2 * pi / float(16)),
                  depth1 ])
        point['initial'].append(
                [pr* np.cos(j * 2 * pi / float(16)),
                 pr* np.sin(j * 2 * pi / float(16)),
                  depth2 ])
point['initial'].append([0,  0,  depth2 ])        
        
#### element        
element['pontoonSurface']=[]
# tri
# side pontoon
# triangular element
for i in range(8):
    for j  in range(16-1):    
        element['pontoonSurface'].append([66+34*i+0,70+2*j+34*i+0,68+2*j+34*i+0])# top
        element['pontoonSurface'].append([66+34*i+1,68+2*j+34*i+1,70+2*j+34*i+1])# bottom
    element['pontoonSurface'].append([66+34*i+0,68+34*i+0,    70+2*j+34*i+0])    # top
    element['pontoonSurface'].append([66+34*i+1,70+2*j+34*i+1,68+34*i+1])    # bottom
# quad
for i in range(8):
    for j  in range(16-1):    
        element['pontoonSurface'].append([68+2*j+34*i,69+2*j+34*i, 70+2*j+34*i,71+2*j+34*i])
    element['pontoonSurface'].append([    70+2*j+34*i,71+2*j+34*i, 68+34*i,    69+34*i    ])    
        
# central pontoon
# triangular element
for j  in range(16-1):    
    element['pontoonSurface'].append([64, 340+2*j, 338+2*j])    # top
    element['pontoonSurface'].append([370,339+2*j, 341+2*j])   # bottom
element['pontoonSurface'].append([    64, 338    ,340+2*j])    # top
element['pontoonSurface'].append([    370,341+2*j, 339    ]) # bottom   

# quad
for j  in range(16-1):    
    element['pontoonSurface'].append([338+2*j,339+2*j,340+2*j,341+2*j])
element['pontoonSurface'].append([    340+2*j,341+2*j,338    ,339    ])    
        

# used for remove free nodes, 
# pure structure function
element['Stru_rigid']=[]
# side pontoon
for i in range(8):
    element['Stru_rigid'].append([67+i*34,68+i*34])
    for j in range(16):
        element['Stru_rigid'].append([67+i*34,69+j*2+i*34])
        element['Stru_rigid'].append([68+i*34,70+j*2+i*34])
        element['Stru_rigid'].append([69+j*2+i*34,70+j*2+i*34])
    for j in range(16-1):
        element['Stru_rigid'].append([69+j*2+i*34,71+j*2+i*34])
        element['Stru_rigid'].append([70+j*2+i*34,72+j*2+i*34])
    element['Stru_rigid'].append([69+i*34,71+j*2+i*34])
    element['Stru_rigid'].append([70+i*34,72+j*2+i*34])
    
# central pontoon

element['Stru_rigid'].append([65,371])
for i in range(16):
    element['Stru_rigid'].append([65, 339+i*2])
    element['Stru_rigid'].append([371,340+i*2])
    element['Stru_rigid'].append([339+i*2,340+i*2])
for i in range(16-1):       
    element['Stru_rigid'].append([339+i*2,341+i*2])
    element['Stru_rigid'].append([340+i*2,342+i*2])
element['Stru_rigid'].append([339,341+i*2])
element['Stru_rigid'].append([340,342+i*2])    
        
element['Stru_rigid']=(np.array(element['Stru_rigid'])-1).tolist()




    
#Step 4 mooring system
# 30 segments for each line
cage_draft=parameters['Environment']['draft']
P_draft=np.array(point['initial'])+[0,0,0.41-cage_draft]
point['initial']=P_draft.tolist()

mooring_point=[[3.4,1,-1],
               [-3.4,1,-1],
               [3.4,-1,-1],
               [-3.4,-1,-1]]
element['mooringline']=[]
N_point_cage=len(point['initial'])
mooring_pair=[[0,11],[1,27],[2,59],[3,43]]
num_per_line=5
for pair in mooring_pair:
    dxyz=(np.array(mooring_point[pair[0]])-np.array(point['initial'][pair[1]-1]))/float(num_per_line)
    for i in range(num_per_line):
        point['initial'].append((np.array(point['initial'][pair[1]-1])+(i+1)*dxyz).tolist())


for j, pair in enumerate(mooring_pair):
    element['mooringline'].append([pair[1],N_point_cage+1+j*num_per_line])
    for i in range(num_per_line-1):
        element['mooringline'].append([N_point_cage+1+i+j*num_per_line,N_point_cage+2+i+j*num_per_line])
        
element['mooringline']=(np.array(element['mooringline'])-1).tolist()



# enlarge?
lam=float(parameters['Environment']['model_ratio'])
P_draft=np.array(point['initial'])*lam
point['initial']=P_draft.tolist()

# move towards x+ for wave theory
P_draft=np.array(point['initial'])+[4.0*lam,0,0]
point['initial']=P_draft.tolist()


# ###############salome########################
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# the below is the command in the Mesh, Salome.
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# ###############salome########################

salome.salome_init()
theStudy = salome.myStudy

smesh = smeshBuilder.New(theStudy)
Mesh_1 = smesh.Mesh()

# add the pints into geometry
for each_node in point['initial']:
    nodeID = Mesh_1.AddNode(float(each_node[0]), float(each_node[1]), float(each_node[2]))


for each_type in element.keys():
    for each_ele in element[each_type]:
        if len(each_ele)<=2:
            edge = Mesh_1.AddEdge([int(each_ele[0]+1), int(each_ele[1]+1)])

isDone = Mesh_1.Compute()

# # naming  the group
# # naming the node
# GROUP_NO
allnodes = Mesh_1.CreateEmptyGroup(SMESH.NODE, 'allnodes')
nbAdd = allnodes.AddFrom(Mesh_1.GetMesh())
smesh.SetName(allnodes, 'allnodes')

# generate the name for each node to assign the hydrodynamic forces.
for i in range(1, len(point['initial']) + 1):
    node1 = Mesh_1.CreateEmptyGroup(SMESH.NODE, 'node%s' % i)
    nbAdd = node1.Add([i])
    smesh.SetName(node1, 'node{}'.format(str(i)))
    
fixed = Mesh_1.CreateEmptyGroup(SMESH.NODE, 'fixed')
# nbAdd = fixed.Add([374,377,380,383])
nbAdd = fixed.Add([(i+1)*num_per_line+N_point_cage for i in range(4)])
smesh.SetName(fixed, 'fixed')

balance = Mesh_1.CreateEmptyGroup(SMESH.NODE, 'balance')
nbAdd = balance.Add([68,102,136,170,204,238,272,306,371])
smesh.SetName(balance, 'balance')

# # GROUP_MA

index_start=0
for each_type in element.keys():
    if len(element[each_type][0])<=2: # the first element is line type
        group=Mesh_1.CreateEmptyGroup(SMESH.EDGE, each_type)
        nbAdd = group.Add([1+i+index_start for i in range(len(element[each_type]))])
        smesh.SetName(group, each_type)
        index_start+=len(element[each_type])


# # give a name to the mesh
meshname = 'Offshore1.med'

try:
    Mesh_1.ExportMED(os.path.join(os.getcwd(), meshname))
    pass
except:
    print('ExportMED() failed. Invalid file name?')

killMyPort(os.getenv('NSPORT'))

for k in element.keys():
    print(k)


with open(os.path.join(os.getcwd(), 'element.json'), 'w') as json_file:
    json.dump(element, json_file,indent=4)
json_file.close()

with open(os.path.join(os.getcwd(), 'Nodes_position.json'), 'w') as json_file:
    json.dump(point, json_file,indent=4)
json_file.close()
if __name__ == '__main__':
    pass

