import os
import json
import numpy as np
import vtk

class netting:
    """
    For Screen hydrodynamic models, the forces on netting are calculated based on individual a panel section of netting.
    The twines and knots in the net panel are considered as an integrated structure. In this module, the net panel is defined by
    three nodes because three (non-collinear) points can determine a unique plane in Euclidean geometry.
    In practice, the force is usually decomposed into two components: drag force F_D and lift force F_L (Cheng et al., 2020).
    """

    def __init__(self, model_index, hydro_element, solidity:float, dw0=0.0):
        """
        :param model_index: [string] Unit: [-]. To indicate the model function, e.g.: 'S1', 'S2', 'S3'.
        :param hydro_element: [[list]] Unit: [-]. A python list to indicate how the net panel are connected. e.g.:[[p1,p2,p3][p2,p3,p4,p5]...]. If the input net panel contains 4 nodes, it will automaticly decomposed to 3 node net panel.
        :param solidity: [float] Unit: [-]. The solidity of netting.
        :param dw0: [float] Unit: [m]. The diameter of the physical net twines. It is used for the hydrodynamic coefficients.
        """
        self.modelIndex = str(model_index)
        self.dw0 = dw0
        self.sn = solidity
        self.element=hydro_element
        # 
        self.area=0.0
        self.vectorN=np.array([0,0,0], dtype=float)
        self.center=np.array([0,0,0], dtype=float)
        self.cd=0.0
        self.cl=0.0
        self.hydro_dynamic_forces = np.zeros_like(self.center)
        self.hydro_static_forces = np.zeros_like(self.center)
        self.hydro_total_forces = np.zeros_like(self.center)

    def __str__(self):
        """Print information of the present object."""
        s0 = "Screen model"
        s1 = "The model index is " + str(self.modelIndex) + "\n"
        return s0 + s1 
    
    def update_state(self,position:np.array):
        self.get_area(position)
        self.get_center(position)
        self.get_normal_vector(position)
    
    
    def get_area(self,position:np.array):        
        for i in range(len(self.element)-2):
            a1 = position[self.element[0]] - position[self.element[i+1]]
            a2 = position[self.element[0]] - position[self.element[i+2]]
            self.area += 0.5 * np.linalg.norm(np.cross(a1, a2))
    
    def get_center(self,position:np.array):
        for i in self.element:
            self.center+=position[i]/len(self.element)
            
    def get_normal_vector(self,position:np.array):
        if len(self.element)==3:
            a1 = position[self.element[0]] - position[self.element[1]]
            a2 = position[self.element[-1]] - position[self.element[-2]]
            self.vectorN=np.cross( a1, a2) / np.linalg.norm(np.cross(a1, a2))
        elif len(self.element)==4: # need to find an advanced way to calculate the normal vector of quad.
            a1 = position[self.element[0]] - position[self.element[1]]
            a2 = position[self.element[0]] - position[self.element[2]]
            self.vectorN=np.cross( a1, a2) / np.linalg.norm(np.cross(a1, a2))
        



def make_trigangle(point_list,triangle_list):
    # ug is the container to store the grid object
    ug = vtk.vtkUnstructuredGrid()
    points = vtk.vtkPoints()
    triangles = vtk.vtkTriangle()  # 3-point elements
    pixel=vtk.vtkPixel() # 4-point elements
    
    # add points   
    for i in range(len(point_list)):
        points.InsertNextPoint(point_list[i])
    ug.SetPoints(points)
    
    # add elements
    for j in range(len(triangle_list)):
        if len(triangle_list[j])==3: # add 3-point elements
            for i in range(3): 
                triangles.GetPointIds().SetId(i, triangle_list[j][i])
            ug.InsertNextCell(triangles.GetCellType(), triangles.GetPointIds())    
        elif len(triangle_list[j])==4:  # add 4-point elements
            for i in range(4): 
                pixel.GetPointIds().SetId(i,triangle_list[j][i])
            ug.InsertNextCell(pixel.GetCellType(), pixel.GetPointIds())    
    return ug

def make_surf_vector(point_list,line_list):
    # ug is the container to store the grid object
    ug = vtk.vtkUnstructuredGrid()
    points = vtk.vtkPoints()
    line = vtk.vtkLine()   # 2-point elements
    
    pn=len(point_list)
   # add points
    for i in range(pn):
        points.InsertNextPoint(point_list[i])
    ug.SetPoints(points)
    
    # add elements
    for j in range(len(line_list)):
        for i in range(2):
            line.GetPointIds().SetId(i, line_list[j][i])
        ug.InsertNextCell(line.GetCellType(), line.GetPointIds())    
    return ug



def save_surf(point_dic,elements,component:str):
    i=0
    writer = vtk.vtkUnstructuredGridWriter()  
    if os.path.exists(os.path.join(os.getcwd(), 'VTKsurfs_'+component)):
        print('\nVTKsurfs_'+component+' exists in the current folder\nPlease remove it and try again.')
        exit()
    else:
        os.mkdir(os.path.join(os.getcwd(), 'VTKsurfs_'+component))
    
    
    for k in point_dic.keys():
        uGrids=make_trigangle(point_dic[k],elements[component])
        writer.SetFileName(os.path.join('VTKsurfs_'+component,'surfTK_'+component+str(i)+'.vtk'))
        writer.SetInputData(uGrids)
        writer.Write()
        i+=1  
    
    
    # to plot the normal vector of netting
    # list_net=[]
    # for each in element[component]:
    #     list_net.append(netting('S3',each,0.2,0.01))
    
    # for each in list_net:
    #     each.update_state(np.array(point_dic[k]))
    # p_center=[]
    # l_center=[]
    # for index, each in enumerate(list_net):
    #     p_center.append(each.center)
    #     p_center.append(each.vectorN+each.center)
    #     l_center.append([2*index,2*index+1])

    # for k in point_dic.keys():
    #     uGrids=make_surf_vector(p_center,l_center)
    #     writer.SetFileName(os.path.join('VTKsurfs_'+component,'lineVTK_'+component+str(i)+'.vtk'))
    #     writer.SetInputData(uGrids)
    #     writer.Write()
    #     i+=1

    
    

if __name__ == '__main__':
    # read nodes' positions
    with open(os.path.join('Nodes_position.json'), 'r') as f:
       point_position = json.load(f)
    f.close()
    # read element 
    with open(os.path.join('element.json'), 'r') as f:
       element = json.load(f)
    f.close()

    save_surf(point_position,element,'netting')
    save_surf(point_position,element,'pontoonSurface')
    # save_surf(point_position,element,'surf_quad')
    # save_surf(point_position,element,'surf_tri')
    