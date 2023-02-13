import os
import json
import numpy as np
import vtk


def make_sea(point_list,element_list):
    # ug is the container to store the grid object
    ug = vtk.vtkUnstructuredGrid()
    points = vtk.vtkPoints()
    pixel=vtk.vtkPixel() # 4-point elements
    
    # add points   
    for i in range(len(point_list)):
        points.InsertNextPoint(point_list[i])
    ug.SetPoints(points)
    
    # add sea surface
    for j in range(len(element_list)):
        for i in range(4): 
            pixel.GetPointIds().SetId(i,element_list[j][i])
        ug.InsertNextCell(pixel.GetCellType(), pixel.GetPointIds())    
       
    return ug

def prepare_pe(elevation:list,x_range:list,y_range:list):
    """prepare the points and elements

    Args:
        elevation (list): [a list of elevation] 
        x_range (list): [x_min, x_max]
        y_range (list): [y_min, y_max]
        
    """
    number_of_elevation=len(elevation)
    point=[]
    dx=(x_range[1]-x_range[0])/number_of_elevation
    
    for index,item in enumerate(elevation):
        point.append([x_range[0]+dx*index,
                      y_range[0],
                      item])
        point.append([x_range[0]+dx*index,
                      y_range[1],
                      item])
    element=[]
    for i in range(number_of_elevation-1):
        element.append([2*i,2*i+1,2*(i+1),2*(i+1)+1])
        
    return point, element
    
    
    

def save_surf(point_dic,x_range:list,y_range:list):
    i=0
    writer = vtk.vtkUnstructuredGridWriter()  
    if os.path.exists(os.path.join(os.getcwd(), 'VTKelevation')):
        print('\nVTKelevation exists in the current folder\nPlease remove it and try again.')
        exit()
    else:
        os.mkdir(os.path.join(os.getcwd(), 'VTKelevation'))
        
    
    
    for k in point_dic.keys():
        point,element=prepare_pe(point_dic[k],x_range,y_range)
        uGrids=make_sea(point,element)
        writer.SetFileName(os.path.join('VTKelevation','elevationVTK'+str(i)+'.vtk'))
        writer.SetInputData(uGrids)
        writer.Write()
        i+=1  
     

if __name__ == '__main__':
    # read nodes' positions
    with open(os.path.join('Elevation.json'), 'r') as f:
       point_position = json.load(f)
    f.close()

    save_surf(point_position,[0,10],[-1,1])
