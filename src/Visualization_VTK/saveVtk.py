
import vtk


point=[[0,0,0],
       [1,1,1],
       [1,0,0],
       [0,1,0],
       [0,0,1]]

line=[[0,1],
      [0,2],
      [0,3],
      [0,4]]

face=[[0,2,3],
      [0,2,4],
      [0,3,4]]


def __MakeMultiPoint(point_list:list):
    ug = vtk.vtkUnstructuredGrid()
    # A polyvertex is a cell represents a set of 0D vertices
    numberOfVertices = len(point_list)

    points = vtk.vtkPoints()
    polyVertex = vtk.vtkPolyVertex()
    polyVertex.GetPointIds().SetNumberOfIds(numberOfVertices)

    for index, item in enumerate(point_list):
        points.InsertNextPoint(item)
        polyVertex.GetPointIds().SetId(index, index)

    
    ug.SetPoints(points)
    ug.InsertNextCell(polyVertex.GetCellType(), polyVertex.GetPointIds())

    return ug



def __MakeMultiLine(point_list:list,line_list:list):
    ug = vtk.vtkUnstructuredGrid()

    points = vtk.vtkPoints()
    multiline = vtk.vtkLine()
    
    for item in point_list:
        points.InsertNextPoint(item)
    ug.SetPoints(points)
        
    for item in line_list:
        for i in range(2):
            multiline.GetPointIds().SetId(i, item[i])
        ug.InsertNextCell(multiline.GetCellType(), multiline.GetPointIds())
    return ug


def __MakeMultiFace(point_list:list,triangle_list:list):
    ug = vtk.vtkUnstructuredGrid()
    points = vtk.vtkPoints()
    triangles = vtk.vtkTriangle()  # 3-point elements
    pixel=vtk.vtkPixel() # 4-point elements
    
    # add points   
    for i in range(len(point_list)):
        points.InsertNextPoint(point_list[i])
    ug.SetPoints(points)
    
    # add elements # no index 
    for item in triangle_list:
    
        if len(item)==3: # add 3-point elements
            for i in range(3): 
                triangles.GetPointIds().SetId(i, item[i])
            ug.InsertNextCell(triangles.GetCellType(), triangles.GetPointIds())    
        elif len(item)==4:  # add 4-point elements
            for i in range(4): 
                pixel.GetPointIds().SetId(i,item[i])
            ug.InsertNextCell(pixel.GetCellType(), pixel.GetPointIds())    
    # add elements # with index 
    # for item in triangle_list:
    #     if len(item)==4: # add 3-point elements
    #         for i in range(1,4): 
    #             triangles.GetPointIds().SetId(i-1, item[i])
    #         ug.InsertNextCell(triangles.GetCellType(), triangles.GetPointIds())    
    #     elif len(item)==5:  # add 4-point elements
    #         for i in range(1,5): 
    #             pixel.GetPointIds().SetId(i-1,item[i])
    #         ug.InsertNextCell(pixel.GetCellType(), pixel.GetPointIds())                
    return ug

def write_vtk(file_name:str,**content):
    """a function to save vtk point, line and face file.
    Args:
        file_name (str): the name for the file
        content (arbitrary keyword arguments):
        point, line, face: python list
        
    """
    writer = vtk.vtkXMLDataSetWriter()
    #print(content)
    
    try:
        u1=__MakeMultiPoint(content['point'])
        writer.SetInputData(u1)
        writer.SetFileName(file_name+'.point.vtu')
        writer.Write()
    except:
        print('No point data is given')
    
    try:
        u2=__MakeMultiLine(content['point'],content['line'])
        writer.SetInputData(u2)
        writer.SetFileName(file_name+'.line.vtu')
        writer.Write()
    except:
        print('No line data is given')
    
    try:
        u3=__MakeMultiFace(content['point'],content['face'])
        writer.SetInputData(u3)
        writer.SetFileName(file_name+'.face.vtu')
        writer.Write()    
    except:
        print('No face data is given')

def write_point_vtk(file_name:str,**content):        
    writer = vtk.vtkXMLDataSetWriter()
    #print(content)
    try:
        u1=__MakeMultiPoint(content['point'])
        writer.SetInputData(u1)
        writer.SetFileName(file_name+'.point.vtu')
        writer.Write()
    except:
        print('No point data is given')
            
def write_line_vtk(file_name:str,**content):
    writer = vtk.vtkXMLDataSetWriter()
    
    u2=__MakeMultiLine(content['point'],content['line'])
    writer.SetInputData(u2)
    writer.SetFileName(file_name+'.line.vtu')
    writer.Write()

        
        
if __name__ == '__main__':
    write_vtk('test',point=point,line=line,face=face)    
        
    