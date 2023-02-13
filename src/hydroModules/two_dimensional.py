"""
-------------------------------------------\n
-         University of Stavanger          \n
-         Hui Cheng (PhD student)          \n
-          Lin Li (Medveileder)            \n
-      Prof. Muk Chen Ong (Supervisor)     \n
-------------------------------------------\n
Any questions about this code,
please email: hui.cheng@uis.no \n

"""
import numpy as np
import sys

np.set_printoptions(threshold=sys.maxsize)
row_air = 1.225  # [kg/m3]   air density
row_water = 1025.0  # [kg/m3]   sea water density
kinematic_viscosity = 1.004e-6  # when the water temperature is 20 degree.
dynamic_viscosity = 1.002e-3  # when the water temperature is 20 degree.
gravity = 9.81

class cage:
    def __init__(self, net_panels:list):
        self.netOBJ=net_panels

        self.element=[]
        for each in net_panels:
            self.element.append(each.element)

        pointID=[]
        for each_element in self.element:
            for i in each_element:
                pointID.append(i)

        self.pointID=list(set(pointID)) # set of points in a cage
        self.center=np.array([0.0,0.0,0.0])
        self.volume=0.0
        self.hydro_dynamic_forces = np.array([0,0,0], dtype=float)
        self.hydro_static_forces =np.array([0,0,0], dtype=float)
        self.hydro_total_forces = np.array([0,0,0], dtype=float)
        self.reduction=[]

    def get_cage_center(self,position:np.array):
        self.center=position[self.pointID,:].mean(axis=0)
    
    def get_cage_volume(self):
        volume=0
        for each_element in self.netOBJ:
            volume+=each_element.area*np.dot(each_element.center/3.0,each_element.vectorN)
        self.volume=volume
   

    def get_hydro_force(self):
        self.hydro_total_forces = np.array([0,0,0], dtype=float)
        for each_element in self.netOBJ:
            self.hydro_total_forces+=each_element.hydro_total_forces
    
    def get_r(self,position:np.array):
        self.reduction=[1.0]*len(position)
        self.get_cage_center(position)
        for i in self.pointID:
            o_p=position[i]-self.center
            if o_p[0]>0:
                self.reduction[i]=0.8
            else:
                self.reduction[i]=1.0
            
    
    
                      
    # def cal_volume_tetrahedron(self, point1, point2, point3):
    #     vector_a = point1-self.center
    #     vector_b = point2-self.center
    #     vector_c = point3-self.center
    #     return abs(np.dot(vector_a, (np.cross(vector_b, vector_c))))/6.0
        
    # def get_cage_volume(self,position:np.array):
    #     self.get_cage_center(position)
    #     volume=0
    #     for each_element in self.element:
    #         if len(each_element)==3:
    #             volume+=cal_volume_tetrahedron(position[each_element[0]],
    #                                            position[each_element[1]]
    #                                            position[each_element[2]])
    #         if len(each_element)==4: # quad = 2 triangle
    #             volume+=cal_volume_tetrahedron(position[each_element[0]],
    #                                            position[each_element[1]]
    #                                            position[each_element[2]])
    #             volume+=cal_volume_tetrahedron(position[each_element[0]],
    #                                            position[each_element[3]]
    #                                            position[each_element[2]])          
    
class netting:
    """
    For Screen hydrodynamic models, the forces on netting are calculated based on individual a panel section of netting.
    The twines and knots in the net panel are considered as an integrated structure. In this module, the net panel is defined by
    three nodes because three (non-collinear) points can determine a unique plane in Euclidean geometry.
    In practice, the force is usually decomposed into two components: drag force F_D and lift force F_L (Cheng et al., 2020).
    """

    def __init__(self, model_index, hydro_element:list, solidity:float, dw0=0.0):
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
        self.v=np.array([0,0,0], dtype=float) # structural velocity of element
        self.u=np.array([0,0,0], dtype=float) # fluid velocity on element
        self.vectorN=np.array([0,0,0], dtype=float)
        self.center=np.array([0,0,0], dtype=float)
        self.theta=0.0 #rad
        self.ratio=0.0 # submerged ratio, 1 -> fully submerged, 0 -> all in air
        self.cd=0.0
        self.cl=0.0
        self.hydro_dynamic_forces = np.array([0,0,0], dtype=float)
        self.hydro_static_forces =np.array([0,0,0], dtype=float)
        self.hydro_total_forces = np.array([0,0,0], dtype=float)

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
        self.area =0.0        
        for i in range(len(self.element)-2):
            a1 = position[self.element[0]] - position[self.element[i+1]]
            a2 = position[self.element[0]] - position[self.element[i+2]]
            self.area += 0.5 * np.linalg.norm(np.cross(a1, a2))
    
    def get_center(self,position:np.array):
        self.center=np.array([0.0,0.0,0.0])
        for i in self.element:
            self.center+=position[i]/len(self.element)
            
    def get_normal_vector(self,position:np.array): # all the normal vector is point out
        if len(self.element)==3:
            a1 = position[self.element[0]] - position[self.element[1]]
            a2 = position[self.element[-1]] - position[self.element[-2]]
            self.vectorN=np.cross( a1, a2) / np.linalg.norm(np.cross(a1, a2))
            
        elif len(self.element)==4: # need to find an advanced way to calculate the normal vector of quad.
            a1 = position[self.element[0]] - position[self.element[1]]
            a2 = position[self.element[0]] - position[self.element[2]]
            self.vectorN=np.cross( a1, a2) / np.linalg.norm(np.cross(a1, a2))
    
    def get_v(self, velocity:np.array):
        """get the structural velocity of the element

        Args:
            velocity (np.array): structural velocity of N points
        """
        self.v=np.array([0,0,0], dtype=float) # structural velocity of element
        for i in self.element:
            self.v+=velocity[i]/len(self.element)
            
    def get_u(self, velocity:np.array):
        """get the fluid velocity on the element

        Args:
            velocity (np.array): fluid velocity of N points.
        """
        self.u=np.array([0,0,0], dtype=float) # fluid velocity on element
        for i in self.element:
            self.u+=velocity[i]/len(self.element)
        
        
        
    def get_inflow_angle(self, position:np.array,u):
        self.get_normal_vector(position)
        self.get_u(u)
        self.theta=np.arccos(abs(np.dot(self.vectorN,self.u))/np.linalg.norm(np.dot(self.vectorN,self.u)))
    
    def get_force_coefficient(self):
        if self.modelIndex == 'S1':  # aarsnes 1990
            self.cd = 0.04 + (-0.04 + self.sn - 1.24 * pow(self.sn, 2) + 13.7 * pow(self.sn, 3)) * np.cos(
                self.theta)
            self.cl = (0.57 * self.sn - 3.54 * pow(self.sn, 2) + 10.1 * pow(self.sn, 3)) * np.sin(
                2 * self.theta)
            
        elif  self.modelIndex=='S2':# Loland 1991
            self.cd = 0.04 + (-0.04 + 0.33 * self.sn + 6.54 * pow(self.sn, 2) - 4.88 * pow(self.sn, 3)) * np.cos(
                self.theta)
            self.cl  = (-0.05 * self.sn + 2.3 * pow(self.sn, 2) - 1.76 * pow(self.sn, 3)) * np.sin(
                2 * self.theta)            
        
        elif self.modelIndex=='S3':  # Kristiansen 2012
            a1 = 0.9
            a2 = 0.1
            b1 = 1.0
            b2 = 0.1
            reynolds_number = row_water * self.dw0 * np.linalg.norm(self.u) / dynamic_viscosity / (1 - self.sn)  # Re
            cd_cylinder = -78.46675 \
                         + 254.73873 * np.log10(reynolds_number) \
                         - 327.88640 * pow(np.log10(reynolds_number), 2) \
                         + 223.64577 * pow(np.log10(reynolds_number), 3) \
                         - 87.922340 * pow(np.log10(reynolds_number), 4) \
                         + 20.007690 * pow(np.log10(reynolds_number), 5) \
                         - 2.4489400 * pow(np.log10(reynolds_number), 6) \
                         + 0.1247900 * pow(np.log10(reynolds_number), 7)
            cn_45 = cd_cylinder * self.sn / (2.0 * pow((1 - self.sn), 2))
            cl_45=np.pi*cn_45/(8.0+cn_45)

            cd_0 = cd_cylinder * (self.sn * (2 - self.sn)) / (2.0 * pow((1 - self.sn), 2))
            cl_0=(0.5 * cd_0 -cl_45)/pow(2,0.5)

            self.cd = cd_0 *(a1 * np.cos(self.theta) + a2 * np.cos(3 * self.theta))
            self.cl = cl_0 *(b1 * np.sin(2 * self.theta) + b2 * np.sin(4 * self.theta))
        else:
            self.cd=0
            self.cl=0
    
    def cal_submerge_ratio(self,position:np.array,elevation:list):
        """
        :param list_z: [ae,be,ce] | from node to water surface, + means under water,- means above water 
        :param elevation: np.array[n,1] |Unit [m]| elevation of the sea, n is the number of nodes
        :return: np.array[m,3] |Unit [N]| buoyancy force on elements, m is the number of elements
        """
        z_list=[]
        for i in self.element:
            z_list.append(elevation[i]-position[i,2])
        sign_set = set(np.sign(z_list))
        ratio_water = 1.0
        
        if sign_set == {-1} or sign_set == {0, -1}:  # ---, 0--, all in air
            ratio_water = 0.0

        elif sign_set == {1} or sign_set == {0, 1}:  # all in water
            ratio_water = 1.0

        elif sign_set == {0}:  # exact half
            ratio_water = 0.5
        else: 
            if np.sign(z_list).sum() > 0:
                ratio_water=0.75
            elif np.sign(z_list).sum() == 0:
                ratio_water=0.5      
            else:    
                ratio_water=0.25      
        self.ratio=min(ratio_water, 1.0)

    def get_hydrostatic_force(self,position:np.array,elevation:list):
        self.hydro_static_forces =np.array([0,0,0], dtype=float)
        # self.cal_submerge_ratio(position,elevation)
        self.ratio=1.0
        element_volume = self.area * self.sn * self.dw0 * 0.25 * np.pi
        # print('sumberged ratio ====='+str(self.ratio))
        # print('element_volume ====='+str(element_volume))
        self.hydro_static_forces=np.array([0,0,self.ratio*gravity*element_volume*row_water])
        
                
        
    def get_hydrodynamic_force(self,position:np.array,v:np.array,u:np.array):
        self.get_inflow_angle(position,u)
        self.get_v(v)
        self.get_force_coefficient()
        
        drag_e=(self.u-self.v)/np.linalg.norm(self.u-self.v+ np.finfo(np.float64).eps)
        lift_e=np.cross(np.cross(drag_e,self.vectorN),drag_e)
        
        row = row_water*self.ratio+row_air*(1-self.ratio)

        fd = 0.5 * row * self.area * self.cd * \
            pow(np.linalg.norm(self.u-self.v), 2) * drag_e
        fl = 0.5 * row * self.area  * self.cd * \
            pow(np.linalg.norm(self.u-self.v), 2) * lift_e
        self.hydro_dynamic_forces=fd+fl
    
    def get_hydro_force(self,position:np.array,elevation:list,v:np.array,u:np.array):
        self.hydro_total_forces = np.array([0,0,0], dtype=float)
        
        self.get_hydrostatic_force(position,elevation)
        self.get_hydrodynamic_force(position,v,u)
        self.hydro_total_forces=self.hydro_dynamic_forces+self.hydro_static_forces
    
    

if __name__ == "__main__":
    pass
