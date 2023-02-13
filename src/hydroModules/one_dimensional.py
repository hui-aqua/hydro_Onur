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


class line:
    """
    For Morison hydrodynamic models, the forces on netting are calculated based on individual twines.
    The twines are taken as cylindrical elements. In practice, the force is usually decomposed into two components:
    normal drag force F_n and tangential drag force F_t (Cheng et al., 2020)
    """

    def __init__(self, model_index, hydro_element:list, dw0=0.0):
        """
        :param model_index: [string] Unit: [-]. To indicate the model function, e.g.: 'M1', 'M2', 'M3'.
        :param hydro_element: [list] Unit: [-]. A python list to indicate how the lines are connected.
        :param solidity: [float] Unit: [-]. The solidity of netting.
        :param dw0: [float] Unit: [m]. The diameter of the physical net twines. It is used for the hydrodynamic coefficients.
        :param dwh: [float] Unit: [m]. The hydrodynamic diameter of the numerical net twines. It is used for the force calculation (reference area)
        """
        self.model_index = str(model_index)
        self.dw0 = dw0  # used for the hydrodynamic coefficients
        self.element=hydro_element
        #
        self.length=0.0
        self.v=np.array([0,0,0], dtype=float) # structural velocity of element
        self.u=np.array([0,0,0], dtype=float) # fluid velocity on element
        self.vectorE=np.array([0,0,0], dtype=float)
        self.center=np.array([0,0,0], dtype=float)
        self.ratio=0.0 # submerged ratio, 1 -> fully submerged, 0 -> all in air
        self.alpha=0.0 # rad
        self.cn=0.0
        self.ct=0.0
        self.hydro_dynamic_forces = np.array([0,0,0], dtype=float)
        self.hydro_static_forces =np.array([0,0,0], dtype=float)
        self.hydro_total_forces = np.array([0,0,0], dtype=float)

    def __str__(self):
        """Print information of the present object."""
        s0 = "Morison model\n"
        s1 = "The model index is " + str(self.modelIndex) + "\n"
        return s0 + s1 
    
    
    def update_state(self,position:np.array):
        self.get_length(position)
        self.get_element_vector(position)
        
    def get_element_vector(self,position:np.array): # all the normal vector is point out
        self.vectorE=position[self.element[0]]-position[self.element[1]]
        
    def get_length(self,position:np.array):
        self.get_element_vector(position)
        self.length=np.linalg.norm(self.vectorE)
        # self.length=0.1
        
    def get_center(self,position:np.array):
        self.center=np.array([0,0,0], dtype=float)
        for i in self.element:
            self.center+=position[i]/len(self.element)
 
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

    def get_attack_angle(self, position:np.array,u):
        self.get_element_vector(position)
        self.get_u(u)
        # (-pi/2,pi/2)
        self.alpha=float(np.arccos(np.dot(self.vectorE,self.u)/(np.finfo(float).eps+np.linalg.norm(np.dot(self.vectorE,self.u)))) )
            
    def get_force_coefficient(self):
        re_num = row_water * self.dw0 * np.linalg.norm(self.u) / dynamic_viscosity
        if self.model_index == 'M3':  # Takagi 2004
            if re_num < 200:
                self.cn = pow(10, 0.7) * pow(re_num, -0.3)
            else:
                self.cn = 1.2
            self.ct = 0.1
            
        elif self.model_index == 'M4':  # choo 1971
            self.ct = np.pi * dynamic_viscosity * (0.55 * np.sqrt(re_num) + 0.084 * pow(re_num, 2.0 / 3.0))
            s = -0.07721565 + np.log(8.0 / re_num)
            if 0 < re_num < 1:
                self.cn = 8 * np.pi * (1 - 0.87 * pow(s, -2)) / (s * re_num)
            elif re_num < 30:
                self.cn = 1.45 + 8.55 * pow(re_num, -0.9)
            elif re_num < 2.33e5:
                self.cn = 1.1 + 4 * pow(re_num, -0.5)
            elif re_num < 4.92e5:
                self.cn = (-3.41e-6) * (re_num - 5.78e5)
            elif re_num < 1e7:
                self.cn = 0.401 * (1 - np.exp(-re_num / 5.99 * 1e5))
            else:
                self.cn = 1.2
    
    def cal_submerge_ratio(self,position:np.array,elevation:list):
        z_list=[]
        for i in self.element:
            z_list.append(elevation[i]-position[i][2])
        sign_set = set(np.sign(z_list))
        ratio_water = 1.0
        if sign_set == {-1} or sign_set == {0, -1}:  # --, 0-, all in air
            ratio_water = 0.0

        elif sign_set == {1} or sign_set == {0, 1}:  # all in water
            ratio_water = 1.0

        elif sign_set == {0}:  # exact half
            ratio_water = 0.5
        else:
            ratio_water=np.array(z_list).sum()/abs(np.array(z_list)).sum()/2.0+0.5
        self.ratio=min(ratio_water,1.0)
        
    def get_hydrostatic_force(self,position:np.array,elevation:list):
        # self.cal_submerge_ratio(position,elevation)
        self.ratio=1
        element_volume = self.length * 0.25*np.pi*pow(self.dw0,2)
        self.hydro_static_forces=np.array([0,0,self.ratio*gravity*element_volume*row_water])
        
    def get_hydrodynamic_force(self,position:np.array,v:np.array,u:np.array):
        self.get_attack_angle(position,u)
        self.get_v(v)
        self.get_force_coefficient()
                
        row = row_water*self.ratio+row_air*(1-self.ratio)
        ut=self.vectorE/self.length*np.cos(self.alpha)*np.linalg.norm(self.u-self.v)
        
        un=(self.u-self.v)-ut
        
        ft = 0.5 * row * self.dw0 * self.length * self.ct * ut*np.linalg.norm(ut)
        
        fn = 0.5 * row * self.dw0 * self.length * self.cn * un*np.linalg.norm(un)
        self.hydro_dynamic_forces=fn+ft               

    def get_hydro_force(self,position:np.array,elevation:list,v:np.array,u:np.array):
        self.hydro_total_forces = np.array([0,0,0], dtype=float)
        
        self.get_hydrostatic_force(position,elevation)
        self.get_hydrodynamic_force(position,v,u)
        self.hydro_total_forces=self.hydro_dynamic_forces+self.hydro_static_forces
    

class pipe(line):
    """
    For Morison hydrodynamic models, the forces on netting are calculated based on individual pipe.
    The twines are taken as cylindrical elements. In practice, the force is usually decomposed into two components:
    normal drag force F_n and tangential drag force F_t (Cheng et al., 2020)
    """
    def __init__(self, model_index, hydro_element, section_diameter,thickness,permeability=False):
        """[summary]

        Args:
            model_index ([type]): [description]
            hydro_element ([type]): [description]
            section_diameter ([type]): [description]
            thickness ([type]): [description]
            permeability (bool, optional): [description]. Defaults to False.
        """
        self.model_index = str(model_index)
        self.element= hydro_element
        self.dw0 = section_diameter  # used for the hydrodynamic coefficients & hydrodynamic forces
        self.t = thickness  # used for the force calculation (reference area)
                
        if permeability: # if the pipe is permeable
            area_out=0.25*np.pi*pow(section_diameter,2)
            area_ini=0.25*np.pi*pow(section_diameter-2*thickness,2)
            cross_area=area_out-area_ini
            self.dwh = pow(4.0*cross_area/np.pi,0.5)
        else:
            self.dwh = section_diameter  # used for the buoyancy force
        # print('dwh is '+str(self.dwh))
        
        self.hydro_dynamic_forces = np.array([0,0,0], dtype=float)
        self.hydro_static_forces =np.array([0,0,0], dtype=float)
        self.hydro_total_forces = np.array([0,0,0], dtype=float)

    def get_hydrostatic_force(self,position:np.array,elevation:list):
        # self.cal_submerge_ratio(position,elevation)
        self.ratio=1
        element_volume = self.length * 0.25*np.pi*pow(self.dwh,2)
        self.hydro_static_forces=np.array([0,0,self.ratio*gravity*element_volume*row_water])



if __name__ == "__main__":
    pass
