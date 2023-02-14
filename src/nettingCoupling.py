"""
-------------------------------------------\n
-         University of Stavanger          \n
-         Hui Cheng (PhD student)          \n
-          Lin Li (Medveileder)            \n
-      Prof. Muk Chen Ong (Supervisor)     \n
-------------------------------------------\n
Any questions about this code,
please email: hui.cheng@uis.no \n
In order to use this module, we recommend ``import nettingFSI as fsi`` in the beginning of your code.
To use this module, one should be import this into the input file for calculations.

"""
import os
import time
import sys
import numpy as np

np.set_printoptions(threshold=sys.maxsize)

# here we assume Code_aster is much faster than OpenFoam, thus OpenFOAM do not need to wait .

def start_flag(cwd, flag):
    """
    :param cwd: work path
    :param flag: file name [string]
    :return:  create a empty file to tell openfoam the file starts writing
    """
    if os.path.isfile(os.path.join(cwd,str(flag))):
        print("path is :"+str(os.path.join(cwd,str(flag))))
        os.remove(os.path.join(cwd,str(flag)))
    else:
        pass


def finish_flag(cwd, flag):
    """
    :param cwd: work path
    :param flag: file name [string]
    :return:  create a empty file to tell openfoam the file finishes writing
    """
    # print("in finish_flag, the cwd is "+cwd)
    os.mknod(os.path.join(cwd,str(flag)))


def write_position(position, cwd):
    """
    :param position: A numpy array of all the nodes' position
    :param cwd: work path,
    :return: write the nodes' positions to "constant" folder as a file named "posi"
    """
    print(cwd)
    start_flag(cwd, "position.flag")
    head_file = ["// Input for the nets in openfoam.",
                 "// Author: Hui Cheng",
                 "// Email: hui.cheng@uis.no",
                 "FoamFile",
                 "{",
                 "    version     2.0;",
                 "    format      ascii;",
                 "    class       dictionary;",
                 "    location    \"constant\";",
                 "    object      posi;",
                 "}",
                 "// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> //",
                 "{",
                 "numOfPoint   " + str(len(position)) + ";"]
    with open(cwd + '/posi.tmp', 'w') as output_file:
        for item in head_file:
            output_file.write(item + "\n")
        for index, item in enumerate(position):
            output_file.write(
                "p" + str(index) + "\t( " + str(item[0]) + "\t" + str(item[1]) + "\t" + str(item[2]) + " );\n")
        output_file.write("}\n")
        output_file.write("// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< //\n")
    output_file.close()
    os.rename(cwd + '/posi.tmp', cwd + '/posi')
    finish_flag(cwd, "position.flag")


def genFSInet(netElement:list):
    fsiNet=[]
    for item in netElement:
        if len(item)==3:
            fsiNet.append(item)
        elif len(item)==4:
            fsiNet.append([item[0],item[1],item[2]])
            fsiNet.append([item[0],item[2],item[3]])
    return fsiNet

def genFh(netElement:list,fsiVolume:list,fh:list):
    netFh=[]
    for index, item in enumerate(netElement):
        if len(item)==3:
            netFh.append(fh[index])
        elif len(item)==4:
            netFh.append((np.array(fh[index])*0.5).tolist())
            netFh.append((np.array(fh[index])*0.5).tolist())
    return np.array(netFh)/np.array(fsiVolume).reshape(len(fsiVolume),1)

def genFSIvolume(netElement:list,thickness:float,position:np.array):
    fsiVolume=[]
    for item in netElement:
        a1 = position[item[0]] - position[item[1]]
        a2 = position[item[0]] - position[item[2]]
        area = 0.5 * np.linalg.norm(np.cross(a1, a2))
        fsiVolume.append(thickness*area)
    return fsiVolume
    
    

def write_element(hydro_element, cwd):
    """
    :param hydro_element: A numpy array of all the net panel element
    :param cwd: work path,
    :return: write the net panels to "constant" folder as a file named "surf"
    """
    head_file = ["// Input for the nets in openfoam.",
                 "// Author: Hui Cheng",
                 "// Email: hui.cheng@uis.no",
                 "FoamFile",
                 "{",
                 "    version     2.0;",
                 "    format      ascii;",
                 "    class       dictionary;",
                 "    location    \"constant\";",
                 "    object      surc;",
                 "}",
                 "// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> //",
                 "{",
                 "numOfSurf   " + str(len(hydro_element)) + ";"]
    with open(cwd + '/surf', 'w') as output_file:
        for item in head_file:
            output_file.write(item + "\n")
        for index, item in enumerate(hydro_element):
            output_file.write(
                "e" + str(index) + "\t( " + str(item[0]) + "\t" + str(item[1]) + "\t" + str(item[2]) + " );\n")
        output_file.write("}\n")
        output_file.write("// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< //\n")
    output_file.close()


def write_fh(hydro_force, time_fe, cwd):
    """
    :param time_fe: simulation in FEM
    :param hydro_force: A numpy array of all the hydrodynamic forces on net panel
    :param cwd: work path,
    :return: write the hydrodynamic forces to "constant" folder as a file named "Fh" and save the total hydrodynamic forces on netting to "forceOnNetting.txt"
    """
    print("Here>>>>>>>>>>>>>Fh>>>  " + str(time_fe))
    start_flag(cwd, "fh.flag")
    head_file = ["// Input for the nets in openfoam.",
                 "// Author: Hui Cheng",
                 "// Email: hui.cheng@uis.no",
                 "FoamFile",
                 "{",
                 "    version     2.0;",
                 "    format      ascii;",
                 "    class       dictionary;",
                 "    location    \"constant\";",
                 "    object      Fh;",
                 "}",
                 "// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> //",
                 "{",
                 "timeInFE   " + str(round(time_fe, 3)) + ";",
                 "numOfFh    " + str(len(hydro_force)) + ";"]

    with open(cwd + '/Fh.tmp', 'w') as output_file:
        for item in head_file:
            output_file.write(item + "\n")
        for index, item in enumerate(hydro_force):
            output_file.write(
                "fh" + str(index) + "\t( " + str(item[0]) + "\t" + str(item[1]) + "\t" + str(item[2]) + " );\n")
        output_file.write("}\n")
        output_file.write("// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< //\n")
    output_file.close()
    os.rename(cwd + '/Fh.tmp', cwd + '/Fh')
    finish_flag(cwd, "fh.flag")

velocity_dict = {'time_record': ['0']}


def get_velocity(cwd, length_velocity, time_aster):
    """
    :param cwd: working path for code aster
    :param length_velocity: the length of the list of element(or vlocity)
    :param time_aster: the time in Code_Aster
    :return: a numpy array of velocity on elements
    """
    print("Here>>>>>>>>>>>>>velo>>>  " + str(time_aster))
    cwd_foam_root = "/".join(cwd.split("/")[0:-1])
    velocity_dict, time_foam = read_velocity(cwd_foam_root, length_velocity)
    time_foam = velocity_dict['time_record'][-1]
    while float(time_aster) > float(time_foam) or str(time_foam) == str(0):
        time.sleep(0.1)
        velocity_dict, time_foam = read_velocity(cwd_foam_root, length_velocity)
    else:
        # velocity_dict, time_foam = read_velocity(cwd_foam_root, length_velocity)
        # print("velocity from get_velocity i s "+str(np.array(velocity_dict["velocities_at_" + str(time_foam)])))
        f = open(cwd + "/velocity_on_element.txt", "w")
        f.write(str(velocity_dict))
        f.close()
        return np.array(velocity_dict["velocities_at_" + str(time_foam)])


def read_velocity(cwd_foam_root, length_velocity):
    """
    :param cwd_foam_root: work path
    :param length_velocity: The length of velocity file, normally equal to the length of elements
    :return: return the velocity dictionary and time in openfoam
    """
    while not os.path.exists(os.path.join(cwd_foam_root,"velocity_on_elements.txt")):
        # print("path= "+ str(os.path.join(cwd_foam_root,"velocity_on_elements.txt")))
        print("Wait for velocity from OpenFoam......")
        time.sleep(0.1)

    else:
        f = open(cwd_foam_root + "/velocity_on_elements.txt", "r")
        lines = f.readlines()
        time_foam = str(0)
        for line in lines[-length_velocity * 5:]:
            if "The velocities at" in line:
                time_foam = line.split(" ")[3]
                if len(line.split()) == 6:
                    start_line = lines.index(line) + 3
                    end_line = lines.index(line) + 3 + length_velocity
                    if (end_line < len(lines)) and (str(time_foam) not in velocity_dict['time_record']):
                        velocity_dict['time_record'].append(str(time_foam))
                        velocity_dict["velocities_at_" + str(time_foam)] = []
                        for item in range(length_velocity):
                            try:
                                velocity_dict["velocities_at_" + str(time_foam)].append(
                                    [float(lines[item + start_line].split()[0][1:]),
                                     float(lines[item + start_line].split()[1]),
                                     float(lines[item + start_line].split()[2][:-1])])
                            except:
                                velocity_dict["velocities_at_" + str(time_foam)].append(
                                    [float(lines[start_line].split()[0][1:]),
                                     float(lines[start_line].split()[1]),
                                     float(lines[start_line].split()[2][:-1])])
    return velocity_dict, velocity_dict['time_record'][-1]

if __name__ == "__main__":
    pass