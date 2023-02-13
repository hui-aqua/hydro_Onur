"""
-------------------------------------------\n
-         University of Stavanger          \n
-         Hui Cheng (PhD student)          \n
-          Lin Li (Medveileder)            \n
-      Prof. Muk Chen Ong (Supervisor)     \n
-------------------------------------------\n
Any questions about this code,
please email: hui.cheng@uis.no \n
wave direction is X+
water lever Z = 0
"""
import numpy as np
from numpy import pi
gravity = 9.81


class Airywave:
    """
    Using Airy wave theory      \n
    Ref. DNV GL-RP205 Ver. 2008:P45
    Linear wave theory (or sinusoidal wave theory).

    """

    def __init__(self, wave_height=1.0, wave_period=10.0, water_depth=60.0, direction=0.0, phase=0.0):
        """
        :param wave_height: [float] | Unit: [m]. wave height.
        :param wave_period: [float] | Unit: [s]. wave period.
        :param water_depth: [float] | Unit: [m]. wave depth.
        :param direction: [float] | Unit: [degree]. direction of propagation. Measured from the X+, 90 degree is Y+.
        :param phase: [float] | Unit: [degree]. phase, usually is 0, but can be random in irregular wave theory
        """

        self.wave_Height = wave_height
        self.wave_Period = wave_period
        self.water_Depth = water_depth
        self.wave_beta = pi * direction / 180.0
        self.wave_phase = pi * phase / 180.0

        # 0 non-dimensional parameter
        self.wave_steepness = 2*pi*wave_height/(gravity*pow(wave_period, 2))
        self.shallow_water = 2*pi*water_depth/(gravity*pow(wave_period, 2))
        self.ursell_number = self.wave_steepness/pow(self.shallow_water, 3)

        # 1 Calculation
        # wave length using an accurate approximation
        alpha = [1, 0.666, 0.445, -0.105, 0.272]
        omega_ba = 4.0 * pow(pi, 2) * water_depth / \
            (gravity*pow(wave_period, 2))
        f_omega = 0.0
        for index, item in enumerate(alpha):
            f_omega += item * pow(omega_ba, index)

        self.wave_Length = wave_period * \
            pow(gravity * water_depth, 0.5) * \
            pow(f_omega / (1 + omega_ba * f_omega), 0.5)

        # wave number
        self.wave_k = 2 * pi / self.wave_Length

        # angular frequency
        self.omega = pow(gravity * self.wave_k * np.tanh(self.wave_k * water_depth), 0.5)
        # phase velocity
        self.wave_phase_velocity = pow(gravity / self.wave_k * np.tanh(self.wave_k * water_depth), 0.5)

        # 2 for easy calculation
        self.theta = 0
        self.pi_h_t = pi * wave_height / wave_period
        self.pi_h_t_2 = wave_height * pow(pi / wave_period, 2)
        self.pi_h_l = pi*wave_height/self.wave_Length
        self.kd = self.wave_k*water_depth
        

    def __str__(self):
        """ Print the information of object. """
        s0 = 'The parameters for wave are:\n'
        s1 = 'water Depth= ' + str(self.water_Depth) + ' m\n'
        s2 = 'wave Period= ' + str(self.wave_Period) + ' s\n'
        s2_1 = 'wave number= ' + str(self.wave_k) + ' 1/m\n'
        s3 = 'wave Length= ' + str(self.wave_Length) + ' m\n'
        s4 = 'wave Height= ' + str(self.wave_Height) + ' m\n'
        s5 = 'wave phase velocity= ' + str(self.wave_phase_velocity) + ' m/s\n'
        s6 = 'wave direction: ' + str(self.wave_beta) + ' degree\n'
        s7 = 'wave steepness = ' +str(self.wave_steepness)+ '\n'
        s8 = 'wave shallowness = ' +str(self.shallow_water)+ '\n'
        s9 = 'wave Ursell number = ' +str(self.ursell_number)+ '\n'
        return s0 + s1 + s2 + s2_1 + s3 + s4 + s5 + s6 +s7+s8+s9

    def calc_theta(self, position, global_time):
        """
        A private function. \n
        :param position: [np.array].shape=(n,3) point coordinates Unit: [m]. \n
        The coordinate of the point which you want to know the wave surface elevation. can be [x,y] or [x,y,z]
        :param global_time: [float] Unit: [s].
        :return: [float] Unit: [m]./theta the targeted position.
        """
        if len(position.shape) == 1:
            # only one theta
            self.theta = self.wave_k * (position[0] * np.cos(self.wave_beta) + position[1] * np.sin(
                self.wave_beta)) - self.omega * global_time + self.wave_phase

        elif len(position.shape) == 2:
            # a list of theta
            self.theta = self.wave_k * (position[:, 0] * np.cos(self.wave_beta) + position[:, 1] * np.sin(
                self.wave_beta)) - self.omega * global_time + self.wave_phase

    def calc_elevation(self, position, global_time):
        """
        A private function. \n
        :param position: [np.array].shape=(n,3) coordinates Unit: [m]. The position of the point which you want to know the wave surface elevation.
        :param global_time: [float] Unit: [s].
        :return: scale or [np.array].shape=(n,) Unit: [m]. The sea surface level in Z direction.
        """
        self.calc_theta(position, global_time)
        return self.wave_Height / 2 * np.cos(self.theta)

    def get_velocity_at_one_node(self, position, global_time, irregularwaves=False):
        """
        :param position: [np.array].shape=(1,3) | Unit: [m]. The position of the point which you want to know the wave velocity
        :param global_time: [float] Unit: [s]. The time which you want to know the wave velocity
        :return:  [np.array].shape=(1,3) Unit: [m/s]. A numpy array of the velocity at the targeted point.
        """
        self.calc_theta(position, global_time)
        eta = self.calc_elevation(position, global_time)

        position_z = position[2]
        velocity = np.zeros((len(global_time), 3))
        velocity[:, 0] = np.cos(self.wave_beta) * self.pi_h_t * np.cosh(self.wave_k * (
            position_z + self.water_Depth)) * np.cos(self.theta) / np.sinh(self.kd)
        velocity[:, 1] = np.sin(self.wave_beta) * self.pi_h_t * np.cosh(self.wave_k * (
            position_z + self.water_Depth)) * np.cos(self.theta) / np.sinh(self.kd)
        velocity[:, 2] = self.pi_h_t * np.sinh(self.wave_k * (
            position_z + self.water_Depth)) * np.sin(self.theta) / np.sinh(self.kd)
        if not irregularwaves:
            velocity[position[2] > eta] = 0.0
        return velocity

    def get_acceleration_at_one_node(self, position, global_time, irregularwaves=False):
        """
        :param position: [np.array].shape=(1,3) | Unit: [m]. The position of the point which you want to know the wave acceleration.
        :param global_time: [float] Unit: [s]. The time which you want to know the wave velocity
        :return: [np.array].shape=(1,3) Unit: [m/s]. A numpy array of the acceleration at the targeted point.
        """
        self.calc_theta(position, global_time)
        eta = self.calc_elevation(position, global_time)
        position_z = position[2]
        acceleration = np.zeros((len(global_time), 3))
        acceleration[:, 0] = np.cos(self.wave_beta) * 2*self.pi_h_t_2 * np.cosh(
            self.wave_k * (position_z + self.water_Depth)) * np.sin(self.theta) / np.sinh(self.kd)
        acceleration[:, 1] = np.sin(self.wave_beta) * 2*self.pi_h_t_2 * np.cosh(
            self.wave_k * (position_z + self.water_Depth)) * np.sin(self.theta) / np.sinh(self.kd)
        acceleration[:, 2] = -2*self.pi_h_t_2 * np.sinh(self.wave_k * (position_z + self.water_Depth)) * np.cos(
            self.theta) / np.sinh(self.kd)
        if not irregularwaves:
            acceleration[position[2] > eta] = 0.0
        return acceleration

    def get_elevations_at_one_node(self, position, time_list):
        """
        Public function.\n
        :param position: [np.array].shape=(n,3) Unit: [m]. The position of one node
        :param time_list: [np.array].shape=(n,1) | Uint: [s]. The time sequence for getting the elevations
        :return: Get a list of elevations at one position with a time sequence \n
        """
        return self.calc_elevation(position, time_list)

    def get_elevation_at_nodes(self, list_of_point, global_time):
        """
        Public function.\n
        :param list_of_point: [np.array].shape=(n,3) Unit: [m]. A list of node positions
        :param global_time: time [s] \n
        :return: Get a list of elevation at a list of point \n
        """
        return self.calc_elevation(list_of_point, global_time)

    def get_velocity_at_nodes(self, list_of_point, global_time, irregularwaves=False):
        """
        Public function.\n
        :param list_of_point:  [np.array].shape=(n,3) Unit: [m]. A list of points's positions
        :param global_time: [float] Unit: [s]. Physical time.
        :return: Get a list of velocity at a list of point\n
        """
        self.calc_theta(list_of_point, global_time)
        eta = self.calc_elevation(list_of_point, global_time)

        positions_z = list_of_point[:, 2]
        velocities = np.zeros_like(list_of_point)
        u1 = self.pi_h_t * np.cosh(self.wave_k * (positions_z +
                                   self.water_Depth)) / np.sinh(self.kd) * np.cos(self.theta)

        velocities[:, 0] = np.cos(self.wave_beta) * u1
        velocities[:, 1] = np.sin(self.wave_beta) * u1

        velocities[:, 2] = self.pi_h_t * np.sinh(self.wave_k * (
            positions_z + self.water_Depth)) * np.sin(self.theta) / np.sinh(self.kd)
        if not irregularwaves:
            velocities[list_of_point[:, 2] > eta] = 0.0
        return velocities

    def get_acceleration_at_nodes(self, list_of_point, global_time, irregularwaves=False):
        """
        Public function.\n
        :param list_of_point: [np.array].shape=(n,3) Unit: [m]. A list of points's positions
        :param global_time: time [s] \n
        :return: Get a list of acceleration at a list of point \n
        """
        self.calc_theta(list_of_point, global_time)
        eta = self.calc_elevation(list_of_point, global_time)

        positions_z = list_of_point[:, 2]
        accelerations = np.zeros_like(list_of_point)

        accelerations[:, 0] = np.cos(self.wave_beta) * 2*self.pi_h_t_2 * np.cosh(
            self.wave_k * (positions_z + self.water_Depth)) * np.sin(self.theta) / np.sinh(self.kd)
        accelerations[:, 1] = np.sin(self.wave_beta) *2*self.pi_h_t_2 * np.cosh(
            self.wave_k * (positions_z + self.water_Depth)) * np.sin(self.theta) / np.sinh(self.kd)
        accelerations[:, 2] = -2*self.pi_h_t_2 * np.sinh(self.wave_k * (positions_z + self.water_Depth)) * np.cos(self.theta) / np.sinh(
            self.kd)
        if not irregularwaves:
            accelerations[list_of_point[:, 2] > eta] = 0.0
        return accelerations


class Stokes2ndwave(Airywave):

    def calc_elevation(self, position, global_time):
        """
        A private function. \n
        :param position: [np.array].shape=(n,3) coordinates Unit: [m]. The position of the point which you want to know the wave surface elevation.
        :param global_time: [float] Unit: [s].
        :return: scale or [np.array].shape=(n,) Unit: [m]. The sea surface level in Z direction.
        """
        self.calc_theta(position, global_time)
        return self.wave_Height / 2 * np.cos(self.theta) + self.wave_Height/8.0*self.pi_h_l*np.cosh(self.kd)/pow(np.sinh(self.kd), 3)*(2+np.cosh(2*self.kd)) * np.cos(2*self.theta)

    def get_elevations_at_one_node(self, position, time_list):
        """
        Public function.\n
        :param position: [np.array].shape=(n,3) Unit: [m]. The position of one node
        :param time_list: [np.array].shape=(n,1) | Uint: [s]. The time sequence for getting the elevations
        :return: Get a list of elevations at one position with a time sequence \n
        """
        self.calc_theta(position, time_list)
        
        return self.wave_Height / 2 * np.cos(self.theta) + self.wave_Height/8.0*self.pi_h_l*np.cosh(self.kd)/pow(np.sinh(self.kd), 3)*(2+np.cosh(2*self.kd)) * np.cos(2*self.theta)


    def get_velocity_at_one_node(self, position, global_time, irregularwaves=False):
        """
        :param position: [np.array].shape=(1,3) | Unit: [m]. The position of the point which you want to know the wave velocity
        :param global_time: [float] Unit: [s]. The time which you want to know the wave velocity
        :return:  [np.array].shape=(1,3) Unit: [m/s]. A numpy array of the velocity at the targeted point.
        """
        return 0

    def get_acceleration_at_one_node(self, position, global_time, irregularwaves=False):
        """
        :param position: [np.array].shape=(1,3) | Unit: [m]. The position of the point which you want to know the wave acceleration.
        :param global_time: [float] Unit: [s]. The time which you want to know the wave velocity
        :return: [np.array].shape=(1,3) Unit: [m/s]. A numpy array of the acceleration at the targeted point.
        """
        return 0

    def get_elevation_at_nodes(self, list_of_point, global_time):
        """
        Public function.\n
        :param list_of_point: [np.array].shape=(n,3) Unit: [m]. A list of node positions
        :param global_time: time [s] \n
        :return: Get a list of elevation at a list of point \n
        """
        return self.calc_elevation(list_of_point, global_time)

    def get_velocity_at_nodes(self, list_of_point, global_time, irregularwaves=False):
        """
        Public function.\n
        :param list_of_point:  [np.array].shape=(n,3) Unit: [m]. A list of points's positions
        :param global_time: [float] Unit: [s]. Physical time.
        :return: Get a list of velocity at a list of point\n
        """
        self.calc_theta(list_of_point, global_time)
        eta = self.calc_elevation(list_of_point, global_time)

        positions_z = list_of_point[:, 2]
        velocities = np.zeros_like(list_of_point)
        a1 = 0.75*self.pi_h_t*self.pi_h_l
        a2 = self.wave_k * (positions_z + self.water_Depth)

        u1 = self.pi_h_t * np.cosh(a2) / np.sinh(self.kd) * np.cos(self.theta)
        u2 = a1*np.cosh(2*a2) / pow(np.sinh(self.kd), 4) * np.cos(2*self.theta)

        velocities[:, 0] = np.cos(self.wave_beta) * (u1+u2)
        velocities[:, 1] = np.sin(self.wave_beta) * (u1+u2)

        u1 = self.pi_h_t * np.sinh(a2) * np.sin(self.theta) / np.sinh(self.kd)
        u2 = a1*np.sinh(2*a2) / pow(np.sinh(self.kd), 4) * np.sin(2*self.theta)
        velocities[:, 2] = u1 + u2
        if not irregularwaves:
            velocities[list_of_point[:, 2] > eta] = 0.0
        return velocities

    def get_acceleration_at_nodes(self, list_of_point, global_time, irregularwaves=False):
        """
        Public function.\n
        :param list_of_point: [np.array].shape=(n,3) Unit: [m]. A list of points's positions
        :param global_time: time [s] \n
        :return: Get a list of acceleration at a list of point \n
        """
        self.calc_theta(list_of_point, global_time)
        eta = self.calc_elevation(list_of_point, global_time)

        positions_z = list_of_point[:, 2]
        accelerations = np.zeros_like(list_of_point)

        a1 = 3*self.pi_h_t_2*self.pi_h_l
        kdz = self.wave_k * (positions_z + self.water_Depth)
        u1 = 2 * self.pi_h_t_2 * np.cosh(kdz) / np.sinh(self.kd) * np.sin(self.theta)
        u2 = a1*np.cosh(2*kdz) / pow(np.sinh(self.kd), 4) * np.sin(2*self.theta)

        accelerations[:, 0] = np.cos(self.wave_beta) * (u1+u2)
        accelerations[:, 1] = np.sin(self.wave_beta) * (u1+u2)

        u1 = -2 * self.pi_h_t_2 * np.sinh(kdz) / np.sinh(self.kd) * np.cos(self.theta)
        u2 = -a1*np.sinh(2*kdz) / pow(np.sinh(self.kd), 4) * np.cos(2*self.theta)
        accelerations[:, 2] = u1+u2

        if not irregularwaves:
            accelerations[list_of_point[:, 2] > eta] = 0.0
        return accelerations


if __name__ == "__main__":
    pass
