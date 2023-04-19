import numpy as np
import matplotlib.pyplot as plt

class Rocket():
    def __init__(self, launch_ang, tmax, wind_speed, wind_ang, dry_mass,
                 wet_mass, length, cd, cl, critical_angle, hcm, hcp, radius,
                 thrustforce, burntime, dt=0.01, rail_length=5.0):
        """Initialize class.
        Args:
            launch_ang (touple): Launch angle (deg [elevation, azimuth]).
            tmax (int): Length of simulation (s).
            wind_speed (float): Wind speed (m/s).
            wind_ang (float): Wind direction (deg).
            dry_mass (float): Dry mass (kg).
            wet_mass (float): Wet mass (kg).
            length (float): Length of the rocket (mm).
            cd (float): Drag coefficient.
            cl (float): Lift coefficient.
            hcm (float): Height of the centre of mass (mm).
            hcp (float): Height of the centre of pressure (mm).
            radius (float): Radius of nose cone (mm).
            thrustforce (float): Force of thrust (N).
            burntime (float): Burntime (s).
            dt (float): Time step (s). Default to 0.01.
            rail_length (float): Lenght of launch rail (m). Default to 5.0.
        """
        self._dt = dt
        self.t = np.linspace(0, tmax, tmax*int(1/self._dt)+1)
        self.r, self.v, self.a = np.zeros((3, tmax*int(1/self._dt)+1, 5))

        #Set initial launch angle in Euler angles pitch and yaw.
        self.r[0][3] = -((np.pi/2)-launch_ang[0]*np.pi/180)*np.cos(
            launch_ang[1]*np.pi/180)
        self.r[0][4] = ((np.pi/2)-launch_ang[0]*np.pi/180)*np.sin(
            launch_ang[1]*np.pi/180)
        for i in range(2):
            if abs(self.r[0][i+3]) < 10E-14:
                self.r[0][i+3] = 0

        self._rodlenght = rail_length
        self._burntime = burntime
        self._thrustforce = thrustforce
        self._mt = wet_mass
        self._mf = wet_mass-dry_mass
        self._m = [wet_mass]
        self._cd = cd
        self._cl = cl
        self._hcm = hcm/1000
        self._hcp = hcp/1000
        self._wind = -wind_speed*np.array([np.sin(wind_ang*np.pi/180),
            np.cos(wind_ang*np.pi/180), 0])
        self._side_area = (length/1000)*(radius/1000)*2
        self._front_area = np.pi*((radius/1000)**2)
        self._critical_angle = critical_angle*np.pi/180
        self._radius = radius/1000
        self._length = length/1000
        self.angle_attack = [np.zeros(3)]
        self.fthrust = [np.zeros(3)]
        self.fdrag = [np.zeros(3)]
        self.flift = [np.zeros(3)]

        self._g = 9.80665 #Gravitational constant (m/s^2).
        self._p0 = 101325 #standard absolute atmospheric pressure (Pa).
        self._air_molar = 0.0289654 #Molar mass of dry air (kg/mol).
        self._gas_constant = 8.314463 #Ideal gas constant (J/(mol*k)).
        self._t0 = 288.15 #Absolute temperature at sea level (K).
        self._t_laps = 0.0065 #Temperature laps rate (K/m).
        self._rho = [1.225] #Air density (kg/m^3).

        #First set of calculations for the air density.
        self._c = (self._p0*self._air_molar/(self._gas_constant*self._t0))
        self._exp = ((self._g*self._air_molar/(
            self._gas_constant*self._t_laps))-1)
        return

    def rotate(self, rot, vector, inverse=False):
        """
        Rotate a vector counter clockwise by its x- and y-axis.

        Args:
            rot (np.array): Rotation (pitch, yaw).
            vector (np.array): input vector  (East, North, Altitude).
            inverse (bool): Invert the rotation if 'True'. Default 'False'.
        Returns:
            np.array: Rotated output vector.
        """
        rotation = np.array([[np.cos(rot[1]), np.sin(rot[1])*np.sin(rot[0]),
            np.sin(rot[1])*np.cos(rot[0])], [0, np.cos(rot[0]),
            -np.sin(rot[0])], [-np.sin(rot[1]), np.cos(rot[1])*np.sin(rot[0]),
            np.cos(rot[1])*np.cos(rot[0])]])

        if inverse:
            output = rotation.T @ vector
        else:
            output = rotation @ vector
        return output

    def update(self, t, r, v):
        """Update mass (kg) as a function of time (s), air density (kg/m^3) as a
        function of height (m) and calculate the angle of attack (rad).
        Args:
            t (float): Time since initialization in seconds.
            r (np.array): Positional vector (East, North, altitude, theta, phi).
            v (np.array): Velocity vector (East, North, altitude, theta, phi).
        """
        if t >= self._burntime:
            self._m.append(self._mt - self._mf)
        else:
            self._m.append(self._mt - (self._mf*t/self._burntime))

        self._rho.append(self._c*(1-(self._t_laps*r[2]/self._t0))**self._exp)

        rel_v = v[:3] - self._wind
        motion_theta = 0
        if abs(rel_v.dot(rel_v)) > 0:
            motion_theta = np.arccos(rel_v[2]/np.sqrt(rel_v.dot(rel_v)))
        rocket_theta = np.arccos(np.cos(r[3])*np.cos(r[4]))
        self.angle_attack.append(np.array([motion_theta, rocket_theta,
            motion_theta - rocket_theta]))
        return

    def weight(self, r):
        """Weight as a vector pointing from the center of mass of the rocket
         towards negative z in the main coordinate system.
        Returns:
            np.array: Force from weight as a vector
                (East, North, altitude, theta, phi).
        """
        weight = np.array([0, 0, -self._g*self._m[-1]])
        if np.sqrt(r.dot(r)) < self._rodlenght:
            l_weight = np.array([0, 0,
                self.rotate(r[3:], weight, inverse=True)[2]])
            weight = self.rotate(r[3:], l_weight)
        return np.array([weight[0], weight[1], weight[2], 0, 0])

    def thrust(self, t, r):
        """Thrust as a vector pointing from the back of the rocket towards
         positive z in the local reference frame, converted into the global
         reference frame.
        Args:
            t (float): Time since initialization in seconds.
            r (np.array): Positional vector (East, North, altitude, theta, phi).
        Returns:
            np.array: Force from thrust as a vector
                (East, North, altitude, theta, phi).
        """
        if t < self._burntime:
            #Thrust is now a vector in the local reference frame.
            l_thrust = self._thrustforce*np.array([0, 0, 1])
            #Convert the vector to the global reference frame.
            g_thrust = self.rotate(r[3:], l_thrust)
            return np.array((g_thrust[0], g_thrust[1], g_thrust[2], 0, 0))
        else:
            return np.zeros(5)

    def drag(self, r, v):
        """Drag as a vector pointing from the center of pressure opposite the
         direction of motion in a global reference frame.
        Args:
            r (np.array): Positional vector (East, North, altitude, theta, phi).
            v (np.array): Velocity vector (East, North, altitude, theta, phi).
        Returns:
            np.array: Force from drag as a vector
                (East, North, altitude, theta, phi).
        """
        rel_v = v[:3] - self._wind
        #Drag in global reference frame.
        g_drag = -0.5*self._cd*self._front_area*self._rho[-1]*rel_v*abs(rel_v)

        return np.array([g_drag[0], g_drag[1], g_drag[2], 0, 0])

    def lift(self, r, v):
        """Lift as a vector pointing from the center of pressure perpendicular
         to the direction of motion in a global reference frame.
        Args:
            r (np.array): Positional vector (East, North, altitude, theta, phi).
            v (np.array): Velocity vector (East, North, altitude, theta, phi).
        Returns:
            np.array: Force from lift as a vector
                (East, North, altitude, theta, phi).
        """
        rel_v = v[:3] - self._wind
        orientation = self.rotate(r[3:], np.array([0, 0, 1]))
        dir_lift = np.cross(rel_v, np.cross(orientation, rel_v))

        cl = 0
        if abs(self.angle_attack[-1][2]) < self._critical_angle:
            cl = self._cl*abs(1 - (abs(self.angle_attack[-1][2] -
                 self._critical_angle)/self._critical_angle))

        lift = 0.5*cl*self._side_area*self._rho[-1]*rel_v**2
        mag_lift = np.sqrt(lift.dot(lift))

        if dir_lift.dot(dir_lift) != 0:
            g_lift = mag_lift*dir_lift/np.sqrt(dir_lift.dot(dir_lift))
        else:
            g_lift = np.array([0, 0, 0])

        #Convert the vector to the local reference frame.
        l_lift = self.rotate(r[3:], g_lift, inverse=True)

        f_lift = np.array([g_lift[0], g_lift[1], g_lift[2],
            -l_lift[1]*(self._hcp-self._hcm), l_lift[0]*(self._hcp-self._hcm)])
        if np.sqrt(r.dot(r)) < self._rodlenght:
            g_lift = self.rotate(r[3:], np.array([0, 0, l_lift[2]]))
            f_lift = np.array([g_lift[0], g_lift[1], g_lift[2], 0, 0])
        return f_lift

    def acceleration(self, t, r, v):
        """Calculate the acceleration using Newtons 2. law.
        Args:
            t (float): Time since initialization in seconds.
            r (np.array): Positional vector (East, North, altitude, theta, phi).
            v (np.array): Velocity vector (East, North, altitude, theta, phi).
        Returns:
            np.array: New acceleration vector
                (East, North, altitude, theta, phi).
        """
        f_thrust = self.thrust(t, r)
        f_drag = self.drag(r, v)
        f_weight = self.weight(r)
        f_lift = self.lift(r, v)
        acc = np.array([(f_thrust[:3] + f_drag[:3] + f_lift[:3] +
            f_weight[:3])/self._m[-1]])
        ang_acc = np.array([(f_thrust[3:] + f_lift[3:])/
            (self._m[-1]*(3*self._radius**2+self._length**2)/12)])
        for i in range(2):
            if abs(ang_acc[-1][i]) < 10E-14:
                ang_acc[-1][i] = 0
        self.fthrust.append(f_thrust[:3])
        self.fdrag.append(f_drag[:3])
        self.flift.append(f_lift[:3])
        return np.concatenate((acc, ang_acc), axis=None)

    def launch(self):
        """Solves the differential equation using the step function."""
        for i in range(len(self.t)-1):
            self.a[i+1] = self.acceleration(self.t[i], self.r[i], self.v[i])
            self.r[i+1], self.v[i+1] = self.step(self.t[i], self.r[i],
                                                 self.v[i], self.a[i+1])
        self.angle_attack = np.array(self.angle_attack, dtype=object)
        self.fdrag = np.array(self.fdrag, dtype=object)
        self.fthrust = np.array(self.fthrust, dtype=object)
        self.flift = np.array(self.flift, dtype=object)
        return

    def step(self, t, r, v, a):
        """A modified Forward Euler step function.
        Args:
            t (float): Time since initialization in seconds.
            r (np.array): Positional vector (East, North, altitude, theta, phi).
            v (np.array): Velocity vector (East, North, altitude, theta, phi).
            a (np.array): Acceleration vector
                (East, North, altitude, theta, phi).
        Returns:
            np.array: New positional vector (East, North, altitude, theta, phi).
            np.array: New velocity vector (East, North, altitude, theta, phi).
         """
        v_nxt = np.zeros(5)
        if r[2] >= 0:
            v_nxt = v + self._dt*a
        r_nxt = r + self._dt*v_nxt
        self.update(t, r_nxt, v_nxt)
        return r_nxt, v_nxt

if __name__ == "__main__":
    inputs = {"tmax":90, "wind_speed":0, "wind_ang":0, "dry_mass":9.85,
        "wet_mass":18.554, "length":2710, "cd":0.75, "cl":0.15,
        "critical_angle":20, "hcm":710, "hcp":510, "radius":51.5,
        "thrustforce":2529, "burntime":6.04}

    mk1 = Rocket(launch_ang=(80, 90), **inputs)
    mk1.launch()

    plt.suptitle(
        f"Rocket trajectory with an apogee of {int(max(mk1.r.T[2]))} m.")
    plt.title(
        f"Touchdown at {int(mk1.r[-1, 0])} m East and {int(mk1.r[-1, 1])} m North.")
    plt.plot(mk1.t, mk1.r[:, 0], 'b', label="Position East")
    plt.plot(mk1.t, mk1.r[:, 1], 'g', label="Position North")
    plt.plot(mk1.t, mk1.r[:, 2], 'r', label="Altitude")
    plt.xlabel("time (s)")
    plt.ylabel("distance (m)")
    plt.legend()
    plt.grid()
    plt.show()
