import numpy as np
import matplotlib.pyplot as plt
from launchsim import Rocket

class Rocket_TVC(Rocket):

    def __init__(self, kp, ki, kd, launch_ang, max_angle=5, **kwargs):
        """Initialize class.
        Args:
            kp (float): proportional coefficient.
            ki (float): integral coefficient.
            kd (float): derivative coefficient.
            launch_ang (touple): Launch angle (deg [theta, phi]).
            max_angle (float): Maximum angle of the thrust vector in any
                direction (deg).
            **kwargs: input for super class.
        """
        super().__init__(launch_ang, **kwargs)
        self._kp = kp
        self._ki = ki
        self._kd = kd
        self._max_angle = max_angle
        self._prev_error = np.zeros(2)
        self._total_error = 0
        self._desired_ang = np.array([self.r[0][3], self.r[0][4]])
        return

    def thrust(self, t, r):
        """Thrust as a vector pointing from the back of the rocket with an
         angle depending on the TVC from positive z in the local frame of
         reference, converted into the global frame of reference.
        Args:
            t (float): Time since initialization (s).
            r (np.array): Positional vector
                (East, North, altitude, pitch, yaw).
        Returns:
            np.array: Force from thrust as a vector
                (East, North, altitude, pitch, yaw).
        """
        if t < self._burntime:
            #Thrust is now a vector in the local frame of reference.
            thrust = self._thrustforce*np.array([0, 0, 1])
            l_thrust = self.rotate(self.pid(r), thrust)
            #Convert the vector to the global frame of reference.
            g_thrust = self.rotate(r[3:], l_thrust)
            return np.array((g_thrust[0], g_thrust[1], g_thrust[2],
                             -l_thrust[1]*self._hcm, l_thrust[0]*self._hcm))
        else:
            return np.zeros(5)

    def pid(self, r):
        """Calculate the angle of the thrust vector using PID.
        Args:
            r (np.array): Positional vector (East, North, altitude, pitch, yaw).
        Returns:
            np.array: Angle of thrust vector (rad [pitch, yaw]).
        """
        error = self._desired_ang - r[3:]
        self._total_error += error*self._dt
        proportional = self._kp*error
        integral = self._ki*self._total_error
        derivative = self._kd*(error - self._prev_error)
        u = proportional + integral + derivative

        for i in range(len(u)):
            if u[i] > self._max_angle*(np.pi/180):
                u[i] = self._max_angle*(np.pi/180)
            elif u[i] < -self._max_angle*(np.pi/180):
                u[i] = -self._max_angle*(np.pi/180)
        self._prev_error = error
        return u

if __name__ == "__main__":
    inputs = {"tmax":90, "wind_speed":0, "wind_ang":0, "dry_mass":9.85,
        "wet_mass":18.554, "length":2710, "cd":0.75, "cl":0.15,
        "critical_angle":20, "hcm":710, "hcp":510, "radius":51.5,
        "thrustforce":2529, "burntime":6.04}

    kp, ki, kd = 80, 600, 30
    mk1 = Rocket_TVC(kp, ki, kd, launch_ang=(80, 90), **inputs)
    mk1.launch()

    plt.figure(figsize=[12.8, 9.6])
    plt.plot(mk1.t[1:605], -mk1.angle_attack.T[1][1:605]+np.pi/2,
        c='tab:orange', label=f"KP={kp}, KI={ki}, KD={kd}")
    plt.plot(mk1.t[1:605], 604*[-mk1.angle_attack.T[1][2]+np.pi/2], '-.',
        c='tab:red', label="Set value")
    plt.xlabel("time (s)", fontsize=20)
    plt.ylabel("vertical orientation (rad)", fontsize=20)
    plt.legend(bbox_to_anchor=(0,1.02,1,0.2), loc="lower left", mode="expand",
        ncol=3, fontsize=16)
    plt.grid()
    plt.tight_layout()
    plt.show()
