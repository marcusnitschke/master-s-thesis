import numpy as np
import matplotlib.pyplot as plt
from launchsim import Rocket
from TVC import Rocket_TVC

def trajectory(mk, tvc):
    plt.figure(figsize=[12.8, 9.6])
    plt.subplot(221)
    plt.plot(mk.t, mk.r.T[2], c='tab:red', label="Altitude")
    plt.plot(mk.t, tvc.r.T[2], ':', c='tab:red', label="Altitude (TVC)")
    plt.plot(mk.t, mk.r.T[0], c='tab:blue', label="East")
    plt.plot(mk.t, tvc.r.T[0], ':', c='tab:blue', label="East (TVC)")
    plt.plot(mk.t, mk.r.T[1], c='tab:green', label="North")
    plt.plot(mk.t, tvc.r.T[1], ':', c='tab:green', label="North (TVC)")
    plt.xlabel("time (s)", fontsize=14)
    plt.ylabel("distance (m)", fontsize=14)
    plt.legend(bbox_to_anchor=(0,1.02,1,0.2), loc="lower left", mode="expand",
        ncol=3, fontsize=12)
    plt.grid()
    plt.subplot(222)
    plt.plot(mk.t, abs(mk.v.T[0])+abs(mk.v.T[1])+abs(mk.v.T[2]),
        c='tab:orange', label="Velocity")
    plt.plot(mk.t, abs(tvc.v.T[0])+abs(tvc.v.T[1])+abs(tvc.v.T[2]), ':',
        c='tab:orange', label="Velocity (TVC)")
    plt.plot(mk.t, abs(mk.a.T[0])+abs(mk.a.T[1])+abs(mk.a.T[2]),
        c='tab:purple', label="Acceleration")
    plt.plot(mk.t, abs(tvc.a.T[0])+abs(tvc.a.T[1])+abs(tvc.a.T[2]), ':',
        c='tab:purple', label="Acceleration (TVC)")
    plt.xlabel("time (s)", fontsize=14)
    plt.ylabel("velocity (m/s) and acceleration (m/s^2)", fontsize=14)
    plt.legend(bbox_to_anchor=(0,1.02,1,0.2), loc="lower left", mode="expand",
        ncol=2, fontsize=12)
    plt.grid()
    plt.subplot(223)
    plt.plot(mk.t[1:], -mk.angle_attack.T[1][1:]+np.pi/2, c='tab:cyan',
        label="Without TVC")
    plt.plot(mk.t[1:], -tvc.angle_attack.T[1][1:]+np.pi/2, ':', c='tab:cyan',
        label="With TVC")
    plt.xlabel("time (s)", fontsize=14)
    plt.ylabel("vertical orientation (rad)", fontsize=14)
    plt.legend(bbox_to_anchor=(0,-0.35,1,0.2), loc="upper left", mode="expand",
        ncol=3, fontsize=12)
    plt.grid()
    plt.subplot(224)
    plt.plot(mk.t[1:605], 604*[-mk.angle_attack.T[1][2]+np.pi/2], '-.',
        c='tab:gray', label="Angle at launch")
    plt.plot(mk.t[1:605], -mk.angle_attack.T[1][1:605]+np.pi/2, c='tab:cyan',
        label="Without TVC")
    plt.plot(mk.t[1:605], -tvc.angle_attack.T[1][1:605]+np.pi/2, ':',
        c='tab:cyan', label="With TVC")
    plt.xlabel("time (s)", fontsize=14)
    plt.ylabel("vertical orientation (rad)", fontsize=14)
    plt.legend(bbox_to_anchor=(0,-0.35,1,0.2), loc="upper left", mode="expand",
        ncol=3, fontsize=12)
    plt.grid()
    plt.tight_layout()

def extreme(mk, tvc, tvce1, tvce2, tvce3):
    plt.figure(figsize=[12.8, 9.6])
    plt.subplot(211)
    plt.plot(mk.t, tvc.r.T[2], c='tab:red', label="Altitude (0 m/s)")
    plt.plot(mk.t, tvc.r.T[0], c='tab:blue', label="East (0 m/s)")
    plt.plot(mk.t, tvc.r.T[1], c='tab:green', label="North (0 m/s)")
    plt.plot(mk.t, tvce1.r.T[2], '--', c='tab:red', label="Altitude (17 m/s)")
    plt.plot(mk.t, tvce1.r.T[0], '--', c='tab:blue', label="East (17 m/s)")
    plt.plot(mk.t, tvce1.r.T[1], '--', c='tab:green', label="North (17 m/s)")
    plt.plot(mk.t, tvce2.r.T[2], ':', c='tab:red', label="Altitude (21 m/s)")
    plt.plot(mk.t, tvce2.r.T[0], ':', c='tab:blue', label="East (21 m/s)")
    plt.plot(mk.t, tvce2.r.T[1], ':', c='tab:green', label="North (21 m/s)")
    plt.plot(mk.t, tvce3.r.T[2], '-.', c='tab:red', label="Altitude (25 m/s)")
    plt.plot(mk.t, tvce3.r.T[0], '-.', c='tab:blue', label="East (25 m/s)")
    plt.plot(mk.t, tvce3.r.T[1], '-.', c='tab:green', label="North (25 m/s)")
    plt.xlabel("time (s)", fontsize=14)
    plt.ylabel("distance (m)", fontsize=14)
    plt.legend(bbox_to_anchor=(0,1.02,1,0.2), loc="lower left", mode="expand",
        ncol=4, fontsize=12)
    plt.grid()
    plt.subplot(223)
    plt.plot(mk.t[1:], -tvc.angle_attack.T[1][1:]+np.pi/2, c='tab:blue',
        label="0 m/s")
    plt.plot(mk.t[1:], -tvce3.angle_attack.T[1][1:]+np.pi/2, '-.',
        c='tab:purple', label="25 m/s")
    plt.plot(mk.t[1:], -tvce1.angle_attack.T[1][1:]+np.pi/2, '--',
        c='tab:green', label="17 m/s")
    plt.plot(mk.t[1:], -tvce2.angle_attack.T[1][1:]+np.pi/2, ':',
        c='tab:red', label="21 m/s")
    plt.xlabel("time (s)", fontsize=14)
    plt.ylabel("vertical orientation (rad)", fontsize=14)
    plt.legend(bbox_to_anchor=(0,-0.35,1,0.2), loc="upper left", mode="expand",
        ncol=3, fontsize=12)
    plt.grid()
    plt.subplot(224)
    plt.plot(mk.t[1:605], -tvc.angle_attack.T[1][1:605]+np.pi/2, c='tab:blue',
        label="0 m/s")
    plt.plot(mk.t[1:605], -tvce3.angle_attack.T[1][1:605]+np.pi/2, '-.',
        c='tab:purple', label="25 m/s")
    plt.plot(mk.t[1:605], -tvce1.angle_attack.T[1][1:605]+np.pi/2, '--',
        c='tab:green', label="17 m/s")
    plt.plot(mk.t[1:605], 604*[-mk.angle_attack.T[1][2]+np.pi/2], c='tab:gray',
        label="Angle at launch")
    plt.plot(mk.t[1:605], -tvce2.angle_attack.T[1][1:605]+np.pi/2, ':',
        c='tab:red', label="21 m/s")
    plt.xlabel("time (s)", fontsize=14)
    plt.ylabel("vertical orientation (rad)", fontsize=14)
    plt.legend(bbox_to_anchor=(0,-0.35,1,0.2), loc="upper left", mode="expand",
        ncol=3, fontsize=12)
    plt.grid()
    plt.tight_layout()

inputs = {"tmax":90, "wind_speed":0, "wind_ang":0, "dry_mass":9.85,
    "wet_mass":18.554, "length":2710, "cd":0.75, "cl":0.15,
    "critical_angle":20, "hcm":710, "hcp":510, "radius":51.5,
    "thrustforce":2529, "burntime":6.04}

#No wind.
mk1 = Rocket(launch_ang=(80, 90), **inputs)
mk1.launch()
tvc1 = Rocket_TVC(80, 600, 30, launch_ang=(80, 90), **inputs)
tvc1.launch()

#Some wind.
inputs.update(wind_speed=3, wind_ang=90)
mk2 = Rocket(launch_ang=(80, 90), **inputs)
mk2.launch()
tvc2 = Rocket_TVC(80, 600, 30, launch_ang=(80, 90), **inputs)
tvc2.launch()

#Strong wind.
inputs.update(wind_speed=8, wind_ang=90)
mk3 = Rocket(launch_ang=(80, 90), **inputs)
mk3.launch()
tvc3 = Rocket_TVC(80, 600, 30, launch_ang=(80, 90), **inputs)
tvc3.launch()

#Very strong wind.
inputs.update(wind_speed=14, wind_ang=90)
mk4 = Rocket(launch_ang=(80, 90), **inputs)
mk4.launch()
tvc4 = Rocket_TVC(80, 600, 30, launch_ang=(80, 90), **inputs)
tvc4.launch()

#Crosswind.
inputs.update(wind_speed=8, wind_ang=0)
mk5 = Rocket(launch_ang=(80, 90), **inputs)
mk5.launch()
tvc5 = Rocket_TVC(80, 600, 30, launch_ang=(80, 90), **inputs)
tvc5.launch()

#Tailwind.
inputs.update(wind_speed=8, wind_ang=270)
mk6 = Rocket(launch_ang=(80, 90), **inputs)
mk6.launch()
tvc6 = Rocket_TVC(80, 600, 30, launch_ang=(80, 90), **inputs)
tvc6.launch()

#Extreme conditions in headwind
inputs.update(wind_speed=17, wind_ang=90)
tvc6 = Rocket_TVC(80, 600, 30, launch_ang=(80, 90), **inputs)
tvc6.launch()
inputs.update(wind_speed=21, wind_ang=90)
tvc7 = Rocket_TVC(80, 600, 30, launch_ang=(80, 90), **inputs)
tvc7.launch()
inputs.update(wind_speed=25, wind_ang=90)
tvc8 = Rocket_TVC(80, 600, 30, launch_ang=(80, 90), **inputs)
tvc8.launch()

#Extreme conditions in crosswind
inputs.update(wind_speed=17, wind_ang=0)
tvc9 = Rocket_TVC(80, 600, 30, launch_ang=(80, 90), **inputs)
tvc9.launch()
inputs.update(wind_speed=21, wind_ang=0)
tvc10 = Rocket_TVC(80, 600, 30, launch_ang=(80, 90), **inputs)
tvc10.launch()
inputs.update(wind_speed=25, wind_ang=0)
tvc11 = Rocket_TVC(80, 600, 30, launch_ang=(80, 90), **inputs)
tvc11.launch()

#Extreme conditions in tailwind
inputs.update(wind_speed=17, wind_ang=270)
tvc12 = Rocket_TVC(80, 600, 30, launch_ang=(80, 90), **inputs)
tvc12.launch()
inputs.update(wind_speed=21, wind_ang=270)
tvc13 = Rocket_TVC(80, 600, 30, launch_ang=(80, 90), **inputs)
tvc13.launch()
inputs.update(wind_speed=25, wind_ang=270)
tvc14 = Rocket_TVC(80, 600, 30, launch_ang=(80, 90), **inputs)
tvc14.launch()

trajectory(mk1, tvc1)
plt.savefig("tvc_no_wind")
plt.clf()

trajectory(mk2, tvc2)
plt.savefig("tvc_3mps_headwind")
plt.clf()
trajectory(mk3, tvc3)
plt.savefig("tvc_8mps_headwind")
plt.clf()
trajectory(mk4, tvc4)
plt.savefig("tvc_14mps_headwind")
plt.clf()
trajectory(mk5, tvc5)
plt.savefig("tvc_8mps_crosswind")
plt.clf()
trajectory(mk6, tvc6)
plt.savefig("tvc_8mps_tailwind")
plt.clf()

extreme(mk1, tvc1, tvc6, tvc7, tvc8)
plt.savefig("tvc_extreme_headwind")
plt.clf()
extreme(mk1, tvc1, tvc9, tvc10, tvc11)
plt.savefig("tvc_extreme_crosswind")
plt.clf()
extreme(mk1, tvc1, tvc12, tvc13, tvc14)
plt.savefig("tvc_extreme_tailwind")
plt.clf()

plt.figure(figsize=[9.6, 7.2])
plt.subplot(211)
plt.plot(mk2.t[:1200], mk2.fthrust[:1200, 0], '-', c='tab:blue', label="thrust")
plt.plot(tvc2.t[:1200], tvc2.fthrust[:1200, 0], '--', c='tab:blue', label="thrust (TVC)")
plt.plot(mk2.t[:1200], mk2.flift[:1200, 0], '-', c='tab:green', label="lift")
plt.plot(tvc2.t[:1200], tvc2.flift[:1200, 0], '--', c='tab:green', label="lift (TVC)")
plt.plot(mk2.t[:1200], mk2.fdrag[:1200, 0], '-', c='tab:red', label="drag")
plt.plot(tvc2.t[:1200], tvc2.fdrag[:1200, 0], '--', c='tab:red', label="drag (TVC)")
plt.ylabel("force (N) in x-axis", fontsize=14)
plt.legend(fontsize=12, bbox_to_anchor=(0,1.02,1,0.2), ncol=6, loc="lower left")
plt.grid()
plt.subplot(212)
plt.plot(mk2.t[:1200], mk2.fthrust[:1200, 2], '-', c='tab:blue')
plt.plot(mk2.t[:1200], mk2.fdrag[:1200, 2], '-', c='tab:red')
plt.plot(mk2.t[:1200], mk2.flift[:1200, 2], '-', c='tab:green')
plt.plot(tvc2.t[:1200], tvc2.fthrust[:1200, 2], '--', c='tab:blue')
plt.plot(tvc2.t[:1200], tvc2.fdrag[:1200, 2], '--', c='tab:red')
plt.plot(tvc2.t[:1200], tvc2.flift[:1200, 2], '--', c='tab:green')
plt.xlabel("time (s)", fontsize=14)
plt.ylabel("force (N) in z-axis", fontsize=14)
plt.grid()
plt.tight_layout()
plt.savefig("forces_3mps")
plt.clf()
