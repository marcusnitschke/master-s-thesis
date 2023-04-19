import numpy as np
import matplotlib.pyplot as plt
from launchsim import Rocket
import pandas as pd

def comparison(df, mk, dfr):
    plt.figure(figsize=[12.8, 9.6])
    plt.subplot(221)
    plt.plot(df['Time (s)'], df['Altitude (m)'], '--', c='tab:red',
        label="Altitude OR")
    plt.plot(dfr['Time (s)'], dfr['Altitude (m)'], ':', c='tab:red',
        label="Altitude OR (roll)")
    plt.plot(mk.t, mk.r.T[2], c='tab:red', label="Altitude Python")
    plt.plot(df['Time (s)'], df['Position East of launch (m)'], '--',
        c='tab:blue', label="East OR")
    plt.plot(dfr['Time (s)'], dfr['Position East of launch (m)'], ':',
        c='tab:blue', label="East OR (roll)")
    plt.plot(mk.t, mk.r.T[0], c='tab:blue', label="East Python")
    plt.plot(df['Time (s)'], df['Position North of launch (m)'], '--',
        c='tab:green', label="North OR")
    plt.plot(dfr['Time (s)'], dfr['Position North of launch (m)'], ':',
        c='tab:green', label="North OR (roll)")
    plt.plot(mk.t, mk.r.T[1], c='tab:green', label="North Python")
    plt.xlabel("time (s)", fontsize=14)
    plt.ylabel("distance (m)", fontsize=14)
    plt.legend(bbox_to_anchor=(0,1.02,1,0.2), loc="lower left", mode="expand",
        ncol=3, fontsize=12)
    plt.grid()
    plt.subplot(222)
    plt.plot(df['Time (s)'], df['Total velocity (m/s)'], '--', c='tab:orange',
        label="Velocity OR")
    plt.plot(dfr['Time (s)'], dfr['Total velocity (m/s)'], ':', c='tab:orange',
        label="Velocity OR (roll)")
    plt.plot(mk.t, abs(mk.v.T[0])+abs(mk.v.T[1])+abs(mk.v.T[2]),
        c='tab:orange', label="Velocity Python")
    plt.plot(df['Time (s)'], df['Total acceleration (m/s)'], '--',
        c='tab:purple', label="Acceleration OR")
    plt.plot(dfr['Time (s)'], dfr['Total acceleration (m/s)'], ':',
        c='tab:purple', label="Acceleration OR (roll)")
    plt.plot(mk.t, abs(mk.a.T[0])+abs(mk.a.T[1])+abs(mk.a.T[2]),
        c='tab:purple', label="Acceleration Python")
    plt.xlabel("time (s)", fontsize=14)
    plt.ylabel("velocity (m/s) and acceleration (m/s^2)", fontsize=14)
    plt.legend(bbox_to_anchor=(0,1.02,1,0.2), loc="lower left", mode="expand",
        ncol=2, fontsize=12)
    plt.grid()
    plt.subplot(223)
    plt.plot(df['Time (s)'], df['Vertical orientation (zenith) (rad)'], '--',
        c='tab:cyan', label="OR")
    plt.plot(dfr['Time (s)'], dfr['Vertical orientation (zenith) (rad)'], ':',
        c='tab:cyan', label="OR (roll)")
    plt.plot(mk.t[1:], -mk.angle_attack.T[1][1:]+np.pi/2, c='tab:cyan',
        label="Python")
    plt.xlabel("time (s)", fontsize=14)
    plt.ylabel("vertical orientation (rad)", fontsize=14)
    plt.legend(bbox_to_anchor=(0,-0.35,1,0.2), loc="upper left", mode="expand",
        ncol=3, fontsize=12)
    plt.grid()
    plt.subplot(224)
    plt.plot(df['Time (s)'][:654],
        df['Vertical orientation (zenith) (rad)'][:654], '--', c='tab:cyan',
        label="OR")
    plt.plot(mk.t[1:605], 604*[-mk.angle_attack.T[1][2]+np.pi/2], '-.',
        c='tab:gray', label="Angle at launch")
    plt.plot(dfr['Time (s)'][:5730],
        dfr['Vertical orientation (zenith) (rad)'][:5730], ':', c='tab:cyan',
        label="OR (roll)")
    plt.plot(mk.t[1:605], -mk.angle_attack.T[1][1:605]+np.pi/2, c='tab:cyan',
        label="Python")
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

#Results from open rocket without roll.
df1 = pd.read_csv('csv/openrocket_0mps.csv', sep=',')
df2 = pd.read_csv('csv/openrocket_3mps90deg.csv', sep=',')
df3 = pd.read_csv('csv/openrocket_8mps90deg.csv', sep=',')
df4 = pd.read_csv('csv/openrocket_14mps90deg.csv', sep=',')
df5 = pd.read_csv('csv/openrocket_8mps0deg.csv', sep=',')
df6 = pd.read_csv('csv/openrocket_8mps270deg.csv', sep=',')

#Results from open rocket with roll.
dfr1 = pd.read_csv('csv/openrocket_0mps_roll.csv', sep=',')
dfr2 = pd.read_csv('csv/openrocket_3mps90deg_roll.csv', sep=',')
dfr3 = pd.read_csv('csv/openrocket_8mps90deg_roll.csv', sep=',')
dfr4 = pd.read_csv('csv/openrocket_14mps90deg_roll.csv', sep=',')
dfr5 = pd.read_csv('csv/openrocket_8mps0deg_roll.csv', sep=',')
dfr6 = pd.read_csv('csv/openrocket_8mps270deg_roll.csv', sep=',')

#Python, no wind.
mk1 = Rocket(launch_ang=(80, 90), **inputs)
mk1.launch()

#Python, some wind.
inputs.update(wind_speed=3, wind_ang=90)
mk2 = Rocket(launch_ang=(80, 90), **inputs)
mk2.launch()

#Python, strong wind.
inputs.update(wind_speed=8, wind_ang=90)
mk3 = Rocket(launch_ang=(80, 90), **inputs)
mk3.launch()

#Python, very strong wind.
inputs.update(wind_speed=14, wind_ang=90)
mk4 = Rocket(launch_ang=(80, 90), **inputs)
mk4.launch()

#Python, crosswind.
inputs.update(wind_speed=8, wind_ang=0)
mk5 = Rocket(launch_ang=(80, 90), **inputs)
mk5.launch()

#Python, tailwind.
inputs.update(wind_speed=8, wind_ang=270)
mk6 = Rocket(launch_ang=(80, 90), **inputs)
mk6.launch()

comparison(df1, mk1, dfr1)
plt.savefig("Comparison_no_wind")
plt.clf()
comparison(df2, mk2, dfr2)
plt.savefig("Comparison_3mps_headwind")
plt.clf()
comparison(df3, mk3, dfr3)
plt.savefig("Comparison_8mps_headwind")
plt.clf()
comparison(df4, mk4, dfr4)
plt.savefig("Comparison_14mps_headwind")
plt.clf()
comparison(df5, mk5, dfr5)
plt.savefig("Comparison_8mps_crosswind")
plt.clf()
comparison(df5, mk6, dfr6)
plt.savefig("Comparison_8mps_tailwind")
plt.clf()

plt.figure(figsize=[12.8, 9.6])
plt.subplot(331)
plt.plot(df1['Time (s)'], df1['Roll rate (rpm)'], c='tab:blue',
    label="OpenRocket (0 deg fin cant)")
plt.plot(dfr1['Time (s)'], dfr1['Roll rate (rpm)'], c='tab:green')
plt.xlabel("time (s)", fontsize=14)
plt.ylabel("roll rate (rpm)", fontsize=14)
plt.legend(fontsize=12, bbox_to_anchor=(0,1.02,1,0.2), loc="lower left")
plt.grid()
plt.subplot(332)
plt.plot(df1['Time (s)'], df1['Air pressure (Pa)'], c='tab:blue')
plt.plot(dfr1['Time (s)'], dfr1['Air pressure (Pa)'], c='tab:green',
    label="OpenRocket (3 deg fin cant)")
plt.plot(mk1.t, np.array(mk1._rho)*(8.3145*(288.15-0.0065*mk1.r.T[2])/0.029),
    c="tab:red")
plt.xlabel("time (s)", fontsize=14)
plt.ylabel("air pressure (Pa)", fontsize=14)
plt.legend(fontsize=12, bbox_to_anchor=(0,1.02,1,0.2), loc="lower left")
plt.grid()
plt.subplot(333)
plt.plot(df1['Time (s)'], df1['Longitudinal moment of inertia (kgm)'],
    c='tab:blue')
plt.plot(dfr1['Time (s)'], dfr1['Longitudinal moment of inertia (kgm)'],
    c='tab:green')
plt.plot(mk1.t, (np.array(mk1._m)*(3*mk1._radius**2+mk1._length**2)/12),
    c="tab:red", label="Python")
plt.xlabel("time (s)", fontsize=14)
plt.ylabel("moment of inertia (kg m^2)", fontsize=14)
plt.legend(fontsize=12, bbox_to_anchor=(0,1.02,1,0.2), loc="lower left")
plt.grid()
plt.subplot(334)
plt.plot(df1['Time (s)'], df1['CP location (cm)'], c='tab:blue')
plt.plot(dfr1['Time (s)'], dfr1['CP location (cm)'], c='tab:green')
plt.plot(mk1.t, len(mk1.t)*[(mk1._length-mk1._hcp)*100], c="tab:red")
plt.xlabel("time (s)", fontsize=14)
plt.ylabel("CP location (cm)", fontsize=14)
plt.grid()
plt.subplot(335)
plt.plot(df1['Time (s)'], df1['CG location (cm)'], c='tab:blue')
plt.plot(dfr1['Time (s)'], dfr1['CG location (cm)'], c='tab:green')
plt.plot(mk1.t, len(mk1.t)*[(mk1._length-mk1._hcm)*100], c="tab:red")
plt.xlabel("time (s)", fontsize=14)
plt.ylabel("CM location (cm)", fontsize=14)
plt.grid()
plt.subplot(336)
plt.plot(df1['Time (s)'], df1['Gravitational acceleration (m/s)'],
    c='tab:blue', label="0 deg fin cant")
plt.plot(dfr1['Time (s)'], dfr1['Gravitational acceleration (m/s)'],
    c='tab:green')
plt.plot(mk1.t, len(mk1.t)*[mk1._g], c="tab:red")
plt.xlabel("time (s)", fontsize=14)
plt.ylabel("g", fontsize=14)
plt.grid()
plt.subplot(337)
plt.plot(df1['Time (s)'], df1['Pressure drag coefficient (?)'], c='tab:blue')
plt.plot(dfr1['Time (s)'], dfr1['Pressure drag coefficient (?)'],
    c='tab:green')
plt.plot(mk1.t, len(mk1.t)*[mk1._cl], c="tab:red")
plt.xlabel("time (s)", fontsize=14)
plt.ylabel("pressure drag coefficient", fontsize=14)
plt.grid()
plt.subplot(338)
plt.plot(df1['Time (s)'], df1['Drag coefficient (?)'], c='tab:blue')
plt.plot(dfr1['Time (s)'], dfr1['Drag coefficient (?)'], c='tab:green')
plt.plot(mk1.t, len(mk1.t)*[mk1._cd], c="tab:red")
plt.xlabel("time (s)", fontsize=14)
plt.ylabel("drag coefficient", fontsize=14)
plt.grid()
plt.tight_layout()
plt.savefig("compare_parameters")
plt.clf()
