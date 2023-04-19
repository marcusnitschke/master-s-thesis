# Marcus Steen Nitschke
## Python simulations for my master thesis
All .py files require Numpy and Matplotlib to be installed, and can be run directly from the terminal. The .csv files are exported from OpenRocket. "rocket verification.ork" is the design file in OpenRocket.


#### launchsim.py
Implements the base model. This is mainly used as import, but can be run to plot the trajectory of a rocket without TVC.

#### TVC.py
Imports the base model, and implement TVC. This is mainly used as import, but can be run to plot vertical orientation.

#### comparison.py
Imports the base model along with the csv-files from OpenRocket, and plot several figures comparing the model to OpenRocket.

#### results.py
Imports the base model and the implementation of TVC, and plot several figures showing the effect of TVC in different wind conditions.
