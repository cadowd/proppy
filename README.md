# proppy
PropPy is an interface and series of scripts to calculate the performance and optimise the selection of RC propeller and motor combinations. It is written in python3 and uses a pyQT interface. Propeller performance is calculated using experimental data from the University of Illinois propeller database and the theoretical results supplied on the APC website. Motor performance is calculated using motor constants.

## Installation:
To install, simply download the codes and run main.py. You will need python3, I use the Anaconda distribution. The code has been developed on windows but should run on linux as well.

## First Use:
The main window shows the power consumption of selected combinations. To plot the power consumption (or efficiency), simply select a motor from the list (new motors can be added if you can't find the one you want), a propeller (currently the only propellers selectable are onces for which the data is available) and an aircraft, then click 'Plot'. Multiple motors, propellers and aircraft can be selected for comparison, the plot is cleared by clicking clear plot in the lower right. The aircraft is defined by its drag polar, which can be calculated from exported XFLR5 results. To use XFLR5 results, a type 2 analysis must be run (fixed lift) and the resulting polar exported to a csv file. This file can then be read by PropPy.

Battery and atmosphere settings are changed in the settings drop down menu.

To find an optimal combination, select a number of propellers and motors and only ONE aircraft, then click 'Optimise selected system'. This opens the optimisation window. Here every combination is analysed for the given options, with the results filtered by the selected constraints. Once you click find optimal combinations, the analysis will begin. This can take some time if it is the first time you run the program, afterwards the results are saved so that future analysis can run faster. The plot showed directly after the analysis is a map showing regions where each combination is the most efficient, with the drag polar of the selected aircraft overlain. You can also plot the efficiency and power consumption of specific combinations, which helps to ensure the combination has the performance you want (eg, high acceleration at high speeds, high launch acceleration etc). Please note that the reported flight times and ranges can be between kind of to very optimistic (in the range of approximately 20%), so take them with a grain of salt. Also many other things can affect the power system performance, such as the interference of the craft's wake for pusher configurations or additional drag for tractor configurations.

The hovering craft and Kv chooser windows are not yet fully functional.

## The optimisation window:
Firstly, you will want to check the conditions of optimisation. Are you looking for the motor/prop combination that gives you the longest range or the longest flight time? Do you need a fixed flight speed or should the program select the optimal speed? Do you have any particular constraints on the combinations (the most important constraint is the percentage by which the motor can be overloaded, this should be kept close to 100% to avoid the risk of the motor burning out). Once this has been selected, click the 'find optimal combinations' button.

The first time a combination is calculated may take some time, but the results are saved to make future calculations speedier (the so called 'power surfaces' to be explained in the theory section). The graph that appears after calculation should resemble the one shown below.

![optimisation_plot](/img/example3.png)

This gives an overview of the regions where each combination is the most efficient, with the steady flight speed on the x axis and the thrust on the y axis. The black line is the drag polar of the selected aircraft. To analyse the plot above, we can see that the best choice propeller considered for slow speeds is the 9x4.5, especially for high thrusts, however at the speeds that this plane is likely to fly (around 15m/s) the 8x8 and 8x6 props are much better choices. The 8x6 would be better if you plan to fly at low speeds but with lots of acceleration (as it is more efficient at higher thrusts in this range), however the 8x8 is a better choice for higher speeds. A good combination will dominate in a wide range. Not all the combinations shown in this plot may meet the constraints however, the successful combinations (the top 5 by default) are shown in the bottom window as shown below:

![optimisation_plot](/img/example_5.png)

It should be noted that the reported values of the max flight time and max range assume a constant, level flight the entire time (you will never get this) and depend strongly on the battery capacity defined (set in the main window). However, as long as the battery voltage doesn't change, the best motor/prop combination is unaffected by the final battery capacity used. Now you have a shortlist of interesting options for your power system. To analyse them further, you can click the buttons to plot the combinations maximum efficiency or power consumption.

### Efficiency plots
This plot shows a coloured contour map of the thrust and speed range for the selected motor/prop combination. This can be very useful for estimating the performance of the combination. For example, the plot below shows the results for a 9x6 prop.
![efficiency_plot](/img/efficiency.png)
Here we can see that the combination is most efficient at relatively high speeds due to it's medium pitch. Note that no combination will show high efficiencies at low speeds as the definition of propulsive efficiency depends on a velocity component (see the  [UI prop database definitions](http://m-selig.ae.illinois.edu/props/propDB.html)). This combination seems to be a good choice, as we have a nice region of high efficiency right in the drag bucket of the aircraft (the minimum of the plane's drag polar, and where the most efficient flight happens). Changes of speed in this region don't result in rapid losses of efficiency. We can also see how much acceleration the craft has available. If the plane is flying at 15m/s (with around 2N of drag) and you give it full throttle, with this combination you'll have another 12N of thrust available for acceleration (a bit more than a kilo of force, the actual acceleration depends on the plane's weight but for a 1kg craft this is just over 1m/s/s). At full throttle at 15m/s you'll also be getting about 35% total efficiency out of this system, which isn't too shabby (for sudden acceleration it's rare to hit high efficiencies).

As a comparison, this is the same motor with a 9x4.7 prop.
![efficiency_plot](/img/efficiency2.png)

At 15m/s, you only have about 8N of available accelerative force, or about two thirds of what you got with the 9x6. This will feel much less responsive in the air, and will take longer to change speed as well as requiring more throttle inputs for manouvres at these speeds. On the other hand though, the 9x4.7 prop is slightly more efficient at lower speeds and gives a slightly higher thrust at take off, so this could be a better choice for a slower flying plane (or even a multicopter, but more on that later).

### Performance plots
![Performance_plot](/img/power.png)
This plot shows the combinations power consumption. Here we see that I was actually pretty lazy for this example case and didn't pick a good motor combination, as a lot of this plot is red with black hatches. This indicates regions where the motor is operating over it's power or current limits, and where there is a risk of burn out. If there's only a little bit of black and red on the plot, you're probably safe, but in this case full throttle will almost certainly burn that motor throughout the entire plotted flight envelope. This plot is also useful however to predict max power draws and currents for your battery selection.
