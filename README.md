# PropPy
PropPy is an interface and series of scripts to calculate the performance and optimise the selection of RC propeller and motor combinations. It is written in python3 and uses a pyQT interface. Propeller performance is calculated using experimental data from the University of Illinois propeller database and the theoretical results supplied on the APC website. Motor performance is calculated using motor constants. It can choose the optimal motor/prop combination from a list, optimising either for range or flight time, as well as giving plots that help to characterise the performance of each motor/prop combination so combinations can also be compared for their accelerations and max speed.

## Installation:
To install, simply download the codes and run main.py. You will need python3, I use the Anaconda distribution. The code has been developed on windows but should run on linux as well.

If this is your first time with python, I recommend downloading the newest version of Anaconda with python 3.5 from here [here](https://www.continuum.io/downloads). Unfortantely the newest version doesn't include pyqt4, so you will have to install that by hand. To do this, open the conda console in administrator mode and type:
```
conda install -c anaconda pyqt=4.11.4
```

If you don't use Anaconda you will have to install all of the required libraries by hand, which may be annoying.


It requires the following libraries:

pyqt4

scipy

numpy (ensure it is the most recent version)

pickle

matplotlib (some back end functions are used for finding contour intersections which may cause issues in early versions).

## First Use:
The main window shows the power consumption of selected combinations. To plot the power consumption (or efficiency), simply select a motor from the list (new motors can be added if you can't find the one you want), a propeller (currently the only propellers selectable are onces for which the data is available) and an aircraft, then click 'Plot'. Multiple motors, propellers and aircraft can be selected for comparison, the plot is cleared by clicking clear plot in the lower right. The aircraft is defined by its drag polar, which can be calculated from exported XFLR5 results. To use XFLR5 results, a type 2 analysis must be run (fixed lift) and the resulting polar exported to a csv file. This file can then be read by PropPy.

Battery and atmosphere settings are changed in the settings drop down menu.

To find an optimal combination, select a number of propellers and motors and only ONE aircraft, then click 'Optimise selected system'. This opens the optimisation window. Here every combination is analysed for the given options, with the results filtered by the selected constraints. Once you click find optimal combinations, the analysis will begin. This can take some time if it is the first time you run the program, afterwards the results are saved so that future analysis can run faster. The plot showed directly after the analysis is a map showing regions where each combination is the most efficient, with the drag polar of the selected aircraft overlain. You can also plot the efficiency and power consumption of specific combinations, which helps to ensure the combination has the performance you want (eg, high acceleration at high speeds, high launch acceleration etc). Please note that the reported flight times and ranges can be between kind of to very optimistic (in the range of approximately 20%), so take them with a grain of salt. Also many other things can affect the power system performance, such as the interference of the craft's wake for pusher configurations or additional drag for tractor configurations.

The hovering craft and Kv chooser windows are not yet fully functional.

# The main window
The main window shows the lists of motors, propellers and planes that you might want to compare, as well as flight speed vs consumption plots. For this introduction, we'll look at comparing several choices of propeller for a single motor choice, however you can also compare motors or aircraft themselves (for this you will need the xflr5 models however). The image below shows the power consumption for a selection of APC props for a turnigy motor.

![optimisation_plot](/img/example1.png)

What you will quickly notice as you start to do more comparisions, is that the difference in optimal propeller choices generally makes minimal impact around the optimal flight speeds of the craft. Where the difference is bigger however, is in maximum speed and in general performance. The combination performance is better visualised in the optimisation window, however from these plots we can see for example that the smaller, high pitch propeller is the best choice for high speed flight (highest max speed, low power consumption near the drag bucket), however for low speed manouvres the larger, 9x4.7 prop may be a better choice. The most informative plots of propeller/motor performance however are in the optimisation window, which will also give estimates for range and flight times.

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

Here we can see that the combination is most efficient at relatively high speeds due to it's medium pitch. Note that no combination will show high efficiencies at low speeds as the definition of propulsive efficiency depends on a velocity component (see the  [UI prop database definitions](http://m-selig.ae.illinois.edu/props/propDB.html)). This combination seems to be a good choice, as we have a nice region of high efficiency right in the drag bucket of the aircraft (the minimum of the plane's drag polar, and where the most efficient flight happens). Changes of speed in this region don't result in rapid losses of efficiency. We can also see how much acceleration the craft has available. If the plane is flying at 15m/s (with around 2N of drag) and you give it full throttle, with this combination you'll have another 12N of thrust available for acceleration (a bit more than a kilo of force, the actual acceleration depends on the plane's weight but for a 1kg craft this is just over 1m/s/s). At full throttle at 15m/s you'll also be getting about 35% total efficiency out of this system, which isn't too shabby (for sudden acceleration it's rare to hit high efficiencies). The exact throttle response will depend on your ESC and whether it has linear end points or not, but generally you can imagine a vertical line on these plots with 0 thrust at 0 throttle and full thrust at 100% throttle for any flight speed.

As a comparison, this is the same motor with a 9x4.7 prop.
![efficiency_plot](/img/efficiency2.png)

At 15m/s, you only have about 8N of available accelerative force, or about two thirds of what you got with the 9x6. This will feel much less responsive in the air, and will take longer to change speed as well as requiring more throttle inputs for manouvres at these speeds. On the other hand though, the 9x4.7 prop is slightly more efficient at lower speeds and gives a slightly higher thrust at take off, so this could be a better choice for a slower flying plane (or even a multicopter, but more on that later).

### Performance plots
![Performance_plot](/img/power.png)

This plot shows the combinations power consumption. Here we see that I was actually pretty lazy for this example case and didn't pick a good motor combination, as a lot of this plot is red with black hatches. This indicates regions where the motor is operating over its power or current limits, and where there is a risk of burn out. If there's only a little bit of black and red on the plot, you're probably safe, but in this case full throttle will almost certainly burn that motor throughout the entire plotted flight envelope. This plot is also useful however to predict max power draws and currents for your battery selection.

## Adding new motors

Adding new motors is a fairly easy task, with accurate estimations of the motor's performance obtainable with only the motor's Kv, internal resistance, no load current and maximum allowed current.

**Display name:**Motor mass: The mass of the motor in grams. This value is only used if the motor mass is a criterion of the optimisation, so can be left as zero if it is of no interest or unknown.
**Kv:** This value is also referred to as the motor constant, speed constant or back emf constant. It is the inverse of the motor torque constant. Generally manufacturers give this constant in RPM/V, if not then the units need to be converted. For the theoretical motor calculations this is the most important value.
**No load current:** This is the current the motor pulls at a nominal voltage with no load attached. Typically this is given at 10V, or a voltage close to the nominal operating conditions. This is a value that helps to characterise some of the losses of the motor.
**Maximum current:** This is the maximum allowed current of the motor for 60 seconds of continuous operation. Surpassing this value will lead to burning out the motor, and is a serious risk of improperly matched motor-propeller combinations. This value is used mainly to eliminate motor-propeller combinations during optimisation that risk over heating the motor.
**Internal resistance:** This value is the electrical resistance of the windings in Ohms, and characterises the iron losses of the motor. It is also used to calculate the heat generated, as electrical resistance is the main heat source within the motor.
**Motor constant reference temperature:** This is the temperature at which these theoretical values of the motor are taken. The motor constants and resistance are affected by temperature, and thus will change in hotter or colder environments. Note that this value is not the ambient temperature, rather the reference temperature of the motor values so should just be set to 20 or 25 degrees. If the ambient temperature is changed in the atmosphere conditions settings, the motor performance will be slightly affected, even without the thermal model enabled.
**Loss constants:** The motor uses a three constant model to characterise its losses. K3 represents the lowest order losses, k2 is the second order losses and k1 is the highest order losses. Typically, most of the losses of the motor can be characterised with a first order model, and so only k3 is needed. The value of k3 is approximately the no load current, which should be used in the absence of real data.
**Thermal characteristics:** The theoretical model of the motor includes a full three component thermal model, a full explanation of the thermal values can be found in the theory section. In general however, the thermal constants are more difficult to calculate and as the thermal effects are generally quite small if the motor is not operating near its limits, the thermal model is not too significant.
