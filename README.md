# proppy
PropPy is an interface and series of scripts to calculate the performance and optimise the selection of RC propeller and motor combinations. It is written in python3 and uses a pyQT interface. Propeller performance is calculated using experimental data from the University of Illinois propeller database and the theoretical results supplied on the APC website. Motor performance is calculated using motor constants.

## Installation:
To install, simply download the codes and run main.py. You will need python3, I use the Anaconda distribution. The code has been developed on windows but should run on linux as well.

## First Use:
The main window shows the power consumption of selected combinations. To plot the power consumption (or efficiency), simply select a motor from the list (new motors can be added if you can't find the one you want), a propeller (currently the only propellers selectable are onces for which the data is available) and an aircraft, then click 'Plot'. Multiple motors, propellers and aircraft can be selected for comparison, the plot is cleared by clicking clear plot in the lower right. The aircraft is defined by its drag polar, which can be calculated from exported XFLR5 results. To use XFLR5 results, a type 2 analysis must be run (fixed lift) and the resulting polar exported to a csv file. This file can then be read by PropPy.

Battery and atmosphere settings are changed in the settings drop down menu.

To find an optimal combination, select a number of propellers and motors and only ONE aircraft, then click 'Optimise selected system'. This opens the optimisation window. Here every combination is analysed for the given options, with the results filtered by the selected constraints. Once you click find optimal combinations, the analysis will begin. This can take some time if it is the first time you run the program, afterwards the results are saved so that future analysis can run faster. The plot showed directly after the analysis is a map showing regions where each combination is the most efficient, with the drag polar of the selected aircraft overlain. You can also plot the efficiency and power consumption of specific combinations, which helps to ensure the combination has the performance you want (eg, high acceleration at high speeds, high launch acceleration etc). Please note that the reported flight times and ranges can be between kind of to very optimistic (in the range of approximately 20%), so take them with a grain of salt.

The hovering craft and Kv chooser windows are not yet fully functional.

## The optimisation window:
Firstly, you will want to check the conditions of optimisation. Are you looking for the motor/prop combination that gives you the longest range or the longest flight time? Do you need a fixed flight speed or should the program select the optimal speed? Do you have any particular constraints on the combinations (the most important constraint is the percentage by which the motor can be overloaded, this should be kept close to 100% to avoid the risk of the motor burning out). Once this has been selected, click the 'find optimal combinations' button.

The first time a combination is calculated may take some time, but the results are saved to make future calculations speedier (the so called 'power surfaces' to be explained in the theory section). The graph that appears after calculation should resemble the one shown below.
![alt tag]( proppy/img/example3.png )

This gives an overview of the regions where each combination is the most efficient. Not all the combinations shown in this plot may meet the constraints however, the successful combinations are shown in the bottom window.
