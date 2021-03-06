# DropletCode
Droplet going through a pore.

Code by E.Benet in F.Vernerey group.

## Input paramters
The code takes 9 input parameters in real values, and returns a non-dimensional solution of the pressure energy and position of the droplet inside the pore. The values are defined as follows:
### Pore geometry
The code defines an axisymmetric pore whose profile follows the following equation:
$$r(z) = a\left(1+\frac{m}{b}z\right)\left(\left(1-\left|\frac{z}{b}\right|\right)^n\right)^{1/n}-\frac{d_{real}}{2}$$

| Variable | Description |
| --- | ----------- |
| ```a_real```|Geometry Width|
| ```b_real```|Geometry Height|
| ```n```|Geometry edge sharpness, large values (above 100) might not give a good solution.|
|```p```|Geometry Slope (0 to 1)|
|```d_real```|Pore smallest diameter|

Ideally, keep the value of ```d_real = 1``` and work with normalized units.

### Droplet properties
The physical properties of the droplet are then defined with these varibles

| Variable | Description |
| --- | ----------- |
|```R0_real```|Vesicle Radius. If ```d_real=1```, a value of 1 means that the droplet is twice the pore size|
|```gamma```|Surface Tension (used only in post-processing)|
|```ContactAngle```|A function handle with the angle as a function of depth. For example: ```@(y) pi```|

Additional plotting paramters

| Variable | Description |
| --- | ----------- |
|```plott```|Plot (1) or not(0)|

## Running the code
To run the code, simply open and run [`MAINCODE.m`](https://github.com/EduardBenet/DropletCode/blob/main/MAINCODE.m). The code has a top section tha allows changing the input paramteres above.

Once the code starts running it performs the following operations:
* The code then normalizes all values based on ```d_real``` and ```gamma```
* Solves the equilibrium problem for each position of the droplet
* Returns a normalized value for Pressure, Energy, and Drop position inside the pore

The problem starts (and ends) with a vesicle outside the pore, which means that it might take a while until it founds a plausible solution, and it might give non-physical solutions ones the vesicle gets out.
Since the problem works controlling the bottom edge of the vesicle the second part of the results will be less accurate than the first part.

## Results
The results produced by the code are normalized. You can interpret them as follows:
| Variable | Description |
| --- | ----------- |
|Normalized Presure difference|$CritPress = \frac{\Delta P}{2\gamma}$|
|Normalized Eenergy difference|$E = \frac{\Delta E}{\gamma d_{real}^2}$|
|Center of mass of the drop inside the pore|$Position$|

