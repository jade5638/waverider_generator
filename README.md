# Waverider Generator

Waverider Generator is a python package which can be used to generate hypersonic waveriders based on an efficient parametrisation described in [1].
The method makes use of the oscultating cone inverse design method and four design parameters `X1`, `X2`, `X3` and `X4`.
## Required Inputs
- Design parameters `X1`, `X2`, `X3` and `X4`. Note this is entered as a list of four elements where the parameters are organised in the order listed here. Refer to the examples.
- Freestream Mach number `M_inf`.
- Shock angle in degrees `beta`.
- Height of the waverider at the base plane `height`.
- Width of the waverider `width`. Note that this refers to half of the total width of the waverider due to the symmetry.
## Optional Inputs
- Number of osculating planes used in the generation of the geometry `n_planes`. Note that this doesn't include the symmetry plane and the tip. A minimum number of planes is set at 10 to preserve quality and this is the default value.
- Number of points in the streamwise direction for the generation of the upper surface as well as the flat part of the shockwave on the lower surface `n_streamwise`. A minimum is set at 10 to preserve quality and this is the default value.
- Maximum step to be used in the tracing of the streamlines for the lower surface `delta_streamwise`. Note this is entered as a percentage of the length of the waverider and a maximum is set at 20% to preserve quality. The default value is 5%.

## Summary of inputs
|Input|Type|Conditions|
|:-------------:|:--------------:|:--------------:|
| `X1`, `X2`, `X3` and `X4` | `list` | $0\leq X1,X2,X3,X4 \leq 1$ <br> $\frac{X2}{\left(1-X1\right)^4}\leq\ \frac{7}{64}\left(\frac{width}{height}\right)^4$<br>(see [1] for more information on this design condition) |
| `M_inf` | `float`, `int` | `M_inf>0` |
|`beta`| `float`, `int`| `0<beta<90`|
|`height`| `float`, `int`| `height>0`|
|`width`| `float`, `int`| `width>0`|
|`n_planes`|`int`|`n_planes>=10`|
|`n_streamise`|`int`|`n_streamwise>=10`|
|`delta_streamise`|`float`|`0<delta_streamwise<=0.2`|

# Usage and functionality

A user can create a `waverider` instance by importing the `waverider` class from `waverider_generator.generator`. Everything takes place inside the class constructor so the user only needs to initialise the instance to create the geometry. An example is shown below: <br>

`from waverider_generator.generator import waverider as wr` <br>
```python
M_inf=5
beta=15
height=1.34
width=3
dp=[0.11,0.63,0,0.46]
n_planes=20
n_streamwise=10
delta_streamwise=0.1
waverider=wr(M_inf=M_inf,beta=beta,height=height,width=width,dp=dp,n_upper_surface=10000,n_shockwave=10000,n_planes=n_planes,n_streamwise=n_streamwise,delta_streamise=delta_streamwise)



-----------
# References

[1] Jiwon Son, Chankyu Son, and Kwanjung Yee. 
'A Novel Direct Optimization Framework for Hypersonic Waverider Inverse Design Methods'.
In: Aerospace 9.7 (June 2022), p. 348. issn: 2226-4310. doi: 10.3390/aerospace9070348.
