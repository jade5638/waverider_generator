# 1. Introduction

Waverider Generator is a python package which can be used to generate hypersonic waveriders based on an efficient parametrisation described in [1].
The method makes use of the oscultating cone inverse design method and four design parameters `X1`, `X2`, `X3` and `X4`. 
Please consult the paper referenced for any further clarification on the geometrical meaning of these inputs. 

<p align="center">
  <img src="https://github.com/jade5638/waverider_generator/blob/refactor/images/waverider_cad_example.png?raw=true" alt="Example of Waverider CAD Export" width="500"/>
</p>
<p align="center">Example of Waverider CAD Export</p>

## 1.1 Coordinate System
The coordinate system is defined such that:
- $$x$$ is the streamwise direction
- $$y$$ is the transverse direction
- $$z$$ is the spanwise direction

with the origin at the tip of the waverider.

<p align="center">
  <img src="https://github.com/jade5638/waverider_generator/blob/refactor/images/waverider_axes.png?raw=true" alt="Coordinate System" width="500"/>
</p>
<p align="center">Coordinate System</p>

## 1.2 Terminology

The terminology for the various curves and parameters referenced in this module can be found in more detail in [2]. The image below summarises the main terms :

<p align="center">
  <img src="https://github.com/jade5638/waverider_generator/blob/refactor/images/waverider_convention.svg?raw=true" alt="Waverider Terminology (taken from [2])" width="500"/>
</p>
<p align="center">Waverider Terminology (taken from [2])</p>

Where :
- $h$ is the height of the waverider at the symmetry plane
- $w$ is the half-width of the waverider
- $\beta$ is the shock angle

## 1.3 Flow Field

### 1.3.1 Oblique Shock

In the flat region of the shockwave, the lower surface is determined via the $\theta$- $\beta$ - $M_{\infty}$ equation, which relates the deflection angle $\theta$ to the shock angle $\beta$ in an oblique shock. 

$$
\tan(\theta) = \frac{2 \cot\left(\beta\right) \left(M_{\infty}^2 \sin^2\left(\beta \right) - 1\right)}{M_{\infty}^2 (\gamma + \cos\left(2 \beta\right)) + 2}
$$

Where $\gamma=1.4$ is the ratio of specific heats for air and $M_{\infty}$ is the freestream Mach number.

### 1.3.2 Osculating Cone Theory

In the curved region of the shockwave, the osculating cone theory is used whereby conical flow is locally applied at each osculating plane. 
The Taylor-Maccoll ODE, which describes conical flow, is solved and the resulting streamlines are propagated until the back of the waverider is reached to produce the lower surface.
See [3] for more information on the Osculating Cone Theory.

# 2. Example Script

An example script is provided [here](https://github.com/jade5638/waverider_generator/blob/refactor/examples/test.py). It goes through the following steps :

0. Import required modules
   ```python
   import matplotlib.pyplot as plt 
   from matplotlib.gridspec import GridSpec 
   from waverider_generator import *
   ```
1. Define the waverider's geometric parameters using the GeometricParameters dataclass
   ```python
   geo_params = GeometricParameters(X1 = 0.2, 
                                    X2 = 0.6, 
                                    X3 = 0.2, 
                                    X4 = 0.1, 
                                    w = 4.2, 
                                    h = 1.876)
   ```
   Where :
   - `X1`, `X2`, `X3` and `X4` have the meaning specified in [1]
   - `w` is the half-width of the waverider 
   - `h` is the height of the waverider at the symmetry plane
    
2. Define the waverider's flow parameters using the FlowParameters dataclass
    ```python
    flow_params = FlowParameters(M_design = 5.561, 
                                 betaDeg = 15, 
                                 gamma = 1.4)
   ```
   Where :
   - `M_design` is the Design Mach Number (see [2])
   - `betaDeg` is the Shock Angle in degrees
   - `gamma` is the ratio of specific heats (1.4 for air)
     
3. Instantiate the waverider
   ```python
   waverider = Waverider(geo_params=geo_params, flow_params=flow_params)
   ```
   
4. Build the waverider
   ```python
   streams, curves = waverider.Build(n_planes=50, n_streamline_pts=20, dx_LS_streamline=0.05)
   ```
   The user is advised to look at the description of the `Build` method to understand its I/O structure
   
6. Create and export a CAD Model of the waverider using `streams` and the `to_CAD` method
   ```python
   solid = waverider.to_CAD(streams=streams, 
                           sides = 'both', 
                           export_filename='waverider_example.step')
   ```  
   This will export the model as waverider_example.step in the user's working directory.
   
   The user is advised to look at the description of the `to_CAD` method to understand its I/O structure

7. Plot waverider from different angles
   
   See script for this step. The following figure is then produced :
<p align="center">
  <img src="https://github.com/jade5638/waverider_generator/blob/refactor/images/waverider_plots.svg?raw=true" alt="Waverider Example Plots" width="750"/>
</p>
<p align="center">Waverider Example Plots</p>

# 3. License
MIT License

Copyright (c) 2024 jade5638

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

# 4. Citing

Please use the reference.bib file included to reference this work. 

# 5. References

[1] Jiwon Son, Chankyu Son, and Kwanjung Yee. 
'A Novel Direct Optimization Framework for Hypersonic Waverider Inverse Design Methods'.
In: Aerospace 9.7 (June 2022), p. 348. issn: 2226-4310. doi: 10.3390/aerospace9070348.

[2] Jade Nassif.
'Multi-objective Multi-point Aerodynamic Optimisation of a Hypersonic Waverider'.
Cranfield University, 2024. 
Available: [https://github.com/jade5638/jade_nassif_thesis/blob/main/s422385_Thesis.pdf](https://github.com/jade5638/jade_nassif_thesis/blob/main/s422385_Thesis.pdf?raw=1)

[3] Helmut Sobieczky, F Dougherty, and Kevin Jones.
“Hypersonic Waverider Design from Given Shock Waves”. In: May 1990
