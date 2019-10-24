# Plot BATS-R-US output in MATLAB

## Example of output:

![The neutral sheet is shown by asking to plot an isosurface of bx=0. One this
surface, we use colors to show the plasma velocity along X (ux). We also show ux at y=0. And, for illustration we show a magnetic field line.](snapshots/sample.png)

## Warning
Working but in development.
Documentation and error handling is minimalist/non-existent.
Only the plot function has documentation and examples.

## What is implemented:

- [x] Conversion from the adaptive grid to uniform grid using interpolation
- [x] All the quantities from the momentum equation (see eq. 7 and 10 from Powell
et al.\ 1999).
- [x] A large number of derived parameters (find exhaustive list below).
- [x] Powerful and user-friendly plot function
- [x] Plot Slice: simple straight cut in the domain
- [x] Plot Quiver: vector field shown using vectors in 3D in a defined rectangular domain
- [x] Plot Contour: level contour of a quantity in the requested plane
- [x] Plot Stream: trace a line along the requested vector field starting at a
certain position
- [x] Plot Surface: Allows *any* cuts in the domain provided by the user. The
user must give a mesh for the surface.
- [x] Plot Isosurface: Find the surface for which a variable has a given value.
The color on the surface may be requested to be that of a field.

## What is *not* implemented:
- [ ] Quiver plot on a surface given by the user
- [ ] Quiver plot on the isosurface
- [ ] Lighting to make the figure sick af!
- [ ] Documentation
- [ ] Good handling of the colorbars positions
- [ ] All for reading inputs that are directly in mesh format (this would extend
the use of this code beyond BATS-R-US).

## Quick start:
This example will reproduce the plot shown on this page.

1. Load the output data:
```matlab
data = bats('file',string_to_cdf_file);
```

2. Convert (interpolate) some variables (var) to a grid of a certain cell-size,
in a certain range (the range is optional)
```matlab
uni = data.toBatsUni(0.125,{'bx','by','bz','ux','uy','uz','jx','jy','jz','rho','p'}, ...
                'xrange',[-40 0],'yrange',[-15 15],'zrange',[-15 15]);
```

Now looking at uni.Output, you will see that it is a mesh rather than a list.

3. Plot using uni.plot(ARGIN) (see 'help uni.plot' for more details)
```matlab
uni.plot('newfigure', 'isosurface','isovariable','bx','level',0,'xrange',[-40 -5],'alpha',0.7,'colorposition','right','variable','ux','color','jet');
uni.plot('slice','variable','ux','yslice',0,'color','jet','colorposition','right','alpha',0.7);
uni.plot('stream','variable','b','start',[-25 -7 0],'color',[1,0,1],'linewidth',2);
```

Other plots are possible by calling uni.plot('plottype') with plottype:
slice/quiver/contour/stream/surface/isosurface

## List-form class: bats
The only thing this class can do is read a cdf file, calculate a bunch of physical quantities (see table below) and create a batsUni object (see below).
No plot or anything else is provided for this class (because I don't know how to
nicely handle the leaves of the adaptive grid so that plots could easily be made).

Basic knowledge about the simulation is stored in:
```matlab
obj.Global
```
It contains the units for the variables, the range, the name of the loaded file.

The output quantities from the simulation are in:
```matlab
obj.Output
```
All the quantities that can be calculated are in:
```matlab
obj.Derived
```
and may be calculated using:
```matlab
obj.calc_*
```
where the * is the quantity to calculate.
Note that some of the calculations depend on each other and the user should call
them in order!

Inputs                  | Calculates                                                                             | Necessary parameters
------------------------|----------------------------------------------------------------------------------------|----------------------
*calc_b*                | magnitude of the magnetic field                                                        | bx,by,bz
*calc_b1*               | magnitude of the magnetic field deviation to the dipole                                | b1x, b1y,b1z
*calc_u*                | magnitude of the bulk velocity                                                         | ux,uy,uz
*calc_rhoU*             | momentum                                                                               | rho,ux,uy,uz
*calc_j*                | current density                                                                        | jx,jy,jz
*calc_jxb*              | $$\mathbf{J}\times\mathbf{B}$$ force                                                   | jx,jy,jz,bx,by,bz
*calc_E*                | Electric field as $$-\mathbf{u}\times\mathbf{B}$$                                      | ux,uy,uz,bx,by,bz
*calc_temp*             | Temperature as: $$T = \frac{P}{n}$$ (in eV)                                            | p, rho
*calc_pb*               | Magnetic pressure                                                                      | b (i.e. must call calc_b first)
*calc_beta*             | Plasma beta                                                                            | p,pb
*calc_alfven*           | Alfven speed                                                                           | b,rho
*calc_vth*              | Thermal speed as: $$v_{th} = \sqrt{2*P/(n*m)}$$                                        | temp
*calc_gyroradius*       | Ion (proton) gyroradius                                                                | u,vth,b
*calc_plasmafreq*       | Ion plasma frequency                                                                   | rho
*calc_inertiallength*   | Ion inertial length (scale below which MHD does not hold)                              | alfven,plasmafreq
*calc_electronVelocity* | Electron velocity as: $$\mathbf{V}_e = \mathbf{u} - (m_p/(m_p+m_e)) * \mathbf{J}/en$$  | rho,jx,jy,jz,ux,uy,uz
*calc_protonVelocity*   | Proton velocity as: $$\mathbf{V}_i = \mathbf{u} + (m_e/(m_p+m_e)) * \mathbf{J}/en$$    | rho,jx,jy,jz,ux,uy,uz

Note that because the data are in list form, we do not calculate any
derivative-based quantities.

### Example:
Say we want to calculate the JxB-force:
1. Load the data
```matlab
data = bats('file',string_to_file,'var',{'x','y','z','jx','jy','jz','bx','by','bz'});
```

2. Calculate it
```matlab
data.calc_jxb;
```

3. See that it is now calculated:
```matlab
data.Derived
```

jxb is no longer empty.

## Mesh-form class: batsUni (subclass of bats)
The second class is called batsUni and is a subclass of bats.
This means that it inherits all the calculation of the derived parameters in the
table above.
This class is made to handle meshed simulations (only BATS-R-US atm) so that plots and derivatives
can easily be evaluated.

Additional calculations are the followings:

Inputs                  | Calculates                                                                                                              | Necessary parameters
------------------------|-------------------------------------------------------------------------------------------------------------------------|----------------------
*calc_gradPb*           | Gradient of magnetic pressure                                                                                           | pb
*calc_gradP*            | Gradient of the plasma pressure                                                                                         | p
*calc_vorticity*        | Vorticity as: $$\mathbf{\omega}= \nabla\times\mathbf{u}$$                                                               | ux,uy,uz
*calc_divBB*            | Divergence of the magnetic energy tensor as: $$\text{Output} = \nabla\cdot\left(\frac{\mathbf{BB}}{\mu_0}}\right)$$     | bx,by,bz
*calc_divRhoUU*         | Divergence of the magnetic energy tensor as: $$\text{Output} = \nabla\cdot\left(\rho\mathbf{uu}\right)$$                | rhoUx,rhoUy,rhoUz,ux,uy,uz

### Example:
1. Load data into a bats class
```matlab
data = bats('file',string_to_file,'var',{'x','y','z','jx','jy','jz','bx','by','bz'});
```

2. Convert the data into a uniform grid
Note that this does interpolation:

```matlab
uni = data.toBatsUni(0.5,{'ux','uy','uz','rho','e','jx','jy','jz','bx','by','bz'},'xrange',[-40 -10],'yrange',[-10 10],'zrange',[-5 5]);
```
Inputs            | Description
------------------|-------------
value (in R units)| COMPULSORY: cell size
{'field',...}     | COMPULSORY: variables to convert to the uniform grid
('xrange',value)  | To limit the domain on which we do the interpolation
('yrange',value)  | To limit the domain on which we do the interpolation
('zrange',value)  | To limit the domain on which we do the interpolation

3. Plot

```matlab
uni.plot('newfigure', 'isosurface','isovariable','bx','level',0,'xrange',[-40 -5],'alpha',0.7,'colorposition','right','variable','ux','color','jet);
```
