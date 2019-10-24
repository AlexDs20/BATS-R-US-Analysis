Working but in development

2 Classes to handle the BATS-R-US part from the SWMF model (single ion species).

One classs is called bats, the only thing it can do is read a cdf file and
calculate a bunch of physical quantities.
The initial quantities are in: obj.Output
All of the quantities that can be calculated are in: obj.Derived
Those quantities may be calculated using: obj.calc_*
where the * is one of the quantity (note that I messed up and the names are not
always just the name of the variable, just auto-complete to find the names).

No plot or anything else is provided for this class (indeed I don't know how
nicely handle the leaves so that plots could easily be made)


The second class is called batsuni and is a subclass of bats.
This class is made to handle meshed simulations so that plots and derivatives
can easily be evaluated.

An object of this class can be obtained from a bats object by calling:
uniformData = bats.toBatsUni(ARGIN);

How to typically use:
% Load the output data:
data = bats('file',string_to_cdf_file);

% Convert (interpolate) some variables (var) to a grid of a certain cell-size,
in a certain range (the range is optional)
uni = data.toBatsUni(0.125,{'bx','by','bz','ux','uy','uz','jx','jy','jz','rho','p'}, ...
                'xrange',[-40 0],'yrange',[-15 15],'zrange',[-15 15]);

% Now looking at uni.Output, you will see that it is a mesh rather than a list.

% You can now do plots using uni.plot(ARGIN) (see 'help uni.plot')
