These are utilities obtained from MathWorks File Exchange which are necessary
to make the plots smoothly.

1)  *colormap and colorbar utilities*.
    From this package, I use the cmapping function. It allows to convert the
    color data of a color plot (e.g. a surface) into "True Color" (Data in RGB)
    This is then used to put all the plots on the same axes so that, even
    though the plots are done at different times, they are on the same axes
    and therefore look nice like in 3D and not superposed on top of each other.

    File Exchange link:
    https://se.mathworks.com/matlabcentral/fileexchange/24371-colormap-and-colorbar-utilities-jul-2014

2)  *Quiver Color and Length management*
    It allows to give the quiver plot colors depending on the magnitude.
    The length function is used to handle the shape of the quiver.
    It allows to make e.g. all vector of the same length, or, limit the
    maximum size of a vector, fix the shape of the head, ...

    File Exchange link:
    https://se.mathworks.com/matlabcentral/fileexchange/71837-quiver-color-and-length-management?s_tid=prof_contriblnk

    github:
    https://github.com/AlexDs20/Quiver

3)  *export_fig*
    Arguably the best user developed function
    The default matlab way of saving figures alters them.
    With *export_fig*, the figure is saved as it looks on the screen.
    I use this in the example that makes the video.

    File Exchange link:
    https://se.mathworks.com/matlabcentral/fileexchange/23629-export_fig

    github:
    https://github.com/altmany/export_fig

4)  *multigradient*
    This allows to create colorbars by giving colors.
    A map of the gradient in color is then created.
    I use this in the plotting examples

    File Exchange:
    https://se.mathworks.com/matlabcentral/fileexchange/68217-multigradient

    github:
    https://github.com/lrkrol/multigradient


