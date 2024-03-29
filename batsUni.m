classdef batsUni < bats

  properties(Hidden)
    %------------------------------
    %     COMPUTED FIELDS
    Vorticityx, Vorticityy, Vorticityz, Vorticity
    GradPx, GradPy, GradPz, GradP
    GradPbx, GradPby, GradPbz, GradPb
    DivBBx, DivBBy, DivBBz, DivBB
    DivRhoUUx, DivRhoUUy, DivRhoUUz, DivRhoUU
  end
  properties (Hidden, Access = protected)
    % Global
    GlobalCellSize
    GlobalXRange
    GlobalYRange
    GlobalZRange
    GlobalInterpolation = false;

  end

  methods
    %----------------------------------------
    %   Constructor
    %----------------------------------------
    function obj = batsUni(varargin)
      if nargin == 0
        return
      end
      obj = varargin{1}.copyObject(obj);
      if find(strcmpi('xrange',varargin))
        xr = varargin{ find(strcmpi('xrange',varargin))+1 };
        obj.GlobalXRange = [min(xr,[],'all') max(xr,[],'all')];
      else
        obj.GlobalXRange = [min(obj.x,[],'all') max(obj.x,[],'all')];
      end
      if find(strcmpi('yrange',varargin))
        yr = varargin{ find(strcmpi('yrange',varargin))+1 };
        obj.GlobalYRange = [min(yr,[],'all') max(yr,[],'all')];
      else
        obj.GlobalYRange = [min(obj.y,[],'all') max(obj.y,[],'all')];
      end
      if find(strcmpi('zrange',varargin))
        zr = varargin{ find(strcmpi('zrange',varargin))+1 };
        obj.GlobalZRange = [min(zr,[],'all') max(zr,[],'all')];
      else
        obj.GlobalZRange = [min(obj.z,[],'all') max(obj.z,[],'all')];
      end
      if find(strcmpi('cellsize',varargin))
        obj.GlobalCellSize = varargin{ find(strcmpi('cellsize',varargin))+1 };
      else
        d = obj.x(2)-obj.x(1);
        if d == 0 d = obj.y(2)-obj.y(1); end
        if d == 0 d = obj.z(2)-obj.z(1); end
      end
      if find(strcmpi('interpolate',varargin))
        obj.GlobalInterpolation = varargin{ find(strcmpi('interpolate',varargin))+1 };
      end
      if find(strcmpi('variables',varargin))
        var = varargin{ find(strcmpi('variables',varargin))+1 };
      end

      if isa(varargin{1},'bats') & ~isempty(obj.GlobalCellSize) ...
                  & find(strcmpi('variables',varargin)) & obj.GlobalInterpolation

        [x,y,z] = obj.toUniformGrid(varargin{1},obj.GlobalCellSize, var, ...
                      'xrange',obj.GlobalXRange,'yrange',obj.GlobalYRange, ...
                      'zrange',obj.GlobalZRange);
        obj.x = x;
        obj.y = y;
        obj.z = z;
      end
    end
    %----------------------------------------
    %   Simple calculations
    %----------------------------------------
    function calc_j(obj)
      % We can calculate the current from the B field
      % This allows us to not have to interpolate the J field
      % and instead only do it for B and calculate it from curlB/mu0

      % (nT/m) / mu0
      % (nT/m) A^2 / (4*pi*1e-7 * N)
      % (1e-9 * N / (A*m^2) ) A^2 / (4*pi*1e-7 * N)
      % (1e-9 * 1 / (1*m^2) ) A^1 / (4*pi*1e-7 * 1)
      % (1e-9 / m^2 ) A^1 / (4*pi*1e-7)
      % (1e-9 / 4*pi*1e-7) * A/m^2
      % (1e-3 / 4*pi*1e-7) * muA/m^2

      mu0 = 4*pi*1e-7;

      % [vecx,vecy,vecz,vec] = obj.calc_curl('b');

      % obj.jx = 1e-3 * vecx / mu0;
      % obj.jy = 1e-3 * vecy / mu0;
      % obj.jz = 1e-3 * vecz / mu0;
      obj.j = sqrt( obj.jx.^2 + obj.jy.^2 + obj.jz.^2 );

      % obj.GlobalUnits.jx = 'muA/m^2';
      % obj.GlobalUnits.jy = 'muA/m^2';
      % obj.GlobalUnits.jz = 'muA/m^2';
      obj.GlobalUnits.j  = 'muA/m^2';
    end
    %----------------------------------------
    function objOut = getData(obj,position)
      % position: [Nx3]
      objOut = batsUni();
      var = obj.listNonEmptyFields;
      for i = 1 : size(position,1)
        dx = abs(obj.x - position(i,1));
        dy = abs(obj.y - position(i,2));
        dz = abs(obj.z - position(i,3));
        idxx = find( dx == min(dx,[],'all') );     % Each give a plane
        idxy = find( dy == min(dy,[],'all') );
        idxz = find( dz == min(dz,[],'all') );
        idx = intersect(idxx,intersect(idxy,idxz));
        [idx_x,idx_y,idx_z] = ind2sub(size(obj.x),idx);
        for j = 1 : numel(var)
          objOut.(var{j})(i,1) = obj.(var{j})(idx_x(1),idx_y(1),idx_z(1));      % In case there 1 point is exactly in the middle...
        end
      end
    end
    %----------------------------------------
    %   Calculating Derivative based things
    %----------------------------------------
    function obj = calc_gradPb(obj)
      [vecx,vecy,vecz,vec] = obj.calc_grad('Pb');

      obj.GradPbx = 1e6 * vecx;
      obj.GradPby = 1e6 * vecy;
      obj.GradPbz = 1e6 * vecz;
      obj.GradPb  = 1e6 * vec;

      obj.GlobalUnits.GradPbx = 'fN/m^3';
      obj.GlobalUnits.GradPby = 'fN/m^3';
      obj.GlobalUnits.GradPbz = 'fN/m^3';
      obj.GlobalUnits.GradPb  = 'fN/m^3';
    end
    %----------------------------------------
    function obj = calc_gradP(obj)
      [vecx,vecy,vecz,vec] = obj.calc_grad('p');

      obj.GradPx = 1e6 * vecx;
      obj.GradPy = 1e6 * vecy;
      obj.GradPz = 1e6 * vecz;
      obj.GradP  = 1e6 * vec;

      obj.GlobalUnits.GradPx = 'fN/m^3';
      obj.GlobalUnits.GradPy = 'fN/m^3';
      obj.GlobalUnits.GradPz = 'fN/m^3';
      obj.GlobalUnits.GradP  = 'fN/m^3';
    end
    %----------------------------------------
    function obj = calc_vorticity(obj)
      [vecx,vecy,vecz,vec] = obj.calc_curl('u');

      obj.Vorticityx = 1e3 * vecx;
      obj.Vorticityy = 1e3 * vecy;
      obj.Vorticityz = 1e3 * vecz;
      obj.Vorticity  = 1e3 * vec;

      obj.GlobalUnits.Vorticityx = '1/s';
      obj.GlobalUnits.Vorticityy = '1/s';
      obj.GlobalUnits.Vorticityz = '1/s';
      obj.GlobalUnits.Vorticity  = '1/s';
    end
    %----------------------------------------
    function obj = calc_divBB(obj)
      %   calculates the divergence term:
      %         div(BB/mu0)
      %   This is stored in the variable divBBx, divBBy, divBBz, divBB

      [vecx,vecy,vecz,vec] = obj.calc_divTensor('b','b');
      mu0 = 4*pi*1e-7;

      % div(BB/mu0)
      % (1/mu0) * nT * nT /m
      % (1/(4*pi*1e-7)) * 1e-9 * nN/m^3

      obj.DivBBx = vecx * 1e-3 / mu0;
      obj.DivBBy = vecy * 1e-3 / mu0;
      obj.DivBBz = vecz * 1e-3 / mu0;
      obj.DivBB  = vec  * 1e-3 / mu0;

      obj.GlobalUnits.DivBBx = 'fN/m^3';
      obj.GlobalUnits.DivBBy = 'fN/m^3';
      obj.GlobalUnits.DivBBz = 'fN/m^3';
      obj.GlobalUnits.DivBB  = 'fN/m^3';
    end
    %----------------------------------------
    function obj = calc_divRhoUU(obj)
      [vecx,vecy,vecz,vec] = obj.calc_divTensor('rhoU','u');

      %    rhoU        u       div
      % amu mum-2 s-1  km s-1    m-1

      mp = 1.67262*1e3;

      obj.DivRhoUUx = mp * vecx ;
      obj.DivRhoUUy = mp * vecy ;
      obj.DivRhoUUz = mp * vecz ;
      obj.DivRhoUU  = mp * vec  ;

      obj.GlobalUnits.DivRhoUUx = 'fN/m^3';
      obj.GlobalUnits.DivRhoUUy = 'fN/m^3';
      obj.GlobalUnits.DivRhoUUz = 'fN/m^3';
      obj.GlobalUnits.DivRhoUU  = 'fN/m^3';
    end
    %--------------------------------------------------
    %     Calc all
    function obj = calc_all(obj)
      calc_all@bats(obj);
      obj.calc_gradPb;
      obj.calc_gradP;
      obj.calc_vorticity;
      obj.calc_divBB;
      obj.calc_divRhoUU;
    end
    %----------------------------------------
    %     Reduce Domain
    %----------------------------------------
    function obj = reduceDomain(obj,varargin)
      % KWARGS: 'xrange', 'yrange' and 'zrange', otherwise
      if find(strcmpi('xrange',varargin))
        xrange = varargin{ find(strcmpi('xrange',varargin))+1 };
      else
        xrange = [min(obj.x,[],'all') max(obj.x,[],'all')];
      end
      if find(strcmpi('yrange',varargin))
        yrange = varargin{ find(strcmpi('yrange',varargin))+1 };
      else
        yrange = [min(obj.y,[],'all') max(obj.y,[],'all')];
      end
      if find(strcmpi('zrange',varargin))
        zrange = varargin{ find(strcmpi('zrange',varargin))+1 };
      else
        zrange = [min(obj.z,[],'all') max(obj.z,[],'all')];
      end

      ix = find(obj.x >= xrange(1) & obj.x <= xrange(end));
      iy = find(obj.y >= yrange(1) & obj.y <= yrange(end));
      iz = find(obj.z >= zrange(1) & obj.z <= zrange(end));
      indices = intersect(ix,intersect(iy,iz));

      sx = numel(intersect(obj.x(ix),obj.x(ix)));
      sy = numel(intersect(obj.y(iy),obj.y(iy)));
      sz = numel(intersect(obj.z(iz),obj.z(iz)));

      % List non empty fields
      var = obj.listNonEmptyFields;

      % Reduce elements in it.
      for i = 1 : numel(var)
        obj.(var{i}) = reshape(obj.(var{i})(indices),[sx, sy, sz]);
      end
    end

    %----------------------------------------
    %     Plotting shit
    %----------------------------------------
    function h = plot(obj,varargin)
      % DESCRPITION:
      % -----------
      %     Create flexible plot for BATS-R-US simulation output once converted into a uniform grid.
      %     Additional calls of the plot function allow to add information on an already made plot.
      %     Possibility to make silces, quiver, contour, streamlines, surface and isosurface plots in the
      %     3D domain.
      %     Each of the plot type have a number of possible input parameters which are usually given as
      %     ('keyword', value) pair.
      %
      % INPUTS:
      % ------
      %   We first list all the possible inputs and what they do.
      %   Then, for each type of plot, we list the inputs that work.
      %
      %   All inputs:
      %   ----------
      %     - 'newfigure'                               Creates a new figure on which the plot is made
      %     - 'sizex',val                               Size of the new figure in units of the screen
      %     - 'sizey',val                               Size of the new figure in units of the screen
      %     - 'slice'                                   plot slices of the domain
      %     - 'quiver'                                  plot vector arrows in the domain
      %     - 'contour'                                 plot contour level in a slice of the domain
      %     - 'stream'                                  plot a streamline for a vector field
      %     - 'surface'                                 plot any surface cutting the domain (user input a mesh)
      %     - 'isosurface'                              plot surface level of a certain quantity (e.g. bx = 0 to identify magnetotail neutral sheet)
      %     - 'variable', field                         the variable we want the plot/color of
      %     - 'colorvariable', field                    variable for which we want to give the color to the isosurface or the streamline
      %     - 'xlim', [min, max]                        Range of the x axis shown on the plot
      %     - 'ylim', [min, max]                        Range of the y axis shown on the plot
      %     - 'zlim', [min, max]                        Range of the z axis shown on the plot
      %     - 'xrange', [min max]                       Range of data to consider for the plot (this allows to plot reduced domain)
      %     - 'yrange', [min max]                       Range of data to consider for the plot (this allows to plot reduced domain)
      %     - 'zrange', [min max]                       Range of data to consider for the plot (this allows to plot reduced domain)
      %     - 'log'                                     If the colorbar should be in log scale
      %     - 'position', values                        Plot position in the figure
      %     - 'color',value                             color of the shown data. Either [R G B] for a single color plot (stream,contour,isosurface). Or a colormap (slice,quiver,stream,surface,isosurface)
      %     - 'alpha',value                             transparency: between [0,1]
      %%TBD - 'colorposition',value                     Position of the colorbar/label of the ploted quantity: either the location or the position (see colorbar doc.)
      %     - 'colorrange', [min max]                   range to use for the colorbar
      %     - 'xslice', value(s)                        value at which we want to do the cut(s)
      %     - 'yslice', value(s)                        value at which we want to do the cut(s)
      %     - 'zslice', value(s)                        value at which we want to do the cut(s)
      %     - 'linewidth',value                         width of the line to show
      %     - 'level', value                            value at which we want equal value of the variable (contour,isosurface)
      %     - 'start', [x,y,z]                          starting position for the streams
      %     - 'increment',value                         (integer) Increment used for quiver plot to not plot all the vectors
      %     - 'HeadAngle',value                         Angle of the head of quiver plot [deg]
      %     - 'HeadLength',value                        Length of the Head, (in units of the XYZ axis i.e. R)
      %     - 'length',value                            Length of the quiver Tail. Only possible value atm: 'equal'
      %     - 'MagUnitLength',value                     Magnitude of the field (quiver plot) corresponding to 1R on the plot, hence, a vector with the given magnitude would have length 1 (in xyz) on the plot
      %%TBD - 'fancylook'                               Make the plot sick as fuck!
      %
      %
      %   Which inputs are allowed?
      %   ------------------------
      %     General simple rules:
      %     --------------------
      %       - The first plot MUST include 'newfigure'
      %       - A type of plot MUST always be given: 'slice', 'quiver', 'contour', 'stream', 'surface', 'isosurface'
      %       - The ('variable',value) pair MUST always be given.
      %
      %     Slice:
      %     -----
      %       MUST:
      %         - 'slice'
      %         - 'variable',value
      %         - 'xslice',values or 'yslice',values or 'zslice',values
      %             several values may be given and will create several cuts showing the variable values in the cut
      %
      %       Optional:
      %         - 'color',value : can be a colormap string or [Nx3]
      %         - 'log'
      %         - 'alpha',value
      %         - 'colorrange',values
      %         - 'xrange', 'yrange', 'zrange': to reduce the domain in which the slice is done
      %         - 'colorposition'
      %
      %     Quiver:
      %     ------
      %       MUST:
      %         - 'quiver'
      %         - 'variable',value : variable must be a vector field quantity e.g. 'u', 'j', 'E', 'b', ...
      %                               it should not include the x,y,z at the end!
      %
      %       Optional:
      %         - 'color',value : colormap string or [Nx3] used to color the vector with their magnitudes
      %         - 'colorrange',values
      %         - 'linewidth',value
      %         - 'xrange','yrange','zrange' to only show the vectors in this sub-domain defined by these inputs
      %         - 'increment',value     : integer number
      %         - 'HeadAngle',value     : Angle of the Head of the arrows, given in degrees
      %         - 'HeadLength',value    : Length of the Head in units of the XYZ axis
      %         - 'RotHead', value      : Angle [deg] with which the head of the quiver will be rotated with the tail of vector as rotation axis
      %         - 'Length',value        : only current possible input: 'equal'
      %                                     if the 'equal' value is passed, all vectors will have the same length in the xyz space
      %         - 'MagUnitLength',value : to manipulate the length of the quiver tail.
      %                                   A vector with the given 'MagUnitLength' value will have length 1 on the plot.
      %           If both ('Length','equal') and ('MagUnitLength',value) are given, all vectors will have the length 'MagUnitLength' (in units of the xyz). Note the big difference of the value to give 'MagUnitLength' when combine to ('Length','equal')
      %         - 'colorposition'
      %
      %     Contour:
      %     ------
      %       MUST:
      %         - 'contour'
      %         - 'variable',value
      %         - 'level',value     : value at which we want to show the variable (we show: variable=level)
      %%TBD     - 'xslice', 'yslice' or 'zslice': only one value for the moment
      %
      %       Optional:
      %         - 'color',value   : color of the line shown [R G B]
      %         - 'linewidth',value
      %         - 'alpha', value
      %         - 'xrange','yrange','zrange' with values [min max] to limit the domain in which to show the plot
      %%TBD     - 'colorposition'
      %
      %     Stream:
      %     ------
      %       MUST:
      %         - 'stream'
      %         - 'variable',value      : of a vector field e.g. 'u', 'j', ...
      %         - 'start',values        : [x,y,z] position at which we start the streamline
      %
      %       Optional:
      %         - 'colorvariable'
      %         - 'colorrange'        (if 'colorvariable')
      %         - 'color' : [R G B]   (or a colormap when using 'colorvariable')
      %         - 'log'               (if colorvariable)
      %         - 'linewidth',value
      %         - 'xrange','yrange','zrange'  : to only show part of the domain
      %         - 'colorposition'
      %
      %     Surface:
      %     -------
      %       MUST:
      %         - 'surface',X,Y,Z       : where X, Y, Z is the mesh of the surface (see example below)
      %         - 'variable', value
      %
      %       Optional:
      %         - 'color',value         : colormap name or [Nx3]
      %         - 'log'
      %         - 'alpha',value
      %         - 'colorrange', [min max]
      %%TBD     - 'colorposition',value
      %
      %     Isosurface:
      %     ----------
      %       MUST:
      %         - 'isosurface'
      %         - 'variable', value       : field for which we want the surface
      %         - 'level', value          : value at which we want to show the fields' surface
      %
      %       Optional:
      %         - 'colorvariable', value       : Field for which we show the color on the isosurface
      %         - 'color',value
      %         - 'log'
      %         - 'alpha',value
      %         - 'colorrange',[min max]
      %%TBD     - 'colorposition',value
      %
      %     Trajectory:
      %     ----------
      %       MUST:
      %         - To use this plot, you MUST first trace particles using the "traceParticles" routine
      %         - 'trajectory'
      %         - 'variable', [i]       : e.g. [i], for ploting the ith particles: Trajectories([i]).
      %                                   ([i] could be [1,3,4] to trace 3 particles)
      %         - ''
      %
      % EXAMPLES:
      % --------
      %
      %   Slice:
      %     uni.plot('newfigure','slice','variable','ux','zslice',0,'color','jet','colorposition','eastoutside');
      %     uni.plot('slice','variable','ux','yslice',0,'color','parula', ...
      %               'colorposition','westoutside','alpha',0.8,'colorrange',[-200 800]);
      %
      %   Quiver:
      %     uni.plot('quiver','variable','u','zrange',0,'color','parula','colorrange',[0 800]);
      %
      %   Contour:
      %     uni.plot('contour','variable','bz','zslice',0,'color',[1 0 0]);
      %
      %   Stream:
      %     uni.plot('stream','variable','b','start',[-25 0 0],'color',[0 1 0], ...
      %             'colorvariable','b','linewidth',2,'xrange',[-30, -10]);
      %
      %   Surface:
      %     % Get the surface you are interested in
      %     [Y,Z] = meshgrid([-10:0.15:10],[-10:0.15:10]);
      %     myX = @(Y,Z)(- sqrt((Y.^2 + Z.^2)) - 10 );
      %     X = myX(Y,Z);
      %     uni.plot('surface',X,Y,Z,'variable','ux','color','cool', 'colorrange',[-400 400],'xlim',[-60 20]);
      %
      %   Isosurface:
      %     uni.plot('isosurface','variable','bx','level',0, ...
      %               'xrange',[-30 -10],'yrange',[-10 10],'alpha',0.7, ...
      %               'colorposition','eastoutside','variable','ux','color','jet');
      %
      %   Trajectory:
      %
      %     Add lighting stuff, so that it looks sick af!
      %

      %----------------------------------------
      %     TREAT INPUT
      [plotType, variable, islog, ...                    % GENERAL
        sizex, sizey, ...
        xl, yl, zl, position, colorposition, ...
        colorName, alpha, colorrange, ...          % Several Types
        xslice, yslice, zslice, ...                % Slice and Contour
        ix, iy, iz, increment, ...                  % Quiver
        level, LineWidth, ...               % Contour
        HeadAngle, HeadLength, RotHead ...
        QuiverLength, MagUnitLength, ...
        start, ColorVariable, ...                                 % Stream
        X, Y, Z ...
      ] ...
        = obj.setInputValues(varargin(:));

      % Check if we want a new figure
      if find(strcmpi('newfigure',varargin))
        h = figure('Units','Normalized','OuterPosition',[0 0 sizex sizey],'Color',[1 1 1]);
        visible = 'on';
        ax = axes(h);
      else
        h = gcf;
        visible = 'off';
        ax_prev = findall(h,'type','axes');
        % Create a new axes, copy of the previous one but I remove the plots that are on it.
        ax = copyobj(ax_prev(end),h);
        delete(get(ax,'Children'));
      end

      if plotType == 1  % Slice
        [x,y,z] = meshgrid(unique(obj.x(ix,iy,iz)),unique(obj.y(ix,iy,iz)),unique(obj.z(ix,iy,iz)));
        field = permute(obj.(variable)(ix,iy,iz),[2 1 3]);

        hp = slice( ax, x,y,z, field, single(xslice),single(yslice),single(zslice) );
        set(hp, ...
            'EdgeColor', 'none', ...
            'FaceAlpha', alpha ...
            );
        colormap(ax,colorName);
        cb = colorbar(ax);
        %%TBD if find(strcmpi('log',varargin))
        %%TBD   % set(hp,'CData',log10(hp.CData));
        %%TBD   set(ax,'colorscale','log');
        %%TBD   %cb.Ruler.Scale = 'log';
        %%TBD   %cb.Ruler.MinorTick = 'on';
        %%TBD end
        cl = [variable,' [',obj.GlobalUnits.(variable),']'];
      elseif plotType == 2  % Quiver
        x = obj.x([ix(1):increment:ix(end)],[iy(1):increment:iy(end)],[iz(1):increment:iz(end)]);
        y = obj.y([ix(1):increment:ix(end)],[iy(1):increment:iy(end)],[iz(1):increment:iz(end)]);
        z = obj.z([ix(1):increment:ix(end)],[iy(1):increment:iy(end)],[iz(1):increment:iz(end)]);
        fieldx = obj.([variable,'x'])([ix(1):increment:ix(end)],[iy(1):increment:iy(end)],[iz(1):increment:iz(end)]);
        fieldy = obj.([variable,'y'])([ix(1):increment:ix(end)],[iy(1):increment:iy(end)],[iz(1):increment:iz(end)]);
        fieldz = obj.([variable,'z'])([ix(1):increment:ix(end)],[iy(1):increment:iy(end)],[iz(1):increment:iz(end)]);

        hp = quiver3(ax, x,y,z, fieldx,fieldy,fieldz );

        mag = sqrt(fieldx.^2+fieldy.^2+fieldz.^2);
        % Colors
        SetQuiverColor(hp,colormap(ax,colorName), ...
                        'range', colorrange, ...
                        'mags',  mag );
        set(hp,'LineWidth',LineWidth);
        cb = colorbar(ax);
        cl = [variable,' [',obj.GlobalUnits.([variable,'x']),']'];

        % Quiver Length
        if isempty(MagUnitLength)
          MagUnitLength = max(mag,[],'all');
          if strcmpi(QuiverLength,'equal')
            Length = 0.9*increment*obj.GlobalCellSize*(mag./mag);
          else
            Length = 1.2*increment*obj.GlobalCellSize*(mag./MagUnitLength);
          end
        else
          if strcmpi(QuiverLength,'equal')
            Length = MagUnitLength * ones(size(mag));
          else
            Length = mag./MagUnitLength;
          end
        end
        drawnow;
        SetQuiverLength(hp,Length,'HeadLength',HeadLength,'HeadAngle',HeadAngle,'RotHead',RotHead);

        % SetLength?
      elseif plotType == 3  % Contour
        [x,y,z] = meshgrid(unique(obj.x(ix,iy,iz)),unique(obj.y(ix,iy,iz)),unique(obj.z(ix,iy,iz)));
        field = permute(obj.(variable)(ix,iy,iz),[2 1 3]);
        hp = contourslice(ax, x,y,z, field, single(xslice),single(yslice),single(zslice), ...
                          [level level]);
        colormap(ax,colorName)
        set(hp(:),'LineWidth',LineWidth,'EdgeAlpha',alpha);

        cb = [];
        cl = [];
      elseif plotType == 4  % Stream
        [x,y,z] = meshgrid( double(unique(obj.x(ix,iy,iz))), ...
                  double(unique(obj.y(ix,iy,iz))),double(unique(obj.z(ix,iy,iz))));
        fieldx = double(permute(obj.([variable,'x'])(ix,iy,iz),[2 1 3]));
        fieldy = double(permute(obj.([variable,'y'])(ix,iy,iz),[2 1 3]));
        fieldz = double(permute(obj.([variable,'z'])(ix,iy,iz),[2 1 3]));

        hp = streamline(ax, stream3(x,y,z, fieldx,fieldy, fieldz, ...
                                start(:,1), start(:,2), start(:,3)) );
        hp2 = streamline(ax, stream3(x,y,z, -fieldx,-fieldy,-fieldz, ...
                                start(:,1), start(:,2), start(:,3)) );
        hp = [hp;hp2]; clear hp2;

        if ~isempty(ColorVariable)
          % Replace the streamplot by a surface plot. Interpolat the field along the field line
          F = griddedInterpolant(obj.x(ix,iy,iz),...
                      obj.y(ix,iy,iz),obj.z(ix,iy,iz), ...
                      double(obj.(ColorVariable)(ix,iy,iz)));
          for j = numel(hp) : -1 : 1
            xsl = hp(j).XData;
            ysl = hp(j).YData;
            zsl = hp(j).ZData;
            fieldinterp = F(xsl,ysl,zsl);
            hp(j) = surface(ax,[xsl;xsl],[ysl;ysl],[zsl;zsl],[fieldinterp;fieldinterp], ...
                            'FaceColor','none','EdgeColor','interp','LineWidth',LineWidth);
          end
          colormap(ax,colorName);
          cb = colorbar(ax);
          cl = [ColorVariable,' [',obj.GlobalUnits.(ColorVariable),']'];
          delete(findall(ax,'type','Line'));
        else
          set(hp,'color',colorName(1,:),'LineWidth',LineWidth);
          cb = [];
          cl = [];
        end
      elseif plotType == 5  % Surface
        F = griddedInterpolant(obj.x,obj.y,obj.z,obj.(variable),'nearest');
        field = F(X,Y,Z);

        hp = surf(ax, X,Y,Z,field);
        set(hp,'EdgeColor','none','FaceAlpha',alpha);

        colormap(ax,colorName);
        cb = colorbar(ax);
        cl = [variable,' [',obj.GlobalUnits.(variable),']'];
      elseif plotType == 6  % Isosurface
        [x,y,z] = meshgrid( double(unique(obj.x(ix,iy,iz))), ...
                  double(unique(obj.y(ix,iy,iz))),double(unique(obj.z(ix,iy,iz))));
        IsoField = double(permute(obj.(variable)(ix,iy,iz),[2 1 3]));

        if ~isempty(ColorVariable)
          colorField = double(permute(obj.(ColorVariable)(ix,iy,iz),[2 1 3]));
          [faces,verts,colors] = isosurface(x,y,z,IsoField,level,colorField);
          hp = patch(ax, 'Vertices',verts,'Faces',faces,'FaceVertexCData',colors,...
                        'FaceColor','flat', 'EdgeColor','flat', 'EdgeColor','none','FaceAlpha',alpha);
          isonormals(x,y,z,IsoField,hp);
          colormap(ax, colorName);
          cb = colorbar(ax);
          cl = [ColorVariable,' [',obj.GlobalUnits.(ColorVariable),']'];
        else
          hp = patch(ax, isosurface(x,y,z,IsoField,level));
          isonormals(x,y,z,IsoField,hp);
          set(hp, 'FaceColor',colorName, 'EdgeColor','none','FaceAlpha',alpha);
          cb = []; cl = [];
        end
        %camlight
        %lighting gouraud
      elseif plotType == 7  % Trajectory
        if ~isempty(ColorVariable) & ~strcmpi(ColorVariable,'t')
          F = griddedInterpolant(obj.x(ix,iy,iz),...
                      obj.y(ix,iy,iz),obj.z(ix,iy,iz), ...
                      double(obj.(ColorVariable)(ix,iy,iz)));
        end
        for j = 1 : numel(variable)
          xd = obj.Trajectories{variable(j)}.pos(:,1)';
          yd = obj.Trajectories{variable(j)}.pos(:,2)';
          zd = obj.Trajectories{variable(j)}.pos(:,3)';
          if isempty(ColorVariable) | ~strcmpi(ColorVariable,'t')
            hp(j) = plot3(ax,xd,yd,zd,...
                        'color',colorName(1,:),'LineWidth',LineWidth);
            cb = [];
            cl = [];
          else
            if strcmpi(ColorVariable,'t')
              fieldinterp = obj.Trajectories{variable(j)}.t';
              cl = [ColorVariable,' [s]'];
            else
              fieldinterp = F(xd,yd,zd);
              cl = [ColorVariable,' [',obj.GlobalUnits.(ColorVariable),']'];
              if j == 1, colormap(ax,colorName); end
            end
            hp(j) = surface(ax,[xd;xd],[yd;yd],[zd;zd],[fieldinterp;fieldinterp], ...
                                'FaceColor','none','EdgeColor','interp','LineWidth',LineWidth);
            cb = colorbar(ax);
          end
        end
      else
        cb = []; cl = [];
      end

      obj.setProperties(ax,cb,xl,yl,zl,position,visible, alpha, ...
              islog,cl,colorposition,colorrange);

    % Link axes, put plots on same axes (needs true colors), delete axes (if no colorbar), ...
      if isempty(find(strcmpi('newfigure',varargin)))
        % Linkaxes:
        Link = linkprop([ax; ax_prev],{'CameraUpVector', 'CameraPosition', 'CameraTarget', 'XLim', 'YLim', 'ZLim'});
        setappdata(h, 'StoreTheLink', Link);
        if plotType == 1
          for j = 1 : numel(hp)
            if islog
              RGB = cmapping(log10(hp(j).CData),colorName,log10(colorrange));
            else
              RGB = cmapping(hp(j).CData,colorName,colorrange);
            end
            hpCopy = copyobj(hp(j),ax_prev(end));
            set(hpCopy,'CData',RGB);
            delete(hp(j));
          end
        elseif plotType == 2
          hpCopy = copyobj(hp,ax_prev(end));
          drawnow; pause(0.1);
          set(hpCopy.Head,'VertexData',get(hp.Head,'VertexData'));    % Somehow this is not conserved when copying
          set(hpCopy.Tail,'VertexData',get(hp.Tail,'VertexData'));
          delete(hp)
        elseif plotType == 3
          for j = 1 : numel(hp)
            C = ones([size(hp(j).CData),3]);
            C(:,:,1) = colorName(1).*C(:,:,1);
            C(:,:,2) = colorName(2).*C(:,:,2);
            C(:,:,3) = colorName(3).*C(:,:,3);
            C(3,:,:) = NaN;
            set(hp(j),'CData',C);
            hpCopy = copyobj(hp(j),ax_prev(end));
          end
          delete(ax);
        elseif plotType == 4
          for j = 1 : numel(hp)
            if ~isempty(ColorVariable)
              if islog
                RGB = cmapping(log10(hp(j).CData),colorName,log10(colorrange));
              else
                RGB = cmapping(hp(j).CData,colorName,colorrange);
              end
              hpCopy = copyobj(hp(j),ax_prev(end));
              set(hpCopy,'CData',RGB);
              delete(hp(j));
            else
              hpCopy = copyobj(hp(j),ax_prev(end));
              delete(hp(j));
            end
          end
          if isempty(ColorVariable)
            delete(ax);
          end
        elseif plotType == 5
          for j = 1 : numel(hp)
            if islog
              RGB = cmapping(log10(hp(j).CData),colorName,log10(colorrange));
            else
              RGB = cmapping(hp(j).CData,colorName,colorrange);
            end
            hpCopy = copyobj(hp(j),ax_prev(end));
            set(hpCopy,'CData',RGB);
            delete(hp(j));
          end
        elseif plotType == 6
          if ~isempty(ColorVariable)
            if islog
              RGB = cmapping(log10(hp.FaceVertexCData),colorName,log10(colorrange));
            else
              RGB = cmapping(hp.FaceVertexCData,colorName,colorrange);
            end
            RGB = cmapping(hp.FaceVertexCData,colorName,colorrange);
            hpCopy = copyobj(hp,ax_prev(end));
            set(hpCopy,'FaceVertexCData',RGB);
            delete(hp);
          else
            hpCopy = copyobj(hp,ax_prev(end));
            delete(ax);
          end
        elseif plotType == 7
          for j = 1 : numel(hp)
            if ~isempty(ColorVariable)
              if islog
                RGB = cmapping(log10(hp(j).CData),colorName,log10(colorrange));
              else
                RGB = cmapping(hp(j).CData,colorName,colorrange);
              end
              hpCopy = copyobj(hp(j),ax_prev(end));
              set(hpCopy,'CData',RGB);
              delete(hp(j));
            else
              hpCopy = copyobj(hp(j),ax_prev(end));
              delete(hp(j));
            end
          end
          if isempty(ColorVariable)
            delete(ax);
          end
        end
      else
        % Change to real colors
        if plotType == 1 || plotType == 5 || (plotType == 4 || plotType == 7 & ~isempty(ColorVariable))
          for j = 1 : numel(hp)
            if islog
              RGB = cmapping(log10(hp(j).CData),colorName,log10(colorrange));
            else
              RGB = cmapping(hp(j).CData,colorName,colorrange);
            end
            set(hp(j),'CData',RGB);
          end
        elseif plotType == 3
          C = ones([size(hp(j).CData),3]);
          C(:,:,1) = colorName(1).*C(:,:,1);
          C(:,:,2) = colorName(2).*C(:,:,2);
          C(:,:,3) = colorName(3).*C(:,:,3);
          C(3,:,:) = NaN;
          set(hp(j),'CData',C);
        elseif plotType == 6
          if ~isempty(variable)
            if islog
              RGB = cmapping(log10(hp.FaceVertexCData),colorName,log10(colorrange));
            else
              RGB = cmapping(hp.FaceVertexCData,colorName,colorrange);
            end
            set(hp,'FaceVertexCData',RGB);
          end
        end
        %axis vis3d
        %h = light;
        %for az = -50:10:50
        %  lightangle(h,az,30)
        %  pause(.1)
        %end
      end
    end
    %----------------------------------------
    %     Plotting spacecraft
    function h = plotSC(obj,positions,var,varargin)
      % h = plotSC(positions,var,varargin)
      %
      % Description:
      %       Plot as if a spacecraft was moving through the simulation
      %
      % Output: h : handle to the figure
      %
      % Input:  positions: [Nx3] array of positions
      %         var: cell arrays of variables to plot, each entry of the cell array is one panel
      %              if an entry of the cell array is a cell array, the plot those variables on the same subplot.
      %         varargin:
      %             'xlinc',value         xlable increment used (so that we do not show too much text)
      %
      % Example:
      %       %Create position array:
      %       pos(:,1) = [-18.00:0.125:-15.00];
      %       pos(:,2) = 0.000 + zeros(size(pos,1),1);
      %       pos(:,3) = 0.000 + zeros(size(pos,1),1);
      %
      %       % Create the cell array of variables to plot (5 subplots).
      %       %          subplot1           subplot2           sp3            sp4   sp5
      %       var = { {'ux','uy','uz'},  {'bx','by','bz'}, {'jx','jy','jz'}, 'rho', 'p' };
      %
      %       obj.plotSC(pos,var);

      data = obj.getData(positions);
      % Check for xlabels increments
      if find(strcmpi('xlinc',varargin))
        xlinc = varargin{ find(strcmpi('xlinc',varargin))+1 };
      else
        xlinc = 1;
      end

      h = figure;
      set(h,'Units','Normalized','OuterPosition',[0 0 1 1],'color','w');
      color = 'kbrcmg';
      L = numel(var);
      M = 1;
      ym = 0.10; yM = 0.98;
      ddy = 0.001;
      dy = (yM-ym-(L-1)*ddy)/L;
      xm = 0.1; xM = 0.9;

      % get the distance between the points so that one can plot them in the given order
      distance = zeros(size(positions,1)-1,1);
      for i = size(positions,1)-1 : -1 : 1
        distance(i) = sqrt( sum((positions(i+1,:)-positions(i,:)).^2) );
      end
      distance = [0; cumsum(distance)];

      Rlabels = sqrt(sum(data.x.^2 + data.y.^2 + data.z.^2,2));
      Xlabels = data.x;
      Ylabels = data.y;
      Zlabels = data.z;

      for i = 1 : L
        VAR = var{i};
        ax(i) = axes;
        axpos = [xm yM-i*dy-(i-1)*ddy xM-xm dy];
        set(ax(i),'Units','Normalized','Position',axpos);
        hold(ax(i),'on');
        if iscell(VAR)
          yl = [];
          for j = 1 : numel(VAR)
            plot(ax(i),distance,data.(VAR{j}),color(j));
            yl = [yl,VAR{j},newline];
          end
          legend(ax(i),VAR);
          ylabel([yl, '[',obj.GlobalUnits.(VAR{1}),']']);
        else
          plot(ax(i),distance,data.(VAR),color(1));
          ylabel([char(VAR), ' [',obj.GlobalUnits.(VAR),']']);
        end
        grid(ax(i),'on');
        box(ax(i),'on');
        axis(ax(i),'tight');
        xticks(ax(i),distance(1:xlinc:end));
        line(ax(i),ax(i).XLim,[0 0],'color',[0.3 0.3 0.3],'LineStyle','--','LineWidth',0.5,'HandleVisibility','off');
        if i ~= L
          set(ax(i),'XTickLabel',[]);
        elseif i == L
          linkaxes(ax(:),'x');
          labels = compose('% 3.2f\\newline% 3.2f\\newline% 3.2f\\newline% 3.2f', ...
                          [Rlabels(1:xlinc:end),Xlabels(1:xlinc:end),Ylabels(1:xlinc:end),Zlabels(1:xlinc:end)]);
          xticklabels( labels );
          yt = get(ax(i),'YTick');
          a = annotation(h,'textbox', [axpos(1)/2 axpos(2) 0 0], 'String', 'R [R]\newlineX [R]\newlineY [R]\newlineZ [R]', 'FitBoxToText', true,'EdgeColor','none');
        end
      end
    end

    %--------------------------------------------------
    function output = copyObject(input, output)
       C = metaclass(input);
       P = C.Properties;
       for k = 1:length(P)
         if ~P{k}.Dependent
           output.(P{k}.Name) = input.(P{k}.Name);
         end
       end
    end
  end
  %--------------------------------------------------

  methods (Static)
    function h = plotEarth(varargin)
      %
      % Tool to plot an Earth
      %
      % Inputs: 'newfigure'

      if find(strcmpi('newfigure',varargin))
        h = figure('Units','Normalized','OuterPosition',[0 0 1 1],'Color',[1 1 1]);
        ax = axes(h);
      else
        h = gcf;
        ax_prev = findall(h,'type','axes');
        ax = copyobj(ax_prev(end),h);
        delete(get(ax,'Children'));
      end

      [y,x,z] = sphere(36);
      c = zeros([size(x),3]);
      IDay1 = find(x>=0&y>=0);
      IDay2 = find(x>0&y<=0);
      IDay = [IDay1(:);IDay2(:)];
      [ix,iy] = ind2sub(size(x),IDay);
      c(IDay) = 1;
      c(numel(y)+IDay) = 1;

      hp = surface(ax,x,y,z,c,'EdgeColor','none');

      if isempty(find(strcmpi('newfigure',varargin)))
        hpCopy = copyobj(hp,ax_prev(end));
        delete(ax);
      end
    end
  end

  %--------------------------------------------------
  methods (Hidden, Access = private)
    %--------------------------------------------------
    %       GENERAL CALCULATION FUNCTIONS
    %       FOR DERIVATIVE
    %--------------------------------------------------
    function [vecx, vecy, vecz, vec] = calc_divTensor(obj,field1,field2)
      %   calculates divergence of a tensor:
      %             \partial_i (field1_i field2_j)
      %   by computing: \partial_i(field1_i) field2_j + field1_i \partial_i (field2_j)
      %   that is:
      %                     div(field1) field2_j        +     field1 dot grad(field2_j)
      %
      %   Units: [field1][field2]/m

      if ndims(obj.x)~=3
        error('The fields must be in xyz in order to calculate the divergence...');
      end
      field1x = [field1,'x']; field1y = [field1,'y']; field1z = [field1,'z'];
      field2x = [field2,'x']; field2y = [field2,'y']; field2z = [field2,'z'];
      Re = 6371.2e3;

      % This needs to be done because matlab permute xyz into yxz for some reasons...
      x = permute(obj.x,[2 1 3]).*Re;
      y = permute(obj.y,[2 1 3]).*Re;
      z = permute(obj.z,[2 1 3]).*Re;
      var1x= permute(obj.(field1x),[2 1 3]);
      var1y= permute(obj.(field1y),[2 1 3]);
      var1z= permute(obj.(field1z),[2 1 3]);
      var2x= permute(obj.(field2x),[2 1 3]);
      var2y= permute(obj.(field2y),[2 1 3]);
      var2z= permute(obj.(field2z),[2 1 3]);


      div = divergence(x,y,z,var1x,var1y,var1z);

      [gradXx,gradXy,gradXz] = gradient(var2x,obj.GlobalCellSize*Re);
      [gradYx,gradYy,gradYz] = gradient(var2y,obj.GlobalCellSize*Re);
      [gradZx,gradZy,gradZz] = gradient(var2z,obj.GlobalCellSize*Re);

      vecx = div .* var2x ...
            + var1x.*gradXx + var1y.*gradXy + var1z.*gradXz;

      vecy = div .* var2y ...
            + var1x.*gradYx + var1y.*gradYy + var1z.*gradYz;

      vecz = div .* var2z ...
            + var1x.*gradZx + var1y.*gradZy + var1z.*gradZz;

      vec = sqrt(vecx.^2 + vecy.^2 + vecz.^2);

      vecx = permute(vecx,[2 1 3]);
      vecy = permute(vecy,[2 1 3]);
      vecz = permute(vecz,[2 1 3]);
      vec  = permute(vec ,[2 1 3]);
    end
    %--------------------------------------------------
    function [vecx, vecy, vecz, vec] = calc_grad(obj,var)
      % Output Units: [inputs]/m
      Re = 6371.2e3;
      [vecx,vecy,vecz] = gradient(obj.(var),obj.GlobalCellSize*Re);

      vec  = sqrt(vecx.^2 + vecy.^2 + vecz.^2);
    end
    %--------------------------------------------------
    function [curlx, curly, curlz, curlMag] = calc_curl(obj,var)
      % Output: Units: [inputs/m];
      varx = [var,'x'];
      vary = [var,'y'];
      varz = [var,'z'];

      % Need to change because the mesh is created using ndgrid
      % while curl wants the grid to be made with meshgrid
      % (which is dumb and swaps dim 1 and 2)
      if ndims(obj.x) == 2    % Means we did a cut
        if all(obj.x == obj.x(1))
          dim1 = 1; dim2 = 2;
        elseif all(obj.y == obj.y(1))
          dim1 = 1; dim2 = 2;
        elseif all(obj.z == obj.z(1))
          dim1 = 2; dim2 = 1;
        end
        x = permute(obj.x,[dim1 dim2]);
        y = permute(obj.y,[dim1 dim2]);
        z = permute(obj.z,[dim1 dim2]);
        varx= permute(obj.(varx),[dim1 dim2]);
        vary= permute(obj.(vary),[dim1 dim2]);
        varz= permute(obj.(varz),[dim1 dim2]);
      elseif ndims(obj.x) == 3
        x = permute(obj.x,[2 1 3]);
        y = permute(obj.y,[2 1 3]);
        z = permute(obj.z,[2 1 3]);
        varx= permute(obj.(varx),[2 1 3]);
        vary= permute(obj.(vary),[2 1 3]);
        varz= permute(obj.(varz),[2 1 3]);
      end

      [curlx,curly,curlz,~] = curl(x,y,z,varx,vary,varz);

      if ndims(obj.x) == 2
        curlx = permute(curlx,[dim1 dim2]);
        curly = permute(curly,[dim1 dim2]);
        curlz = permute(curlz,[dim1 dim2]);
      elseif ndims(obj.x) == 3
        curlx = permute(curlx,[2 1 3]);
        curly = permute(curly,[2 1 3]);
        curlz = permute(curlz,[2 1 3]);
      end
      Re = 6371.2e3;

      curlx = curlx/Re;
      curly = curly/Re;
      curlz = curlz/Re;
      curlMag  = sqrt( curlx.^2 + curly.^2 + curlz.^2 );
    end
    %--------------------------------------------------
    %     Plot stuff
    %--------------------------------------------------
    function [plotType, variable, islog, ...                   % GENERAL
              sizex, sizey, ...
              xl, yl, zl, position, colorposition, ...
              colorName, alpha, colorrange, ...         % Several Types
              xslice, yslice, zslice, ...               % Slice
              ix, iy, iz, increment, ...                % Quiver
              level, LineWidth, ...                     % Contour
              HeadAngle, HeadLength, RotHead, ...
              QuiverLength, MagUnitLength, ...
              start, ColorVariable, ...                                % Stream
              X, Y, Z ...                               % Surface
              ] ...
            = setInputValues(obj,var)
      %----------------------------------------
      %           GENERAL
        if find(strcmpi('slice',var))
          plotType = 1;
        elseif find(strcmpi('quiver',var))
          plotType = 2;
        elseif find(strcmpi('contour',var))
          plotType = 3;
        elseif find(strcmpi('stream',var))
          plotType = 4;
        elseif find(strcmpi('surface',var))
          plotType = 5;
        elseif find(strcmpi('isosurface',var))
          plotType = 6;
        elseif find(strcmpi('trajectory',var))
          plotType = 7;
        end

        if find(strcmpi('variable',var))
          variable = var{ find(strcmpi('variable',var))+1 };
        else
          variable = [];
        end

        ca = findall(gcf,'type','axes');
        if find(strcmpi('xlim',var))
          xl = var{ find(strcmpi('xlim',var))+1 };
        else
          if isempty(find(strcmpi('newfigure',var)))
            xl = ca.XLim;
          else
            xl = obj.GlobalXRange;
          end
        end
        if find(strcmpi('ylim',var))
          yl = var{ find(strcmpi('ylim',var))+1 };
        else
          if isempty(find(strcmpi('newfigure',var)))
            yl = ca.YLim;
          else
            yl = obj.GlobalYRange;
          end
        end
        if find(strcmpi('zlim',var))
          zl = var{ find(strcmpi('zlim',var))+1 };
        else
          if isempty(find(strcmpi('newfigure',var)))
            zl = ca.ZLim;
          else
            zl = obj.GlobalZRange;
          end
        end

        if find(strcmpi('position',var))
          position = var{ find(strcmpi('position',var))+1 };
        else
          position = [0.126 0.11 0.7513 0.815];
        end

        if find(strcmpi('log',var))
          islog = 1;
        else
          islog = 0;
        end
        % Size of Figure in units of screen size
        if find(strcmpi('sizex',var))
          sizex =  var{ find(strcmpi('sizex',var))+1 };
        else
          sizex = 1;
        end
        if find(strcmpi('sizey',var))
          sizey =  var{ find(strcmpi('sizey',var))+1 };
        else
          sizey = 1;
        end

      %----------------------------------------
      %   May be common to different plots
        if find(strcmpi('alpha',var))
          alpha = var{ find(strcmpi('alpha',var))+1 };
        else
          alpha = 1;
        end
      %----------------------------------------
      %       Slice
        if find(strcmpi('xslice',var)), xslice = var{ find(strcmpi('xslice',var))+1 };
        else, xslice = [];
        end
        if find(strcmpi('yslice',var)), yslice = var{ find(strcmpi('yslice',var))+1 };
        else, yslice = [];
        end
        if find(strcmpi('zslice',var)), zslice = var{ find(strcmpi('zslice',var))+1 };
        else, zslice = [];
        end
      %----------------------------------------
      %      Domain, linear and dimensional indices
        if find(strcmpi('xrange',var)), xrange = var{ find(strcmpi('xrange',var))+1 };
        else, xrange = [min(obj.x,[],'all') max(obj.x,[],'all')];
        end
        if find(strcmpi('yrange',var)), yrange = var{ find(strcmpi('yrange',var))+1 };
        else, yrange = [min(obj.y,[],'all') max(obj.y,[],'all')];
        end
        if find(strcmpi('zrange',var)), zrange = var{ find(strcmpi('zrange',var))+1 };
        else, zrange = [min(obj.z,[],'all') max(obj.z,[],'all')];
        end
        ix = find( obj.x(:,1,1)>=xrange(1) & obj.x(:,1,1)<=xrange(end) );
        iy = find( obj.y(1,:,1)>=yrange(1) & obj.y(1,:,1)<=yrange(end) );
        iz = find( obj.z(1,1,:)>=zrange(1) & obj.z(1,1,:)<=zrange(end) );
      %----------------------------------------
      %       Color
        if find(strcmpi('color',var))
          colorName = var{ find(strcmpi('color',var))+1 };
        else
          if (plotType == 1 | plotType == 2 | plotType == 5 | plotType == 6 | ...
             ( (plotType == 4 | plotType == 7) & ~isempty(find(strcmpi('colorvariable',var)))) )
            colorName = parula;
          elseif plotType == 3
            colorName = [1 0 1];
          elseif plotType == 4 | plotType == 7
            colorName = [0.66, 0.66, 0.66];
          end
        end

        if find(strcmpi('colorvariable',var))
          ColorVariable = var{ find(strcmpi('colorvariable',var))+1 };
        else
          ColorVariable = [];
        end

        if find(strcmpi('colorrange',var))
          colorrange = var{ find(strcmpi('colorrange',var))+1 };
        else
          if plotType == 1
            colorrange = [min(obj.(variable)(ix,iy,iz),[],'all') max(obj.(variable)(ix,iy,iz),[],'all')];
          elseif plotType == 2
            mag = sqrt( obj.([variable,'x']).^2 + obj.([variable,'y']).^2 + obj.([variable,'z']).^2 );
            colorrange = [min(mag(ix,iy,iz),[],'all') max(mag(ix,iy,iz),[],'all')];
          elseif plotType == 4 | plotType == 6
            if find(strcmpi('colorvariable',var))
              colorrange = [min(obj.(variable)(ix,iy,iz),[],'all') max(obj.(variable)(ix,iy,iz),[],'all')];
            else
              colorrange = [0 1];
            end
          elseif plotType == 7
            if find(strcmpi('colorvariable',var)) & ~strcmpi(ColorVariable,'t')
              colorrange = [min(obj.(ColorVariable)(ix,iy,iz),[],'all') max(obj.(ColorVariable)(ix,iy,iz),[],'all')];
            elseif find(strcmpi('colorvariable',var)) & strcmpi(ColorVariable,'t')
              colorrange = [ min(cellfun(@(x) x.t(1),obj.Trajectories,'UniformOutput', true)),...
                              max(cellfun(@(x) x.t(end),obj.Trajectories,'UniformOutput', true)) ];
            else
              colorrange = [0 1];
            end
          else
            colorrange = [0 1];  % This does not matter as for streams and contours, I only have 1 color in the colormap
          end
        end

        if find(strcmpi('colorposition',var))
          colorposition = var{ find(strcmpi('colorposition',var))+1 };
        else
          colorposition = [0.92 0.1105 0.0112 0.8143];
          %colorposition = [0.07 0.1105 0.0112 0.8143];
        end

      %----------------------------------------
      %      Line Prop
        if find(strcmpi('linewidth',var)), LineWidth = var{ find(strcmpi('linewidth',var))+1 };
        else, LineWidth = 0.5;
        end
      %----------------------------------------
      %      Levels (contour plots and isosurface)
        if find(strcmpi('level',var)), level = var{ find(strcmpi('level',var))+1 };
        else, level = [];
        end
      %----------------------------------------
      %      Stream
        if find(strcmpi('start',var)), start = var{ find(strcmpi('start',var))+1 };
        else, start = [-20 0 0];
        end
      %----------------------------------------
      %      Grid
        if find(strcmpi('surface',var))
          X = var{ find(strcmpi('surface',var))+1 };
          Y = var{ find(strcmpi('surface',var))+2 };
          Z = var{ find(strcmpi('surface',var))+3 };
        else
          X = []; Y = []; Z = [];
        end
      %----------------------------------------
      %      Quiver
        if find(strcmpi('increment',var))
          increment = var{ find(strcmpi('increment',var))+1 };
        else
          increment = 1;
        end
        if find(strcmpi('HeadAngle',var))
          HeadAngle = var{ find(strcmpi('HeadAngle',var))+1 };
        else
          HeadAngle = 30;
        end
        if find(strcmpi('HeadLength',var))
          HeadLength = var{ find(strcmpi('HeadLength',var))+1 };
        else
          HeadLength = (1/3)*obj.GlobalCellSize;
        end
        if find(strcmpi('RotHead',var))
          RotHead = var{ find(strcmpi('RotHead',var))+1 };
        else
          RotHead = 0;
        end
        if find(strcmpi('length',var))
          QuiverLength = var{ find(strcmpi('length',var))+1 };
        else
          QuiverLength = [];
        end
        if find(strcmpi('MagUnitLength',var))
          MagUnitLength = var{ find(strcmpi('MagUnitLength',var))+1 };
        else
          MagUnitLength = [];
        end
    end
    %--------------------------------------------------
    function setProperties(obj,ax,cb, ...
              xl, yl, zl, position, ...
              visible, alpha, ...
              islog, cl,colorposition,colorrange)
      axis(ax,'equal');
      set(ax, ...
          'XLim', xl, 'YLim', yl, 'ZLim', zl, ...
          'Units','Normalized', ...
          'Position', position, ...
          'Visible', visible ...
          );

      % Labels
      xlabel(ax, ['X [', obj.GlobalUnits.x,']'] );
      ylabel(ax, ['Y [', obj.GlobalUnits.y,']'] );
      zlabel(ax, ['Z [', obj.GlobalUnits.z,']'] );

      if islog
        set(ax,'colorscale','log');
      else
        set(ax,'colorscale','lin');
      end

      % Colorbar stuff
      if ~isempty(cb)
        if ischar(colorposition)
          set(cb,'Location',colorposition);
        else
          set(cb,'Position',colorposition);
        end
        ylabel(cb,cl);
        caxis(ax,colorrange)
        % Colorbar transparency: (found at: https://stackoverflow.com/questions/37423603/setting-alpha-of-colorbar-in-matlab-r2015b )
        % Somehow it does not work when run here but it works if use keyboard (however it does not stick! fuck this shit)
        if 0
          drawnow
          cdata = cb.Face.Texture.CData;
          cdata(end,:) = uint8(alpha * cdata(end,:));
          cb.Face.Texture.CData = cdata;
          cb.Face.Texture.ColorType = 'truecoloralpha';
          cb.Face.ColorBinding = 'discrete';
        end
      end
      %%POSITION
    end
    %--------------------------------------------------
    function indices = domainIndex(obj,varargin)
      % KWARGS: 'xrange', 'yrange' and 'zrange', otherwise
      if find(strcmpi('xrange',varargin))
        xrange = varargin{ find(strcmpi('xrange',varargin))+1 };
      else
        xrange = [min(obj.x,[],'all') max(obj.x,[],'all')];
      end
      if find(strcmpi('yrange',varargin))
        yrange = varargin{ find(strcmpi('yrange',varargin))+1 };
      else
        yrange = [min(obj.y,[],'all') max(obj.y,[],'all')];
      end
      if find(strcmpi('zrange',varargin))
        zrange = varargin{ find(strcmpi('zrange',varargin))+1 };
      else
        zrange = [min(obj.z,[],'all') max(obj.z,[],'all')];
      end

      ix = find(obj.x >= xrange(1) & obj.x <= xrange(end));
      iy = find(obj.y >= yrange(1) & obj.y <= yrange(end));
      iz = find(obj.z >= zrange(1) & obj.z <= zrange(end));
      indices = intersect(ix,intersect(iy,iz));
    end
  end
  %--------------------------------------------------

  methods (Hidden)
  %--------------------------------------------------
  %     Getters
    function Global = getGlobal(obj)
      Global = getGlobal@bats(obj);
      Global.CellSize = obj.GlobalCellSize;
      Global.XRange = obj.GlobalXRange;
      Global.YRange = obj.GlobalYRange;
      Global.ZRange = obj.GlobalZRange;
      Global.Interpolation = obj.GlobalInterpolation;
    end
  %--------------------------------------------------
    function Derived = getDerived(obj)
      Derived = getDerived@bats(obj);

      Derived.Vorticityx = obj.Vorticityx;
      Derived.Vorticityy = obj.Vorticityy;
      Derived.Vorticityz = obj.Vorticityz;
      Derived.Vorticity  = obj.Vorticity;

      Derived.GradPx = obj.GradPx;
      Derived.GradPy = obj.GradPy;
      Derived.GradPz = obj.GradPz;
      Derived.GradP  = obj.GradP;

      Derived.GradPbx = obj.GradPbx;
      Derived.GradPby = obj.GradPby;
      Derived.GradPbz = obj.GradPbz;
      Derived.GradPb  = obj.GradPb;

      Derived.DivBBx = obj.DivBBx;
      Derived.DivBBy = obj.DivBBy;
      Derived.DivBBz = obj.DivBBz;
      Derived.DivBB  = obj.DivBB;

      Derived.DivRhoUUx = obj.DivRhoUUx;
      Derived.DivRhoUUy = obj.DivRhoUUy;
      Derived.DivRhoUUz = obj.DivRhoUUz;
      Derived.DivRhoUU  = obj.DivRhoUU;

    end
  %--------------------------------------------------
  end
end
