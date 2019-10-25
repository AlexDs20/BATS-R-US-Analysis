classdef batsUni < bats
  properties (Hidden, Access = protected)
    % Global
    GlobalCellSize
    GlobalXRange
    GlobalYRange
    GlobalZRange
    GlobalInterpolation = false;

    %------------------------------
    %     COMPUTED FIELDS
    Vorticityx, Vorticityy, Vorticityz, Vorticity
    GradPx, GradPy, GradPz, GradP
    GradPbx, GradPby, GradPbz, GradPb
    DivBBx, DivBBy, DivBBz, DivBB
    DivRhoUUx, DivRhoUUy, DivRhoUUz, DivRhoUU
  end

  methods
    %----------------------------------------
    %   Constructor
    %----------------------------------------
    function obj = batsUni(varargin)
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
    function obj = getData(obj,position)
      dx = obj.x - position(1);
      dy = obj.y - position(2);
      dz = obj.z - position(3);
      idxx = find( dx == min(dx) );
      idxy = find( dy == min(dy) );
      idxz = find( dz == min(dz) );

      idx = intersect(idxx,intersect(idxy,idxz));

      [idx_x,idx_y,idx_z] = ind2sub(size(obj.x),idx);

      var = obj.listNonEmptyFields;

      for i = 1 : numel(var)
        obj.(var{i}) = obj.(var{i})(idx_x,idx_y,idx_z);
      end
    end
    %----------------------------------------
    function obj = clearField(obj,var)
      for i = 1 : numel(var)
        obj.(var{i}) = [];
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
      %     - 'slice'                                   plot slices of the domain
      %     - 'quiver'                                  plot vector arrows in the domain
      %     - 'contour'                                 plot contour level in a slice of the domain
      %     - 'stream'                                  plot a streamline for a vector field
      %     - 'surface'                                 plot any surface cutting the domain (user input a mesh)
      %     - 'isosurface'                              plot surface level of a certain quantity (e.g. bx = 0 to identify magnetotail neutral sheet)
      %     - 'variable', field                         the variable we want the plot/color of
      %     - 'isovariable', field                      variable for which we want to find the isosurface
      %     - 'xlim', [min, max]                        Range of the x axis shown on the plot
      %     - 'ylim', [min, max]                        Range of the y axis shown on the plot
      %     - 'zlim', [min, max]                        Range of the z axis shown on the plot
      %     - 'xrange', [min max]                       Range of data to consider for the plot (this allows to plot reduced domain)
      %     - 'yrange', [min max]                       Range of data to consider for the plot (this allows to plot reduced domain)
      %     - 'zrange', [min max]                       Range of data to consider for the plot (this allows to plot reduced domain)
      %     - 'position', values                        Plot position in the figure
      %     - 'color',value                             color of the shown data. Either [R G B] for a single color plot (stream,contour,isosurface). Or a colormap (slice,quiver,surface,isosurface)
      %     - 'alpha',value                             transparency: between [0,1]
      %%TBD - 'colorposition',value                     Position of the colorbar/label of the ploted quantity
      %     - 'colorrange', [min max]                   range to use for the colorbar
      %     - 'xslice', value(s)                        value at which we want to do the cut(s)
      %     - 'yslice', value(s)                        value at which we want to do the cut(s)
      %     - 'zslice', value(s)                        value at which we want to do the cut(s)
      %     - 'linewidth',value                         width of the line to show
      %     - 'level', value                            value at which we want equal value of the variable (contour,isosurface)
      %     - 'start', [x,y,z]                          starting position for the stream
      %     - 'forward'                                 only follow the stream along the field
      %     - 'backward'                                only follow the stream in the direction opposite to the field (if none of 'forward','backward' are given, stream in both directions)
      %%TBD - 'interp'
      %%TBD - 'fancylook'                               Make the plot sick as fuck!
      %%TBD - 'increment',value                         Increment used for quiver plot to not plot all the vectors
      %
      %
      %   Which inputs are allowed?
      %   ------------------------
      %     General simple rules:
      %     --------------------
      %
      %       - The first plot MUST include 'newfigure'
      %       - A type of plot MUST always be given: 'slice', 'quiver', 'contour', 'stream', 'surface', 'isosurface'
      %       - The ('variable',value) pair MUST always be given except for isosurface where it is optional.
      %         For isosurface, the 'isovariable' is a MUST.
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
      %         - 'contour'
      %         - 'variable',value      : of a vector field e.g. 'u', 'j', ...
      %         - 'start',values        : [x,y,z] position at which we start the streamline
      %
      %       Optional:
      %         - 'color',value : [R G B]
      %         - 'linewidth',value
      %         - 'forward'             : to only show along the field
      %         - 'backward'            : to only show anti-parallel to the field
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
      %         - 'alpha',value
      %         - 'colorrange', [min max]
      %%TBD     - 'colorposition',value
      %
      %     Isosurface:
      %     ----------
      %       MUST:
      %         - 'isosurface'
      %         - 'isovariable', value    : field for which we want the surface
      %         - 'level', value          : value at which we want to show the fields' surface
      %
      %       Optional:
      %         - 'variable', value       : Field for which we show the color on the isosurface
      %         - 'color',value
      %         - 'alpha',value
      %         - 'colorrange',[min max]
      %%TBD     - 'colorposition',value
      %
      % EXAMPLES:
      % --------
      %
      %   Slice:
      %     uni.plot('newfigure','slice','variable','ux','zslice',0,'color','jet','colorposition','right');
      %     uni.plot('slice','variable','ux','yslice',0,'color','parula', ...
      %               'colorposition','left','alpha',0.8,'colorrange',[-200 800]);
      %
      %   Quiver:
      %     uni.plot('quiver','variable','u','zrange',0,'color','parula','colorrange',[0 800]);
      %
      %   Contour:
      %     uni.plot('contour','variable','bz','zslice',0,'color',[1 0 0]);
      %
      %   Stream:
      %     uni.plot('stream','variable','b','start',[-25 0 0],'color',[0 1 0],'linewidth',2,'forward', ...
      %               'xrange',[-30, -10]);
      %
      %   Surface:
      %     % Get the surface you are interested in
      %     [Y,Z] = meshgrid([-10:0.15:10],[-10:0.15:10]);
      %     myX = @(Y,Z)(- sqrt((5*Y.^2 + 10*Z.^2)) - 10 );
      %     X = myX(Y,Z);
      %     uni.plot('surface',X,Y,Z,'variable','ux','color','cool', 'colorrange',[-400 400],'xlim',[-60 20]);
      %
      %   Isosurface:
      %     uni.plot('isosurface','isovariable','bx','level',0, ...
      %               'xrange',[-30 -10],'yrange',[-10 10],'alpha',0.7, ...
      %               'colorposition','right','variable','ux','color','jet');
      %
      %
      %     Add lighting stuff, so that it looks sick af!
      %
      %

      % Check if we want a new figure
      if find(strcmpi('newfigure',varargin))
        h = figure('Units','Normalized','OuterPosition',[0 0 1 1],'Color',[1 1 1]);
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

      %----------------------------------------
      %     TREAT INPUT
      [plotType, variable, ...                    % GENERAL
       xl, yl, zl, position, colorposition, ...
       colorName, alpha, colorrange, ...          % Several Types
       xslice, yslice, zslice, ...                % Slice and Contour
       ix, iy, iz, increment, ...                  % Quiver
       level, LineWidth, ...               % Contour
       start, ...                                 % Stream
       X, Y, Z ...
       ] ...
        = obj.setInputValues(varargin(:));


      if plotType == 1
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
      elseif plotType == 2
        x = obj.x([ix(1):increment:ix(end)],[iy(1):increment:iy(end)],[iz(1):increment:iz(end)]);
        y = obj.y([ix(1):increment:ix(end)],[iy(1):increment:iy(end)],[iz(1):increment:iz(end)]);
        z = obj.z([ix(1):increment:ix(end)],[iy(1):increment:iy(end)],[iz(1):increment:iz(end)]);
        fieldx = obj.([variable,'x'])([ix(1):increment:ix(end)],[iy(1):increment:iy(end)],[iz(1):increment:iz(end)]);
        fieldy = obj.([variable,'y'])([ix(1):increment:ix(end)],[iy(1):increment:iy(end)],[iz(1):increment:iz(end)]);
        fieldz = obj.([variable,'z'])([ix(1):increment:ix(end)],[iy(1):increment:iy(end)],[iz(1):increment:iz(end)]);

        hp = quiver3(ax, x,y,z, fieldx,fieldy,fieldz );

        % Colors
        SetQuiverColor(hp,colormap(ax,colorName), ...
                        'range', colorrange, ...
                        'mags',  sqrt(fieldx.^2+fieldy.^2+fieldz.^2) );
        set(hp,'LineWidth',LineWidth);
        cb = colorbar(ax);
        cl = [variable,' [',obj.GlobalUnits.([variable,'x']),']'];

        % SetLength?
      elseif plotType == 3
        [x,y,z] = meshgrid(unique(obj.x(ix,iy,iz)),unique(obj.y(ix,iy,iz)),unique(obj.z(ix,iy,iz)));
        field = permute(obj.(variable)(ix,iy,iz),[2 1 3]);
        hp = contourslice(ax, x,y,z, field, single(xslice),single(yslice),single(zslice), ...
                          [level level]);
        colormap(ax,colorName)
        set(hp(:),'LineWidth',LineWidth,'EdgeAlpha',alpha);

        cb = [];
        cl = [];
      elseif plotType == 4
        [x,y,z] = meshgrid( double(unique(obj.x(ix,iy,iz))), ...
                  double(unique(obj.y(ix,iy,iz))),double(unique(obj.z(ix,iy,iz))));
        fieldx = double(permute(obj.([variable,'x'])(ix,iy,iz),[2 1 3]));
        fieldy = double(permute(obj.([variable,'y'])(ix,iy,iz),[2 1 3]));
        fieldz = double(permute(obj.([variable,'z'])(ix,iy,iz),[2 1 3]));

        if isempty(find(strcmpi('backward',varargin)))
          hp = streamline(ax, stream3(x,y,z, fieldx,fieldy, fieldz, ...
                                  start(:,1), start(:,2), start(:,3)) );
          set(hp,'color',colorName,'LineWidth',LineWidth);
        end
        if isempty(find(strcmpi('forward',varargin)))
          hp2 = streamline(ax, stream3(x,y,z, -fieldx,-fieldy,-fieldz, ...
                                  start(:,1), start(:,2), start(:,3)) );
          set(hp2,'color',colorName,'LineWidth',LineWidth);
        end

        cb = [];
        cl = [];
      elseif plotType == 5
        F = griddedInterpolant(obj.x,obj.y,obj.z,obj.(variable),'nearest');
        field = F(X,Y,Z);

        hp = surf(ax, X,Y,Z,field);
        set(hp,'EdgeColor','none','FaceAlpha',alpha);

        colormap(ax,colorName);
        cb = colorbar(ax);
        cl = [variable,' [',obj.GlobalUnits.(variable),']'];
      elseif plotType == 6
        IsoVariable = varargin{ find(strcmpi('isovariable',varargin))+1 };

        [x,y,z] = meshgrid( double(unique(obj.x(ix,iy,iz))), ...
                  double(unique(obj.y(ix,iy,iz))),double(unique(obj.z(ix,iy,iz))));
        IsoField = double(permute(obj.(IsoVariable)(ix,iy,iz),[2 1 3]));
        if ~isempty(variable)
          colorField = double(permute(obj.(variable)(ix,iy,iz),[2 1 3]));
          [faces,verts,colors] = isosurface(x,y,z,IsoField,level,colorField);
          hp = patch(ax, 'Vertices',verts,'Faces',faces,'FaceVertexCData',colors,...
                        'FaceColor','flat', 'EdgeColor','flat', 'EdgeColor','none','FaceAlpha',alpha);
                      %'FaceColor','interp','EdgeColor','interp');
          isonormals(x,y,z,IsoField,hp);
          colormap(ax, colorName);
          cb = colorbar(ax);
          cl = [variable,' [',obj.GlobalUnits.(variable),']'];
        else
          hp = patch(ax, isosurface(x,y,z,IsoField,level));
          isonormals(x,y,z,IsoField,hp);
          set(hp, 'FaceColor',colorName, 'EdgeColor','none','FaceAlpha',alpha);

          cb = []; cl = [];
        end
        %camlight
        %lighting gouraud

      else
        cb = []; cl = [];
      end

      obj.setProperties(ax,cb,xl,yl,zl,position,visible, ...
              cl,colorposition,colorrange);

    % Link axes, put plots on same axes (needs true colors), delete axes (if no colorbar), ...
      if isempty(find(strcmpi('newfigure',varargin)))
        % Linkaxes:
        Link = linkprop([ax; ax_prev],{'CameraUpVector', 'CameraPosition', 'CameraTarget', 'XLim', 'YLim', 'ZLim'});
        setappdata(h, 'StoreTheLink', Link);
        if plotType == 1
          for j = 1 : numel(hp)
            RGB = cmapping(hp(j).CData,colorName,colorrange);
            hpCopy = copyobj(hp(j),ax_prev(end));
            set(hpCopy,'CData',RGB);
            delete(hp(j));
          end
        elseif plotType == 2
          hpCopy = copyobj(hp,ax_prev(end));
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
          if exist('hp')
            hpCopy = copyobj(hp,ax_prev(end));
            delete(hp)
          end
          if exist('hp2')
            hpCopy = copyobj(hp2,ax_prev(end));
            delete(hp2)
          end
          delete(ax);
        elseif plotType == 5
          for j = 1 : numel(hp)
            RGB = cmapping(hp(j).CData,colorName,colorrange);
            hpCopy = copyobj(hp(j),ax_prev(end));
            set(hpCopy,'CData',RGB);
            delete(hp(j));
          end
        elseif plotType == 6
          if ~isempty(variable)
            RGB = cmapping(hp.FaceVertexCData,colorName,colorrange);
            hpCopy = copyobj(hp,ax_prev(end));
            set(hpCopy,'FaceVertexCData',RGB);
            delete(hp);
          else
            hpCopy = copyobj(hp,ax_prev(end));
            delete(ax);
          end
        end
      else
        % Change to real colors
        if plotType == 1 || plotType == 5
          for j = 1 : numel(hp)
            RGB = cmapping(hp(j).CData,colorName,colorrange);
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
            RGB = cmapping(hp.FaceVertexCData,colorName,colorrange);
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
  end
  %--------------------------------------------------

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
    function [plotType, variable, ...                   % GENERAL
              xl, yl, zl, position, colorposition, ...
              colorName, alpha, colorrange, ...         % Several Types
              xslice, yslice, zslice, ...               % Slice
              ix, iy, iz, increment ...       % Quiver
              level, LineWidth, ...                     % Contour
              start, ...                                % Stream
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
          if ~find(strcmpi('newfigure',var))
            xl = ca.XLim;
          else
            xl = obj.GlobalXRange;
          end
        end
        if find(strcmpi('ylim',var))
          yl = var{ find(strcmpi('ylim',var))+1 };
        else
          if ~find(strcmpi('newfigure',var))
            yl = ca.YLim;
          else
            yl = obj.GlobalYRange;
          end
        end
        if find(strcmpi('zlim',var))
          zl = var{ find(strcmpi('zlim',var))+1 };
        else
          if ~find(strcmpi('newfigure',var))
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

      %----------------------------------------
      %   May be common to different plots
        if find(strcmpi('alpha',var))
          alpha = var{ find(strcmpi('alpha',var))+1 };
        else
          alpha = 1;
        end

        % INTERP
      %----------------------------------------
      %       Color
        if find(strcmpi('color',var))
          colorName = var{ find(strcmpi('color',var))+1 };
        else
          if (plotType == 1 | plotType == 2 | plotType == 5 | plotType == 6)
            colorName = 'parula';
          elseif plotType == 3
            colorName = [1 0 1];
          elseif plotType == 4
            colorName = [0.66, 0.66, 0.66];
          end
        end

        if find(strcmpi('colorrange',var))
          colorrange = var{ find(strcmpi('colorrange',var))+1 };
        else
          if plotType == 1
            colorrange = [min(obj.(variable),[],'all') max(obj.(variable),[],'all')];
          elseif plotType == 2
            mag = sqrt( obj.([variable,'x']).^2 + obj.([variable,'y']).^2 + obj.([variable,'z']).^2 );
            colorrange = [min(mag,[],'all') max(mag,[],'all')];
          elseif plotType == 6
            if find(strcmpi('variable',var))
              colorrange = [min(obj.(variable),[],'all') max(obj.(variable),[],'all')];
            else
              colorrange = [0 1];
            end
          else
            colorrange = [0 1];  % This does not matter as for streams and contours, I only have 1 color in the colormap
          end
        end

        if find(strcmpi('colorposition',var))
          colorposition = var{ find(strcmpi('colorposition',var))+1 };
          if ischar(colorposition) & strcmpi(colorposition,'left')
            colorposition = [0.07 0.1105 0.0112 0.8143];
          elseif ischar(colorposition) & strcmpi(colorposition,'right')
            colorposition = [0.92 0.1105 0.0112 0.8143];
          end
        else
          colorposition = [0.07 0.1105 0.0112 0.8143];
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
      %      Quiver increment
        if find(strcmpi('increment',var))
          increment = var{ find(strcmpi('increment',var))+1 };
        else
          increment = 1;
        end
    end
    %--------------------------------------------------
    function setProperties(obj,ax,cb, ...
              xl, yl, zl, position, ...
              visible, ...
              cl,colorposition,colorrange)
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

      % Colorbar stuff
      if ~isempty(cb)
        set(cb, ...
            'Position',colorposition ...
            );
        ylabel(cb,cl);
        caxis(ax,colorrange)
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
