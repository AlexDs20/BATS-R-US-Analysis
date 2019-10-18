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
      if find(strcmp('xrange',varargin))
        xr = varargin{ find(strcmp('xrange',varargin))+1 };
        obj.GlobalXRange = [min(xr,[],'all') max(xr,[],'all')];
      else
        obj.GlobalXRange = [min(obj.x,[],'all') max(obj.x,[],'all')];
      end
      if find(strcmp('yrange',varargin))
        yr = varargin{ find(strcmp('yrange',varargin))+1 };
        obj.GlobalYRange = [min(yr,[],'all') max(yr,[],'all')];
      else
        obj.GlobalYRange = [min(obj.y,[],'all') max(obj.y,[],'all')];
      end
      if find(strcmp('zrange',varargin))
        zr = varargin{ find(strcmp('zrange',varargin))+1 };
        obj.GlobalZRange = [min(zr,[],'all') max(zr,[],'all')];
      else
        obj.GlobalZRange = [min(obj.z,[],'all') max(obj.z,[],'all')];
      end
      if find(strcmp('cellsize',varargin))
        obj.GlobalCellSize = varargin{ find(strcmp('cellsize',varargin))+1 };
      else
        d = obj.x(2)-obj.x(1);
        if d == 0 d = obj.y(2)-obj.y(1); end
        if d == 0 d = obj.z(2)-obj.z(1); end
      end
      if find(strcmp('interpolate',varargin))
        obj.GlobalInterpolation = varargin{ find(strcmp('interpolate',varargin))+1 };
      end
      if find(strcmp('variables',varargin))
        var = varargin{ find(strcmp('variables',varargin))+1 };
      end

      if isa(varargin{1},'bats') & ~isempty(obj.GlobalCellSize) ...
                  & find(strcmp('variables',varargin)) & obj.GlobalInterpolation

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

      [vecx,vecy,vecz,vec] = obj.calc_curl('b');

      obj.jx = 1e-3 * vecx / mu0;
      obj.jy = 1e-3 * vecy / mu0;
      obj.jz = 1e-3 * vecz / mu0;
      obj.j = sqrt( obj.jx.^2 + obj.jy.^2 + obj.jz.^2 );

      obj.GlobalUnits.jx = 'muA/m^2';
      obj.GlobalUnits.jy = 'muA/m^2';
      obj.GlobalUnits.jz = 'muA/m^2';
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
    %----------------------------------------
    function obj = calc_gradPb(obj)
      [obj.GradPbx,obj.GradPby,obj.GradPbz,obj.GradPb] = ...
                   obj.calc_grad('Pb');

      obj.GlobalUnits.GradPbx = 'nN/m^3';
      obj.GlobalUnits.GradPby = 'nN/m^3';
      obj.GlobalUnits.GradPbz = 'nN/m^3';
      obj.GlobalUnits.GradPb  = 'nN/m^3';
    end
    %----------------------------------------
    function obj = calc_gradP(obj)
      [obj.GradPx,obj.GradPy,obj.GradPz,obj.GradP] = ...
                   obj.calc_grad('p');

      obj.GlobalUnits.GradPx = 'nN/m^3';
      obj.GlobalUnits.GradPy = 'nN/m^3';
      obj.GlobalUnits.GradPz = 'nN/m^3';
      obj.GlobalUnits.GradP  = 'nN/m^3';
    end
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

      obj.DivBBx = vecx * 1e-9 / mu0;
      obj.DivBBy = vecy * 1e-9 / mu0;
      obj.DivBBz = vecz * 1e-9 / mu0;
      obj.DivBB  = vec  * 1e-9 / mu0;

      obj.GlobalUnits.DivBBx = 'nN/m^3';
      obj.GlobalUnits.DivBBy = 'nN/m^3';
      obj.GlobalUnits.DivBBz = 'nN/m^3';
      obj.GlobalUnits.DivBB  = 'nN/m^3';
    end
    %----------------------------------------
    function obj = calc_divRhoUU(obj)
      [vecx,vecy,vecz,vec] = obj.calc_divTensor('rhoU','u');

      %    rhoU        u       div
      % amu cm-2 s-1  km s-1    m-1
      % mp(kg) *1e12 m^-2 s^-2
      % mp(kg) * 1e16 * nN/m^3

      mp = 1.67262*1e-27;

      obj.DivRhoUUx = mp * vecx * 1e16;
      obj.DivRhoUUy = mp * vecy * 1e16;
      obj.DivRhoUUz = mp * vecz * 1e16;
      obj.DivRhoUU  = mp * vec  * 1e16;

      obj.GlobalUnits.DivRhoUUx = 'nN/m^3';
      obj.GlobalUnits.DivRhoUUy = 'nN/m^3';
      obj.GlobalUnits.DivRhoUUz = 'nN/m^3';
      obj.GlobalUnits.DivRhoUU  = 'nN/m^3';
    end
    %----------------------------------------
    %     Reduce Domain
    %----------------------------------------
    function obj = reduceDomain(obj,varargin)
      % KWARGS: 'xrange', 'yrange' and 'zrange', otherwise
      if find(strcmp('xrange',varargin))
        xrange = varargin{ find(strcmp('xrange',varargin))+1 };
      else
        xrange = [min(obj.x,[],'all') max(obj.x,[],'all')];
      end
      if find(strcmp('yrange',varargin))
        yrange = varargin{ find(strcmp('yrange',varargin))+1 };
      else
        yrange = [min(obj.y,[],'all') max(obj.y,[],'all')];
      end
      if find(strcmp('zrange',varargin))
        zrange = varargin{ find(strcmp('zrange',varargin))+1 };
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
    function plot(obj,varargin)
      %
      % INPUTS:
      %       GLOBALS:
      %         - 'slice', 'quiver', 'contour', 'stream'
      %         - 'variable', value
      %         - 'xlim', value ([min max])
      %         - 'ylim', value ([min max])
      %         - 'zlim', value ([min max])
      %         - 'position', value ([x0 y0 xsize ysize]) (in Normalized units)
      %         - 'color', colorName    : a colormap or a [r g b] value
      %
      %       Slice and Quiver
      %         - 'alpha', value
      %         - 'colorposition', value ('left' or 'right')
      %         - 'colorrange', value ([min max])
      %
      %       Slice and Contour:
      %         - 'xslice', 'yslice', 'zslice', value (in Re)
      %
      %       Quiver and Contour:
      %         - 'xrange', 'yrange', 'zrange', value ([min max])
      %
      %       Contour and Stream:
      %         - 'LineWidth'
      %
      %       Contour:
      %         - 'level'
      %
      %       Stream:
      %         - 'start', values ( [Nx3] array )
      %         - 'forward', 'backward' : to only trace forward/backward, otherwise both

      %----------------------------------------
      %     TREAT INPUT

      [plotType, variable, ...                    % GENERAL
       xl, yl, zl, position, colorposition, ...
       colorName, alpha, colorrange, ...          % Several Types
       xslice, yslice, zslice, ...                % Slice
       qIndices,  ...                             % Quiver
       contourLevel, LineWidth, ...               % Contour
       start ...                                  % Stream
       ] ...
        = obj.setInputValues(varargin(:));

      %------------------------------
      % To apply on all axes

      % Check if we want a new figure
      if find(strcmp('newfigure',varargin))
        h = figure; set(h,'Units','Normalized','OuterPosition',[0 0 1 1],'Color',[1 1 1]);
        visible = 'on';
        ax = axes;
      else
        h = gcf;
        visible = 'off';
        ax_prev = findall(h,'type','axes');
        % Create a new axes, copy of the previous one
        ax = copyobj(ax_prev(1),h);
      end
        % GET ALL INFO OF AXES

      if plotType == 1
        [x,y,z] = meshgrid(unique(obj.x),unique(obj.y),unique(obj.z));
        field = permute(obj.(variable),[2 1 3]);

        hp = slice( ax, x,y,z, field, single(xslice),single(yslice),single(zslice) );
        set(hp, ...
            'EdgeColor', 'none', ...
            'FaceAlpha', alpha ...
            );
        colormap(ax,colorName);
        cb = colorbar(ax);
        cl = [variable,' [',obj.GlobalUnits.(variable),']'];

      elseif plotType == 2
        x = obj.x(qIndices);
        y = obj.y(qIndices);
        z = obj.z(qIndices);
        fieldx = obj.([variable,'x'])(qIndices);
        fieldy = obj.([variable,'y'])(qIndices);
        fieldz = obj.([variable,'z'])(qIndices);

        hp = quiver3(ax, x,y,z, fieldx,fieldy,fieldz );

        % Colors
        SetQuiverColor(hp,colormap(ax,colorName), ...
                        'range', colorrange, ...
                        'mags',  sqrt(fieldx.^2+fieldy.^2+fieldz.^2) );
        cb = colorbar(ax);
        cl = [variable,' [',obj.GlobalUnits.([variable,'x']),']'];

        % SetLength?
      elseif plotType == 3
        [x,y,z] = meshgrid(unique(obj.x),unique(obj.y),unique(obj.z));
        field = permute(obj.(variable),[2 1 3]);
        hp = contourslice(ax, x,y,z, field, single(xslice),single(yslice),single(zslice), ...
                          [contourLevel contourLevel]);
        colormap(ax,colorName)
        set(hp(:),'LineWidth',LineWidth);

        cb = [];
        cl = [];
      elseif plotType == 4
        [x,y,z] = meshgrid( double(unique(obj.x)),double(unique(obj.y)),double(unique(obj.z)));
        fieldx = double(permute(obj.([variable,'x']),[2 1 3]));
        fieldy = double(permute(obj.([variable,'y']),[2 1 3]));
        fieldz = double(permute(obj.([variable,'z']),[2 1 3]));

        if isempty(find(strcmp('backward',varargin)))
          hp = streamline(ax, stream3(x,y,z, fieldx,fieldy, fieldz, ...
                                  start(:,1), start(:,2), start(:,3)) );
        end

        if isempty(find(strcmp('forward',varargin)))
          hp2= streamline(ax, stream3(x,y,z, -fieldx,-fieldy,-fieldz, ...
                                  start(:,1), start(:,2), start(:,3)) );
        end

        set(hp,'color',colorName,'LineWidth',LineWidth);
        set(hp2,'color',colorName,'LineWidth',LineWidth);

        cb = [];
        cl = [];
      else
        cb = []; cl = [];
      end

      obj.setProperties(ax,cb,xl,yl,zl,position,visible, ...
              cl,colorposition,colorrange);

      if (plotType == 1 | plotType == 2 | plotType == 3) ...
           & isempty(find(strcmp('newfigure',varargin)))
        % Linkaxes:
        Link = linkprop([ax; ax_prev],{'CameraUpVector', 'CameraPosition', 'CameraTarget', 'XLim', 'YLim', 'ZLim'});
        setappdata(gcf, 'StoreTheLink', Link);

        %scopy = copyobj(s,ax_prev);
        %delete(s);
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
              qIndices, ...                             % Quiver
              contourLevel, LineWidth, ...              % Contour
              start ...                                 % Stream
              ] ...
            = setInputValues(obj,var)

      %----------------------------------------
      %           GENERAL
      %----------------------------------------
        if find(strcmp('slice',var))
          plotType = 1;
        elseif find(strcmp('quiver',var))
          plotType = 2;
        elseif find(strcmp('contour',var))
          plotType = 3;
        elseif find(strcmp('stream',var))
          plotType = 4;
        end

        if find(strcmp('variable',var))
          variable = var{ find(strcmp('variable',var))+1 };
        end


        if find(strcmp('xlim',var))
          xl = var{ find(strcmp('xlim',var))+1 };
        else
          xl = obj.GlobalXRange;
        end
        if find(strcmp('ylim',var))
          yl = var{ find(strcmp('ylim',var))+1 };
        else
          yl = obj.GlobalYRange;
        end
        if find(strcmp('zlim',var))
          zl = var{ find(strcmp('zlim',var))+1 };
        else
          zl = obj.GlobalZRange;
        end

        if find(strcmp('position',var))
          position = var{ find(strcmp('position',var))+1 };
        else
          position = [0.126 0.11 0.7513 0.815];
        end

      %----------------------------------------
      %   May be common to different plots
      %----------------------------------------
        if find(strcmp('alpha',var))
          alpha = var{ find(strcmp('alpha',var))+1 };
        else
          alpha = 1;
        end

        if find(strcmp('color',var))
          colorName = var{ find(strcmp('color',var))+1 };
        else
          if (plotType == 1 | plotType == 2)
            colorName = 'parula';
          elseif plotType == 3
            colorName = [1 0 1];
          elseif plotType == 4
            colorName = [0.7, 0.7, 0.7];
          end
        end

      %----------------------------------------
      %       Slice and Quiver
      %----------------------------------------
        if find(strcmp('colorrange',var))
          colorrange = var{ find(strcmp('colorrange',var))+1 };
        else
          if plotType == 1
            colorrange = [min(obj.(variable),[],'all') max(obj.(variable),[],'all')];
          elseif plotType == 2
            mag = sqrt( obj.([variable,'x']).^2 + obj.([variable,'y']).^2 + obj.([variable,'z']).^2 );
            colorrange = [min(mag,[],'all') max(mag,[],'all')];
          else
            colorrange = [];
          end
        end

        if find(strcmp('colorposition',var))
          colorposition = var{ find(strcmp('colorposition',var))+1 };
          if ischar(colorposition) & strcmp(colorposition,'left')
            colorposition = [0.07 0.1105 0.0112 0.8143];
          elseif ischar(colorposition) & strcmp(colorposition,'right')
            colorposition = [0.92 0.1105 0.0112 0.8143];
          end
        else
          colorposition = [0.07 0.1105 0.0112 0.8143];
        end
      %----------------------------------------
      %       Slice and Contour
      %----------------------------------------
        if plotType == 1 | plotType == 3
          if find(strcmp('xslice',var)), xslice = var{ find(strcmp('xslice',var))+1 };
          else, xslice = [];
          end
          if find(strcmp('yslice',var)), yslice = var{ find(strcmp('yslice',var))+1 };
          else, yslice = [];
          end
          if find(strcmp('zslice',var)), zslice = var{ find(strcmp('zslice',var))+1 };
          else, zslice = [];
          end
        else
          xslice = [];
          yslice = [];
          zslice = [];
        end
      %----------------------------------------
      %      Quiver and Contour
      %----------------------------------------
        if (plotType == 2 | plotType == 3)
          % Domain
          if find(strcmp('xrange',var)), xrange = var{ find(strcmp('xrange',var))+1 };
          else, xrange = [min(obj.x,[],'all') max(obj.x,[],'all')];
          end
          if find(strcmp('yrange',var)), yrange = var{ find(strcmp('yrange',var))+1 };
          else, yrange = [min(obj.y,[],'all') max(obj.y,[],'all')];
          end
          if find(strcmp('zrange',var)), zrange = var{ find(strcmp('zrange',var))+1 };
          else, zrange = [min(obj.z,[],'all') max(obj.z,[],'all')];
          end
          qIndices = obj.domainIndex('xrange',xrange,'yrange',yrange,'zrange',zrange);
        else
          qIndices = [];
        end
      %----------------------------------------
      %      Contour and Stream
      %----------------------------------------
        if (plotType == 3 | plotType == 4)
          if find(strcmp('linewidth',var)), LineWidth = var{ find(strcmp('linewidth',var))+1 };
          else, LineWidth = [];
          end
        else
          LineWidth = [];
        end
      %----------------------------------------
      %      Contour
      %----------------------------------------
        if plotType == 3
          if find(strcmp('level',var)), contourLevel = var{ find(strcmp('level',var))+1 };
          else, contourLevel = [];
          end
        else, contourLevel = [];
        end
      %----------------------------------------
      %      Contour
      %----------------------------------------
        if plotType == 4
          if find(strcmp('start',var)), start = var{ find(strcmp('start',var))+1 };
          else, start = [-20 0 0];
          end
        else, start = [];
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
      if find(strcmp('xrange',varargin))
        xrange = varargin{ find(strcmp('xrange',varargin))+1 };
      else
        xrange = [min(obj.x,[],'all') max(obj.x,[],'all')];
      end
      if find(strcmp('yrange',varargin))
        yrange = varargin{ find(strcmp('yrange',varargin))+1 };
      else
        yrange = [min(obj.y,[],'all') max(obj.y,[],'all')];
      end
      if find(strcmp('zrange',varargin))
        zrange = varargin{ find(strcmp('zrange',varargin))+1 };
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
