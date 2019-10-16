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
    function plot(varargin)

      if find(strcmp('pcolor',varargin))
        plotType = 1;
      elseif find(strcmp('quiver',varargin))
        plotType = 2;
      elseif find(strcmp('quiver3',varargin))
        plotType = 3;
      elseif find(strcmp('contour',varargin))
        plotType = 4;
      end

      if find(strcmp('variable',varargin))
        var = varargin{ find(strcmp('variable',varargin))+1 };
      end
      %------------------------------
      % To apply on all axes
      if find(strcmp('xlim',varargin))
        xl = varargin{ find(strcmp('xlim',varargin))+1 };
      else
        xl = obj.XRange;
      end
      if find(strcmp('ylim',varargin))
        yl = varargin{ find(strcmp('ylim',varargin))+1 };
      else
        yl = obj.YRange;
      end
      if find(strcmp('zlim',varargin))
        zl = varargin{ find(strcmp('zlim',varargin))+1 };
      else
        zl = obj.ZRange;
      end
      if find(strcmp('facealpha',varargin))
        facealpha = varargin{ find(strcmp('facealpha',varargin))+1 };
      else
        facealpha = 1.0;
      end
      if find(strcmp('position',varargin))
        position = varargin{ find(strcmp('position',varargin))+1 };
      else
        position = [0.126 0.11 0.7513 0.815];
      end
      if find(strcmp('colorposition',varargin))
        colorposition = varargin{ find(strcmp('colorposition',varargin))+1 };
        if colorposition == 'left'
          colorposition = [0.08 0.1105 0.0112 0.8143];
        elseif colorposition == 'right'
        end
      else
        colorposition = [0.08 0.1105 0.0112 0.8143];
      end

      % Check if we want a new figure
      if find(strcmp('newfigure',varargin))
        h = figure; set(h,'Units','Normalized','OuterPosition',[0 0 1 1],'Color',[1 1 1]);
        ax = axes; set(ax,'XLim',xl,'YLim',yl,'ZLim',zl,'Units','Normalized','Position',position);
      else
        h = gcf; set(h,'Units','Normalized','OuterPosition',[0 0 1 1],'Color',[1 1 1]);
        % If we already have a fig. I want a copy of one of the axes.
        ax = copyobj()
        set(ax,'Visible','off','XTick',[],'YTick',[],'ZTick',[]);
      end

      if plotType == 1
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
