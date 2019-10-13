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

        [x,y,z] = obj.toUniformGrid(varargin{1},obj.GlobalCellSize, var);
        obj.x = x;
        obj.y = y;
        obj.z = z;
      end
    end

    %----------------------------------------
    %   Overloading some simple calc
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

    units = 'muA/m^2'
    end
    %----------------------------------------
    %   Some additional variables
    %----------------------------------------
    function calc_speciesVelocity(obj)
    % uses the expression for the bulk velocity
    % and the current to recover the electron
    % and ion velocities
    %
    %   V = (mi*Vi + me*Ve)/(mi+me)
    %   J = e*n*(Vi-Ve);
    %
    %   Vi = J/(e*n) + Ve
    %   V = (mi*(J/en + Ve) + me*Ve)/(mi+me)
    %   V = (mi*J/en + (mi+me)*Ve)/(mi+me)
    %   V = (mi*J/en)/(mi+me) + Ve
    %
    %   Ve = V - (mi*J/en)/(mi+me)
    %   Ve = V - (mi/(mi+me))*J/en
    %
    %   Vi = V - (mi/(mi+me))*J/en + J/en
    %   Vi = V + (me/(mi+me))*J/en
    %
    %   Ve = V - (mi/(mi+me))*J/en
    %   Vi = V + (me/(mi+me))*J/en
      mp = 1.6726219*1e-27;
      me = 9.1093835*1e-31;
      q  = 1.6021766*1e-19;

      Vex = ux - (mp/(mp+me))*jx./(q*rho);
      Vey = uy - (mp/(mp+me))*jy./(q*rho);
      Vez = uz - (mp/(mp+me))*jz./(q*rho);

      Vix = ux + (me/(mp+me))*jx./(q*rho);
      Viy = uy + (me/(mp+me))*jy./(q*rho);
      Viz = uz + (me/(mp+me))*jz./(q*rho);

    end

    %----------------------------------------
    %   Calculating Derivative based things
    %----------------------------------------
    function obj = calc_vorticity(obj)
      [Vortx, Vorty, Vortz, Vort] = obj.calc_curl('u');

      obj.Vorticityx = Vortx*1e3;
      obj.Vorticityy = Vorty*1e3;
      obj.Vorticityz = Vortz*1e3;
      obj.Vorticity  = Vort *1e3;

      obj.GlobalUnits.Vorticityx = '1/s';
      obj.GlobalUnits.Vorticityy = '1/s';
      obj.GlobalUnits.Vorticityz = '1/s';
      obj.GlobalUnits.Vorticity  = '1/s';
    end

    %----------------------------------------
    function obj = calc_gradPb(obj)
      [obj.gradPbx,obj.gradPby,obj.gradPbz,obj.gradPb] = ...
                   obj.calc_grad('Pb');

      obj.GlobalUnits.gradPbx = 'nN/m^3';
      obj.GlobalUnits.gradPby = 'nN/m^3';
      obj.GlobalUnits.gradPbz = 'nN/m^3';
      obj.GlobalUnits.gradPb  = 'nN/m^3';
    end
    %----------------------------------------
    function obj = calc_gradP(obj)
      [obj.gradPx,obj.gradPy,obj.gradPz,obj.gradP] = ...
                   obj.calc_grad('P');

      obj.GlobalUnits.gradPx = 'nN/m^3';
      obj.GlobalUnits.gradPy = 'nN/m^3';
      obj.GlobalUnits.gradPz = 'nN/m^3';
      obj.GlobalUnits.gradP  = 'nN/m^3';
    end
    %----------------------------------------
    function obj = calc_divBB(obj)
    %   calculates divergence of a tensor:
    %             \partial_i (field1_i field2_j)
    %   by computing: \partial_i(field1_i) field2_j + field1_i \partial_i (field2_j)
    %   that is:
    %                     div(field1) field2        +     field1 dot grad(field2)

      [vecx,vecy,vecz,vec] = obj.calc_divTensor('b','b');
      mu0 = 4*pi*1e-7;

      % div(BB/mu0)
      % (1/mu0) * nT * nT /m
      % (1/(4*pi*1e-7)) * 1e-9 * nN/m^3

      vecx = vecx * 1e-9 / mu0;
      vecy = vecy * 1e-9 / mu0;
      vecz = vecz * 1e-9 / mu0;
      vec  = vec  * 1e-9 / mu0;

      units = 'nN/m^3';

    end
    %----------------------------------------
    function obj = calc_divRhoUU(obj)
      [vecx,vecy,vecz,vec] = obj.calc_divTensor('rhoU','u');

      % kg m-2 s-1 km s-1 m-1
      % 1e3 kg m s-2 m-3
      % 1e3 1e9 nN/m-3

      vecx = vecx * 1e12;
      vecy = vecy * 1e12;
      vecz = vecz * 1e12;
      vec  = vec  * 1e12;

      units = 'nN/m^3'

    end
    %----------------------------------------

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
        if colorposition = 'left'
          colorposition = [0.08 0.1105 0.0112 0.8143];
        elseif colorposition = 'right'
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

      if plotType = 1
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
    %                     div(field1) field2        +     field1 dot grad(field2)
    %
    %   Units: [field1][field2]/m

      field1x = [field1,'x']; field1y = [field1,'y']; field1z = [field1,'z'];
      field2x = [field2,'x']; field2y = [field2,'y']; field2z = [field2,'z'];

      Re = 6371.2e3;

      div = divergence(obj.x,obj.y,obj.z, ...
                  obj.(field1x),obj.(field1y),obj.(field1z));

      [gradxX,gradyX,gradzX] = gradient(obj.(field2x),obj.GlobalCellSize*Re);
      [gradxY,gradyY,gradzY] = gradient(obj.(field2y),obj.GlobalCellSize*Re);
      [gradxZ,gradyZ,gradzZ] = gradient(obj.(field2z),obj.GlobalCellSize*Re);

      vecx = div .* obj.(field2x) ...
            + obj.(field1x).*gradxX + obj.(field1y).*gradyX + obj.(field1z).*gradzX;

      vecy = div .* obj.(field2y) ...
            + obj.(field1x).*gradxY + obj.(field1y).*gradyY + obj.(field1z).*gradzY;

      vecz = div .* obj.(field2z) ...
            + obj.(field1x).*gradxZ + obj.(field1y).*gradyZ + obj.(field1z).*gradzZ;

      vec = sqrt(vecx.^2 + vecy.^2 + vecz.^2);
    end
    %--------------------------------------------------

    function [vecx, vecy, vecz, vec] = calc_grad(obj,var)
    % Output Units: [inputs]/m
      Re = 6371.2e3;
      [vecx,vecy,vecz] = gradient(obj.(var),obj.GlobalCellSize*Re);

      vec  = sqrt(vecx.^2 + vecy.^2 + vecz.^2);
    end
    %--------------------------------------------------

    function [curlx, curly, curlz, curl] = calc_curl(obj,var)
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
        varx= permute(obj.varx,[dim1 dim2]);
        vary= permute(obj.vary,[dim1 dim2]);
        varz= permute(obj.varz,[dim1 dim2]);
      elseif ndims(obj.x) == 3
        x = permute(obj.x,[2 1 3]);
        y = permute(obj.y,[2 1 3]);
        z = permute(obj.z,[2 1 3]);
        varx= permute(obj.varx,[2 1 3]);
        vary= permute(obj.vary,[2 1 3]);
        varz= permute(obj.varz,[2 1 3]);
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
      curl  = sqrt( curlx.^2 + curly.^2 + curlz.^2 );
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
    function Derived = getDerived(obj)
      Derived = getDerived@bats(obj);

      Derived.Vorticity  = obj.Vorticity;
      Derived.Vorticityx = obj.Vorticityx;
      Derived.Vorticityy = obj.Vorticityy;
      Derived.Vorticityz = obj.Vorticityz;

      Derived.GradP  = obj.GradP;
      Derived.GradPx = obj.GradPx;
      Derived.GradPy = obj.GradPy;
      Derived.GradPz = obj.GradPz;

      Derived.GradPb  = obj.GradPb;
      Derived.GradPbx = obj.GradPbx;
      Derived.GradPby = obj.GradPby;
      Derived.GradPbz = obj.GradPbz;
    end
  end

end
