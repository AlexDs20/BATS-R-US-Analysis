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
  properties (Hidden, Access = private)
    %------------------------------
    %     Derived Fields used
    %     for div(tensor) computation
    BxMu0, ByMu0, BzMu0
    rhoUx, rhoUy, rhoUz

    %------------------------------
    %  Units of additional fields

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
    %   Calculating things (Derivative based)
    %----------------------------------------
    function obj = calc_vorticity(obj)

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
        ux= permute(obj.ux,[dim1 dim2]);
        uy= permute(obj.uy,[dim1 dim2]);
        uz= permute(obj.uz,[dim1 dim2]);
      elseif ndims(obj.x) == 3
        x = permute(obj.x,[2 1 3]);
        y = permute(obj.y,[2 1 3]);
        z = permute(obj.z,[2 1 3]);
        ux= permute(obj.ux,[2 1 3]);
        uy= permute(obj.uy,[2 1 3]);
        uz= permute(obj.uz,[2 1 3]);
      end

      [curlx,curly,curlz,~] = curl(x,y,z,ux,uy,uz);

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

      obj.Vorticityx = curlx*1e3/Re;
      obj.Vorticityy = curly*1e3/Re;
      obj.Vorticityz = curlz*1e3/Re;
      obj.Vorticity =  sqrt( curlx.^2 + curly.^2 + curlz.^2 )*1e3/Re;

      obj.GlobalUnits.Vorticityx = '1/s';
      obj.GlobalUnits.Vorticityy = '1/s';
      obj.GlobalUnits.Vorticityz = '1/s';
      obj.GlobalUnits.Vorticity  = '1/s';
    end

    %----------------------------------------
    function obj = calc_grad(obj,var)
      d = obj.x(2)-obj.x(1);
      if d == 0 d = obj.y(2)-obj.y(1); end
      if d == 0 d = obj.z(2)-obj.z(1); end
      Re = 6371.2e3;
      [gx,gy,gz] = gradient(obj.(var),d*Re);

      newvarx = ['grad',var,'x'];
      newvary = ['grad',var,'y'];
      newvarz = ['grad',var,'z'];
      newvar  = ['grad',var    ];

      obj.newvarx = gx;
      obj.newvary = gy;
      obj.newvarz = gz;
      obj.newvar  = sqrt(gx.^2 + gy.^2 + gz.^2);

      if strcmp(var,'p') | strcmp(var,'Pb')
        obj.GlobalUnits.newvarx = 'nN/m^3';
        obj.GlobalUnits.newvary = 'nN/m^3';
        obj.GlobalUnits.newvarz = 'nN/m^3';
        obj.GlobalUnits.newvar  = 'nN/m^3';
      else
        obj.GlobalUnits.newvarx = [obj.GlobalUnits.var,'/m'];
        obj.GlobalUnits.newvary = [obj.GlobalUnits.var,'/m'];
        obj.GlobalUnits.newvarz = [obj.GlobalUnits.var,'/m'];
        obj.GlobalUnits.newvar  = [obj.GlobalUnits.var,'/m'];
      end
    end
    %----------------------------------------
    function obj = calc_divTensor(obj,field1,field2)
    %   calculates divergence of a tensor:
    %             \partial_i (field1_i field2_j)
    %   by computing: \partial_i(field1_i) field2_j + field1_i \partial_i (field2_j)
    %   that is:
    %                     div(field1) field2        +     field1 dot grad(field2)

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
    end


%     function obj = calc_gradPb(obj)
%       d = obj.xmesh(2)-obj.xmesh(1);
%       if d == 0 d = obj.ymesh(2)-obj.ymesh(1); end
%       if d == 0 d = obj.zmesh(2)-obj.zmesh(1); end
%       Re = 6371.2e3;
%       [gpx,gpy,gpz] = gradient(obj.Pb.data,d*Re);
%
%       obj.gradPbx.data = gpx;
%       obj.gradPby.data = gpy;
%       obj.gradPbz.data = gpz;
%       obj.gradPb.data = sqrt(gpx.^2 + gpy.^2 + gpz.^2);
%       obj.gradPbx.units= 'nN/m^3';
%       obj.gradPby.units= 'nN/m^3';
%       obj.gradPbz.units= 'nN/m^3';
%       obj.gradPb.units = 'nN/m^3';
%       obj.gradPbx.name = 'Magnetic Pressure gradient in x';
%       obj.gradPby.name = 'Magnetic Pressure gradient in y';
%       obj.gradPbz.name = 'Magnetic Pressure gradient in z';
%       obj.gradPb.name  = 'Magnetic Pressure gradient mag';
%     end

  end

  methods (Hidden, Access = private)
    function obj = calc_divTensor(obj,field1,field2)
    end

    function obj = calc_grad(obj,var)
    end

    function obj = calc_curl(obj,var)
    end
  end

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
