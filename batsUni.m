classdef batsUni < bats

  properties (Hidden, GetAccess = protected, SetAccess = protected)

    % Global
    GlobalCellSize
    GlobalXRange
    GlobalYRange
    GlobalZRange
    GlobalInterpolation = false;

    %------------------------------
    %     COMPUTED FIELDS
    Xmesh, Ymesh, Zmesh

    Vorticityx, Vorticityy, Vorticityz, Vorticity
    GradPx, GradPy, GradPz, GradP
    GradPbx, GradPby, GradPbz, GradPb

  end

  methods
    function obj = batsUni(varargin)
      if isa(varargin{1},'bats')
        obj = varargin{1};
      end
      if find(strcmp('xrange',varargin))
        obj.GlobalXRange = varargin{ find(strcmp('xrange',varargin))+1 };
      else
        obj.GlobalXRange = [min(obj.x,[],'all') max(obj.x,[],'all')];
      end
      if find(strcmp('yrange',varargin))
        obj.GlobalYRange = varargin{ find(strcmp('yrange',varargin))+1 };
      else
        obj.GlobalYRange = [min(obj.y,[],'all') max(obj.y,[],'all')];
      end
      if find(strcmp('zrange',varargin))
        obj.GlobalZRange = varargin{ find(strcmp('zrange',varargin))+1 };
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
        obj.Xmesh = x;
        obj.Ymesh = y;
        obj.Zmesh = z;
      end
    end

    %--------------------------------------------------
    %     Getters
    function Global = get.Global(obj)
      Global.File = obj.GlobalFile;
      Global.Time = obj.GlobalTime;
      Global.Units= obj.GlobalUnits;
      Global.CoordinateSystem = obj.GlobalCoordinateSystem;
      Global.CellSize = obj.GlobalCellSize;
      Global.XRange = obj.GlobalXRange;
      Global.YRange = obj.GlobalYRange;
      Global.ZRange = obj.GlobalZRange;
    end
    %--------------------------------------------------
    function Derived = get.Derived(obj)
      % From bats class
      Derived.b = obj.b;
      Derived.b1 = obj.b1;
      Derived.u = obj.u;
      Derived.j = obj.j;
      Derived.jxbx = obj.jxbx;
      Derived.jxby = obj.jxby;
      Derived.jxbz = obj.jxbz;
      Derived.jxb = obj.jxbz;
      Derived.Ex = obj.Ex;
      Derived.Ey = obj.Ey;
      Derived.Ez = obj.Ez;
      Derived.E = obj.E;
      Derived.Temp = obj.Temp;
      Derived.Pb = obj.Pb;
      Derived.Beta = obj.Beta;
      Derived.Alfven = obj.Alfven;
      Derived.Vth = obj.Vth;
      Derived.Gyroradius = obj.Gyroradius;
      Derived.PlasmaFrequency = obj.PlasmaFreq;
      Derived.InertialLength = obj.InertialLength;

      % Specific for batsUni
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
    %--------------------------------------------------


%     %--------------------------------------------------
%     function obj = calc_vorticity(obj)
%
%       % Need to change because the xmesh are created using ndgrid
%       % while curl wants the grid to be made with meshgrid
%       % (which is dumb and swaps dim 1 and 2)
%       if ndims(obj.xmesh) == 2    % Means we did a cut
%         if all(obj.xmesh == obj.xmesh(1))
%           dim1 = 1; dim2 = 2;
%         elseif all(obj.ymesh == obj.ymesh(1))
%           dim1 = 1; dim2 = 2;
%         elseif all(obj.zmesh == obj.zmesh(1))
%           dim1 = 2; dim2 = 1;
%         end
%         x = permute(obj.xmesh,[dim1 dim2]);
%         y = permute(obj.ymesh,[dim1 dim2]);
%         z = permute(obj.zmesh,[dim1 dim2]);
%         ux= permute(obj.ux.data,[dim1 dim2]);
%         uy= permute(obj.uy.data,[dim1 dim2]);
%         uz= permute(obj.uz.data,[dim1 dim2]);
%       elseif ndims(obj.xmesh) == 3
%         x = permute(obj.xmesh,[2 1 3]);
%         y = permute(obj.ymesh,[2 1 3]);
%         z = permute(obj.zmesh,[2 1 3]);
%         ux= permute(obj.ux.data,[2 1 3]);
%         uy= permute(obj.uy.data,[2 1 3]);
%         uz= permute(obj.uz.data,[2 1 3]);
%       end
%
%       [curlx,curly,curlz,~] = curl(x,y,z,ux,uy,uz);
%
%       if ndims(obj.xmesh) == 2
%         curlx = permute(curlx,[dim1 dim2]);
%         curly = permute(curly,[dim1 dim2]);
%         curlz = permute(curlz,[dim1 dim2]);
%       elseif ndims(obj.xmesh) == 3
%         curlx = permute(curlx,[2 1 3]);
%         curly = permute(curly,[2 1 3]);
%         curlz = permute(curlz,[2 1 3]);
%       end
%       Re = 6371.2e3;
%
%       obj.vorticityx.data = curlx*1e3/Re;
%       obj.vorticityy.data = curly*1e3/Re;
%       obj.vorticityz.data = curlz*1e3/Re;
%       obj.vorticityx.name = 'Vorticity x';
%       obj.vorticityy.name = 'Vorticity y';
%       obj.vorticityz.name = 'Vorticity z';
%       obj.Vorticityx.units= '1/s';
%       obj.Vorticityy.units= '1/s';
%       obj.Vorticityz.units= '1/s';
%     end
%     function obj = calc_gradP(obj)
%       d = obj.xmesh(2)-obj.xmesh(1);
%       if d == 0 d = obj.ymesh(2)-obj.ymesh(1); end
%       if d == 0 d = obj.zmesh(2)-obj.zmesh(1); end
%       Re = 6371.2e3;
%       [gpx,gpy,gpz] = gradient(obj.p.data,d*Re);
%
%       obj.gradPx.data = gpx;
%       obj.gradPy.data = gpy;
%       obj.gradPz.data = gpz;
%       obj.gradP.data = sqrt(gpx.^2 + gpy.^2 + gpz.^2);
%       obj.gradPx.units= 'nN/m^3';
%       obj.gradPy.units= 'nN/m^3';
%       obj.gradPz.units= 'nN/m^3';
%       obj.gradP.units = 'nN/m^3';
%       obj.gradPx.name = 'Pressure gradient in x';
%       obj.gradPy.name = 'Pressure gradient in y';
%       obj.gradPz.name = 'Pressure gradient in z';
%       obj.gradP.name  = 'Pressure gradient mag';
%     end
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
  end   % Methods

end
