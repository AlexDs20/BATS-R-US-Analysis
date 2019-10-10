classdef bats < handle

%--------------------------------------------------
%     PROPERTIES

  properties
    time
    x
    y
    z
    bx
    by
    bz
    b
    b1x
    b1y
    b1z
    ux
    uy
    uz
    u
    jx
    jy
    jz
    j
    rho
    p
    e
    temp
    jxb
    Ex
    Ey
    Ez
    E
    beta
    alfven
    vorticityx
    vorticityy
    vorticityz
    gradPx
    gradPy
    gradPz
    gradP
    JxBx
    JxBy
    JxBz
    Pb
    gradPbx
    gradPby
    gradPbz
    gradPb
    vth
    gyroradius
    plasmafreq
    inertiallength
    coordinateSystem = '';
  end

%--------------------------------------------------
%     METHODS
  methods
    function obj = bats(s,varargin)
      %
      % INPUT: Compulsory:
      %         s: string: path to cdf file
      %        KWARGS:
      %         'Variables', var : with var a cell array of strings with the name of the variables to load.
      %                            If not specified, load all variables
      %                            {'x','y','z','bx','by','bz','b1x','b1y','b1z',
      %                             'ux','uy','uz','jx','jy','jz','rho','p','e'}
      % METHODS:
      %       toUniformGrid()
      %       calc_
      %
      % USE:

      var = {};
      if find(strcmp('Variables',varargin))
        var = varargin{ find(strcmp('Variables',varargin))+1 };
        var{end+1} = 'x'; var{end+1} = 'y'; var{end+1} = 'z';
        var = intersect(var,var);
      else
        % This should be handled differently.
        var = {'x','y','z','bx','by','bz','b1x','b1y','b1z','ux','uy','uz','jx','jy','jz','rho','p','e'};
      end
      % Treat varargin
      [data,info] = cdfread(s,'Variables',var,'ConvertEpochToDatenum',true,'CombineRecords',true);

      obj.time = [];

      for i = 1 : numel(var)
        obj.(var{i}).name = var{i};
        obj.(var{i}).units = char(info.VariableAttributes.units(strcmp(info.VariableAttributes.units(:,1),var{i}),2));
        obj.(var{i}).data = data{i};
        obj.(var{i}).coordinateSystem = info.GlobalAttributes.grid_1_type;
      end
    end % End of bats

    function uniData = toUniformGrid(obj,cellSize,var,varargin)
    % function uniData = toUniformGrid(obj,cellSize,var,varargin)
    %   INPUT:
    %           obj = bats object
    %           cellSize = in R units
    %           var = cellarray of variables e.g. {'bx','by'}
    %     varargin:
    %           'xrange': giving the range of values for the interp
    %           'yrange': giving the range of values for the interp
    %           'zrange': giving the range of values for the interp
    %
    %   OUTPUT:
    %           uniData properties:
    %               xmesh
    %               ymesh
    %               zmesh
    %               (var{i})
    %

      % Warning
      warning('This may take a long time. Make sure you have enough memory!');

      xrange = [];
      yrange = [];
      zrange = [];

      if isempty(var)
        var = {'bx','by','bz','ux','uy','uz','jx','jy','jz','rho','p'};
      end

      if find(strcmp('xrange',varargin)) & find(strcmp('yrange',varargin)) & find(strcmp('zrange',varargin))
        xr = varargin{ find(strcmp('xrange',varargin))+1 };
        yr = varargin{ find(strcmp('yrange',varargin))+1 };
        zr = varargin{ find(strcmp('zrange',varargin))+1 };
        xrange = [xr(1):single(cellSize):xr(end)];
        yrange = [yr(1):single(cellSize):yr(end)];
        zrange = [zr(1):single(cellSize):zr(end)];
      else
        xrange = [min(obj.x.data):single(cellSize):max(obj.x.data)];
        yrange = [min(obj.y.data):single(cellSize):max(obj.y.data)];
        zrange = [min(obj.z.data):single(cellSize):max(obj.z.data)];
      end

      % Constrain the xyz used for interpolation
      shift = 5*cellSize;
      ix = find(obj.x.data >= xrange(1)-shift & obj.x.data <= xrange(end)+shift);
      iy = find(obj.y.data >= yrange(1)-shift & obj.y.data <= yrange(end)+shift);
      iz = find(obj.z.data >= zrange(1)-shift & obj.z.data <= zrange(end)+shift);
      ind = intersect(ix,intersect(iy,iz));
      x = double(obj.x.data(ind));
      y = double(obj.y.data(ind));
      z = double(obj.z.data(ind));

      [xmesh,ymesh,zmesh] = ndgrid(xrange,yrange,zrange);
      uniData.xmesh = squeeze(xmesh);
      uniData.ymesh = squeeze(ymesh);
      uniData.zmesh = squeeze(zmesh);

      % Keep same info:
      for i = 1 : numel(var)
        disp(['Reshaping variable: ',var{i}]);
        uniData.(var{i}).name =  obj.(var{i}).name;
        uniData.(var{i}).units = obj.(var{i}).units;
        uniData.(var{i}).coordinateSystem = obj.(var{i}).coordinateSystem;

        F = scatteredInterpolant(x,y,z,double(obj.(var{i}).data(ind)));
        uniData.(var{i}).data = single(F(double(uniData.xmesh), double(uniData.ymesh), ...
                                         double(uniData.zmesh)));
      end
    end

    %----------------------------------------
    %   Calc functions.
    %   Those with derivatives should only
    %   be computed on uniform grids

    function obj = calc_temp(obj)
    % Calculated through P = nkT
    % it uses rho for n, it is fine as long as its only H+
      eV = 6241.50935;    % nPa/cm^3 --> eV

      obj.temp.data = eV*obj.p.data/obj.rho.data;
      obj.temp.units= 'eV';
      obj.temp.name = 'Temp';
    end
    function obj = calc_b(obj)
      keyboard
      obj.b = obj.bx;
      obj.b.data = sqrt( obj.bx.data.^2 + obj.by.data.^2 + obj.bz.data.^2 );
    end
    function obj = calc_j(obj)
      obj.j = obj.jx;
      obj.j.data = sqrt( obj.jx.data.^2 + obj.jy.data.^2 + obj.jz.data.^2 );
    end
    function obj = calc_u(obj)
      obj.u = obj.ux;
      obj.u.data = sqrt( obj.ux.data.^2 + obj.uy.data.^2 + obj.uz.data.^2 );
    end
    function obj = calc_E(obj)
      ux = obj.ux.data; uy = obj.uy.data; uz = obj.uz.data;
      bx = obj.bx.data; by = obj.by.data; bz = obj.bz.data;

      obj.Ex.data = -1.0*(uy.*bz - uz.*by) / 1000.0
      obj.Ey.data = -1.0*(uz.*bx - ux.*bz) / 1000.0
      obj.Ez.data = -1.0*(ux.*by - uy.*bx) / 1000.0
      obj.E.data  = sqrt( obj.Ex.data.^2 + obj.Ey.data.^2 + obj.Ez.data.^2 );

      obj.Ex.units= 'mV/m'; obj.Ey.units= 'mV/m'; obj.Ez.units= 'mV/m'; obj.E.units = 'mV/m';
      obj.Ex.name = 'Ex'; obj.Ey.name = 'Ey'; obj.Ez.name = 'Ez'; obj.E.name = 'E';
      obj.Ex.coordinateSystem = obj.ux.coordinateSystem;
      obj.Ey.coordinateSystem = obj.uy.coordinateSystem;
      obj.Ez.coordinateSystem = obj.uz.coordinateSystem;
      obj.E.coordinateSystem = {};
    end
    function obj = calc_Pb(obj)
    % Calculates the magnetic pressure
      mu0 = 4*pi*1e-7;
      B = [obj.bx.data obj.by.data obj.bz.data];

      obj.Pb.data = sum(B.^2,2)/(2*mu0)*1e-9;
      obj.Pb.units= 'nPa';
      obj.Pb.name = 'Magnetic Pressure';
    end
    function obj = calc_beta(obj)
    % Calculates beta factor
      obj.beta.data = obj.p.data./obj.Pb.data;
      obj.beta.name = 'beta';
      obj.beta.units= '';
    end
    function obj = calc_alfven(obj)
    % Only works for pure protons
      mu0 = 4*pi*1e-7 * 1.6726*1e-27 * 1e6;   % #/cm^3 to kg/m^3

      obj.alfven.data = obj.b.data./sqrt(mu0*obj.rho.data);
      obj.alfven.name = 'Alfven Speed';
      obj.alfven.units= 'km/s';
    end
    function obj = calc_JxB(obj)

      obj.JxBx.data = (obj.jy.data.*obj.bz.data - obj.jz.data.*obj.by.data)*1e-6;
      obj.JxBy.data = (obj.jz.data.*obj.bx.data - obj.jx.data.*obj.bz.data)*1e-6;
      obj.JxBz.data = (obj.jx.data.*obj.by.data - obj.jy.data.*obj.bx.data)*1e-6;
      obj.JxB.data  = sqrt( obj.JxBx.data.^2 + obj.JxBy.data.^2 + obj.JxBz.data.^2 );

      obj.JxBx.units = 'nN/m^3';
      obj.JxBy.units = 'nN/m^3';
      obj.JxBz.units = 'nN/m^3';
      obj.JxB.units  = 'nN/m^3';
      obj.JxBx.name = 'JxB_x';
      obj.JxBy.name = 'JxB_y';
      obj.JxBz.name = 'JxB_z';
      obj.JxB.name  = 'JxB';
      obj.JxBx.coordinateSystem = obj.bx.coordinateSystem;
      obj.JxBy.coordinateSystem = obj.bx.coordinateSystem;
      obj.JxBz.coordinateSystem = obj.bx.coordinateSystem;
    end
    function obj = calc_vth(obj)
      m = 1.6276e-27; % H+ mass kg
      n = (obj.rho.data/m)*1e6;

      obj.vth.data = sqrt( obj.p.data*1e-9/(n*m) )/ 1e3;
      obj.vth.units= 'km/s';
      obj.vth.name = 'Thermal Velocity';
    end
    function obj = calc_gyroradius(obj)
      v_squared = obj.u.data.^2 + obj.vth.data.^2;
      v = sqrt(v_squared)*1000;
      b = obj.b.data*1e-9;
      m = 1.6276e-27;
      q = 1.6022e-19;
      R = 6371.2e3;

      obj.gyroradius.data = m*v./(q.*B*R);
      obj.gyroradius.units= 'Re';
      obj.gyroradius.name = 'Gyroradius';
    end
    function obj = calc_plasmafreq(obj)
      m = 1.6276e-27;
      n = (obj.rho.data/m)*1e6;
      q = 1.6022e-19;

      obj.plasmafreq.data = sqrt(4*pi*n*q^2/m);
      obj.plasmafreq.units= 'rad/s';
      obj.plasmafreq.name = 'Ion plasma frequency';
    end
    function obj = calc_inertiallength(obj)
      R = 6371.2e3;

      obj.inertiallength.data = obj.alfven.data./(obj.plasmafreq*R);
      obj.inertiallength.units= 'Re';
      obj.inertiallength.name = 'Ion inertial length';
    end
    function obj = calc_vorticity(obj)

      % Need to change because the xmesh are created using ndgrid
      % while curl wants the grid to be made with meshgrid
      % (which is dumb and swaps dim 1 and 2)
      if ndims(obj.xmesh) == 2    % Means we did a cut
        if all(obj.xmesh == obj.xmesh(1))
          dim1 = 1; dim2 = 2;
        elseif all(obj.ymesh == obj.ymesh(1))
          dim1 = 1; dim2 = 2;
        elseif all(obj.zmesh == obj.zmesh(1))
          dim1 = 2; dim2 = 1;
        end
        x = permute(obj.xmesh,[dim1 dim2]);
        y = permute(obj.ymesh,[dim1 dim2]);
        z = permute(obj.zmesh,[dim1 dim2]);
        ux= permute(obj.ux.data,[dim1 dim2]);
        uy= permute(obj.uy.data,[dim1 dim2]);
        uz= permute(obj.uz.data,[dim1 dim2]);
      elseif ndims(obj.xmesh) == 3
        x = permute(obj.xmesh,[2 1 3]);
        y = permute(obj.ymesh,[2 1 3]);
        z = permute(obj.zmesh,[2 1 3]);
        ux= permute(obj.ux.data,[2 1 3]);
        uy= permute(obj.uy.data,[2 1 3]);
        uz= permute(obj.uz.data,[2 1 3]);
      end

      [curlx,curly,curlz,~] = curl(x,y,z,ux,uy,uz);

      if ndims(obj.xmesh) == 2
        curlx = permute(curlx,[dim1 dim2]);
        curly = permute(curly,[dim1 dim2]);
        curlz = permute(curlz,[dim1 dim2]);
      elseif ndims(obj.xmesh) == 3
        curlx = permute(curlx,[2 1 3]);
        curly = permute(curly,[2 1 3]);
        curlz = permute(curlz,[2 1 3]);
      end
      Re = 6371.2e3;

      obj.vorticityx.data = curlx*1e3/Re;
      obj.vorticityy.data = curly*1e3/Re;
      obj.vorticityz.data = curlz*1e3/Re;
      obj.vorticityx.name = 'Vorticity x';
      obj.vorticityy.name = 'Vorticity y';
      obj.vorticityz.name = 'Vorticity z';
      obj.Vorticityx.units= '1/s';
      obj.Vorticityy.units= '1/s';
      obj.Vorticityz.units= '1/s';
    end
    function obj = calc_gradP(obj)
      d = obj.xmesh(2)-obj.xmesh(1);
      if d == 0 d = obj.ymesh(2)-obj.ymesh(1); end
      if d == 0 d = obj.zmesh(2)-obj.zmesh(1); end
      Re = 6371.2e3;
      [gpx,gpy,gpz] = gradient(obj.p.data,d*Re);

      obj.gradPx.data = gpx;
      obj.gradPy.data = gpy;
      obj.gradPz.data = gpz;
      obj.gradP.data = sqrt(gpx.^2 + gpy.^2 + gpz.^2);
      obj.gradPx.units= 'nN/m^3';
      obj.gradPy.units= 'nN/m^3';
      obj.gradPz.units= 'nN/m^3';
      obj.gradP.units = 'nN/m^3';
      obj.gradPx.name = 'Pressure gradient in x';
      obj.gradPy.name = 'Pressure gradient in y';
      obj.gradPz.name = 'Pressure gradient in z';
      obj.gradP.name  = 'Pressure gradient mag';
    end
    function obj = calc_gradPb(obj)
      d = obj.xmesh(2)-obj.xmesh(1);
      if d == 0 d = obj.ymesh(2)-obj.ymesh(1); end
      if d == 0 d = obj.zmesh(2)-obj.zmesh(1); end
      Re = 6371.2e3;
      [gpx,gpy,gpz] = gradient(obj.Pb.data,d*Re);

      obj.gradPbx.data = gpx;
      obj.gradPby.data = gpy;
      obj.gradPbz.data = gpz;
      obj.gradPb.data = sqrt(gpx.^2 + gpy.^2 + gpz.^2);
      obj.gradPbx.units= 'nN/m^3';
      obj.gradPby.units= 'nN/m^3';
      obj.gradPbz.units= 'nN/m^3';
      obj.gradPb.units = 'nN/m^3';
      obj.gradPbx.name = 'Magnetic Pressure gradient in x';
      obj.gradPby.name = 'Magnetic Pressure gradient in y';
      obj.gradPbz.name = 'Magnetic Pressure gradient in z';
      obj.gradPb.name  = 'Magnetic Pressure gradient mag';
    end

    %----------------------------------------
    %   Visualisation methods

  end   % Methods

end     % Class

% classdef bats2d
% %--------------------------------------------------
% %     PROPERTIES
%   properties (Dependent = true)
%     time
%     data
%     coordinateSystem = '';
%     cut = '';
%     dim1
%     dim2
%     name = '';
%     units = '';
%   end
%
% %--------------------------------------------------
% %     METHODS
%   methods
%     function b2d = bats2d(bats,cut,cut_val,varargin)
%       % Takes the bats with uniform grid (used toUniformGrid) and cut through in the cut dimension at cut_val
%       dim1 = [];
%       dim2 = [];
%     end
%     function h = batsPlot()
%       h = [];
%     end
%   end
% end
