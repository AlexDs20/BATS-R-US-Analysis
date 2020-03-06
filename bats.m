classdef bats < dynamicprops

%--------------------------------------------------
%     PROPERTIES
  properties
    % General properties
    Global = struct();

    % BATS-R-US Ouput
    Output = struct();
  end

  properties (Dependent)
    % BATS-R-US Derived
    Derived = struct();
  end

  properties (Hidden, GetAccess = protected, SetAccess = protected)
    %------------------------------
    %   Physical Variables
    % Global
    GlobalFile = '';
    GlobalTime = '';
    GlobalUnits = struct();
    GlobalCoordinateSystem = '';

    % Bats Output
    x, y, z
    bx, by, bz
    b1x, b1y, b1z
    ux, uy, uz
    jx, jy, jz
    rho
    p
    e

    % BATS-R-US Derived quantities
    % Magnitudes
    b, b1, u, j

    % Vectors
    jxbx, jxby, jxbz, jxb
    Ex, Ey, Ez, E
    rhoUx, rhoUy, rhoUz, rhoU

    % Velocities
    Vix, Viy, Viz, Vi
    Vex, Vey, Vez, Ve

    % Scalars
    Pb, Vth, Temp, Beta, Alfven
    Gyroradius, PlasmaFrequency, InertialLength

    % Particles' Trajectories
    Trajectories = {};

    %------------------------------
    %   Internal Variables

    % Reduced domain indices
    IndicesReducedDomain

    % Interpolation functions (when going to batsUni)
    Interpolation = struct();

  end

%--------------------------------------------------
%     METHODS
  methods
    function obj = bats(varargin)
      % function obj = bats(varargin)
      %
      % INPUT:  If no inputs => returns empty object
      %   KWARGS:
      %         'file': string: '/path/to/cdf/file.cdf'
      %         'variables', var : with var a cell array of strings with the name of the variables to load.
      %                            If not specified, load all variables
      %                            {'x','y','z','bx','by','bz','b1x','b1y','b1z',
      %                             'ux','uy','uz','jx','jy','jz','rho','p','e'}
      % METHODS:
      %       toUniformGrid()
      %       calc_
      %
      % USE:
      if nargin == 0
        return
      elseif mod(nargin,2) == 1
        disp('Wrong number of inputs. Inputs go by pair, keyword + value');
        return
      end
      if find(strcmp('file',varargin))
        s = varargin{ find(strcmp('file',varargin))+1 };
      end
      if find(strcmp('variables',varargin))
        var = varargin{ find(strcmp('variables',varargin))+1 };
        var{end+1} = 'x'; var{end+1} = 'y'; var{end+1} = 'z';
        var = intersect(var,var);
      else
        % This should be handled differently.
        var = {'x','y','z','bx','by','bz','b1x','b1y','b1z','ux','uy','uz','jx','jy','jz','rho','p','e'};
      end
      % Treat varargin
      [data,info] = cdfread(s,'Variables',var,'ConvertEpochToDatenum',true,'CombineRecords',true);

      obj.GlobalFile = s;
      obj.GlobalTime = '';
      obj.GlobalCoordinateSystem = info.GlobalAttributes.grid_1_type;

      for i = 1 : numel(var)
        obj.GlobalUnits.(var{i}) = ...
                  char(info.VariableAttributes.units(strcmp(info.VariableAttributes.units(:,1),var{i}),2));
        obj.(var{i}) = data{i};
      end
    end % End of bats
    %------------------------------------------------------------

    %--------------------------------------------------
    %     Getters
    function Global = get.Global(obj)
      Global = obj.getGlobal;
    end
    %------------------------------------------------------------
    function Output = get.Output(obj)
      Output = obj.getOutput;
    end
    %------------------------------------------------------------
    function Derived = get.Derived(obj)
      Derived = obj.getDerived;
    end
    %------------------------------------------------------------


    %------------------------------------------------------------
    %   Modify the data in obj
    function obj = loadFields(obj,var)
      [data,info] = cdfread(obj.GlobalFile,'Variables',var,'ConvertEpochToDatenum',true,'CombineRecords',true);

      for i = 1 : numel(var)
        if numel(var) == 1
          field = var;
          fieldData = data;
        else
          field = var{i};
          fieldData = data{i};
        end

        obj.GlobalUnits.(field) = ...
                  char(info.VariableAttributes.units(strcmp(info.VariableAttributes.units(:,1),field),2));
        obj.(field) = fieldData;

        % Check the length of the data is the same as what we already have
        % In case we reduced the domain already
        if numel(obj.x)~=numel(obj.(field)) & numel(obj.IndicesReducedDomain)==numel(obj.x)
          obj.(field) = obj.(field)(obj.IndicesReducedDomain);
        elseif numel(obj.x)~=numel(obj.(field))
          warning('Problem with the size of the newly loaded data. We suggest you reduce the domain again.');
        end
      end
    end
    %------------------------------------------------------------
    function uniData = toBatsUni(obj,cellSize,var,varargin)
      xrange = [];
      yrange = [];
      zrange = [];

      if isempty(var)
        Fields = obj.listNonEmptyFields;
        var = {'bx','by','bz','ux','uy','uz','jx','jy','jz','rho','p'};
        var = intersect(var,Fields);
      end

      if find(strcmp('xrange',varargin))
        xr = varargin{ find(strcmp('xrange',varargin))+1 };
        xrange = [min(xr,[],'all') : single(cellSize) : max(xr,[],'all')];
      else
        xrange = [min(obj.x,[],'all') : single(cellSize) : max(obj.x,[],'all')];
      end
      if find(strcmp('yrange',varargin))
        yr = varargin{ find(strcmp('yrange',varargin))+1 };
        yrange = [min(yr,[],'all') : single(cellSize) : max(yr,[],'all')];
      else
        yrange = [min(obj.y,[],'all') : single(cellSize) : max(obj.y,[],'all')];
      end
      if find(strcmp('zrange',varargin))
        zr = varargin{ find(strcmp('zrange',varargin))+1 };
        zrange = [min(zr,[],'all') : single(cellSize) : max(zr,[],'all')];
      else
        zrange = [min(obj.z,[],'all') : single(cellSize) : max(obj.z,[],'all')];
      end

      uniData = batsUni(obj,'xrange',xrange,'yrange',yrange,'zrange',zrange, ...
                        'cellsize',cellSize,'interpolate',true,'variables', var);
    end
    %------------------------------------------------------------
    function [xmesh,ymesh,zmesh] = toUniformGrid(obj,obj1,cellSize,var,varargin)
      % function [xmesh,ymesh,zmesh] = toUniformGrid(obj,cellSize,var,varargin)
      %
      %   INPUT:
      %           obj = object that will contain the uniform grid
      %           obj1 = object from which it will take the data
      %           cellSize = in R units
      %           var = cellarray of variables e.g. {'bx','by'}
      %     varargin:
      %           'xrange': giving the range of values for the interp
      %           'yrange': giving the range of values for the interp
      %           'zrange': giving the range of values for the interp
      %
      %   OUTPUT:
      %           xmesh
      %           ymesh
      %           zmesh

      % Warning
      warning('This may take a long time. Make sure you have enough memory!');

      xrange = [];
      yrange = [];
      zrange = [];

      if isempty(var)     % We only care about the original fields. The rest can be quickly recomputed
        Fields = obj1.listNonEmptyFields;
        var = {'bx','by','bz','ux','uy','uz','jx','jy','jz','rho','p'};
        var = intersect(var,Fields);
      end

      if find(strcmp('xrange',varargin))
        xr = varargin{ find(strcmp('xrange',varargin))+1 };
        xrange = [min(xr,[],'all') : single(cellSize) : max(xr,[],'all')];
      else
        xrange = [min(obj1.x,[],'all') : single(cellSize) : max(obj1.x,[],'all')];
      end
      if find(strcmp('yrange',varargin))
        yr = varargin{ find(strcmp('yrange',varargin))+1 };
        yrange = [min(yr,[],'all') : single(cellSize) : max(yr,[],'all')];
      else
        yrange = [min(obj1.y,[],'all') : single(cellSize) : max(obj1.y,[],'all')];
      end
      if find(strcmp('zrange',varargin))
        zr = varargin{ find(strcmp('zrange',varargin))+1 };
        zrange = [min(zr,[],'all') : single(cellSize) : max(zr,[],'all')];
      else
        zrange = [min(obj1.z,[],'all') : single(cellSize) : max(obj1.z,[],'all')];
      end

      % Constrain the xyz used for interpolation
      shift = 5*cellSize;
      ix = find(obj1.x >= xrange(1)-shift & obj1.x <= xrange(end)+shift);
      iy = find(obj1.y >= yrange(1)-shift & obj1.y <= yrange(end)+shift);
      iz = find(obj1.z >= zrange(1)-shift & obj1.z <= zrange(end)+shift);
      ind = intersect(ix,intersect(iy,iz));
      x = double(obj1.x(ind));
      y = double(obj1.y(ind));
      z = double(obj1.z(ind));

      [xmesh,ymesh,zmesh] = ndgrid(xrange,yrange,zrange);
      xmesh = squeeze(xmesh);
      ymesh = squeeze(ymesh);
      zmesh = squeeze(zmesh);

      % Keep same info:
      for i = 1 : numel(var)
        disp(['Reshaping variable: ',var{i}]);

        F = scatteredInterpolant(x,y,z,double(obj1.(var{i})(ind)),'nearest');

        obj.(var{i}) = single(F(double(xmesh), double(ymesh), double(zmesh)));

        % Save the interpolant in case we want to trace particles
        obj.Interpolation.(var{i}) = F;
      end
    end
    %----------------------------------------
    %   Calc functions.
    %------------------------------------------------------------
    function obj = calc_b(obj)
      if isempty(obj.bx) obj.loadFields({'bx','by','bz'}); end

      obj.b = sqrt( obj.bx.^2 + obj.by.^2 + obj.bz.^2 );
      obj.GlobalUnits.b = obj.GlobalUnits.bx;
    end
    %------------------------------------------------------------
    function obj = calc_b1(obj)
      if isempty(obj.b1x) obj.loadFields({'b1x','b1y','b1z'}); end

      obj.b1 = sqrt( obj.b1x.^2 + obj.b1y.^2 + obj.b1z.^2 );
      obj.GlobalUnits.b1 = obj.GlobalUnits.b1x;
    end
    %------------------------------------------------------------
    function obj = calc_u(obj)
      if isempty(obj.ux) obj.loadFields({'ux','uy','uz'}); end

      obj.u = sqrt( obj.ux.^2 + obj.uy.^2 + obj.uz.^2 );
      obj.GlobalUnits.u = obj.GlobalUnits.ux;
    end
    %------------------------------------------------------------
    function obj = calc_rhoU(obj)
      % [rhoU]: amu \mum^-2 s-1
      if isempty(obj.ux)    obj.loadFields({'ux','uy','uz'});   end
      if isempty(obj.rho)   obj.loadFields('rho');              end

      %m = 1.6726*1e-27;
      m = 1;
      obj.rhoUx = obj.rho.* obj.ux * 1e-3 * m;
      obj.rhoUy = obj.rho.* obj.uy * 1e-3 * m;
      obj.rhoUz = obj.rho.* obj.uz * 1e-3 * m;
      obj.rhoU  = sqrt( obj.rhoUx.^2 + obj.rhoUy.^2 + obj.rhoUz.^2 );

      obj.GlobalUnits.rhoUx = 'amu mu m^-2 s^-1';
      obj.GlobalUnits.rhoUy = 'amu mu m^-2 s^-1';
      obj.GlobalUnits.rhoUz = 'amu mu m^-2 s^-1';
      obj.GlobalUnits.rhoU  = 'amu mu m^-2 s^-1';
    end
    %------------------------------------------------------------
    function obj = calc_j(obj)
      if isempty(obj.jx) obj.loadFields({'jx','jy','jz'}); end

      obj.j = sqrt( obj.jx.^2 + obj.jy.^2 + obj.jz.^2 );
      obj.GlobalUnits.j = obj.GlobalUnits.jx;
    end
    %------------------------------------------------------------
    function obj = calc_jxb(obj)
      % units:  [j]: muA/m^2
      %         [b]: nT
      %         [jxb]: fN/m^3
      if (isempty(obj.jx)|isempty(obj.ux))
        obj.loadFields({'jx','jy','jz','ux','uy','uz'});
      end

      obj.jxbx = (obj.jy.*obj.bz - obj.jz.*obj.by);
      obj.jxby = (obj.jz.*obj.bx - obj.jx.*obj.bz);
      obj.jxbz = (obj.jx.*obj.by - obj.jy.*obj.bx);
      obj.jxb  = sqrt( obj.jxbx.^2 + obj.jxby.^2 + obj.jxbz.^2 );

      obj.GlobalUnits.jxbx = 'fN/m^3';
      obj.GlobalUnits.jxby = 'fN/m^3';
      obj.GlobalUnits.jxbz = 'fN/m^3';
      obj.GlobalUnits.jxb  = 'fN/m^3';
    end
    %------------------------------------------------------------
    function obj = calc_Econv(obj)
      % units:  [u]: km/s
      %         [b]: nT
      %         [E]: mV/m
      if (isempty(obj.bx)|isempty(obj.ux))
        obj.loadFields({'bx','by','bz','ux','uy','uz'});
      end

      obj.Ex = -1.0*(obj.uy.*obj.bz - obj.uz.*obj.by) / 1e3;
      obj.Ey = -1.0*(obj.uz.*obj.bx - obj.ux.*obj.bz) / 1e3;
      obj.Ez = -1.0*(obj.ux.*obj.by - obj.uy.*obj.bx) / 1e3;
      obj.E  = sqrt( obj.Ex.^2 + obj.Ey.^2 + obj.Ez.^2 );

      obj.GlobalUnits.Ex = 'mV/m';
      obj.GlobalUnits.Ey = 'mV/m';
      obj.GlobalUnits.Ez = 'mV/m';
      obj.GlobalUnits.E = 'mV/m';
    end
    %------------------------------------------------------------
    function obj = calc_temp(obj)
      % Calculated through P = nkT
      % it uses rho for n, it is fine as long as its only H+
      if (isempty(obj.p)|isempty(obj.rho))
        obj.loadFields({'p','rho'});
      end

      eV = 1.602176634*1e-4;    % Includes the 1e15
      obj.Temp = obj.p./(obj.rho*eV);
      obj.GlobalUnits.Temp = 'eV';
    end
    %------------------------------------------------------------
    function obj = calc_pb(obj)
      % Calculates the magnetic pressure
      if isempty(obj.b) obj.calc_b; end

      mu0 = 4*pi*1e-7;
      obj.Pb = obj.b.^2/(2*mu0)*1e-9;
      obj.GlobalUnits.Pb = 'nPa';
    end
    %------------------------------------------------------------
    function obj = calc_beta(obj)
      % Calculates beta factor
      if isempty(obj.p)   obj.loadFields('p');  end
      if isempty(obj.Pb)  obj.calc_pb;          end

      obj.Beta = obj.p./obj.Pb;
      obj.GlobalUnits.Beta = '';
    end
    %------------------------------------------------------------
    function obj = calc_alfven(obj)
      % Only works for pure protons
      if isempty(obj.b)     obj.calc_b;               end
      if isempty(obj.rho)   obj.loadFields('rho');  end

      mu0 = 4*pi*1e-7;
      m = 1.6726*1e-27;
      factor = mu0*m;
      obj.Alfven = (obj.b./sqrt(factor*obj.rho))*1e-15;
      obj.GlobalUnits.Alfven = 'km/s';
    end
    %------------------------------------------------------------
    function obj = calc_vth(obj)
      % Compute it without the temperature
      % vth = sqrt(2*p/(n*m))   with p: pressure, n: number density, m: mass
      if isempty(obj.Temp) obj.calc_temp; end

      m =  1.6726219;
      eV = 1.6021766;    % Includes the 1e15
      eVm = 1e8 * eV/m;

      obj.Vth = sqrt( 2* eVm * obj.Temp ) * 1e-3;
      obj.GlobalUnits.Vth = 'km/s';
    end
    %------------------------------------------------------------
    function obj = calc_gyroradius(obj)
      if isempty(obj.u)     obj.calc_u;     end
      if isempty(obj.Vth)   obj.calc_vth;   end
      if isempty(obj.b)     obj.calc_b;   end

      v = sqrt(obj.u.^2 + obj.Vth.^2)*1000;
      b = obj.b*1e-9;
      m = 1.6276219e-27;
      q = 1.6021766e-19;
      R = 6371.2e3;

      obj.Gyroradius = m*v./(q*b*R);
      obj.GlobalUnits.Gyroradius = 'Re';
    end
    %------------------------------------------------------------
    function obj = calc_plasmafreq(obj)
      % Ion plasma frequency
      % To get the electron, just change the mass (in kg)
      if isempty(obj.rho)   obj.loadFields('rho');  end

      m = 1.6276e-27;
      n = (obj.rho)*1e6;
      eps0 = 8.8542e-12;
      q = 1.6022e-19;

      obj.PlasmaFrequency = sqrt(n*q^2/(m*eps0));
      obj.GlobalUnits.PlasmaFrequency = 'rad/s';
    end
    %------------------------------------------------------------
    function obj = calc_inertiallength(obj)
      if isempty(obj.Alfven)      obj.calc_alfven;      end
      if isempty(obj.PlasmaFrequency)  obj.calc_plasmafreq;  end

      R = 6371.2;

      obj.InertialLength = obj.Alfven./(obj.PlasmaFrequency*R);
      obj.GlobalUnits.InertialLength = 'Re';
    end
    %------------------------------------------------------------
    function calc_electronVelocity(obj)
      % uses the expression for the bulk velocity
      % and the current to recover the electron
      % and ion velocities
      %
      %   V = (mi*Vi + me*Ve)/(mi+me)
      %   J = e*n*(Vi-Ve);
      %
      %   Ve = V - (mi/(mi+me))*J/en
      %   Vi = V + (me/(mi+me))*J/en
      if ( isempty(obj.rho) | isempty(obj.jx) | isempty(obj.ux) )
        obj.loadFields({'rho','jx','jy','jz','ux','uy','uz'});
      end
      mp = 1.6726219*1e-27;
      me = 9.1093835*1e-31;
      q  = 1.6021766*1e-19;

      obj.Vex = obj.ux - 1e-15 * (mp/(mp+me))*obj.jx./(q*obj.rho);
      obj.Vey = obj.uy - 1e-15 * (mp/(mp+me))*obj.jy./(q*obj.rho);
      obj.Vez = obj.uz - 1e-15 * (mp/(mp+me))*obj.jz./(q*obj.rho);
      obj.Ve = sqrt( obj.Vex.^2 + obj.Vey.^2 + obj.Vez.^2 );

      obj.GlobalUnits.Vex = 'km/s';
      obj.GlobalUnits.Vey = 'km/s';
      obj.GlobalUnits.Vez = 'km/s';
      obj.GlobalUnits.Ve  = 'km/s';
    end
    %------------------------------------------------------------
    function calc_protonVelocity(obj)
      % uses the expression for the bulk velocity
      % and the current to recover the electron
      % and ion velocities
      %
      %   V = (mi*Vi + me*Ve)/(mi+me)
      %   J = e*n*(Vi-Ve);
      %
      %   Ve = V - (mi/(mi+me))*J/en
      %   Vi = V + (me/(mi+me))*J/en
      if ( isempty(obj.rho) | isempty(obj.jx) | isempty(obj.ux) )
        obj.loadFields({'rho','jx','jy','jz','ux','uy','uz'});
      end
      mp = 1.6726219*1e-27;
      me = 9.1093835*1e-31;
      q  = 1.6021766*1e-19;

      obj.Vix = obj.ux + 1e-15 * (me/(mp+me))*obj.jx./(q*obj.rho);
      obj.Viy = obj.uy + 1e-15 * (me/(mp+me))*obj.jy./(q*obj.rho);
      obj.Viz = obj.uz + 1e-15 * (me/(mp+me))*obj.jz./(q*obj.rho);
      obj.Vi = sqrt( obj.Vix.^2 + obj.Viy.^2 + obj.Viz.^2 );

      obj.GlobalUnits.Vix = 'km/s';
      obj.GlobalUnits.Viy = 'km/s';
      obj.GlobalUnits.Viz = 'km/s';
      obj.GlobalUnits.Vi  = 'km/s';
    end
    %------------------------------------------------------------
    %------------------------------------------------------------
    function obj = calc_all(obj)
      obj.calc_b;
      obj.calc_b1;
      obj.calc_u;
      obj.calc_rhoU;
      obj.calc_j;
      obj.calc_jxb;
      obj.calc_Econv;
      obj.calc_temp;
      obj.calc_pb;
      obj.calc_beta;
      obj.calc_alfven;
      obj.calc_vth;
      obj.calc_gyroradius;
      obj.calc_plasmafreq;
      obj.calc_inertiallength;
      obj.calc_electronVelocity;
      obj.calc_protonVelocity;
    end
    %------------------------------------------------------------

    %------------------------------------------------------------
    %   Return Data
    function obj = getData(obj,position)
      [~,idx_x] = min(abs( obj.x - position(1) ));
      [~,idx_y] = min(abs( obj.y - position(2) ));
      [~,idx_z] = min(abs( obj.z - position(3) ));

      idx = intersect(idx_x,intersect(idx_y,idx_z));

      var = obj.listNonEmptyFields;

      for i = 1 : numel(var)
        obj.(var{i}) = obj.(var{i})(idx);
      end
    end
    function obj = clearFields(obj,var)
      for i = 1 : numel(var)
        obj.(var{i}) = [];
      end
    end

    %--------------------------------------------------
    %   Domain
    function obj = reduceDomain(obj,varargin)
      if nargin == 0
        disp('No changes done. Pass arguments to reduce domain.')
        return
      end

      % Set the property to get the indices of the reduced domain
      obj.getIndexDomain(varargin{:});

      % Go through each field and only keep the desired points
      Fields = obj.listNonEmptyFields;

      for i = 1 : numel(Fields)
        obj.(Fields{i}) = obj.(Fields{i})(obj.IndicesReducedDomain);
      end
    end
    %------------------------------------------------------------
    function obj = getIndexDomain(obj,varargin)
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

      ix = find( obj.x>=min(xrange,[],'all') & obj.x<=max(xrange,[],'all') );
      iy = find( obj.y>=min(yrange,[],'all') & obj.y<=max(yrange,[],'all') );
      iz = find( obj.z>=min(zrange,[],'all') & obj.z<=max(zrange,[],'all') );
      obj.IndicesReducedDomain = intersect(ix,intersect(iy,iz));
    end
    %------------------------------------------------------------
    %   list fields
    function Fields = listNonEmptyFields(obj,varargin)
      Fields = {};
      if isempty(varargin) | strcmp(varargin{1},'Output')
        if ~isempty(obj.x)    Fields{end+1} = 'x';      end
        if ~isempty(obj.y)    Fields{end+1} = 'y';      end
        if ~isempty(obj.z)    Fields{end+1} = 'z';      end

        if ~isempty(obj.bx)   Fields{end+1} = 'bx';     end
        if ~isempty(obj.by)   Fields{end+1} = 'by';     end
        if ~isempty(obj.bz)   Fields{end+1} = 'bz';     end

        if ~isempty(obj.b1x)  Fields{end+1} = 'b1x';    end
        if ~isempty(obj.b1y)  Fields{end+1} = 'b1y';    end
        if ~isempty(obj.b1z)  Fields{end+1} = 'b1z';    end

        if ~isempty(obj.ux)   Fields{end+1} = 'ux';     end
        if ~isempty(obj.uy)   Fields{end+1} = 'uy';     end
        if ~isempty(obj.uz)   Fields{end+1} = 'uz';     end

        if ~isempty(obj.jx)   Fields{end+1} = 'jx';     end
        if ~isempty(obj.jy)   Fields{end+1} = 'jy';     end
        if ~isempty(obj.jz)   Fields{end+1} = 'jz';     end

        if ~isempty(obj.rho)  Fields{end+1} = 'rho';    end
        if ~isempty(obj.p)    Fields{end+1} = 'p';      end
        if ~isempty(obj.e)    Fields{end+1} = 'e';      end
      end

      if isempty(varargin) | strcmp(varargin{1},'Derived')
        if ~isempty(obj.b)    Fields{end+1} = 'b';      end
        if ~isempty(obj.b1)   Fields{end+1} = 'b1';     end
        if ~isempty(obj.u)    Fields{end+1} = 'u';      end
        if ~isempty(obj.j)    Fields{end+1} = 'j';      end

        if ~isempty(obj.jxbx) Fields{end+1} = 'jxbx';   end
        if ~isempty(obj.jxby) Fields{end+1} = 'jxby';   end
        if ~isempty(obj.jxbz) Fields{end+1} = 'jxbz';   end
        if ~isempty(obj.jxb)  Fields{end+1} = 'jxb';    end

        if ~isempty(obj.Ex)   Fields{end+1} = 'Ex';     end
        if ~isempty(obj.Ey)   Fields{end+1} = 'Ey';     end
        if ~isempty(obj.Ez)   Fields{end+1} = 'Ez';     end
        if ~isempty(obj.E)    Fields{end+1} = 'E';      end

        if ~isempty(obj.rhoUx) Fields{end+1} = 'rhoUx';   end
        if ~isempty(obj.rhoUy) Fields{end+1} = 'rhoUy';   end
        if ~isempty(obj.rhoUz) Fields{end+1} = 'rhoUz';   end
        if ~isempty(obj.rhoU)  Fields{end+1} = 'rhoU';    end

        if ~isempty(obj.Pb)               Fields{end+1} = 'Pb';               end
        if ~isempty(obj.Vth)              Fields{end+1} = 'Vth';              end
        if ~isempty(obj.Temp)             Fields{end+1} = 'Temp';             end
        if ~isempty(obj.Beta)             Fields{end+1} = 'Beta';             end
        if ~isempty(obj.Alfven)           Fields{end+1} = 'Alfven';           end
        if ~isempty(obj.Gyroradius)       Fields{end+1} = 'Gyroradius';       end
        if ~isempty(obj.InertialLength)   Fields{end+1} = 'InertialLength';   end
        if ~isempty(obj.PlasmaFrequency)  Fields{end+1} = 'PlasmaFrequency';  end
      end
    %------------------------------------------------------------
    end
    %------------------------------------------------------------
    function output = copyObject(input, output)
       C = metaclass(input);
       P = C.Properties;
       for k = 1:length(P)
         if ~P{k}.Dependent
           output.(P{k}.Name) = input.(P{k}.Name);
         end
       end
    end
    %------------------------------------------------------------
    %   Trace Particle
    function obj = traceParticles(obj,particles)
    % obj = traceParticles(particles)
    %
    % Trace particles in an electromagnetic fields using Boris pushing algorithm
    %
    % The fields are the already interpolated fields when gone to uniform grid and are used to get the value of their value at exact particles' position
    %
    % INPUT:
    % -----
    %
    %   fields:   structure fields.x, fields.y, fields.z          [m]
    %                       fields.bx, fields.by, fields.bz       [T]
    %                       fields.ex, fields.ey, fields.ez       [V/m]
    %
    %
    %   particles: cell array of structures
    %               particles{i} contains the information of the ith particle to be traced
    %               particles{i}.m
    %               particles{i}.q
    %               particles{i}.x
    %               particles{i}.y
    %               particles{i}.z
    %               particles{i}.vpar
    %               particles{i}.vperp
    %               particles{i}.gyrophase                     Gyrophase [degrees], angle measured from the ExB direction towards the E direction
    %               particles{i}.RefDir                        To provide for the gyrophase if there is no electric       field. It serves as the reference direction from which the gyrophase is measured
    %               particles{i}.Tmax                          Total tracing time
    %               particles{i}.max_step_length                Max length for Boris and Max distance for saving
    %               particles{i}.max_dt                        Max dt used for Boris algo
    %               particles{i}.save_dt                        Max dt at which we save
    %               particles{i}.gyroperiod_steps               Minimum amount of steps for a single gyration
    %               particles{i}.forward
    %
    %
    % EXAMPLE:
    % -------
    %     Re = 6371.2e3;
    %     res = 0.125;                            % Grid resolution
    %     particles{1}.m = 1.6726219e-27;         % Proton mass
    %     particles{1}.q = 1.602177e-19;          % e
    %     particles{1}.x = -20*Re;
    %     particles{1}.y = 0*Re;
    %     particles{1}.z = 0*Re;
    %     particles{1}.vpar = 3e5;                % 300km/s
    %     particles{1}.vperp = 2e5;               % 200km/s
    %     particles{1}.gyrophase = 0;             % gyrophase in degrees, angle from the ExB_n vector towards E, 0 = along ExB
    %     particles{1}.RefDir = [0 1 0];          % If no electric field, provide RefDir to know the reference direction from which to measure the gyrophase.
    %     particles{1}.Tmax = 10*60;                 % Trace for 10 minutes trajectory
    %     particles{1}.max_step_length = res/2;   % time reso based on speed
    %     particles{1}.max_dt = 10;               % Arbitrary time resolution 10 sec
    %     particles{1}.save_dt = 50;              % Max time at which we save
    %     particles{1}.gyroperiod_steps = 20;     % Time reso: 20 push per gyroperiod
    %     particles{1}.forward = 1;               % Trace forward in time


    % Start by checking if interpolation exists
    x = double(obj.x(:));
    y = double(obj.y(:));
    z = double(obj.z(:));
    if isfield(obj.Interpolation,'bx') | isfield(obj.Interpolation,'by') | isfield(obj.Interpolation,'bz')
    else
      disp('Missing interpolation for the magnetic field. Doing it now:');
      disp('Interpolation Bx');
      obj.Interpolation.bx = scatteredInterpolant(x,y,z,double(obj.bx(:)),'nearest');
      disp('Interpolation By');
      obj.Interpolation.by = scatteredInterpolant(x,y,z,double(obj.by(:)),'nearest');
      disp('Interpolation Bz');
      obj.Interpolation.bz = scatteredInterpolant(x,y,z,double(obj.bz(:)),'nearest');
    end
    if isfield(obj.Interpolation,'Ex') | isfield(obj.Interpolation,'Ey') | isfield(obj.Interpolation,'Ez')
    else
      obj.calc_Econv;
      disp('Missing interpolation for the electric field. Doing it now:');
      warning('Note that this uses the Convective E Field. This is probably not true for your simulation and should NOT be used in an article!');
      disp('Interpolation Ex');
      obj.Interpolation.Ex = scatteredInterpolant(x,y,z,double(obj.Ex(:)),'nearest');
      disp('Interpolation Ey');
      obj.Interpolation.Ey = scatteredInterpolant(x,y,z,double(obj.Ey(:)),'nearest');
      disp('Interpolation Ez');
      obj.Interpolation.Ez = scatteredInterpolant(x,y,z,double(obj.Ez(:)),'nearest');
    end

    obj.Trajectories = particles;

    Bx = obj.Interpolation.bx;
    By = obj.Interpolation.by;
    Bz = obj.Interpolation.bz;
    Ex = obj.Interpolation.Ex;
    Ey = obj.Interpolation.Ey;
    Ez = obj.Interpolation.Ez;

    % Trace each particle:
    for i = 1 : numel(particles)
      % New position variable
      pos(1,:) = [particles{i}.x, particles{i}.y, particles{i}.z];

      % Get the field at the particle's position
      b(1) = Bx(pos(1,1),pos(1,2),pos(1,3));
      b(2) = By(pos(1,1),pos(1,2),pos(1,3));
      b(3) = Bz(pos(1,1),pos(1,2),pos(1,3));
      e(1) = Ex(pos(1,1),pos(1,2),pos(1,3));
      e(2) = Ey(pos(1,1),pos(1,2),pos(1,3));
      e(3) = Ez(pos(1,1),pos(1,2),pos(1,3));

      % Make sure we only have Eperp
      ep = obj.PerpendicularE(e,b);

      % Initialise the velocity: give the vpar and vperp
      vel(1,:) = obj.initialiseVelocity(particles{i},ep,b);

      % Get the time step
      dt_2 = particles{i}.max_dt;
      [dt,dt_2] = obj.calculateTimeStep(particles{i},vel(1,:),b,dt_2);

      % Need to make half a step back first
      vel(1,:) = obj.boris(vel(1,:),ep,b,-sign(particles{i}.forward)*dt/2, ...
                  particles{i}.q, particles{i}.m);

      % Variable to keep track of the time
      t = -sign(particles{i}.forward)*dt/2;

      % These will keep track the positions and velocities between the known position/vel and the one calculated
      pos(2,:) = pos(1,:);
      vel(2,:) = vel(1,:);

      % Variables to save
      save{i}.t(1) = t;
      save{i}.pos(1,:) = pos(1,:);
      save{i}.vel(1,:) = vel(1,:);

      % Loop over the time
      j = 0;
      while t < particles{i}.Tmax
        pos(1,:) = pos(2,:);
        vel(1,:) = vel(2,:);
        % Get local E and B
        b(1) = Bx(pos(1,1),pos(1,2),pos(1,3));
        b(2) = By(pos(1,1),pos(1,2),pos(1,3));
        b(3) = Bz(pos(1,1),pos(1,2),pos(1,3));
        e(1) = Ex(pos(1,1),pos(1,2),pos(1,3));
        e(2) = Ey(pos(1,1),pos(1,2),pos(1,3));
        e(3) = Ez(pos(1,1),pos(1,2),pos(1,3));

        % Make sure we only have Eperp
        ep = obj.PerpendicularE(e,b);

        % Get the time step to make
        [dt,dt_2] = calculateTimeStep(particles{i},vel(1,:),b,dt_2);

        % Get new velocity
        vel(2,:) = obj.boris(vel(1,:),ep,b,sign(particles{i}.forward)*dt, ...
                        particles{i}.q,particles{i}.m);
        % Push the positions
        pos(2,:) = pos(1,:) + vel(2,:).*sign(particles{i}.forward)*dt_2;

        % Increment Time/iteration
        t = t + dt_2;
        j = j+1;

        % Save the previous data if needed
        save{i} = obj.UpdateSaveData(save{i},pos,vel,t,dt_2,particles{i},j);

        disp( sprintf('Simulation time %0.2f s of %0.2f s',t,particles{i}.Tmax) );
      end   % Time loop

      % Save the info in the particles cell array
      obj.Trajectories{i}.t = save{i}.t;
      obj.Trajectories{i}.pos = save{i}.pos;
      obj.Trajectories{i}.vel = save{i}.vel;

    end     % particle loop

    end     % function
    %------------------------------------------------------------

  end   % Methods

  methods (Hidden)
    %------------------------------
    %     ACTUAL GETTERS
    function Global = getGlobal(obj)
      Global.File = obj.GlobalFile;
      Global.Time = obj.GlobalTime;
      Global.Units= obj.GlobalUnits;
      Global.CoordinateSystem = obj.GlobalCoordinateSystem;
    end
    %------------------------------------------------------------
    function Output = getOutput(obj)
      Output.x = obj.x;
      Output.y = obj.y;
      Output.z = obj.z;
      Output.bx = obj.bx;
      Output.by = obj.by;
      Output.bz = obj.bz;
      Output.b  = obj.b;
      Output.b1x = obj.b1x;
      Output.b1y = obj.b1y;
      Output.b1z = obj.b1z;
      Output.b1  = obj.b1;
      Output.ux = obj.ux;
      Output.uy = obj.uy;
      Output.uz = obj.uz;
      Output.u  = obj.u;
      Output.jx = obj.jx;
      Output.jy = obj.jy;
      Output.jz = obj.jz;
      Output.j  = obj.j;
      Output.rho = obj.rho;
      Output.p = obj.p;
      Output.e = obj.e;
    end
    %------------------------------------------------------------
    function Derived = getDerived(obj)
      Derived.b = obj.b;
      Derived.b1 = obj.b1;
      Derived.u = obj.u;
      Derived.j = obj.j;
      Derived.jxbx = obj.jxbx;
      Derived.jxby = obj.jxby;
      Derived.jxbz = obj.jxbz;
      Derived.jxb = obj.jxb;
      Derived.Ex = obj.Ex;
      Derived.Ey = obj.Ey;
      Derived.Ez = obj.Ez;
      Derived.E  = obj.E;
      Derived.rhoUx = obj.rhoUx;
      Derived.rhoUy = obj.rhoUy;
      Derived.rhoUz = obj.rhoUz;
      Derived.rhoU  = obj.rhoU;
      Derived.Vix = obj.Vix;
      Derived.Viy = obj.Viy;
      Derived.Viz = obj.Viz;
      Derived.Vi  = obj.Vi;
      Derived.Vex = obj.Vex;
      Derived.Vey = obj.Vey;
      Derived.Vez = obj.Vez;
      Derived.Ve  = obj.Ve;
      Derived.Temp = obj.Temp;
      Derived.Pb = obj.Pb;
      Derived.Beta = obj.Beta;
      Derived.Alfven = obj.Alfven;
      Derived.Vth = obj.Vth;
      Derived.Gyroradius = obj.Gyroradius;
      Derived.PlasmaFrequency = obj.PlasmaFrequency;
      Derived.InertialLength = obj.InertialLength;
      Derived.Trajectories = obj.Trajectories;
    end
    %------------------------------------------------------------
  end

  methods (Hidden, Access = protected)
    %------------------------------------------------------------
    %   Functions for particle tracing
    function save = UpdateSaveData(obj,save,pos,vel,t,dt,particle,j)
      % Checks if the data needs to be save and does it

      % Spatial step since save
      DR = sqrt(sum((save.pos(end,:) - pos(2,:)).^2));
      % Time step since save
      DT = save.t(end) - t;
      % Check iteration
      if mod(j,nearest(particle.gyroperiod_steps/20)) == 0
        iter = true;
      else
        iter = false;
      end
      if DR > particle.max_step_length | DT > particle.save_dt | iter
        save.t(end+1,1) = t-dt;
        save.pos(end+1,:) = pos(1,:);
        save.vel(end+1,:) = vel(1,:);           % Note this is not fully correct as the pos and vel are staggered
      end
    end
    %------------------------------
    function [dt,dt_2] = calculateTimeStep(obj,particle,vel,b,dt_2)
      % calculates the time steps needed for the current interation given the conditions

      q = particle.q;
      m = particle.m;
      dt0 = particle.max_dt;
      steps = particle.gyroperiod_steps;
      max_length = particle.max_step_length;

      % Cyclotron frequency
      wc = (q/m) * sqrt(dot(b,b));
      dt1 = abs(2*pi/wc)/steps;

      % Step length
      dt2 = max_length / sqrt(dot(vel,vel));

      % mean value
      dt_1 = min([dt0,dt1,dt2]);    % New step
      dt = 0.5*(dt_1 + dt_2);       % Mean step
      dt_2 = dt_1;                  % Update the old stop to hold the new one
    end
    %------------------------------
    function ep = PerpendicularE(obj,e,b)
      % Only get E perpendicular (this should only be done if running ideal MHD!!)
      % B*B
      b2 = dot(b,b);

      % b normalized
      bn = b./sqrt(b2);

      % Make sure E perp to b
      ep = e - dot(e,bn)*bn;
    end
    %------------------------------
    function new_vel = boris(obj,V,e,b,dt,q,m)
      % Gets the updated velocity using Brois algo.

      % Set constant local variables
      eta = q/m;
      etadt2 = eta*dt/2;

      % Calculate t and vmin (first half acceleration)
      t = etadt2*b;
      vmin = V + etadt2*e;

      % Cross product (vmid = vmin + vmin x t)
      vmid = vmin + cross(vmin,t);

      % Calculate s
      s = 2*t./(1+t.^2);

      % Another cross product (vplus = vmin + vmid x s)
      vplus = vmin + cross(vmid,s);

      % Loop and set new velocities (second half acceleration)
      new_vel = vplus + etadt2*e;
    end
    %------------------------------
    function V = initialiseVelocity(obj,particle,ep,b)
      % If no electric field, and not given, the reference direction from which the gyrophase is measured is [0 1 0];

      % B^2
      b2 = dot(b,b);

      % B normalized
      bn = b./sqrt(b2);

      % E and ExB directions
      if ep==[0,0,0]
        % Get a reference direction (not ExB because there is no E field)
        if isfield(particle,'RefDir')
          ref = particle.RefDir;
        else
          ref = [0 1 0];
        end
        ExB = ref - dot(ref,bn).*bn;
        ExB_n = ExB./sqrt(dot(ExB,ExB));
        en = cross(bn,ExB_n);
      else
        en = ep./sqrt(dot(ep,ep));
        % ExB-drift direction
        ExB_n = cross(en,bn);
        % get the en to be perp to bn and ExB_n
        en = cross(bn,ExB_n);
      end

      % Add all the velocities
      % 3 origins:  - the component along the b field, as given by user
      %             - thermal component as given by the user taking into account the gyrophase, part of it is along the ExB component and some along the component of E perp to B

      % ????????    - We add the ExB velocity here and it should therefore not be given by the user
      V = particle.vpar.*bn + ...
          particle.vperp.*cosd(particle.gyrophase).*ExB_n + ...
          particle.vperp.*sind(particle.gyrophase).*en;
          %particle.vperp.*sind(particle.gyrophase).*en + ...
          %cross(ep,b)./b2;
    end
    %------------------------------------------------------------
  end   % Methods (hidden, protected)

end     % Class
