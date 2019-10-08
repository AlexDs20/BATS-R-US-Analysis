classdef Bats
  % Bats.BatsField

  properties (Dependent = true)
    date
  end
end

classdef BatsField

%--------------------------------------------------
%     PROPERTIES

  properties (Dependent = true)
    time
    data
    coordinateSystem = '';
    name = '';
    units = '';
  end

%--------------------------------------------------
%     METHODS
  methods
    function obj = BatsField(s,varargin)
      %
      % INPUT: Compulsory:
      %         s: string: path to cdf file
      %        KWARGS:
      %         'Variables', var : with var a cell array of strings with the name of the variables to load.
      %                            If not specified, load all variables
      % METHODS:
      %       toUniformGrid()
      %       calc_
      %
      % USE:

      var = {};
      % Treat varargin
      [data,info] = cdfread(s,'Variables',var,'ConvertEpochToDatenum',true,'CombineRecords',true);
      for i = 1 : 1 : numel(var)
        obj.(var{i}).time = [];
        obj.(var{i}).name = var{i};
        obj.(var{i}).units = char(info.VariableAttributes.units(strcmp(info.VariableAttributes.units(:,1),var{i}),2));
        obj.(var{i}).data = data{i};
        obj.(var{i}).coordinateSystem = info.GlobalAttributes.grid_1_type;
      end
    end
    function uniData = toUniformGrid(cellSize,var,varargin)
      xrange = [];
      yrange = [];
      zrange = [];

      x = obj.x.data;
      y = obj.y.data;
      z = obj.z.data;
      [xmesh,ymesh,zmesh] = ndgrid(xrange,yrange,zrange);

      uniData.xmesh = squeeze(xmesh);
      uniData.ymesh = squeeze(ymesh);
      uniData.zmesh = squeeze(zmesh);

      % Keep same info:
      for i : numel(var)
        uniData.(var{i}).time =  obj.time;
        uniData.(var{i}).name =  obj.name;
        uniData.(var{i}).units = obj.units;
        uniData.(var{i}).coordinateSystem = obj.coordinateSystem;

        uniData.(var{i}).data = interpn(x,y,z,obj.(var{i}), ...
                                uniData.xmesh, uniData.ymesh, uniData.zmesh);
      end
      function calc_gradP()
      end
      function calc_JxB()
      end
      function calc_Pb()
      end
      function calc_gradPb()
      end
      function calc_E()
      end

    end
  end

end

classdef bats2d
%--------------------------------------------------
%     PROPERTIES
  properties (Dependent = true)
    time
    data
    coordinateSystem = '';
    cut = '';
    dim1
    dim2
    name = '';
    units = '';
  end

%--------------------------------------------------
%     METHODS
  methods
    function b2d = bats2d(BatsField,cut,cut_val,varargin)
      % Takes the BatsField with uniform grid (used toUniformGrid) and cut through in the cut dimension at cut_val
      dim1 = [];
      dim2 = [];
    end
    function h = batsPlot()
      h = [];
    end
  end
end
