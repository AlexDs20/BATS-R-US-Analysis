classdef super < handle

  properties
    Global
  end

  properties (Hidden, GetAccess = protected, SetAccess = protected)
    GlobalFile = '';
  end


  methods
    function obj = super(varargin)
      disp('Created super!');
    end
    function Global = get.Global(obj)
      Global = obj.getGlobalContent;
    end
  end
  methods(Hidden)
    function Global = getGlobalContent(obj)
      Global = struct();
      Global.File = obj.GlobalFile;
    end
  end
end
