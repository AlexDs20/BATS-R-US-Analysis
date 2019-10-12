classdef sub < super
  methods
    function obj = sub(varargin)
      disp('sub Created');
    end
  end
  methods (Hidden)
    function Global = getGlobalContent(obj)
      Global = getGlobalContent@super(obj);
      Global.time = 10;
    end
  end
end
