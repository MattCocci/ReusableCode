function [opt] = FillDefaultSettings(defaults, varargin)

  if nargin == 1
    opt = defaults;
  else
    opt = varargin{1};

    flds = setdiff(fieldnames(defaults), fieldnames(opt));
    for f = 1:length(flds)
      if ~isfield(opt, flds{f})
        opt.(flds{f}) = defaults.(flds{f});
      end
    end
  end

end
