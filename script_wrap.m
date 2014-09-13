% OVERVIEW
% --------
% script_wrap.m: A function that can wrap a script (or a cell of
%                commands) so that you can keep environments/workspaces
%                separate.
%
% INPUTS
% ------
% script: Either a script name (as a string) within the working path or
%         a cell of instructions to step through
%
% OPTIONAL INPUTS
% ---------------
% vararagin   Additional arguments allow you to return variables after
%             running a script or set of commands
% - unpack    The first optional argument, this is an indicator for
%             whether to unpack the variables and return them
%             individually (unpack=1) or to return them in a cell
%             structure (unpack=0)
% - ret_vars  The second optional argument. This is a cell with the
%             names of the variables to be returned
% 
function [varargout] = script_wrap(script, varargin)

% If given a cell, rather than a script, run the elements one by one
if iscell(script)
  for cmd = 1:length(script)
    eval(script{cmd});
  end
else
  eval(script);
end

%% Check if there are items to be returned and return them
if nargin > 1
  unpack = varargin{1};
  ret_vars = varargin{2};

  % If you want to return the variables individually
  if unpack
    for ii = 1:length(ret_vars)
      varargout(ii) = {eval(ret_vars{ii})};
    end

  % If you want to return the variables all in one big struct
  else
    for ii = 1:length(ret_vars)
      to_ret.(ret_vars{ii}) = eval(ret_vars{ii});
    end
    varargout{1} = to_ret;
  end
end

end
