% ERRORPROOF - To wrap functions so that any errors don't kill a program
function [returned, errored_out] = errorproof(fcn, varargin)

try
  returned = feval(fcn, args{:});;
  errored_out = 0;
catch
  returned = lasterror();
  errored_out = 1;
end
