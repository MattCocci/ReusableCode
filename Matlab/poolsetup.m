function [ poolobj ] = poolsetup(parworkers, varargin)
% POOLSETUP - Wrapper and easier way to set up a parallel pool.
%
% First checks if one is open already. If so, do nothing and just return
% the object for that parallel session. If not, open up a pool.
%
% Can also pass in files to add to the session.

  poolcheck = gcp('nocreate');
  if isempty(poolcheck)
    poolobj = parpool(parworkers); % Open up a pool
  else
    poolobj = poolcheck;
  end

  if nargin > 1
    addAttachedFiles(poolobj, varargin{1});
  end

end
