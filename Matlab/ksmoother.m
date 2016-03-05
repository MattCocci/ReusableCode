function [R] = ksmoother(varargin)
% Kalman Smoother progam
%
% Accepts either:
% 1. The sysmats structure (see kfilter.m for what it should contain)
%    and the KF structure that is output from kfilter.m when fullout=1
% 2. All the usual arguments to kfilter.m, in which case this program
%    will call kfilter.m to first filter the data, then smooth it.
%
% Smoothed estimates get appended to KF structure

  % Make sure you have the filtered data in R
  if nargin == 2
    sysmats = varargin{1};
    R       = varargin{2};
  else
    [L, R]  = kfilter(varargin{:});
    sysmats = varargin{2};
  end

  % Extract objects/info that we need for smoothing
  capT    = size(R.s_filt,2);
  T       = sysmats.T;
  s_filt  = R.s_filt;
  s_pred  = R.s_pred;
  ss_filt = R.ss_filt;
  ss_pred = R.ss_pred;

  % Initialize matrices to hold smoothed estimates
  s_smooth  = repmat(s_filt(:,end),  1, capT);
  ss_smooth = repmat(ss_filt(:,end), 1, 1, capT);

  % Smoothe the dataset, starting at T-1 since T is already in there
  for t = capT-1:-1:1
    J = ss_filt(:,t)*T*inv(ss_pred(:,t+1));

    s_smooth(:,t)    = s_filt(:,t)    + J*(s_smooth(:,t+1)    - s_pred(:,t+1));
    ss_smooth(:,:,t) = ss_filt(:,:,t) + J*(ss_smooth(:,:,t+1) - ss_pred(:,:,t+1))*J';
  end


  R.s_smooth  = s_smooth;
  R.ss_smooth = ss_smooth;

end
