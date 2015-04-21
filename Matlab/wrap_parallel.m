function [varargout] = wrap_parallel(parworkers, path_dep, Nloop, loop_args, fcn, args)

  Nargs  = length(args)

  % A cell with functions that take a thing "x" and an index "i_" and
  % either return "x(i_)" or just "x"
  if length(loop_args) < Nargs % loop args is a vector indices
    temp = logical(zeros(1,Nargs));
    temp(loop_args) = 1;
    loop_args = temp;
  else
    loop_args = logical(loop_args)
  end
  indx_args = cell(Nargs,1);
  indx_args(loop_args)  = {@(x,i_) x{i_}};
  indx_args(~loop_args) = {@(x,i_) x};

  % Function to return the arguments on the ith iteration (accounting
  % for the fact that some things should be indexed from on the ith
  % step, others not
  get_ith = @(i_) arrayfun(@(n) indx_args{n}(args{n},i_), 1:Nargs, 'un', 0)
  get_ith(1)
  get_ith(2)
  unpack  = @(x) x{:};

  %if parworkers
    %poolcheck = gcp('nocreate');
    %if isempty(poolcheck)
      %poolobj = parpool(parworkers); % Open up a pool
    %else
      %poolobj = poolcheck;
    %end
    %addAttachedFiles(poolobj, path_dep);

    %parfor i_ = 1:Nloop
      %feval(fcn, unpack(get_ith(i_)));
    %end
    %delete(poolobj);
  %else
    %for i_ = 1:Nloop
      %joe = feval(fcn, unpack(get_ith(i_)));
    %end
  %end


end

