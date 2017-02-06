function [R] = shootingProgram(fcns, beta, delta, lambda, T)

  % Unpack functions
  f     = fcns.f;
  df    = fcns.df;
  dfinv = fcns.dfinv;
  u     = fcns.u;
  du    = fcns.du;
  duinv = fcns.duinv;

  % crit*kss (steady state) is the convergence criterion
  crit = .00000015;

  % (Consumption) = (output) + (nondepreciated k) - (new k)
  c = @(kold,knew) f(kold) + (1-delta)*kold - knew;

  % Steady state capital
  kss   = dfinv((1/beta) - (1-delta));

  % Initial capital level
  k0 = lambda*kss;

  % Solution for paths csol and ksol
  ksol    = nan(T,1);
  ksol(1) = k0;

  % At each iteration, take k_{t-1} as given, want to find k_t
  for t=2:T

    % Store k_{t-1}
    k_t1 = ksol(t-1);

    % Initialize bounds for k_t, the thing we want to guess
    kmin = k_t1;
    kmax = kss;

    % Keep trying different guesses for k_t until you find a good guess
    keepGuessing = 1;
    while keepGuessing

      % A guess for k_t, the midpoint of the interval [kmin,kmax]
      k_t = (kmin+kmax)/2;

      % Initialize k_{s-1} and k_{s-2}
      k_s1 = k_t;
      k_s2 = k_t1;

      % Iterate on that initial guess to find k_s for s >= t+1
      keepIterating = 1;
      while keepIterating

        % Find k_s given k_{s-1} and k_{s-2}
        k_s = f(k_s1) + (1-delta)*k_s1 ...
              - duinv(...
                  du( c(k_s2,k_s1) ) / ...
                  (beta*( df(k_s1) + (1-delta))) ...
              );

        %if c(k_s2,k_s1)<=0
          %[k_s kss]
        %end
        %if t <=2
          %[duinv(du(c(k_s2,k_s1)) / ...
                  %(beta*( df(k_s1) + (1-delta))))]
          %keyboard
        %end

        % Check to see if k path is going to zero (k_t was too low)
        if k_s <= k_s1
          kmin = k_t;
          keepIterating = 0; % We'll exit inner while loop

        % Check to see if kguess is going beyond kss (k_t was too high)
        elseif k_s > kss
          kmax = k_t;
          keepIterating = 0; % We'll exit inner while loop
        end

        % Update k_s1 and k_s2
        k_s2 = k_s1;
        k_s1 = k_s;
      end

      % Check whether interval [kmin,kmax] sufficiently small enough to
      % say that we've found k_t and should exit outer while loop
      if abs(kmax-kmin) < crit*kss
        keepGuessing = 0;
      end

    end

    % If you exited the loop, you have your k_t
    ksol(t) = k_t;
  end

  % Pack output
  R.kss   = kss;
  R.css   = c(kss,kss);
  R.ysol  = f(ksol);
  R.ksol  = ksol;
  R.csol  = [c(ksol(1:(T-1)), ksol(2:T)); NaN];
  R.srate = 1 - R.csol./R.ysol;

end

