function [] = ValueFunctionIteration(F, beta, choice, state, shock)

  thresh = 10e-2;

  nChoices     = length(choice);
  nGridChoices = cellfun(@length, choice);
  nGridState   = length(state);
  nGridShock   = length(state);
  [outputs{1:Nchoices}] = ndgrid(choice{:});

  Vold = 0*ones(nGridState, nGridShock);
  V    = Vold;

  while max(max(abs(V-Vold))) > thresh || niter < 1

    Vold = V;

    %% SLOW BUT MORE INTUITIVE CODE BLOCK %%
    % Loop over possible starting values of capital
    for i = 1:Ngridk

      % Loop over possible values of z
      for j = 1:Ngridz


        % Matrix of output levels under different labor choices;
        % - Each row is a different labor choice;
        % - Each column is the same
        % - All entries are conditional on zdraw z(j) and initial
        %   capital level k(i)
        % - Entry a,b in the matrix is output using L(a) units of labor
        ymat = repmat(f(z(j),k(i),L), 1, Ngridk);

        % Compute c under different next-period capital choices
        cmat = ymat - kmat;

        % Any negative consumption has to go
        cmat(cmat<0) = 0;

        % Form the RHS of the Bellman Equation as a matrix
        BellmanRHS = u(cmat,Lmat) + repmat([beta*(Vold*P(j,:)')]', NgridL, 1);

        % Find the maximum value function
        [V(i,j), maxloc] = max(BellmanRHS(:));
        [nstar(i,j), jstar(i,j)] = ind2sub([NgridL, Ngridk], maxloc);
      end
    end

    niter = niter + 1
  end


end



