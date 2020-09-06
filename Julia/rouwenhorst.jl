  function rouwenhorst(μ,ϕ,σ,N)
  # Inputs: μ, ϕ, σ of AR(1) process y(t) = ϕ y(t-1) + e(t)
  # μ is the unconditional mean of the process
  # ϕ is the AR coefficient
  # σ is the standard deviation of e(t)
  # # Summary: The function discretizes the AR(1) process into an N equally
  # spaced state Markov chain with transition probabilities given by the
  # Rouwenhorst (1995) method. States are in the set [μ - nu , μ + nu]

  # Outputs: Lambda, P
  # Lambda is a Nx1 vector of equally spaced states centered around mu
  # P is the Markov transition matrix

    nu = sqrt( ((N-1)/(1-ϕ^2)) )*σ

    Lambda = linspace(μ-nu,μ+nu,N)'

    p = (1+ϕ)/2
    q = p

    if N > 1
      P1 = [ [p 1-p]; [1-q q] ]
    else
      P1 = [1]
    end

    if N == 2
        P = P1
    else
      for n = 3:N
        zcol = zeros(n-1,1)
        zrow = zcol'

        A = [ [P1 zcol];  [zrow 0] ]
        B = [ [zcol P1];  [0 zrow] ]
        C = [ [zrow 0]; [P1 zcol] ]
        D = [ [0 zrow]; [zcol P1] ]

        P1 = p*A + (1-p)*B + (1-q)*C + q*D
        P1[2:end-1,:] = P1[2:end-1,:]/2
      end
      P = P1
    end
    return collect(Lambda), P
  end

