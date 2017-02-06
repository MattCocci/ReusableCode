#=using PyPlot=#

########################################################################
## Define Parameters and Functions #####################################
########################################################################

  # Production related
  # f     = A*k^α
  # df    = df/dk
  # dfinv = functional inverse of df
  A        = 1
  α        = 0.3
  f(k)     = A*k.^α
  df(k)    = A*α*k.^(α-1)
  dfinv(y) = (y/(A*α)).^(1/(α-1))

  # Preferences/utility related
  # u     = (c^(1-σ))/(1-σ)
  # du    = du/dc
  # duinv = functional inverse of du
  β        = 0.96
  σ        = 0.999
  u(c)     = (c.^(1-σ))/(1-σ)
  du(c)    = complex(c).^(-σ)
  duinv(v) = v.^(-1/σ)

  # Depreciation rate
  δ = 0.08

  # Start capital path at lambda*kss
  λ = 0.1

########################################################################
## Shooting Program Function ###########################################
########################################################################

  # Define the shooting function
  function shootingProgram(f, df, dfinv, u, du, duinv, β, δ, λ, T)

    # crit*kss (steady state) is the convergence criterion
    crit = .00000015

    # (Consumption) = (output) + (nondepreciated k) - (new k)
    c(kold, knew) = f(kold) + (1-δ)*kold - knew

    # Steady state capital
    kss = dfinv( (1/β)-(1-δ) )

    # Initial capital level
    k0 = λ*kss

    # Solution for path ksol
    ksol    = Array(Float64, T)
    ksol[1] = k0

    # At each iteration, take k_{t-1} as given, want to find k_t
    for t = 2:T

      # Store k_{t-1}
      k_t1 = ksol[t-1]

      # Initialize bounds for k_t and k_t (the thing we want to guess)
      kmin = k_t1
      kmax = kss
      k_t  = 0

      # Keep trying different guesses for k_t until you find a good
      # guess. A "good guess" means the interval [kmin,kmax]
      # is sufficiently small enough to say that we've found k_t and
      # should exit outer while loop
      while abs(kmax-kmin) > crit*kss

        # A guess for k_t, the midpoint of the interval [kmin,kmax]
        k_t = (kmin+kmax)/2

        # Initialize k_{s-1} and k_{s-2}
        k_s1 = k_t
        k_s2 = k_t1

        # Iterate on that initial guess to find k_s for s >= t+1
        keepIterating = true
        while keepIterating

          println(
              f(k_s1) + (1-δ)k_s1
                  - duinv(
                      du( c(k_s2,k_s1) ) /
                      ( β*(df(k_s1)+(1-δ)) )
                  )
          )

          #=if c(k_s2, k_s1) >= 0=#
            # Find k_s given k_{s-1} and k_{s-2}
            k_s = f(k_s1) + (1-δ)k_s1
                  - duinv(
                      du( c(k_s2,k_s1) ) /
                      ( β*(df(k_s1)+(1-δ)) )
                  )
          #=else=#
            #=k_s = kss + 1=#
          #=end=#

          # Check to see if k path is going to zero (k_t was too low)
          if k_s <= k_s1
            kmin = k_t
            keepIterating = false    # We'll exit inner while loop

          # Check to see if kguess is going beyond kss (k_t was too high)
          elseif k_s > kss
            kmax = k_t
            keepIterating = false    # We'll exit inner while loop
          end

          # Update k_s1 and k_s2
          k_s2 = k_s1
          k_s1 = k_s
        end
      end

      # If you exited the loop, you have your k_t
      ksol[t] = k_t
    end

    # Construct some stuff to return
    css   = c(kss,kss)
    ysol  = f(ksol)
    csol  = [c(ksol[1:(T-1)], ksol[2:T]); NaN]
    srate = 1 - csol./ysol

    return kss, css, ysol, ksol, csol, srate

  end


########################################################################
## Run the Program and Plot ############################################
########################################################################

kss, css, ysol, ksol, csol, srate =
  shootingProgram(f, df, dfinv, u, du, duinv, β, δ, λ, 200)



