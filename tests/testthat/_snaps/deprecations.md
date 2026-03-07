# synthdid_estimate() is soft-deprecated

    Code
      tau.hat <- synthdid_estimate(setup$Y, setup$N0, setup$T0)
    Condition
      Warning:
      `synthdid_estimate()` was deprecated in synthdid 2.0.0.
      i Please use `synthdid()` instead.

# sc_estimate() is soft-deprecated

    Code
      tau.hat <- sc_estimate(setup$Y, setup$N0, setup$T0)
    Condition
      Warning:
      `sc_estimate()` was deprecated in synthdid 2.0.0.
      i Please use `synthdid()` instead.

# did_estimate() is soft-deprecated

    Code
      tau.hat <- did_estimate(setup$Y, setup$N0, setup$T0)
    Condition
      Warning:
      `did_estimate()` was deprecated in synthdid 2.0.0.
      i Please use `synthdid()` instead.

# synthdid_se() is soft-deprecated

    Code
      se <- synthdid_se(tau.hat, method = "jackknife")
    Condition
      Warning:
      `synthdid_se()` was deprecated in synthdid 2.0.0.
      i Please use `vcov()` instead.

# panel.matrices() is soft-deprecated

    Code
      setup <- panel.matrices(california_prop99)
    Condition
      Warning:
      `panel.matrices()` was deprecated in synthdid 2.0.0.
      i The formula interface `synthdid()` handles panel conversion internally.

# synthdid_effect_curve() is soft-deprecated

    Code
      curve <- synthdid_effect_curve(tau.hat)
    Condition
      Warning:
      `synthdid_effect_curve()` was deprecated in synthdid 2.0.0.
      i Please use `predict()` instead.

# synthdid_converged() is soft-deprecated

    Code
      conv <- synthdid_converged(tau.hat)
    Condition
      Warning:
      `synthdid_converged()` was deprecated in synthdid 2.0.0.
      i Convergence status is shown in the `print()` output.

# synthdid_convergence_info() is soft-deprecated

    Code
      info <- synthdid_convergence_info(tau.hat)
    Condition
      Warning:
      `synthdid_convergence_info()` was deprecated in synthdid 2.0.0.
      i Convergence diagnostics are shown in the `print()` and `summary()` output.

# synthdid_controls() is soft-deprecated

    Code
      ctrls <- synthdid_controls(tau.hat)
    Condition
      Warning:
      `synthdid_controls()` was deprecated in synthdid 2.0.0.
      i Use the formula interface with `synthdid()` instead.

# synthdid_memory_estimate() is soft-deprecated

    Code
      mem <- synthdid_memory_estimate(N = 39, T = 31)
    Condition
      Warning:
      `synthdid_memory_estimate()` was deprecated in synthdid 2.0.0.
      i Use the formula interface with `synthdid()` instead.

# timesteps() is soft-deprecated

    Code
      ts <- timesteps(attr(tau.hat, "setup")$Y)
    Condition
      Warning:
      `timesteps()` was deprecated in synthdid 2.0.0.
      i Use the formula interface with `synthdid()` instead.

# synthdid_plot() is soft-deprecated

    Code
      p <- synthdid_plot(tau.hat)
    Condition
      Warning:
      `synthdid_plot()` was deprecated in synthdid 2.0.0.
      i Please use `plot()` instead.
      Warning:
      `timesteps()` was deprecated in synthdid 2.0.0.
      i Use the formula interface with `synthdid()` instead.

# synthdid_placebo() is soft-deprecated

    Code
      placebo <- synthdid_placebo(tau.hat)
    Condition
      Warning:
      `synthdid_placebo()` was deprecated in synthdid 2.0.0.
      i Use the formula interface with `synthdid()` instead.

