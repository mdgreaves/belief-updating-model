# Static Belief Updating Model

This repository contains MATLAB code for simulating, inverting, and validating a static self-belief updating model. The model captures how participants adjust their self-belief ratings in response to feedback, allowing for asymmetric learning when feedback is favourable versus unfavourable.



## Repository Structure

```
.
├── upd_data_bmc.m    % Empirical analysis: inversion + BPA + Bayesian model comparison (BMC)
├── upd_ffx_bpa.m     % Fixed-effects Bayesian parameter averaging (BPA)
├── upd_fig.m         % Helper for plotting example trials
├── upd_invert.m      % Wrapper for variational Bayesian analysis (VBA) model inversion
├── upd_log.m         % Exact moments for exp-transformed Gaussian parameters
├── upd_model.m       % Observation model
├── upd_sim.m         % Simulator for trial-level data
├── upd_toy.m         % Generate toy synthetic data for testing pipeline
├── upd_verify.m      % Simulation-recovery and model confusion analysis
└── README.md
```



## Analyses Included

### 1. Simulation–Recovery and Model Confusion

* Generates datasets across a grid of parameter values.
* Inverts each dataset and computes RMSE between true and recovered parameters.
* Evaluates model discriminability between full and reduced models.

### 2. Empirical Analysis

* Inverts the model on participant trial-by-trial data.
* Pools posteriors across subjects using fixed-effects BPA.
* Performs posterior inference on group-level parameters (one-tailed test for positivity).
* Reports means, standard deviations, and posterior probabilities.
* Compares the three models at the group level using random-effects Bayesian model comparison (VBA_groupBMC), printing expected model frequencies and exceedance probabilities.

### 3. Toy Data Generation

* Generates synthetic data in the required format to test the full pipeline.
* Useful for validating analysis steps without needing real data.


## Quick Start

1. Generate toy data for testing:

```matlab
upd_toy
```

2. Run the empirical analysis pipeline:

```matlab
upd_data_bmc
```

3. Verify parameter recovery and model discrimination:

```matlab
upd_verify
```


## Prerequisites

* MATLAB R2021a or later
* [VBA toolbox](https://mbb-team.github.io/VBA-toolbox/) installed and on the MATLAB path.

## Contact
Questions? Please feel free to reach out.

