{smcl}
{* *! version 1.0.0  04sep2025}{...}
{title:Title}

{phang}
{bf:ebayes_shrink} — Empirical Bayes shrinkage of one or more fixed-effect estimates using per-observation VCVs


{title:Syntax}

{p 8 15 2}
{cmd:ebayes_shrink} {it:varlist} {ifin}{cmd:,} {opt wcov(prefix)} 
[{opt method(joint|separate)} {opt shprefix(name)} {opt wprefix(name)} {opt tol(#)}]

{synoptset 24 tabbed}{...}
{synopthdr:Options}
{synoptline}
{synopt:{opt wcov(prefix)}}prefix for per-observation lower-triangle VCV entries (required){p_end}
{synopt:{opt method(joint|separate)}}estimate a joint K×K signal covariance (default {bf:joint}) or diagonal-only ({bf:separate}){p_end}
{synopt:{opt shprefix(name)}}prefix for shrunken variables; default {bf:sh}{p_end}
{synopt:{opt wprefix(name)}}prefix for shrink "weight" variables; default {bf:w}{p_end}
{synopt:{opt tol(#)}}PSD safeguard tolerance for covariance clipping; default {bf:1e-10}{p_end}
{synoptline}

{title:Description}

{pstd}
{cmd:ebayes_shrink} implements a general Empirical Bayes shrinkage step for any set of
per-observation fixed-effect estimates. You supply:

{p 10 14 2}
(1) a list of K variables in data {it:fevarlist} containing the estimated effects {it:b_i} per
observation (e.g. worker and job effects: {bf:alpha theta})

{p 10 14 2}
(2) a {it:prefix} giving the per-observation estimated VCV entries {it:V_i} for those effects
as variables named {it:prefix}{bf:ij} for the lower triangle (e.g. {bf:S11 S21 S22} when K=2 and wcov(S).).

{pstd}
The command estimates the (population) signal covariance matrix {it:tau} by method-of-moments
and computes the per-observation posterior means {it: b_i^{EB} = tau*(tau+V_i)^{-1} * b_i }.
It creates shrunken variables {bf:sh_}{it:fevar} (prefix configurable) and per-dimension
shrink "weights" {bf:w_}{it:fevar} (the diagonal of {it:w_i = tau*(tau+V_i)^{-1}}).
It also returns the estimates of {it:tau} in {cmd:e()}.

{pstd}
When {opt method(joint)} (default), {it:tau} is a full K×K covariance estimated jointly.
When {opt method(separate)}, {it:tau} is restricted to be diagonal; off-diagonal {it:V_i}
entries (e.g., {bf:S21}) may be present but are ignored.


{title:Inputs and naming of VCV variables}

{pstd}
Let {it:fevarlist} contain K variables in the order {it:b1 b2 … bK}.
For each observation, provide the lower-triangle of its estimated sampling-VCV {it:V_i}
using variables named:

{p2colset 9 30 30 2}{...}
{p2col :{it:prefix}11}{it:Var}(b1){p_end}
{p2col :{it:prefix}21}{it:Cov}(b2,b1){p_end}
{p2col :{it:prefix}22}{it:Var}(b2){p_end}
{p2col :…}… up to {it:prefix}{bf:KK}{p_end}

{pstd}
For {opt method(joint)}, {bf:all} lower-triangle entries {it:prefix}{bf:ij} with i≥j must exist.
For {opt method(separate)}, only the K diagonal entries {it:prefix}{bf:ii} are required;
if off-diagonals are present, they are ignored (no error).


{title:Options}

{phang}
{opt wcov(prefix)} (required) specifies the prefix used to find the VCV-variable names
{it:prefix}{bf:ij}. The command checks that all required variables exist and are not missing in-sample.

{phang}
{opt method(joint|separate)} chooses whether to estimate the signal covariance {it:tau} jointly
(K×K, default {bf:joint}) or separately by dimension (diagonal only, {bf:separate}).
The joint method uses all {it:V_i} entries; separate uses only diagonals.

{phang}
{opt shprefix(name)} sets the prefix for the shrunken variables. The default is {bf:sh},
producing e.g. {bf:sh_alpha} for input {bf:alpha}.

{phang}
{opt wprefix(name)} sets the prefix for the weight variables. The default is {bf:w},
producing e.g. {bf:w_alpha}. Weights are the diagonal elements of {it:W_i = tau*(tau+V_i)^{-1}}.
Note that shrunken values use the full matrix {it:w_i}, not only the diagonal.

{phang}
{opt tol(#)} sets the tolerance for PSD clipping of covariance matrices. The estimated {it:Tau}
(and any covariance formed internally) is symmetrized and small negative eigenvalues (up to
{it:tol}×{it:max|eigen|}) are floored at zero. Default {bf:1e-10}. Increase slightly if you
encounter tiny negative variances due to numerical noise.


{title:Remarks}

{pstd}
{it:Identification/normalization.} {cmd:ebayes_shrink} treats the input {it:b_i} and their
{it:V_i} as the estimand and estimator variance in the user's chosen parameterization.
If your fixed effects are normalized (mean-zero or baseline-drop), that choice should be
consistent across all observations and already reflected in {it:b_i} and {it:V_i}.

{pstd}
{it:MoM for Tau.} the signal covariance is estimated by
{it:tau} = {bf:mean}(b_i b_i') − {bf:mean}(V_i),
computed over the analysis sample ({cmd:if}/ {cmd:in} apply; observations with any missing
required entry are dropped). For {opt method(separate)}, off-diagonal elements of {it:tau}
are set to zero after ensuring nonnegativity of the diagonal.

{pstd}
{it:Posterior means.} For each observation i, with observed effects vector {it:b_i} (K×1) and
sampling VCV {it:V_i} (K×K), the Empirical Bayes posterior mean is
{it: b_i^{EB} = w_i b_i } with {it: w_i = tau (tau + V_i)^{-1} }.
The reported weights {bf:w_*} are the diagonal entries of {it:w_i}; shrunken estimates
use the full matrix {it:q_i}.

{pstd}
{it:Robustness.} If the MoM {it:tau} is mildly indefinite because of small samples,
overlap patterns, or numerical rounding, the command applies a PSD safeguard (like
{help reghdfe}) controlled by {opt tol()}.


{title:Examples}

{pstd}
Two-way FE example (worker {bf:alpha}, job {bf:theta}) with per-observation V entries {bf:S11 S21 S22}:

{cmd}
. twoway_fe_fill w, pid(pid) jid(jid) bread(exact) center(weighted)
. ebayes_shrink alpha theta, wcovvariables(S)

. ebayes_shrink alpha theta, wcovvariables(S) method(separate)

. ebayes_shrink alpha theta, wcovvariables(S) shprefix(eb) wprefix(ww) tol(1e-9)
{txt}

{pstd}
Using {cmd:if} to restrict the sample:

{cmd}
. ebayes_shrink alpha theta if inrange(year,2022,2025), wcovvariables(S)
{txt}


{title:Stored results}

{pstd}
{cmd:ebayes_shrink} stores the following in {cmd:r()}:

{synoptset 20 tabbed}{...}
{synopt:{cmd:r(Tau)}}K×K matrix: estimated signal covariance {it:Tau}{p_end}
{synopt:{cmd:r(tau_diag)}}1×K rowvector: diagonal of {it:Tau}{p_end}
{synopt:{cmd:r(share)}}1×K rowvector: {it:diag(Tau)}/sum{it:diag(Tau)} (component shares){p_end}
{synopt:{cmd:r(K)}}scalar: number of FE dimensions (K){p_end}
{synopt:{cmd:r(traceTau)}}scalar: trace({it:Tau}){p_end}


{title:Methods and formulas}

{pstd}
Let {it:b_i} be the vector of length K of estimated effects for observation i, and {it:V_i} its
estimated sampling VCV (provided via {opt wcov()}). Assume the hierarchical model
{it:b_i = θ_i + ε_i}, with {it:θ_i} drawn from a population with covariance {it:tau} and
{it:ε_i | data} having covariance {it:V_i}. The method-of-moments estimator is
{it:tau} = E[{it:b_i b_i'}] − E[{it:V_i}], implemented as sample means.
Posterior means are {it:b_i^{EB} = tau (tau + V_i)^{-1} b_i}. When {opt method(separate)},
{it:tau} is restricted to be diagonal.


{title:Diagnostics and troubleshooting}

{pstd}
If you receive "missing WCOV variable …", verify that all required {it:prefix}{bf:ij}
variables exist for the requested method and that the ordering of {it:fevarlist} matches
the {it:prefix}{bf:ij} indexing. Observations with any missing required inputs are
excluded silently.

{pstd}
Large or zero weights can indicate very imprecise {it:b_i} (large {it:V_i}) or a near-zero
estimated {it:tau}. Consider restricting the sample, revisiting normalization, or using
{opt method(separate)} for stability.


{title:Author}

{pstd}
Martin Eckhoff Andresen, Department of Economics, University of Oslo.
Developed for {browse "https://arxiv.org/pdf/2606.02503":Pay Beliefs and the Amenity-Pay Tradeoff}, joint with Manudeep
Bhuller and Alfred Løvgren. 


{title:Also see}

{psee}
{help reghdfe}, {help mixed}, {help boottest}
