capture program drop ebayes_shrink
program define ebayes_shrink, rclass
    version 16.0
    // varlist = FE estimates to shrink
    // WCOV() = prefix of per-observation V entries, e.g. S -> S11 S21 S22 ...
    // Options: method(joint|separate) shprefix() wprefix() tol()
    syntax varlist(min=1 numeric) [if] [in], WCOV(string) ///
        [ METHOD(string) SHPREFIX(string) WPREFIX(string) TOL(real 1e-10) replace]

	quietly{
    // basics
	local wcov = trim("`wcov'")
	if ("`wcov'"=="") {
		di as err "option wcov() is required and cannot be empty"
		exit 198
	}

    local K : word count `varlist'
    local method = lower("`method'")
    if ("`method'"=="") local method "joint"
    if !inlist("`method'","joint","separate") {
        di as err "method() must be joint or separate"
        exit 198
    }
    local shpx = cond("`shprefix'"=="","sh","`shprefix'")
    local wpx  = cond("`wprefix'" =="","w","`wprefix'")
    marksample touse

    // --- verify WCOV variables exist and build the needed list ---
    tempname needS
    local needS
    if ("`method'"=="joint") {
        forvalues i=1/`K' {
            forvalues j=1/`i' {
                local sname `wcov'`i'`j'
                capture confirm variable `sname'
                if _rc {
					if "`replace'"=="" {
						di as err "Missing WCOV variable: `sname'"
						exit 111
					}
				else drop `sname'
                }
                local needS `needS' `sname'
            }
        }
    }
    else {
        forvalues i=1/`K' {
            local sname `wcov'`i'`i'
            capture confirm variable `sname'
            if _rc {
				if "`replace'"=="" {
					di as err "Missing WCOV diagonal: `sname'"
					exit 111
				}
				else drop `sname'
            }
            local needS `needS' `sname'
        }
        // if S21, S31... exist, we simply ignore them
    }

    // --- enforce complete cases in-sample ---
    tempvar mask
    gen byte `mask' = `touse'
    foreach v of local varlist {
        quietly replace `mask' = 0 if missing(`v')
    }
    foreach s of local needS {
        quietly replace `mask' = 0 if missing(`s')
    }
    quietly count if `mask'
    if r(N)==0 {
        di as err "No observations after applying mask and non-missing checks."
        exit 2000
    }

    // --- create outputs ---
    local shouts
    local wouts
    foreach v of local varlist {
        local shv = "`shpx'_`v'"
        local wv  = "`wpx'_`v'"
        capture drop `shv'
        capture drop `wv'
        gen double `shv' = .
        gen double `wv'  = .
        local shouts `shouts' `shv'
        local wouts  `wouts'  `wv'
    }

    // --- temp matrices to collect returns ---
    tempname Tau TauDiag Share
    matrix `Tau'     = J(`K',`K',.)
    matrix `TauDiag' = J(1,`K',.)
    matrix `Share'   = J(1,`K',.)

    // --- call Mata core ---
    mata: ebayes_shrink_core(tokens("`varlist'"), "`wcov'", "`mask'", ///
         tokens("`shouts'"), tokens("`wouts'"), ///
         `= ("`method'"=="joint")', `tol', ///
         "`Tau'", "`TauDiag'", "`Share'")
	
	//Naive VCV
	corr `varlist', cov
	tempname naive
	mat `naive'=r(C)
	
    // --- label matrices and return ---
    matrix colnames `Tau'     = `varlist'
    matrix rownames `Tau'     = `varlist'
    matrix colnames `TauDiag' = `varlist'
	
	return scalar traceTau = trace(`Tau')
    return matrix tau=`Tau'
	return matrix tau_diag = `TauDiag'
	return matrix tau_naive= `naive'
    //return matrix share    = `Share'
    return scalar K        = `K'
    
	
	// Human-friendly method text
	local methlab = cond("`method'"=="joint","joint","separate")

	// Label each shrunken var and its weight
	foreach v of local varlist {
		local shv = "`shpx'_`v'"
		local wv  = "`wpx'_`v'"

		// Use the source variable's label if it exists
		local srcLab : variable label `v'
		local src    = cond("`srcLab'"=="", "`v'", "`srcLab'")

		// Shrunken value label
		capture confirm variable `shv'
		if (!_rc) {
			label var `shv' ///
				"EB shrink of [`src'] (method: `methlab', wcov: `wcov')"
		}

		// Weight label (diag of W = Tau*(Tau+V)^(-1))
		// Note: for method(joint), this is the diagonal element of the multivariate weight matrix.
		capture confirm variable `wv'
		if (!_rc) {
			label var `wv' ///
				"EB weight for [`src'] (diag W; method: `methlab', wcov: `wcov')"
		}
	}

	}
	
    // friendly display
    /*di as txt "Empirical Bayes shrinkage (" as res "`method'" as txt ")"
    mat list `Tau', noheader format(%9.6f)
    di as txt "diag(Tau):" _n as res rownames(`TauDiag') _n as res `TauDiag'
    di as txt "shares (diag(Tau)/sum diag):" _n as res `Share'
    di as txt "trace(Tau) = " as res %9.6f r(traceTau)*/
	

end


cap mata: mata drop psd_clip()
mata
real matrix psd_clip(real matrix M, real scalar tol)
{
    real rowvector wr
    real colvector w
    real matrix    Z, D
    real scalar    s, cutoff, i

    if (rows(M)==0 | cols(M)==0) return(M)
    M = 0.5*(M + M')
    symeigensystem(M, wr, Z)
    w = wr'

    s = max(abs(w))
    if (missing(s) | s<=0) s = 1
    cutoff = tol * s
    for (i=1; i<=rows(w); i++) {
        if (w[i] < 0 & abs(w[i]) <= cutoff) w[i] = 0
    }
    D = diag(w)
    M = Z * D * Z'
    M = 0.5*(M + M')
    return(M)
}
end

cap mata: mata drop ebayes_shrink_core()
mata
void ebayes_shrink_core(
    string rowvector feVars,      // K FE varnames
    string scalar      Sprefix,   // WCOV prefix; expects lower-tri Sij
    string scalar      touseVar,  // mask
    string rowvector   shVars,    // K output varnames (shrunken)
    string rowvector   wVars,     // K output varnames (diag weights)
    real   scalar      doJoint,   // 1 joint; 0 separate
    real   scalar      tol,       // PSD / ridge tol
    string scalar      TauName,   // where to store KxK Tau
    string scalar      TauDiagName, // where to store 1xK diag(Tau)
    string scalar      ShareName    // where to store 1xK shares
)
{
    real colvector    mask
    real matrix       Y
    real scalar       N, K, i, j, n, colidx, vijval
    real matrix       sumBB, sumV, Tau, V, M, Winv, W
    real colvector    b
    real rowvector    shares
    real colvector    diagTau
    real matrix       Sstack
    real scalar       ncols

    // Pull FE matrix & mask
    mask = st_data(., touseVar)
    Y    = st_data(selectindex(mask:==1), feVars)
    N    = rows(Y)
    K    = cols(Y)
    if (N==0 | K==0) {
        st_matrix(TauName,        J(K,K,.))
        st_matrix(TauDiagName,    J(1,K,.))
        st_matrix(ShareName,      J(1,K,.))
        return
    }

    // Stack all required Sij columns in one matrix [N x (K*(K+1)/2)]
    ncols  = K*(K+1)/2
    Sstack = J(N, ncols, .)
    colidx = 0
    for (i=1; i<=K; i++) {
        for (j=1; j<=i; j++) {
            string scalar vij
            vij     = Sprefix + strofreal(i) + strofreal(j)
            colidx  = colidx + 1
            Sstack[, colidx] = st_data(selectindex(mask:==1), vij)
        }
    }

    // Moments: Tau = mean(b b') - mean(V)
    sumBB = J(K,K,0)
    sumV  = J(K,K,0)
    for (n=1; n<=N; n++) {
        b = Y[n,.]'
        sumBB = sumBB + (b * b')

        V = J(K,K,0)
        colidx = 0
        for (i=1; i<=K; i++) {
            for (j=1; j<=i; j++) {
               
                colidx   = colidx + 1
                vijval   = Sstack[n, colidx]
                V[i,j]   = vijval
                V[j,i]   = vijval
            }
        }
        sumV = sumV + V
    }
    Tau = (sumBB - sumV) / N

    if (doJoint==0) {
        // diagonal only
        real matrix Tdiag
        Tdiag = J(K,K,0)
        for (i=1; i<=K; i++) {
            if (Tau[i,i] < 0) Tau[i,i] = 0
            Tdiag[i,i] = Tau[i,i]
        }
        Tau = Tdiag
    } else {
        // PSD clip tiny negatives
        //Tau = psd_clip(Tau, tol)
    }

    // store Tau and diag/share
    diagTau = diagonal(Tau)
    shares  = (sum(diagTau)==0 ? J(1,K,.) : (diagTau' :/ sum(diagTau))')
    st_matrix(TauName,      Tau)
    st_matrix(TauDiagName,  diagTau')
    st_matrix(ShareName,    shares)

    // Posterior means and "weights"
    for (n=1; n<=N; n++) {
        b = Y[n,.]'

        V = J(K,K,0)
        colidx = 0
        for (i=1; i<=K; i++) {
            for (j=1; j<=i; j++) {
            
                colidx   = colidx + 1
                vijval   = Sstack[n, colidx]
                V[i,j]   = vijval
                V[j,i]   = vijval
            }
        }

        // W = Tau * (Tau + V)^(-1)
        M    = Tau + V
        M    = 0.5*(M + M')
        Winv = invsym(M + tol*I(K))
        W    = Tau * Winv

        real colvector sh
        sh = W * b

        // write outputs for this row
        for (i=1; i<=K; i++) {
            st_store(selectindex(mask:==1)[n], shVars[i], sh[i])
            st_store(selectindex(mask:==1)[n],  wVars[i], W[i,i])
        }
    }
}
end
