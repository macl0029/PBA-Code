/**************************************************************************
 Project: Exposure Misclassification adjustment -- SUMMARY LEVEL
 Stata command version — revised 2026-09-05
 This _2 version repairs macro quoting and retains systematic estimates
 independently of later total-error exclusions.
 Parentheses around squared correlation macros preserve negative-rho copulas.

 Calling pba_summary_exp_rr:
   - runs the simulations
   - prints the final results table
   - draws the forest plot unless nograph is specified
   - returns results in r()
**************************************************************************/

set trace off

capture program drop pba_summary_exp_rr
program define pba_summary_exp_rr, rclass
    version 16.0

    syntax, ///
        A(integer) B(integer) C(integer) D(integer) ///
        SENCASEA(real) SENCASEB(real) ///
        SENCTRLA(real) SENCTRLB(real) ///
        SPECASEA(real) SPECASEB(real) ///
        SPECTRLA(real) SPECTRLB(real) ///
        TYPE(string) SIMS(string) ///
        [SEED(integer -1) RHOSEN(real 0.80) RHOSPE(real 0.80) NOGRAPH ///
         DRAWSSAVING(string)]

    /**********************************************************************
    Basic checks
    **********************************************************************/

    local type = lower(strtrim("`type'"))

    if !inlist("`type'", "diff", "nondiff") {
        di as error "type() must be either diff or nondiff"
        exit 198
    }

    capture local sims_n = floor(`sims')
    if _rc {
        di as error "sims() must be a positive integer or numeric expression"
        exit 198
    }

    if missing(`sims_n') | `sims_n' <= 0 | `sims_n' != (`sims') {
        di as error "sims() must evaluate to a positive integer"
        exit 198
    }

    local sims = `sims_n'

    foreach x in a b c d sims {
        if ``x'' <= 0 {
            di as error "`x' must be positive"
            exit 198
        }
    }

    foreach x in sencasea sencaseb senctrla senctrlb ///
               specasea specaseb spectrla spectrlb {
        if ``x'' <= 0 {
            di as error "`x' must be positive"
            exit 198
        }
    }

    if `seed' < -1 {
        di as error "seed() must be nonnegative"
        exit 198
    }

    if abs(`rhosen') >= 1 {
        di as error "rhosen() must be strictly between -1 and 1"
        exit 198
    }

    if abs(`rhospe') >= 1 {
        di as error "rhospe() must be strictly between -1 and 1"
        exit 198
    }

    if `seed' != -1 {
        set seed `seed'
    }

    /**********************************************************************
    Constants
    **********************************************************************/

    local n_case = `a' + `b'
    local n_ctrl = `c' + `d'

    local observed_rr = (`a' / (`a' + `c')) / (`b' / (`b' + `d'))
    local se_re_only  = sqrt(1 / `a' + 1 / `b' - 1 / (`a' + `c') - 1 / (`b' + `d'))

    preserve

        quietly {
            clear
            set obs `sims'
            gen long draw = _n

            /**************************************************************
            Draw correlated sensitivities using a Gaussian copula
            **************************************************************/

            gen double z1 = rnormal()
            gen double z2 = `rhosen' * z1 + sqrt(1 - (`rhosen')^2) * rnormal()

            gen double u1 = normal(z1)
            gen double u2 = normal(z2)

            gen double sen_case = invibeta(`sencasea', `sencaseb', u1)
            gen double sen_ctrl = invibeta(`senctrla', `senctrlb', u2)

            if "`type'" == "nondiff" {
                replace sen_ctrl = sen_case
            }

            drop z1 z2 u1 u2

            /**************************************************************
            Draw correlated specificities using a Gaussian copula
            **************************************************************/

            gen double z1 = rnormal()
            gen double z2 = `rhospe' * z1 + sqrt(1 - (`rhospe')^2) * rnormal()

            gen double u1 = normal(z1)
            gen double u2 = normal(z2)

            gen double spe_case = invibeta(`specasea', `specaseb', u1)
            gen double spe_ctrl = invibeta(`spectrla', `spectrlb', u2)

            if "`type'" == "nondiff" {
                replace spe_ctrl = spe_case
            }

            drop z1 z2 u1 u2

            /**************************************************************
            Bias-adjusted cell frequencies
            **************************************************************/

            gen double ac = (`a' - `n_case' * (1 - spe_case)) / ///
                (sen_case - (1 - spe_case))

            gen double bc = `n_case' - ac

            gen double cc = (`c' - `n_ctrl' * (1 - spe_ctrl)) / ///
                (sen_ctrl - (1 - spe_ctrl))

            gen double dc = `n_ctrl' - cc

            /**************************************************************
            Exposure prevalence among cases and controls,
            accounting for sampling error
            **************************************************************/

            gen double prevE_cases    = rbeta(ac, bc)
            gen double prevE_controls = rbeta(cc, dc)

            /**************************************************************
            Predictive values
            **************************************************************/

            gen double ppv_case = ///
                (sen_case * prevE_cases) / ///
                ((sen_case * prevE_cases) + (1 - spe_case) * (1 - prevE_cases))

            gen double ppv_control = ///
                (sen_ctrl * prevE_controls) / ///
                ((sen_ctrl * prevE_controls) + (1 - spe_ctrl) * (1 - prevE_controls))

            gen double npv_case = ///
                (spe_case * (1 - prevE_cases)) / ///
                ((1 - sen_case) * prevE_cases + spe_case * (1 - prevE_cases))

            gen double npv_control = ///
                (spe_ctrl * (1 - prevE_controls)) / ///
                ((1 - sen_ctrl) * prevE_controls + spe_ctrl * (1 - prevE_controls))

            /**************************************************************
            Bias-adjusted 2 x 2 tables
            **************************************************************/

            gen double ab = rbinomial(`a', ppv_case) + ///
                            rbinomial(`b', 1 - npv_case)

            gen double bb = `n_case' - ab

            gen double cb = rbinomial(`c', ppv_control) + ///
                            rbinomial(`d', 1 - npv_control)

            gen double db = `n_ctrl' - cb

            /**************************************************************
            Risk ratios
            **************************************************************/

            gen double rr_bb = (ab / (ab + cb)) / (bb / (bb + db))

            gen double se_bb = sqrt( ///
                1 / ab + 1 / bb - 1 / (ab + cb) - 1 / (bb + db) ///
            )

            gen double rr_total = exp(log(rr_bb) - rnormal() * se_bb)

            gen double rr_re = exp(log(`observed_rr') - rnormal() * `se_re_only')

            gen double rr_syst = (ac / (ac + cc)) / (bc / (bc + dc))

            /**************************************************************
            Valid draws
            **************************************************************/

            * Diagnostics use separate, sequential rejection stages.
            gen byte corrected_ok = ///
                !missing(ac, bc, cc, dc) & ac > 0 & bc > 0 & cc > 0 & dc > 0
            gen byte predictive_ok = corrected_ok & ///
                !missing(ppv_case, ppv_control, npv_case, npv_control) & ///
                inrange(ppv_case, 0, 1) & inrange(ppv_control, 0, 1) & ///
                inrange(npv_case, 0, 1) & inrange(npv_control, 0, 1)
            gen byte reclassified_ok = predictive_ok & ///
                !missing(ab, bb, cb, db) & ab > 0 & bb > 0 & cb > 0 & db > 0

            * Systematic estimates depend only on the corrected table.
            gen byte valid_systematic = corrected_ok & ///
                !missing(rr_syst) & rr_syst > 0

            gen byte valid = ///
                !missing(ac, bc, cc, dc, ab, bb, cb, db, rr_total, rr_syst) & ///
                ac > 0 & bc > 0 & cc > 0 & dc > 0 & ///
                ab > 0 & bb > 0 & cb > 0 & db > 0

            replace rr_total = . if !valid
            replace rr_syst  = . if !valid_systematic

            if "`drawssaving'" != "" {
                export delimited draw rr_re rr_syst rr_total corrected_ok ///
                    predictive_ok reclassified_ok valid_systematic valid ///
                    using "`drawssaving'", replace
            }
        }

        quietly count if valid_systematic
        local n_valid_systematic = r(N)
        local n_impossible_systematic = `sims' - `n_valid_systematic'
        quietly count if !corrected_ok
        local rejected_corrected = r(N)
        quietly count if corrected_ok & !predictive_ok
        local rejected_predictive = r(N)
        quietly count if predictive_ok & !reclassified_ok
        local rejected_zero_cells = r(N)
        quietly count if reclassified_ok & !valid
        local rejected_effect = r(N)
        assert valid <= reclassified_ok

        quietly count if valid
        local n_valid = r(N)
        local n_impossible = `sims' - `n_valid'

        if `n_valid_systematic' == 0 {
            restore
            di as error "No valid systematic-error simulations."
            exit 498
        }

        if `n_valid' == 0 {
            di as text "No total-error draws remained; valid systematic estimates were retained."
        }

        /******************************************************************
        Percentiles
        ******************************************************************/

        quietly _pctile rr_re, p(2.5 50 97.5)
        local re_l = r(r1)
        local re_m = r(r2)
        local re_u = r(r3)
        local re_w = `re_u' / `re_l'

        quietly _pctile rr_syst if valid_systematic, p(2.5 50 97.5)
        local syst_l = r(r1)
        local syst_m = r(r2)
        local syst_u = r(r3)
        local syst_w = `syst_u' / `syst_l'

        local total_l = .
        local total_m = .
        local total_u = .
        local total_w = .
        if `n_valid' > 0 {
            quietly _pctile rr_total if valid, p(2.5 50 97.5)
            local total_l = r(r1)
            local total_m = r(r2)
            local total_u = r(r3)
            local total_w = `total_u' / `total_l'
        }

        /******************************************************************
        Store results in a matrix
        ******************************************************************/

        tempname results

        matrix `results' = ///
            (`re_m',    `re_l',    `re_u',    `re_w',    0 \ ///
             `syst_m',  `syst_l',  `syst_u',  `syst_w',  `n_impossible_systematic' \ ///
             `total_m', `total_l', `total_u', `total_w', `n_impossible')

        matrix colnames `results' = Median Lower Upper Width N_Impossible
        matrix rownames `results' = Random_error Systematic_error Total_error

        /******************************************************************
        Print final results
        ******************************************************************/

        di as text _n "Probabilistic bias analysis: exposure misclassification"
        di as text "Type:               " as result "`type'"
        di as text "Simulations:        " as result %12.0fc `sims'
        di as text "Systematic draws:   " as result %12.0fc `n_valid_systematic'
        di as text "Total-error draws:  " as result %12.0fc `n_valid'
        di as text "Rejected systematic:" as result %12.0fc `n_impossible_systematic'
        di as text "Rejected total:     " as result %12.0fc `n_impossible'
        di as text "Observed RR:        " as result %9.4f `observed_rr'

        clear
        set obs 3

        gen str20 analysis = ""
        replace analysis = "Random error"     in 1
        replace analysis = "Systematic error" in 2
        replace analysis = "Total error"      in 3

        gen double median = .
        gen double lower  = .
        gen double upper  = .
        gen double width  = .
        gen long n_impossible = .

        replace median = `re_m'    in 1
        replace median = `syst_m'  in 2
        replace median = `total_m' in 3

        replace lower = `re_l'    in 1
        replace lower = `syst_l'  in 2
        replace lower = `total_l' in 3

        replace upper = `re_u'    in 1
        replace upper = `syst_u'  in 2
        replace upper = `total_u' in 3

        replace width = `re_w'    in 1
        replace width = `syst_w'  in 2
        replace width = `total_w' in 3

        replace n_impossible = 0              in 1
        replace n_impossible = `n_impossible_systematic' in 2
        replace n_impossible = `n_impossible' in 3

        format median lower upper width %9.4f
        format n_impossible %12.0fc

        di as text _n "Final results"
        list analysis median lower upper width n_impossible, noobs clean abbreviate(20)

        /******************************************************************
        Graph
        ******************************************************************/

        if "`nograph'" == "" {
            gen byte y = .
            replace y = 2 if analysis == "Random error"
            replace y = 1 if analysis == "Systematic error"
            replace y = 3 if analysis == "Total error"

            quietly summarize lower, meanonly
            local x_min = r(min) * 0.9

            quietly summarize upper, meanonly
            local x_max = r(max) * 1.1

            twoway ///
                (rcap lower upper y, horizontal) ///
                (scatter y median, msymbol(O) msize(medium)), ///
                xscale(log range(`x_min' `x_max')) ///
                xlabel(, format(%9.2f)) ///
                ylabel( ///
                    3 "Total error" ///
                    2 "Random error only" ///
                    1 "Systematic error only", ///
                    angle(0) noticks ///
                ) ///
                yscale(range(0.5 3.5)) ///
                title("Summary of Bias Adjustment Methods") ///
                xtitle("Risk ratio estimate") ///
                ytitle("") ///
                legend(off)
        }

    restore

    /**********************************************************************
    Return results in r()
    **********************************************************************/

    return matrix results = `results'

    return scalar observed_rr  = `observed_rr'
    return scalar sims         = `sims'
    return scalar valid        = `n_valid'
    return scalar impossible   = `n_impossible'
    return scalar valid_systematic = `n_valid_systematic'
    return scalar impossible_systematic = `n_impossible_systematic'
    return scalar rejected_corrected = `rejected_corrected'
    return scalar rejected_predictive = `rejected_predictive'
    return scalar rejected_zero_cells = `rejected_zero_cells'
    return scalar rejected_effect = `rejected_effect'

    return scalar re_median    = `re_m'
    return scalar re_lower     = `re_l'
    return scalar re_upper     = `re_u'
    return scalar re_width     = `re_w'

    return scalar syst_median  = `syst_m'
    return scalar syst_lower   = `syst_l'
    return scalar syst_upper   = `syst_u'
    return scalar syst_width   = `syst_w'

    return scalar total_median = `total_m'
    return scalar total_lower  = `total_l'
    return scalar total_upper  = `total_u'
    return scalar total_width  = `total_w'

end


/**************************************************************************
 Example call
 No seed is set here. Calls within a running Stata session continue its RNG.
 A fresh Stata process starts at Stata's default RNG state; the separate
 validation batch advances it using operating-system entropy before any runs.
 drawssaving("path.csv") optionally saves simulation draws and stage flags.
**************************************************************************/

pba_summary_exp_rr, ///
    a(215) b(1449) c(668) d(4296) ///
    sencasea(50.6) sencaseb(14.3) ///
    senctrla(50.6) senctrlb(14.3) ///
    specasea(70) specaseb(1) ///
    spectrla(70) spectrlb(1) ///
    type(nondiff) ///
    sims(100000)
