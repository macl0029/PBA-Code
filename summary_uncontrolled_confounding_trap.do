/**************************************************************************
 Version 6: RR and OR options, modeled after summary_confounding.r.
 Original input files are preserved. Load this file with do; no example executes.

 Data layout: cases a/b and noncases c/d, exposed/unexposed columns.
 Prevalence and confounder-disease RR priors are independent trapezoids.
 effect(RR) is the default; effect(OR) selects the Mantel-Haenszel odds ratio.
 OR uses full stratum totals and the Robins-Breslow-Greenland log variance.
 Systematic draws are retained independently of later random allocation zeros.
 Invalid expected tables, allocation probabilities and systematic effects are
 discarded at their own stage. Zero-cell allocations affect total error only.
 Random-only output uses all requested simulations, matching the R PBA.
 All three components return separate discarded counts.

 This program does not set a seed or reset the current RNG stream.
 drawssaving(path) optionally saves full-precision CSV diagnostics and refuses
 to overwrite an existing file. No user data are replaced (preserve/restore).
 Numeric arguments have their previous option names. Use, for example:
   pba_conf_summary_mh, a(48) b(378) c(52) d(522) ///
     ponemin(.36) ponemoda(.384) ponemodb(.416) ponemax(.44) ///
     pzeromin(.56) pzeromoda(.584) pzeromodb(.616) pzeromax(.64) ///
     rrmin(1.12) rrmoda(1.138) rrmodb(1.162) rrmax(1.18) ///
     sims(100000) effect(OR) nograph

 r(results): rows Random_error, Systematic_error, Total_error;
 columns Median, Lower, Upper, Width, N_Impossible.
 r(observed_effect), r(crude_se), r(effect) are measure-neutral;
 r(observed_rr) is retained for RR, r(observed_or) is returned for OR.
**************************************************************************/

* Internal vectorized calculator: the same equations are used by both stages.
capture program drop _pba_conf_mh_calc
program define _pba_conf_mh_calc
    version 16.0
    syntax varlist(min=8 max=8 numeric), EFFECT(string) GENerate(name) ///
        [VARiance(name) NOVARiance]
    local effect = upper(strtrim("`effect'"))
    if !inlist("`effect'", "RR", "OR") exit 198
    if "`novariance'" == "" & "`variance'" == "" exit 198
    tokenize `varlist'
    tempvar m0 n0 t0 m1 n1 t1 num den v0 v1 vc1 vc2 vc3
    quietly {
        gen double `m0' = `1' + `3'
        gen double `n0' = `2' + `4'
        gen double `t0' = `m0' + `n0'
        gen double `m1' = `5' + `7'
        gen double `n1' = `6' + `8'
        gen double `t1' = `m1' + `n1'
        if "`effect'" == "RR" {
            gen double `num' = `1'*`n0'/`t0' + `5'*`n1'/`t1'
            gen double `den' = `2'*`m0'/`t0' + `6'*`m1'/`t1'
            gen double `generate' = `num'/`den'
            if "`novariance'" == "" {
                gen double `v0' = (`1'+`2')*`m0'*`n0'/`t0'^2 - `1'*`2'/`t0'
                gen double `v1' = (`5'+`6')*`m1'*`n1'/`t1'^2 - `5'*`6'/`t1'
                gen double `variance' = (`v0'+`v1')/(`num'*`den')
            }
        }
        else {
            gen double `num' = `1'*`4'/`t0' + `5'*`8'/`t1'
            gen double `den' = `2'*`3'/`t0' + `6'*`7'/`t1'
            gen double `generate' = `num'/`den'
            if "`novariance'" == "" {
                gen double `vc1' = ((`1'+`4')*(`1'*`4')/`t0'^2 + ///
                    (`5'+`8')*(`5'*`8')/`t1'^2)/`num'^2
                gen double `vc2' = ((`1'+`4')*(`2'*`3')/`t0'^2 + ///
                    (1-(`1'+`4')/`t0')*(`1'*`4')/`t0' + ///
                    (`5'+`8')*(`6'*`7')/`t1'^2 + ///
                    (1-(`5'+`8')/`t1')*(`5'*`8')/`t1')/(`num'*`den')
                gen double `vc3' = ((1-(`1'+`4')/`t0')*(`2'*`3')/`t0' + ///
                    (1-(`5'+`8')/`t1')*(`6'*`7')/`t1')/`den'^2
                gen double `variance' = (`vc1'+`vc2'+`vc3')/2
            }
        }
        if "`novariance'" == "" replace `variance' = . if `variance'<0
    }
end

capture program drop pba_conf_summary_mh
program define pba_conf_summary_mh, rclass
    version 16.0

    syntax, ///
        A(integer) B(integer) C(integer) D(integer) ///
        PONEMIN(real) PONEMODA(real) PONEMODB(real) PONEMAX(real) ///
        PZEROMIN(real) PZEROMODA(real) PZEROMODB(real) PZEROMAX(real) ///
        RRMIN(real) RRMODA(real) RRMODB(real) RRMAX(real) ///
        SIMS(string) ///
        [EFFECT(string) NOGRAPH DRAWSSAVING(string)]

    /**********************************************************************
    Basic checks
    **********************************************************************/

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
        if missing(``x'') | ``x'' <= 0 {
            di as error "`x' must be positive"
            exit 198
        }
    }

    if "`effect'" == "" local effect "RR"
    local effect = upper(strtrim("`effect'"))
    if !inlist("`effect'", "RR", "OR") {
        di as error "effect() must be RR or OR"
        exit 198
    }

    /**********************************************************************
    Check trapezoid parameters
    **********************************************************************/

    foreach x in ponemin ponemoda ponemodb ponemax pzeromin pzeromoda pzeromodb pzeromax rrmin rrmoda rrmodb rrmax {
        if missing(``x'') {
            di as error "`x' must be finite"
            exit 198
        }
    }

    if !(`ponemin' <= `ponemoda' & `ponemoda' <= `ponemodb' & `ponemodb' <= `ponemax') {
        di as error "p1 trapezoid parameters must satisfy ponemin() <= ponemoda() <= ponemodb() <= ponemax()"
        exit 198
    }

    if !(`pzeromin' <= `pzeromoda' & `pzeromoda' <= `pzeromodb' & `pzeromodb' <= `pzeromax') {
        di as error "p0 trapezoid parameters must satisfy pzeromin() <= pzeromoda() <= pzeromodb() <= pzeromax()"
        exit 198
    }

    if !(`rrmin' <= `rrmoda' & `rrmoda' <= `rrmodb' & `rrmodb' <= `rrmax') {
        di as error "RR trapezoid parameters must satisfy rrmin() <= rrmoda() <= rrmodb() <= rrmax()"
        exit 198
    }

    if `ponemin' < 0 | `ponemax' > 1 {
        di as error "p1 trapezoid parameters must be between 0 and 1"
        exit 198
    }

    if `pzeromin' < 0 | `pzeromax' > 1 {
        di as error "p0 trapezoid parameters must be between 0 and 1"
        exit 198
    }

    if `rrmin' <= 0 {
        di as error "RR trapezoid parameters must be positive"
        exit 198
    }

    local den_pone  = `ponemax'  + `ponemodb'  - `ponemin'  - `ponemoda'
    local den_pzero = `pzeromax' + `pzeromodb' - `pzeromin' - `pzeromoda'
    local den_rr    = `rrmax'    + `rrmodb'    - `rrmin'    - `rrmoda'

    if `den_pone' <= 0 {
        di as error "p1 trapezoid has zero or negative area"
        exit 198
    }

    if `den_pzero' <= 0 {
        di as error "p0 trapezoid has zero or negative area"
        exit 198
    }

    if `den_rr' <= 0 {
        di as error "RR trapezoid has zero or negative area"
        exit 198
    }

    local f1_pone  = (`ponemoda'  - `ponemin')  / `den_pone'
    local f2_pone  = (2 * `ponemodb'  - `ponemoda'  - `ponemin')  / `den_pone'

    local f1_pzero = (`pzeromoda' - `pzeromin') / `den_pzero'
    local f2_pzero = (2 * `pzeromodb' - `pzeromoda' - `pzeromin') / `den_pzero'

    local f1_rr    = (`rrmoda'    - `rrmin')    / `den_rr'
    local f2_rr    = (2 * `rrmodb'    - `rrmoda'    - `rrmin')    / `den_rr'

    /**********************************************************************
    Constants
    **********************************************************************/

    local m = `a' + `c'
    local n = `b' + `d'

    if "`effect'" == "RR" {
        local observed_effect = (`a' / `m') / (`b' / `n')
        local v_re = 1 / `a' + 1 / `b' - 1 / `m' - 1 / `n'
        local axis_title "Risk ratio estimate"
    }
    else {
        local observed_effect = (`a' / `b') / (`c' / `d')
        local v_re = 1 / `a' + 1 / `b' + 1 / `c' + 1 / `d'
        local axis_title "Odds ratio estimate"
    }

    preserve

        quietly {
            clear
            set obs `sims'

            /**************************************************************
            Draw p1 from trapezoidal distribution:
            prevalence of confounder among exposed
            **************************************************************/

            gen double u = runiform()
            gen double pone = .

            replace pone = `ponemin' + ///
                sqrt(u * (`ponemoda' - `ponemin') * `den_pone') ///
                if u <= `f1_pone'

            replace pone = (`ponemin' + `ponemoda' + u * `den_pone') / 2 ///
                if u > `f1_pone' & u <= `f2_pone'

            replace pone = `ponemax' - ///
                sqrt((1 - u) * (`ponemax' - `ponemodb') * `den_pone') ///
                if u > `f2_pone'

            drop u

            /**************************************************************
            Draw p0 from trapezoidal distribution:
            prevalence of confounder among unexposed
            **************************************************************/

            gen double u = runiform()
            gen double pzero = .

            replace pzero = `pzeromin' + ///
                sqrt(u * (`pzeromoda' - `pzeromin') * `den_pzero') ///
                if u <= `f1_pzero'

            replace pzero = (`pzeromin' + `pzeromoda' + u * `den_pzero') / 2 ///
                if u > `f1_pzero' & u <= `f2_pzero'

            replace pzero = `pzeromax' - ///
                sqrt((1 - u) * (`pzeromax' - `pzeromodb') * `den_pzero') ///
                if u > `f2_pzero'

            drop u

            /**************************************************************
            Draw confounder-disease RR from trapezoidal distribution
            **************************************************************/

            gen double u = runiform()
            gen double rrcd = .

            replace rrcd = `rrmin' + ///
                sqrt(u * (`rrmoda' - `rrmin') * `den_rr') ///
                if u <= `f1_rr'

            replace rrcd = (`rrmin' + `rrmoda' + u * `den_rr') / 2 ///
                if u > `f1_rr' & u <= `f2_rr'

            replace rrcd = `rrmax' - ///
                sqrt((1 - u) * (`rrmax' - `rrmodb') * `den_rr') ///
                if u > `f2_rr'

            drop u

            /**************************************************************
            Expected stratified cell frequencies
            **************************************************************/

            gen double mone = pone  * `m'
            gen double none = pzero * `n'

            gen double bone = rrcd * none * `b' / (rrcd * none + `n' - none)
            gen double aone = rrcd * mone * `a' / (rrcd * mone + `m' - mone)

            gen double cone = mone - aone
            gen double done = none - bone

            gen double azero = `a' - aone
            gen double bzero = `b' - bone
            gen double czero = `c' - cone
            gen double dzero = `d' - done

            gen double mzero = `m' - mone
            gen double nzero = `n' - none

            gen byte expected_ok = !missing(azero,bzero,czero,dzero,aone,bone,cone,done) & ///
                azero>=0 & bzero>=0 & czero>=0 & dzero>=0 & ///
                aone>=0 & bone>=0 & cone>=0 & done>=0 & ///
                aone<=`a' & bone<=`b' & cone<=`c' & done<=`d'

            /**************************************************************
            Binomial draws for second source of uncertainty
            **************************************************************/

            gen double prob_aone = aone / `a'
            gen double prob_bone = bone / `b'
            gen double prob_cone = cone / `c'
            gen double prob_done = done / `d'

            gen byte probability_ok = expected_ok & ///
                inrange(prob_aone,0,1) & inrange(prob_bone,0,1) & ///
                inrange(prob_cone,0,1) & inrange(prob_done,0,1)

            * Retain systematic estimates before any binomial zero-cell filtering.
            _pba_conf_mh_calc azero bzero czero dzero aone bone cone done, ///
                effect(`effect') generate(effect_syst) novariance
            replace effect_syst = . if !probability_ok | effect_syst<=0
            gen byte valid_systematic = !missing(effect_syst)

            gen double aone_draw = rbinomial(`a', prob_aone) ///
                if probability_ok

            gen double bone_draw = rbinomial(`b', prob_bone) ///
                if probability_ok

            gen double cone_draw = rbinomial(`c', prob_cone) ///
                if probability_ok

            gen double done_draw = rbinomial(`d', prob_done) ///
                if probability_ok

            gen double azero_draw = `a' - aone_draw
            gen double bzero_draw = `b' - bone_draw
            gen double czero_draw = `c' - cone_draw
            gen double dzero_draw = `d' - done_draw

            gen double mone_draw  = aone_draw  + cone_draw
            gen double mzero_draw = azero_draw + czero_draw

            gen double none_draw  = bone_draw  + done_draw
            gen double nzero_draw = bzero_draw + dzero_draw

            /**************************************************************
            Impossible draws
            Matches R flag:
              any of A0, B0, C0, D0, A1, B1, C1, D1 equals zero
            Missing values are also treated as impossible in Stata.
            **************************************************************/

            gen byte simulated_ok = probability_ok & ///
                !missing(azero_draw, bzero_draw, czero_draw, dzero_draw, ///
                         aone_draw, bone_draw, cone_draw, done_draw) & ///
                azero_draw > 0 & bzero_draw > 0 & ///
                czero_draw > 0 & dzero_draw > 0 & ///
                aone_draw > 0 & bone_draw > 0 & ///
                cone_draw > 0 & done_draw > 0

            * Identical effect and SE formulas to calc_mh_effect() in the R PBA.
            _pba_conf_mh_calc azero_draw bzero_draw czero_draw dzero_draw ///
                aone_draw bone_draw cone_draw done_draw, ///
                effect(`effect') generate(effect_mh) variance(var_mh)
            replace effect_mh = . if !simulated_ok
            replace var_mh = . if !simulated_ok
            gen byte effect_ok = simulated_ok & !missing(effect_mh,var_mh) & ///
                effect_mh>0 & var_mh>=0

            gen double effect_total = exp(log(effect_mh) + rnormal() * sqrt(var_mh)) ///
                if effect_ok
            gen byte valid = effect_ok & !missing(effect_total) & effect_total>0
            replace effect_total = . if !valid

            * Random-only uses all requested iterations, matching the R PBA.
            gen double effect_re = exp(log(`observed_effect') + rnormal() * sqrt(`v_re'))
            gen byte valid_random = !missing(effect_re) & effect_re>0
            replace effect_re = . if !valid_random

            gen long draw = _n
            gen byte rejected_expected = !expected_ok
            gen byte rejected_probability = expected_ok & !probability_ok
            gen byte rejected_zero_cells = probability_ok & !simulated_ok
            gen byte rejected_effect = simulated_ok & !effect_ok
            gen byte rejected_total_nonfinite = effect_ok & !valid
            gen byte rejected_systematic_effect = probability_ok & !valid_systematic
            count if valid_systematic
            local n_systematic = r(N)
            count if valid
            local n_valid = r(N)
            count if valid_random
            local n_random = r(N)
            local n_impossible = `sims' - `n_valid'
            local n_impossible_systematic = `sims' - `n_systematic'
            local n_impossible_random = `sims' - `n_random'
            foreach stage in expected probability zero_cells effect total_nonfinite systematic_effect {
                count if rejected_`stage'
                local n_rejected_`stage' = r(N)
            }
            assert rejected_expected + rejected_probability + rejected_zero_cells + ///
                rejected_effect + rejected_total_nonfinite == !valid
            assert rejected_expected + rejected_probability + rejected_systematic_effect == !valid_systematic
            if "`drawssaving'" != "" {
                format pone pzero rrcd effect_re effect_syst effect_total effect_mh var_mh %24.17g
                export delimited draw pone pzero rrcd effect_re effect_syst effect_total effect_mh var_mh ///
                    expected_ok probability_ok simulated_ok effect_ok valid valid_systematic valid_random ///
                    rejected_expected rejected_probability rejected_zero_cells rejected_effect ///
                    rejected_total_nonfinite rejected_systematic_effect using "`drawssaving'", datafmt
            }
        }

        /******************************************************************
        Percentiles
        ******************************************************************/

        foreach component in re syst total {
            local `component'_l = .
            local `component'_m = .
            local `component'_u = .
            local `component'_w = .
            quietly count if !missing(effect_`component')
            if r(N)>0 {
                quietly _pctile effect_`component', p(2.5 50 97.5)
                local `component'_l = r(r1)
                local `component'_m = r(r2)
                local `component'_u = r(r3)
                local `component'_w = ``component'_u' / ``component'_l'
            }
        }

        /******************************************************************
        Store results in a matrix
        ******************************************************************/

        tempname results

        matrix `results' = ///
            (`re_m',    `re_l',    `re_u',    `re_w',    `n_impossible_random' \ ///
             `syst_m',  `syst_l',  `syst_u',  `syst_w',  `n_impossible_systematic' \ ///
             `total_m', `total_l', `total_u', `total_w', `n_impossible')

        matrix colnames `results' = Median Lower Upper Width N_Impossible
        matrix rownames `results' = Random_error Systematic_error Total_error

        /******************************************************************
        Print final results
        ******************************************************************/

        di as text _n "Probabilistic bias analysis: uncontrolled confounding"
        di as text "Simulations:        " as result %12.0fc `sims'
        di as text "Valid simulations:  " as result %12.0fc `n_valid'
        di as text "Impossible draws:   " as result %12.0fc `n_impossible'
        di as text "Effect measure:     " as result "`effect'"
        di as text "Systematic retained:" as result %12.0fc `n_systematic'
        di as text "Observed `effect':        " as result %9.4f `observed_effect' 

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
        gen long nimpossible = .

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

        replace nimpossible = `n_impossible_random' in 1
        replace nimpossible = `n_impossible_systematic' in 2
        replace nimpossible = `n_impossible' in 3

        format median lower upper width %9.4f
        format nimpossible %12.0fc

        di as text _n "Final results"
        list analysis median lower upper width nimpossible, noobs clean abbreviate(20)

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
                xtitle("`axis_title'") ///
                ytitle("") ///
                legend(off)
        }

    restore

    /**********************************************************************
    Return results in r()
    **********************************************************************/

    return matrix results = `results'

    return local effect "`effect'"
    return scalar observed_effect = `observed_effect'
    return scalar crude_se = sqrt(`v_re')
    if "`effect'" == "RR" return scalar observed_rr = `observed_effect'
    if "`effect'" == "OR" return scalar observed_or = `observed_effect'
    return scalar valid_systematic = `n_systematic'
    return scalar impossible_systematic = `n_impossible_systematic'
    return scalar valid_random = `n_random'
    return scalar impossible_random = `n_impossible_random'
    foreach stage in expected probability zero_cells effect total_nonfinite systematic_effect {
        return scalar rejected_`stage' = `n_rejected_`stage''
    }
    return scalar sims         = `sims'
    return scalar valid        = `n_valid'
    return scalar impossible   = `n_impossible'

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
