******** full estimation ********

clear all
*set tracedepth 2
*set trace on
*program drop _all
*macro drop _all
set seed 2020

cd "/Users/AndrewKao/Documents/Grad/G2/QuantEcon/Modeling/wheat/"
* cd "~/Dropbox/GitHub/wheat/"

cap erase  "results.dta"

*** get datasets

** demand
import delim "raw/population_panel.csv", clear varn(4) 
keep if countryname == "World"
reshape long v, i(countryname) j(year)
replace year = year + 1955
ren v population
keep year population
global pop_2021 = population[62]
global pop_2022 = population[63]
tempfile pop
save `pop'

import delim "raw/wheat_time_series.csv", clear
merge 1:1 year using `pop', keep(3) nogen

gen log_q_demand = log(q_demand)
gen log_pop = log(population)
tempfile demand
save `demand'


** supply shock
import excel "raw/cost_data.xlsx", sheet("data") clear first
tempfile shock
save `shock'


** panel data
insheet using "raw/wheat_production_panel.csv", clear
reshape long y, i(country) j(year)
rename y Q
gen logQ = log(Q)
tempfile wheat_panel
save `wheat_panel'

insheet using "raw/wheat_time_series.csv", clear
merge 1:m year using `wheat_panel' , keep(3) nogen
gen logP = log(p)
drop if inlist(country, "Serbia and Montenegro", "USSR", "Czechoslovakia", "Yugoslav SFR") //dropping countries that don't exist in 2020
// drop if inlist(country, "Botswana", "Chad", "Malta", "New Caledonia", "Palestine", "Qatar", "United Arab Emirates")
drop if missing(Q) //dropping countries with no observations (newer countries before they were formed)
encode country, gen(id)
su id if country == "Ukraine"
local ukraine_id = r(mean)

tempfile panel
save `panel'


********* MAIN RESULT (POINT EST.) *********

* estimate demand
use `demand', clear
reg log_q_demand log_pop, r
local gamma_d = exp(e(b)[1,2])
local delta_d = e(b)[1,1]
local Q_D_2022 = `gamma_d'*$pop_2022^`delta_d'
predict demand_hat
gen log_zeta = log_q_demand - demand_hat
tempfile zeta
save `zeta'

*Expectation and variance of log supply shocks 
use `shock', clear
sum log_nu, d
local nu_mean = r(mean)
local nu_var = r(Var)


* estimate beta
use `panel', clear
reg logP logQ i.id, r 
local beta = e(b)[1,1]
bysort id: egen m_q = mean(logQ)
gen q = logQ - m_q
sum q, d
local q_var = r(Var)

*drop if inlist(country, "Botswana", "Chad", "Malta", "New Caledonia", "Palestine", "Qatar", "United Arab Emirates")
* estimate alpha
local alpha = 1+`beta'*`q_var'/(`q_var' - (1-1/29)*`nu_var')

* estimate gamma_c
gen gamma_help = logP - log(`alpha') + (`alpha' -1)*(`nu_mean' - logQ)
bysort id: egen gamma_help_2 = mean(gamma_help)
gen gamma_c = exp(gamma_help_2)

gen log_nu = (logP - log(`alpha') + (`alpha'-1)*(-logQ) - log(gamma_c))/(`alpha' -1)

tempfile x
save `x'
local ptest = 0
forv i = 1/1000 {
use `x', clear
preserve
gen fake = runiform(0,1)
collapse (sum) fake,by(year)
sort fake
gen year2 = _n + 1991
tempfile z
qui save `z'
restore
qui merge n:1 year using `z', keep(1 3) nogen
drop year
ren year2 year
qui merge n:1 year using `zeta'
qui gen nu = exp(log_nu)
qui gen zeta = exp(log_zeta)
qui reg nu zeta, r
if (2 * ttail(e(df_r), abs(_b[zeta]/_se[zeta])))*4000 < 0.01 {
	local ptest = `ptest' + 1
}
qui reg log_nu log_zeta, r
if (2 * ttail(e(df_r), abs(_b[log_zeta]/_se[log_zeta])))*4000 < 0.01 {
	local ptest = `ptest' + 1
}
qui reg log_nu zeta, r
if (2 * ttail(e(df_r), abs(_b[zeta]/_se[zeta])))*4000 < 0.01 {
	local ptest = `ptest' + 1
}
qui reg nu log_zeta, r
if (2 * ttail(e(df_r), abs(_b[log_zeta]/_se[log_zeta])))*4000 < 0.01 {
	local ptest = `ptest' + 1
}
if mod(`i',100) == 0 {
	di `i'
}
}
di `ptest'
block

keep id gamma_c
duplicates drop

preserve
keep id
save "country_id", replace
restore

* getting quantities that firms plan to produce in 2022
foreach i of numlist 1/124 {
gen gamma_c_alpha_`i' = (gamma_c/gamma_c[`i'])^(1/(1-`alpha'))
egen  sum_gamma_c_alpha_`i' = total(gamma_c_alpha_`i')
gen Q_tilde_`i' = `Q_D_2022'/sum_gamma_c_alpha_`i'
gen gamma_c_`i' = gamma_c[`i']
}

* equilibrium price based on those quantities in 2022
gen P_2022 = Q_tilde_1^(`alpha'-1)*`alpha'*gamma_c[1]

*dropping gamma_ukraine
drop if id == `ukraine_id'

* getting counterfactual quantities that firms plan to produce in 2022
foreach i of numlist 1/113 115/124 {
egen  sum_gamma_c_alpha_`i'_co = total(gamma_c_alpha_`i')
gen Q_tilde_`i'_co = `Q_D_2022'/sum_gamma_c_alpha_`i'_co	
}

* counterfactual equilibrium price based on those quantities in 2022
gen P_2022_co = Q_tilde_1_co^(`alpha'-1)*`alpha'*gamma_c[1]

drop sum_gamma_c_alpha_* gamma_c_alpha_*

***saving the results of the actual simulation
keep if id == 1
drop id gamma_c
gen exp_nu = `nu_mean'
gen var_nu = `nu_var'
gen gamma_d = `gamma_d'
gen delta = `delta_d'
gen Q_D_2022_hat = `Q_D_2022'
gen beta = `beta'
gen var_q = `q_var'
gen alpha = `alpha'

gen draw = 0

cap confirm file "results.dta"
if _rc==601 {
	save "results.dta", replace
} 
else {
	append using "results.dta", force
	save "results.dta", replace
}


********* BOOTSTRAP (CONFIDENCE SET) *********
foreach b of numlist 1/1000 {

* estimate demand
use `demand', clear
bsample
qui reg log_q_demand log_pop, r
local gamma_d = exp(e(b)[1,2])
local delta_d = e(b)[1,1]
local Q_D_2022 = `gamma_d'*$pop_2022^`delta_d'

*Expectation and variance of log supply shocks 
use `shock', clear
bsample
qui sum log_nu, d
local nu_mean = r(mean)
local nu_var = r(Var)


* estimate beta
use `panel', clear
bsample, strata(id)
qui reg logP logQ i.id, r 
local beta = e(b)[1,1]
bysort id: egen m_q = mean(logQ)
gen q = logQ - m_q
qui sum q, d
local q_var = r(Var)

* estimate alpha
local alpha = 1+`beta'*`q_var'/(`q_var' - (1-1/29)*`nu_var')

if (`alpha' > 1) {

* estimate gamma_c
gen gamma_help = logP - log(`alpha') + (`alpha' -1)*(`nu_mean' - logQ)
bysort id: egen gamma_help_2 = mean(gamma_help)
gen gamma_c = exp(gamma_help_2)

keep id gamma_c
qui duplicates drop

* getting quantities that firms plan to produce in 2022
foreach i of numlist 1/124 {
gen gamma_c_alpha_`i' = (gamma_c/gamma_c[`i'])^(1/(1-`alpha'))
egen  sum_gamma_c_alpha_`i' = total(gamma_c_alpha_`i')
gen Q_tilde_`i' = `Q_D_2022'/sum_gamma_c_alpha_`i'
gen gamma_c_`i' = gamma_c[`i']
}

* equilibrium price based on those quantities in 2022
gen P_2022 = Q_tilde_1^(`alpha'-1)*`alpha'*gamma_c[1]

if missing(P_2022) {
	replace P_2022 = Q_tilde_2^(`alpha'-1)*`alpha'*gamma_c[2]
}

if missing(P_2022) {
	replace P_2022 = Q_tilde_3^(`alpha'-1)*`alpha'*gamma_c[3]
}

*dropping gamma_ukraine
drop if id == `ukraine_id'

* getting counterfactual quantities that firms plan to produce in 2022
foreach i of numlist 1/113 115/124 {
egen  sum_gamma_c_alpha_`i'_co = total(gamma_c_alpha_`i')
gen Q_tilde_`i'_co = `Q_D_2022'/sum_gamma_c_alpha_`i'_co	
}

* counterfactual equilibrium price based on those quantities in 2022
gen P_2022_co = Q_tilde_1_co^(`alpha'-1)*`alpha'*gamma_c[1]
if missing(P_2022_co) {
	replace P_2022_co = Q_tilde_2_co^(`alpha'-1)*`alpha'*gamma_c[2]
}
if missing(P_2022_co) {
	replace P_2022_co = Q_tilde_3_co^(`alpha'-1)*`alpha'*gamma_c[3]
}


drop sum_gamma_c_alpha_* gamma_c_alpha_*

***saving the results of the actual simulation
keep if id == 1
drop id gamma_c

}
else {
	clear 
	set obs 1
}

gen exp_nu = `nu_mean'
gen var_nu = `nu_var'
gen gamma_d = `gamma_d'
gen delta = `delta_d'
gen Q_D_2022_hat = `Q_D_2022'
gen beta = `beta'
gen var_q = `q_var'
gen alpha = `alpha'

gen draw = `b'

cap confirm file "results.dta"
if _rc==601 {
	save "results.dta", replace
} 
else {
	append using "results.dta", force
	save "results.dta", replace
}

}

gen P_diff = P_2022_co - P_2022

*GETTING 95% CI FOR P_2022_co and other variables
foreach v in "P_2022_co" "Q_tilde_1_co" "P_diff" {
	
preserve
*drop simulations that are rejected by our model
drop if alpha <= 1


*dropping the point estimate to get CS
drop if draw == 0
drop if missing(`v')

sort `v'
qui sum `v', d
local `v'_l95 = round(`v'[round(r(N)/100*2.5)],0.01)
local `v'_u95 = round(`v'[round(r(N)/100*97.5)],0.01)
local nsim = r(N)
restore 

di "`v' is " round(`v'[1001],0.01) " with 95\% confidence interval [" ``v'_l95' "; " ``v'_u95' "] based on " `nsim' " simulations"

}

*** graphics
block
use "results.dta", clear
keep if draw == 0
keep Q_* draw

preserve
* real data
keep Q_*_co draw
ren (Q_*_co) (Q_*)
reshape long Q_tilde_, i(draw) j(id)
ren Q_tilde_ Q_co
tempfile co
save `co'
restore

drop Q_*_co
reshape long Q_tilde_, i(draw) j(id)
ren Q_tilde_ Q_tilde
merge 1:1 id using `co', keep(1 3) nogen

** TODO: merge id data
insheet using "raw/wheat_time_series.csv", clear
merge 1:m year using `wheat_panel' , keep(3) nogen
gen logP = log(p)
drop if inlist(country, "Serbia and Montenegro", "USSR", "Czechoslovakia", "Yugoslav SFR") //dropping countries that don't exist in 2020
// drop if inlist(country, "Botswana", "Chad", "Malta", "New Caledonia", "Palestine", "Qatar", "United Arab Emirates")
drop if missing(Q) //dropping countries with no observations (newer countries before they were formed)
encode country, gen(id)

** then make graphic...
