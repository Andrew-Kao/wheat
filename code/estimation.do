******** full estimation ********

clear all
set tracedepth 2
set trace on
program drop _all
macro drop _all
set seed 2020

cd "/Users/AndrewKao/Documents/Grad/G2/QuantEcon/Modeling/wheat/"

*** get datasets

** demand
import delim "raw/population_panel.csv", clear varn(4) 
keep if countryname == "World"
reshape long v, i(countryname) j(year)
replace year = year + 1955
ren v population
keep year population
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
merge 1:m year using `wheat_panel', keep(3) nogen
gen logP = log(p)

encode country, gen(id)
tempfile panel
save `panel'



********* MAIN RESULT (POINT EST.) *********




* estimate beta
reg log_q_demand log_pop, r



********* BOOTSTRAP (CONFIDENCE SET) *********




































