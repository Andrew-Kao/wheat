clear all
set tracedepth 2
set trace on
set scheme cleanplots
program drop _all
macro drop _all
set seed 2020


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
tsset id year


** bootstrap at year level

reg logP logQ i.id, r



