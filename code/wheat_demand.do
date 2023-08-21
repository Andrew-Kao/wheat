*** demand 
set more off
import delim "/Users/AndrewKao/Documents/Grad/G2/QuantEcon/Modeling/wheat/raw/population_panel.csv", clear varn(4)
keep if countryname == "World"
reshape long v, i(countryname) j(year)
replace year = year + 1955
ren v population
keep year population
tempfile pop
save `pop'

import delim "/Users/AndrewKao/Documents/Grad/G2/QuantEcon/Modeling/wheat/raw/wheat_time_series.csv", clear

merge 1:1 year using `pop', keep(3) nogen

gen p_bushel = p/36.74

gen overflow = q_production - q_demand 
gen cumulative_stock = sum(overflow) + 150

// twoway line q_demand q_production q_ukraine year
twoway line overflow cumulative_stock q_production year

gen log_q_demand = log(q_demand)
gen log_pop = log(population)
gen log_p = log(p)

reg log_q_demand log_pop log_p, r
local epsilon = _b[p]
reg log_q_demand log_pop, r
// log pop  



* effect of removing Ukraine from supply


* model representative firm in Ukraine
* model global firm: what is marginal cost curve? 

* convex cost? farms are price takers

* demand curve...

* "In recent years, low international wheat prices have often encouraged farmers in the United States to change to more profitable crops. In 1998, the price at harvest of a 60 pounds (27 kg) bushel[142] was $2.68 per.[143] Some information providers, following CBOT practice, quote the wheat market in per ton denomination.[144] A USDA report revealed that in 1998, average operating costs were $1.43 per bushel and total costs were $3.97 per bushel.[143] In that study, farm wheat yields averaged 41.7 bushels per acre (2.2435 metric ton/hectare), and typical total wheat production value was $31,900 per farm, with total farm production value (including other crops) of $173,681 per farm, plus $17,402 in government payments. There were significant profitability differences between low- and high-cost farms, due to crop yield differences, location, and farm size. "

* technology is improving -- what about increases?
* demand is increasing (harder to estimate?)

** see sb974-5_1 for CDF of wheat costs
