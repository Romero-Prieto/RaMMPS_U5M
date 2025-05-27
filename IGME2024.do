tempfile   temp
tempfile   temp2
tempfile   temp3

import     delimited "/Users/lshjr3/Documents/RaMMPS/under-five mortality/UN IGME 2024.csv", clear
keep if    geographicarea  == "Malawi"
keep if    seriesname      == "UN IGME estimate"
keep if    wealthquintile  == "Total" 
keep       indicator observationvalue lowerbound upperbound referencedate sex

generate   name             = ""
replace    name             = "child"                  if indicator == "Child Mortality rate age 1-4"
replace    name             = "infant"                 if indicator == "Infant mortality rate"
replace    name             = "non_neonatal_under5"    if indicator == "Mortality rate 1-59 months"
replace    name             = "post_neonatal"          if indicator == "Mortality rate age 1-11 months"
replace    name             = "neonatal"               if indicator == "Neonatal mortality rate"
replace    name             = "stillbirth"             if indicator == "Stillbirth rate"
replace    name             = "under5"                 if indicator == "Under-five mortality rate"
drop if    name            == ""

drop       indicator
rename     referencedate Year
order      sex Year name
label      variable sex              ""
label      variable Year             ""
label      variable observationvalue ""
label      variable lowerbound       ""
label      variable upperbound       ""
save      `temp'

local      lISt             = "neonatal infant under5 stillbirth post_neonatal child non_neonatal_under5"
drop if    Year            != .
keep       sex Year
save      `temp2'

foreach var of local lISt {
	dis       "`var'"
	use       `temp', clear
	keep if    name         == "`var'"
	rename     observationvalue `var'
	rename     lowerbound `var'_LB
	rename     upperbound `var'_UB
	keep       sex Year `var'*
	save      `temp3', replace
	
	use       `temp2', clear
	merge 1:1  sex Year using `temp3', nogenerate noreport
	save      `temp2', replace
	}

keep if    sex             == "Total"
drop       sex
gsort     -Year
export     delimited using "/Users/lshjr3/Documents/RaMMPS/under-five mortality/IGMEmalawi.csv", replace
