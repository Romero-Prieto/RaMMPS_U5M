local      pATh         = "/Users/lshjr3/Documents/RaMMPS/under-five mortality" /*adjust path*/
tempfile   temp
tempfile   temp2
tempfile   temp3
tempfile   data
tempfile   sample
tempfile   SBH
tempfile   TPH
tempfile   FPH
tempfile   Histories
tempfile   HHDeaths
tempfile   mics
tempfile   dhs

use      "/Users/lshjr3/Desktop/RaMMPS/MW_AnalyticSample.dta", clear
local      lISt         = ""
foreach var of varlist _all {
	local      lISt         = "`lISt'" + " " + "`var'"
	if substr("`: type `var''",1,3) != "str" & "`var'" != "SubmissionDate" & "`var'" != "starttime" & "`var'" != "endtime" { /*to all variables as strings*/
		quiet generate   temp = string(`var',"`: format `var''") if `var' != .
		drop      `var' 
		quiet generate  `var' = temp
		drop       temp
		order     `lISt'
		}
	}
replace    caseid       = "A" + caseid
save      `data', replace

keep       SubmissionDate starttime endtime main_language_label_1 username subscriberid user_language caseid region_1 b6_1 b7_1 b8_1 d1 d2 d3* e3_1 E4 e4* e2_label_1 e4_label_1 hd1_a duration phone_call_duration household_duration
duplicates drop

destring   region_1 b6_1 b7_1 b8_1 d1 d2 e3_1 e4*, replace
rename     hd1_a kinship
rename     region_1 Region
rename     e4_label_1 district
rename     main_language_label_1 language
rename     subscriberid enumerator
replace    enumerator     = "E" + enumerator if enumerator != ""
rename     user_language enumerator_language
rename     username enumerator_name
rename     e3_1 age  
generate   ageG           = 0                if age  < 18
replace    ageG           = 1                if age  < 50 & ageG  == .
replace    ageG           = 2                if age  < 65 & ageG  == .
recode     ageG        (. = 3)
generate   sex            = 1                if e2_label_1 == "MALE"
recode     sex         (. = 2)               if e2_label_1 == "FEMALE"

generate   UR             = 1                
recode     UR          (1 = 2)               if e4b_1 == 2
generate   Electricity    = 1                if b7_1  == 1
recode     Electricity (. = 0)               if b7_1  == 2
generate   Roofing        = 1                if b6_1  == 6  |  b6_1  == 7 /*metal or cement*/
recode     Roofing     (. = 0)               if b6_1  != .  & (b6_1  != 8 | b6_1  != 9)  
generate   Water          = 1                if b8_1  == 1  | b8_1  == 2  | b8_1  == 3 | b8_1  == 4 | b8_1  == 5 | b8_1  == 7 | b8_1  == 9 | b8_1  == 10 | b8_1  == 11
recode     Water       (. = 0)               if b8_1  != .  & b8_1  != 16
recode     Electricity Roofing Water (0 = 1) (1 = 2)

/*Education cleaning*/
replace    d3_specify_1   = strlower(d3_specify_1)
local      lISt           = "advanced college degree diploma diplama certificate bachelor technical tertiary tetiary tetially teaching months short police business level"
foreach name of local lISt {
	replace    d3_specify_1 = "`name'" if strpos(d3_specify_1,"`name'") > 0
	}
replace    d3_specify_1   = "tertiary"       if d3_specify_1 == "tetiary"            | d3_specify_1 == "tetially"
replace    d3_specify_1   = "diploma"        if d3_specify_1 == "diplama"
replace    d3_specify_1   = "months"         if d3_specify_1 == "short"
replace    d3_specify_1   = "undisclosed"    if d3_specify_1 == "she refused to say" | d3_specify_1 == "don't want to disclose" 

replace    d3_1           = "7"              if d3_1         == "t"
destring   d3_1, replace
replace    d3_1           = 7                if d3_specify_1 == "grade 7"
replace    d3_1           = 0                if d3_specify_1 == "no school attended"
replace    d3_1           = 0                if d1_1         == 0                    | d1_1         == 98
tabulate   d3_specify_1 d3_1

/*Classification of Attainment following DHS hv109*/
generate   Education      = 8                if d3_1         == .  /*not reported*/
recode     Education   (. = 0)               if d3_1         == 0  /*no education (not attended)*/
recode     Education   (. = 1)               if d3_1          < 8  /*incomplete primary*/
recode     Education   (. = 2)               if d3_1         == 8  /*complete primary*/
recode     Education   (. = 3)               if d3_1          < 12 /*incomplete secondary*/
recode     Education   (. = 4)               if d3_1         == 12 /*complete secondary*/
recode     Education   (. = 5)               if d3_1          > 12 /*higher education*/
recode     Education   (8 = .)
recode     Education   (0 = 1) (3 = 2) (4 = 3) (5 = 3)

rename     duration SVY_duration
rename     phone_call_duration call_duration
rename     household_duration_1 HHD_duration
keep       SubmissionDate starttime endtime enumerator* language district caseid Region age ageG Education* sex kinship UR Electricity Roofing Water SVY_duration HHD_duration call_duration
order      SubmissionDate starttime endtime enumerator* language district caseid Region age ageG Education* sex kinship UR Electricity Roofing Water SVY_duration HHD_duration call_duration
save      `sample', replace

use       `data', clear
keep       caseid hd* ph* fph* *pregnancy* FPH_TPH_SSH interviewer_gender_1 correct*
rename     fph228_1_*_*_* fph228_days_*_*_* /*to have the actual name in the questionnaire and codebook*/
rename     fph228_2_*_*_* fph228_months_*_*_* 
rename     fph228_3_*_*_* fph228_years_*_*_*

rename     hd8a1* hd8aX* 
rename     hd8b1* hd8bX*
forvalues j = 1(1)9 {
	capture rename     hd8a_`j'_1_* hd8a`j'_1_*
	capture rename     hd8b_`j'_1_* hd8b`j'_1_*
	}

drop       *phone*
duplicates drop

generate   sUFiX          = ""
order      caseid sUFiX 
save      `temp', replace

drop if    sUFiX         == ""
keep       caseid sUFiX
save      `temp2', replace
save      `temp3', replace
forvalues j = 1(1)4 {
	forvalues i = 0(1)20 {
		use       `temp', clear
		local      GO = 0
		foreach var of varlist _all {
			if substr("`var'",-floor(`i'/10) - 2,.) == "_`i'" {
				local  GO = 1
				}
			}
		
		if `GO' == 1 {
			keep       caseid sUFiX *_`i'
			duplicates drop
			generate   tESt           = ""
			foreach var of varlist _all {
				replace    tESt = "OK" if "`var'" != "caseid" & "`var'" != "sUFiX" & "`var'" != "tESt" & `var' != ""
				}
			keep if    tESt          == "OK"
			replace    sUFiX          = "_`i'" + sUFiX
			rename     *_`i' *
			drop       tESt
			append     using `temp2'
			save      `temp2', replace
			
			use       `temp', clear
			drop       *_`i'
			save      `temp', replace
			}
		}
		use       `temp3', clear
		merge 1:1  caseid sUFiX using `temp', nogenerate
		save      `temp3', replace
		
		use       `temp2', clear
		save      `temp', replace
		drop if    sUFiX         == "" | sUFiX != "" 
		keep       caseid sUFiX
		save      `temp2', replace
	}

use       `temp3', clear
generate   tESt           = ""
foreach var of varlist _all {
	replace    tESt = "OK" if "`var'" != "caseid" & "`var'" != "sUFiX" & "`var'" != "tESt" & `var' != ""
	}
keep if    tESt          == "OK"
drop       tESt

order      caseid sUFiX hd* FPH_TPH_SSH ph* fph* *pregnancy* interviewer_gender *duration
sort       caseid sUFiX
save      `temp', replace

keep       caseid sUFiX FPH_TPH_SSH ph1* ph2 ph3* ph4* ph5* ph6* ph7* ph8* fph201 fph202* fph203* fph204* fph205* fph206* fph207* interviewer_gender
drop if    ph2           == "" & fph201 == ""
generate   mother         = 1                            if ph2         == "1" | fph201      == "1"
recode     mother      (. = 0)                           if ph2         == "0" | fph201      == "2"
recode     mother      (. = 99)                          if ph2         == "99"| fph201      == "3"
recode     mother      (. = 98)                          if ph2         != ""  | fph201      != ""

destring   ph4* ph6* ph8* fph203* fph205* fph207*, replace
recode     ph4* ph6* ph8* (. = 0)                        if ph2         == "1" | ph2         == "0"
recode     fph203* fph205* fph207* (. = 0)               if fph201      == "1" | fph201      == "2"
 
generate   sons           = ph4a + ph6a + ph8a
generate   daughters      = ph4b + ph6b + ph8b
generate   sonsD          = ph8a 
generate   daughtersD     = ph8b

replace    sons           = fph203a + fph205a + fph207a  if sons        == .
replace    daughters      = fph203b + fph205b + fph207b  if daughters   == .
replace    sonsD          = fph207a                      if sonsD       == .
replace    daughtersD     = fph207b                      if daughtersD  == .
destring   FPH_TPH_SSH, replace

generate   sample         = ""
replace    sample         = "SBH"                        if FPH_TPH_SSH == .
replace    sample         = "FPH"                        if FPH_TPH_SSH <= .5  & sample      == "" 
replace    sample         = "TPH"                        if FPH_TPH_SSH  > .5  & sample      == ""
drop if    mother        == 99 /*refuses*/
recode     mother     (98 = .)
bysort     caseid: generate   k = _n
bysort     caseid: generate   K = _N
keep       caseid sample sons daughters sonsD daughtersD mother interviewer_gender k K
order      caseid sample
save      `SBH', replace


use       `temp', clear
keep       caseid sUFiX FPH_TPH_SSH ph* *pregnancy* interviewer_gender *duration
generate   t              = 1                            if FPH_TPH_SSH != ""           & ph10        != ""
bysort     caseid: egen       T = max(t)
keep if    T             == 1

generate   pregnant       = 1                            if ph10        == "1"
recode     pregnant    (. = 0)                           if ph10        == "2"
recode     pregnant    (. = 98)                          if ph10        == "98"
recode     pregnant    (. = 99)                          if ph10        == "99"
destring   ph11_months ph11_weeks, replace
generate   weeks2term     = ph11_months*4.34             if ph11        == "1"
replace    weeks2term     = ph11_weeks                   if ph11        == "2"

generate   everPregnant   = 1                            if ph12        == "1"          | ph13        == "1"
recode     everPreg    (. = 0)                           if ph12        == "2"          | ph13        == "2"
drop       ph10* ph11* ph12* ph13*

generate   birth          = "livebirth"                  if ph16        == "1"
generate   flagPQ         = "re-classified"              if ph17        == "1"          & birth       == ""
replace    birth          = "livebirth"                  if ph17        == "1"         	
replace    birth          = "unknown"                    if birth       == ""           & (ph16       == "98" | ph16       == "99")
replace    birth          = "stillbirth"                 if birth       == ""           & ph16        != ""  
generate   name           = "given"                      if ph18        == "1"
replace    name           = "not given"                  if ph18        == "2"
generate   sex            = ph19                         if ph19        == "1"          | ph19        == "2"
destring   ph20* ph21* ph25* ph26* ph24* ph23*, replace
recode     ph20* ph21* ph25* ph26* ph24* (98 = .) (99 = .) (2222 = .)
replace    ph20b          = ph25b                        if ph20b       == .            & ph25b       != .
replace    ph20           = ph25                         if ph20        == .            & ph25        != .
generate   B_min          = mdy(ph20,1,ph20b)
replace    B_min          = mdy(1,1,ph20b)               if B_min       == . 
generate   B_max          = mdy(1 + mod(ph20,12),1,ph20b + floor(ph20/12)) - 1
replace    B_max          = mdy(12,31,ph20b)             if B_max       == . 
format     %tdDD/NN/CCYY B_*

replace    ph21           = ph26                         if ph21        == .            & ph26        != .
replace    ph21b          = ph26c                        if ph21b       == .            & ph26c       != .
replace    ph21a          = ph26b                        if ph21a       == .            & ph26b       != .
generate   gestation      = ph21a                        if ph21        == 1
replace    gestation      = ph21b*4.34                   if ph21        == 2
generate   age_children   = ph23b/365.25                 if ph23        == 1
replace    age            = ph23c/12                     if ph23        == 2
replace    age            = ph23d                        if ph23        == 3
generate   survival       = "alive"                      if ph22        == "1"
replace    survival       = "dead"                       if ph22        == "2"
replace    survival       = "unknown"                    if ph22        == "99"         | year(B_max)  < 2014
replace    survival       = ""                           if birth       != "livebirth"
generate   D_min          = 0                            if birth       == "livebirth"  & survival    == ""
generate   D_max          = 0                            if birth       == "livebirth"  & survival    == ""
replace    survival       = "dead"                       if birth       == "livebirth"  & survival    == ""
replace    D_min          = ph24b                        if ph24        == 1            & (D_min      == .  | D_min   == 0) & survival == "dead"
replace    D_min          = ph24c*365.25/12              if ph24        == 2            & (D_min      == .  | D_min   == 0) & survival == "dead"
replace    D_min          = ph24d*365.25                 if ph24        == 3            & (D_min      == .  | D_min   == 0) & survival == "dead"
replace    D_max          = ph24b + 1                    if ph24        == 1            & (D_max      == .  | D_max   == 0) & survival == "dead"
replace    D_max          = (ph24c + 1)*365.25/12        if ph24        == 2            & (D_max      == .  | D_max   == 0) & survival == "dead"
replace    D_max          = (ph24d + 1)*365.25           if ph24        == 3            & (D_max      == .  | D_max   == 0) & survival == "dead"

generate   temp           = 1                            if birth       == "livebirth"
bysort     caseid: egen       mother = max(temp)
recode     mother      (. = 0)
keep       caseid sUFiX pregnant weeks2term everPregnant mother *duration interviewer_gender flagPQ name birth sex B_* gestation age survival D_*
order      caseid sUFiX pregnant weeks2term everPregnant mother *duration interviewer_gender flagPQ birth name sex B_* gestation age survival D_*
save      `temp2', replace

keep if    pregnant      != .
keep       caseid pregnant weeks2term everPregnant mother *duration interviewer_gender
save      `temp3', replace

use       `temp2', clear 
keep if    birth         != ""
drop       sUFiX
merge m:1  caseid using `temp3', nogenerate update
drop if    pregnant      == 99 /*refuses*/
recode     pregnant   (98 = .) 
generate   sample         = "TPH"
order      caseid sample
sort       caseid B_min
bysort     caseid: generate   k = _n
bysort     caseid: generate   K = _N
save      `TPH', replace


use       `temp', clear
keep       caseid sUFiX FPH_TPH_SSH fph* *pregnancy* interviewer_gender *duration
generate   t              = 1                            if FPH_TPH_SSH != ""           & fph201      != ""
bysort     caseid: egen       T = max(t)
keep if    T             == 1

generate   mother         = 1                            if fph201      == "1"
recode     mother      (. = 0)                           if fph201      == "2"
recode     mother      (. = 99)                          if fph201      == "3"
recode     mother      (. = 98)                          if fph201      != ""

destring   fph211 fph_total, replace
generate   births_To      = fph_total
generate   stillbirths_To = fph211                       if fph210      == "1"
recode     stillbirths (. = 0)                           if fph210      == "2"

generate   birth          = "livebirth"                  if fph216      == "1"
generate   flagPQ         = "re-classified"              if fph217      == "1"          & birth       == ""
replace    birth          = "livebirth"                  if fph217      == "1"         	
replace    birth          = "unknown"                    if birth       == ""           & (fph216     == "98" | fph216      == "99") 
replace    birth          = "stillbirth"                 if birth       == ""           & fph216      != ""  
generate   name           = "given"                      if fph218      == "1"
replace    name           = "not given"                  if fph218      == "2"
generate   sex            = fph219                       if fph219      == "1"          | fph219      == "2"
destring   fph220 fph221 fph222 fph229* fph228* fph230* fph225, replace
recode     fph220 fph221 fph222 fph229* fph228* fph230* fph225 (98 = .) (99 = .) (-98 = .)
replace    fph220         = fph229_years                                            if (birth  == "stillbirth" & fph229   == 1) | birth  == "unknown"
replace    fph221         = fph229b                                                 if (birth  == "stillbirth" & fph229   == 1) | birth  == "unknown"
generate   B_min          = mdy(fph221,fph222,fph220)
replace    B_min          = mdy(fph221,1,fph220)                                    if B_min == . 
replace    B_min          = mdy(1,fph221,fph220)                                    if B_min == .
replace    B_min          = mdy(1,1,fph220)                                         if B_min == .
generate   B_max          = mdy(fph221,fph222,fph220)
replace    B_max          = mdy(1 + mod(fph221,12),1,fph220 + floor(fph221/12)) - 1 if B_max == .
replace    B_max          = mdy(12,fph222,fph220)                                   if B_max == .
replace    B_max          = mdy(12,31,fph220)                                       if B_max == .
format     %tdDD/NN/CCYY B_*

generate   gestation      = fph230b                                                 if fph230      == 2
replace    gestation      = fph230c*4.34                                            if fph230      == 1
generate   age_children   = fph225
generate   survival       = "alive"                                                 if fph224      == "1"
replace    survival       = "dead"                                                  if fph224      == "2"
replace    survival       = "unknown"                                               if fph224      == "3"
generate   D_min          = 0                                                       if birth       == "livebirth"  & survival    == ""
generate   D_max          = 0                                                       if birth       == "livebirth"  & survival    == ""
replace    survival       = "dead"                                                  if birth       == "livebirth"  & survival    == ""

replace    D_min          = fph228_days                                             if fph228      == 1            & D_min       == .
replace    D_min          = fph228b_months*365.25/12                                if fph228b     == 2            & D_min       == .
replace    D_min          = fph228_months*365.25/12                                 if fph228      == 2            & D_min       == .
replace    D_min          = fph228_years*365.25                                     if fph228      == 3            & D_min       == .
replace    D_max          = fph228_days + 1                                         if fph228      == 1            & D_max       == .
replace    D_max          = (fph228b_months + 1)*365.25/12                          if fph228b     == 2            & D_max       == .
replace    D_max          = (fph228_months + 1)*365.25/12                           if fph228      == 2            & D_max       == .
replace    D_max          = (fph228_years + 1)*365.25                               if fph228      == 3            & D_max       == .

keep       caseid sUFiX mother *duration interviewer_gender flagPQ name birth sex B_* gestation age survival D_* 
order      caseid sUFiX mother *duration interviewer_gender flagPQ name birth sex B_* gestation age survival D_*
save      `temp2', replace

keep if    mother        != .
keep       caseid mother *duration interviewer_gender
save      `temp3', replace

use       `temp2', clear 
keep if    birth         != ""
drop       sUFiX
merge m:1  caseid using `temp3', nogenerate update
drop if    mother        == 99 /*refuses*/
recode     mother     (98 = .) 
generate   sample         = "FPH"
order      caseid sample
sort       caseid B_min
bysort     caseid: generate   k = _n
bysort     caseid: generate   K = _N
save      `FPH', replace

use       `TPH', clear
merge 1:1  caseid sample k using `FPH', nogenerate
merge m:1  caseid sample using `SBH', nogenerate
rename     fullpregnancy_duration  FPH_duration 
rename     truncpregnancy_duration TPH_duration
replace    FPH_duration   = ""                           if sample      != "FPH"
replace    TPH_duration   = ""                           if sample      != "TPH"
destring   *_duration, replace
recode     *duration   (0 = .)
recode     D_max       (0 = 0.25)                        if D_min       == 0
replace    gestation      = .                            if gestation    < 24           & survival    != "dead"
save      `Histories', replace

keep       sample caseid birth survival B_* D_* gestation
drop if    B_min         == .
sort       sample caseid B_min birth
bysort     sample caseid B_min: generate   number = _n
keep if    number        == 1

generate   gestation_I    = gestation
recode     gestation_I (. = 37)                                            if birth       == "livebirth"
replace    gestation_I    = max(gestation_I,20)                            if birth       == "livebirth"
recode     gestation_I (. = 19)                                            if birth       != "livebirth"
  
generate   P1             = 365.25*2/12  
generate   P2             = min(max((D_min + D_max)/2 + P1/2,P1),365.25/2) if survival    == "dead"
replace    P2             = 365.25/2                                       if survival    == "alive"
replace    P2             = P1                                             if birth       != "livebirth"
replace    P2             = P1                                             if P2          == .
generate   G              = (B_min + B_max)/2 - gestation_I*7
generate   P              = (B_min + B_max)/2 + P1
generate   P_ext          = (B_min + B_max)/2 + P2
drop       P1 P2 gestation* survival D_* number
format     %tdDD/NN/CCYY G P* 

sort       sample caseid B_min
bysort     sample caseid: generate   n = _n
bysort     sample caseid: generate   N = _N
drop if    N             == 1
save      `temp2', replace

replace    n              = n + 1
rename    (G P* birth) previous_=
keep       sample caseid n previous_*
save      `temp3', replace

use       `temp2', clear
merge 1:1  sample caseid n using `temp3', nogenerate noreport keep(master match)
save      `temp2', replace

replace    n              = n - 1
rename    (G P* birth) next_=
keep       sample caseid n next_*
save      `temp3', replace

use       `temp2', clear
merge 1:1  sample caseid n using `temp3', nogenerate noreport keep(master match)
generate   FlagBI         = .
replace    FlagBI         = 1                                              if P                 > next_G          & next_birth     == "livebirth"
replace    FlagBI         = 1                                              if G                 < previous_P      & previous_birth == "livebirth"
replace    FlagBI         = 1                                              if G                 < previous_P_ext  & previous_birth == "livebirth"
keep if    FlagBI        != .
keep       sample caseid B_min FlagBI
save      `temp2', replace

use       `Histories', clear
merge m:1  sample caseid B_min using `temp2', nogenerate noreport keep(master match)
save      `Histories', replace


use       `temp', clear
keep       caseid sUFiX hd* correct*

generate   tESt           = ""
foreach var of varlist _all {
	replace    tESt = "OK" if "`var'" != "caseid" & "`var'" != "sUFiX" & "`var'" != "tESt" & `var' != ""
	}
keep if    tESt          == "OK"
drop       tESt

destring   hd2* hd3* hd4* correct* hd5b* hd6* hd7* hd10* hd11* hd12*, replace
recode     hd2* hd3* hd4* ( . = 0) if sUFiX == "_1"
generate   under5         = hd2_number + hd2b_number + hd2c_number
replace    under5         = correct1                                                  if correct1   != .      & (hd3c == 2 | hd3c == 4)
generate   above5         = hd3_number + hd3b_number
replace    above5         = correct2                                                  if correct2   != .      & (hd3c == 3 | hd3c == 4)
generate   under5_m3      = min(hd3d2,under5)
generate   above5_m3      = min(hd3d4,above5)
save      `temp2', replace

generate   days_resident  = hd5b2                                                     if hd5b1      == 1
replace    days_resident  = hd5b3/12*365.25                                           if hd5b1      == 2
replace    days_resident  = hd5b4*365.25                                              if hd5b1      == 3
replace    days_resident  = hd10c                                                     if hd10b      == 1
replace    days_resident  = hd10d/12*365.25                                           if hd10b      == 2

generate   age_death      = hd6_days/365.25                                           if hd6        == 1
replace    age_death      = hd6_months/12                                             if hd6        == 2
replace    age_death      = hd6_years                                                 if hd6        == 3
replace    age_death      = hd11_days/365.25                                          if hd11       == 1
replace    age_death      = hd11_months/12                                            if hd11       == 2

generate   days_dead      = hd7_days                                                  if hd7        == 1
replace    days_dead      = hd7_months/12*365.25                                      if hd7        == 2
replace    days_dead      = hd12_days                                                 if hd12       == 1
replace    days_dead      = hd12_months/12*365.25                                     if hd12       == 2
replace    days_dead      = hd7_t2                                                    if days_dead  == .      & hd7_t2     != .


generate   hd8a_Reasons   = ""
generate   hd8b_Reasons   = ""    
forvalues i = 1(1)9 {
	capture replace    hd8a_Reasons = hd8a_Reasons + " " + "`i'" if hd8a`i' == "1"
	capture replace    hd8b_Reasons = hd8b_Reasons + " " + "`i'" if hd8b`i' == "1"
	}

replace    hd8a_Reasons   = hd8a                                                      if hd8a_Reaso	== ""     & hd8a       != ""
replace    hd8b_Reasons   = hd8b                                                      if hd8b_Reaso	== ""     & hd8b       != ""
drop       hd8a hd8b
rename     *_Reasons *
order      hd8 hd8a hd8aX hd8b hd8bX hd13

generate   CRVS           = "Registered"                                              if hd8        == "1"    | hd13       == "1"    
replace    CRVS           = "Not registered"                                          if hd8        == "2"    | hd13       == "2"
replace    CRVS           = "Don't know"                                              if hd8        == "98"   | hd13       == "98"
replace    CRVS           = "Refuse"                                                  if hd8        == "99"   | hd13       == "99"
rename    (hd8 hd8a hd8b hd8*X hd13) CRVS_=

generate   under5_d3      = 1                                                         if age_death   < 5
recode     under5_d3   (. = 0)                                                        if age_death  != .
recode     under5_d3   (. = 1)                                                        if hd11       != .
recode     under5_d3   (. = 0)                                                        if hd6        != .
keep if    CRVS          != ""
keep       caseid days_resident age_death days_dead CRVS* under5_d3 hd7* hd12*    
order      caseid days_resident age_death days_dead CRVS* under5_d3
generate   agegroup       = "A. under5"                                               if under5_d3  == 1
replace    agegroup       = "B. above5"                                               if under5_d3  == 0
generate   status         = "dead"
replace    status         = "dead migrant"                                            if days_resi   < 365.25/4
replace    status         = "outdated " + status                                      if days_dead  >= 365.25/4 & days_dead  != .
drop       under5_d3
save      `temp3', replace

use       `temp2', clear
keep       caseid under5 above5 under5_m3 above5_m3
keep if    under5        != .
generate   total          = under5 + above5
expand     total
drop if    total         == 0
bysort     caseid: generate   k = _n
generate   agegroup       = "A. under5"                                               if k          <= under5
replace    agegroup       = "B. above5"                                               if k           > under5
drop       k
bysort     caseid agegroup: generate   k = _n
generate   status         = "migrant"                                                 if k          <= under5_m3 & agegroup == "A. under5"
replace    status         = "resident"                                                if k           > under5_m3 & agegroup == "A. under5"
replace    status         = "migrant"                                                 if k          <= above5_m3 & agegroup == "B. above5"
replace    status         = "resident"                                                if k           > above5_m3 & agegroup == "B. above5"
keep       caseid agegroup status
append     using `temp3'

merge m:1  caseid using `sample', nogenerate keep(master match) keepusing(starttime)
sort       caseid agegroup status
generate   O              = mdy(month(dofc(starttime)),day(dofc(starttime)),year(dofc(starttime)))
generate   A              = O - 365.25/12*3
format     %tdDD/NN/CCYY A O
replace    O              = O - days_dead                                             if days_dead  != .
replace    A              = max(O - min(days_resident,age_death*365.25),A)            if days_resid != .
drop if    A              > O
drop       starttime
save      `HHDeaths', replace

use       `temp2', clear
keep       caseid under5 above5
keep if    under5        != . /*to identify those responding to the HHD questionnaire*/
generate   HH_size        = under5 + above5 + 1 /*family members plus the proxy*/
keep       caseid HH_size
generate   HH_sample      = 1
save      `temp2', replace

use       `sample', clear
merge 1:1  caseid using `temp2', nogenerate
list       caseid              if HH_size    == .
generate   household      = 1  if HH_size     < 5
replace    household      = 2  if HH_size     < 9  & house == .
recode     household   (. = 3) if HH_size    != .
drop       HH_size
generate   reproductive   = 1  if sex        == 2  & age    < 50
recode     reproductiv (. = 0)
save      `sample', replace


import     spss using "`pATh'/Malawi MICS6 SPSS/hh.sav", clear
keep       HH1 HH2 HH6 HH7 HH12 HH48 WS1 HC5 HC8
keep if    HH12          == 1 /*consented interviews*/

generate   UR             = HH6
generate   Region         = HH7
generate   Roofing        = 0     if HC5   != .
recode     Roofing     (0 = 1)    if HC5   == 31 | HC5   == 33 | HC5   == 34 | HC5   == 35 | HC5   == 36
generate   Electricity    = 0
recode     Electricity (0 = 1)    if HC8   ==  1 | HC8   ==  2
generate   Water          = 0     if WS1   != .
recode     Water       (0 = 1)    if WS1   == 11 | WS1   == 12 | WS1   == 13 | WS1   == 14 | WS1   == 21 | WS1   == 31 | WS1   == 41 | WS1   == 51 | WS1   == 71 | WS1   == 72 | WS1   == 91
generate   household      = 1     if HH48   < 5
replace    household      = 2     if HH48   < 9  & house == .
recode     household   (. = 3)
recode     Electricity Roofing Water (0 = 1) (1 = 2)
keep       HH1 HH2 UR Region Electricity Roofing Water household
save      `mics', replace

import     spss using "`pATh'/Malawi MICS6 SPSS/wm.sav", clear
keep if    WM9           == 1 /*consented interviews*/
keep       HH1 HH2 PSU stratum WDOB WM1 WM2 WM3 WM6D WM6M WM6Y WM9 WM17 WM7H WM7M WM10H WM10M WB3M WB3Y WB4 MT11 MT12 CM1 CM2 CM3 CM4 CM5 CM6 CM7 CM8 CM9 CM10 CM11 CM12 CM15 CM17 HH7 wmweight wscore WB6A WB6B welevel
generate   caseid         = _n

generate   Education      = 1
replace    Education      = 2                                                       if (WB6A == 1 & WB6B == 8) | WB6A == 2 | WB6A == 3
replace    Education      = 3                                                       if (WB6A == 3 & WB6B >= 4 & WB6B < 98)
replace    Education      = 3                                                       if  WB6A == 4 | WB6A == 5
 
rename     WB4 age
generate   interview      = mdy(WM6M,WM6D,WM6Y)
recode     WB3M       (99 = .) (98 = .)
recode     WB3Y     (9999 = .) (9998 = .)
replace    WB3M           = 1 + mod(WDOB - 1,12)                                    if WB3M    == .
replace    WB3Y           = 1900 + floor((WDOB - 1)/12)                             if WB3Y    == .

generate   DOB_min        = mdy(WB3M,1,WB3Y)
generate   DOB_max        = mdy(1 + mod(WB3M,12),1,WB3Y + floor(WB3M/12))
replace    DOB_min        = mdy(1,1,WB3Y)                                           if DOB_min == .
replace    DOB_max        = mdy(1,1,WB3Y + 1)                                       if DOB_max == .
replace    DOB_min        = mdy(1,1,WM6Y - age - 1)                                 if DOB_min == .
replace    DOB_max        = mdy(1,1,WM6Y - age)                                     if DOB_max == .
format     %tdDD/NN/CCYY interview DOB_*

rename     HH7 region
rename     wmweight W
keep if    W              > 0                
generate   mobile         = 0
recode     mobile      (0 = 1)                                                      if MT11 == 1
rename     PSU cluster
keep       caseid HH1 HH2 WM1 WM2 WM3 interview DOB_* Education age W mobile cluster stratum wscore
sort       WM1 WM2 WM3

sort       cluster caseid
bysort     cluster: generate   woman          = _n
bysort     cluster: generate   Women          = _N
generate   w              = Women                                                   if woman   == 1
generate   iNDeX          = sum(w)
replace    iNDeX          = iNDeX - Women
drop       w
egen       GO             = cut(age), at(18,30,40,50,65) icodes
replace    GO             = GO + 1
merge m:1  HH1 HH2 using `mics', nogenerate keep(master match)
drop       HH1 HH2
save      `mics', replace

import     spss using "`pATh'/Malawi MICS6 SPSS/bh.sav", clear
recode     BH4M BH4D  (99 = .) (98 = .)
recode     BH4Y     (9999 = .) (9998 = .)
replace    BH4M           = 1 + mod(BH4C - 1,12)                                    if BH4M    == .
replace    BH4Y           = 1900 + floor((BH4C - 1)/12)                             if BH4Y    == .

generate   B_min          = mdy(BH4M,BH4D,BH4Y)
replace    B_min          = mdy(1 + mod(BH4M,12),1,BH4Y + floor((BH4M + 1)/12)) - 1 if B_min   == . & BH4D    != .
generate   B_max          = B_min
replace    B_min          = mdy(BH4M,1,BH4Y)                                        if B_min   == . 
replace    B_min          = mdy(1,1,BH4Y)                                           if B_min   == .
replace    B_max          = mdy(1 + mod(BH4M,12),1,BH4Y + floor(BH4M/12))           if B_max   == .
replace    B_max          = mdy(1,1,BH4Y + 1)                                       if B_max   == .
format     %tdDD/NN/CCYY B_*

generate   D_min          = BH9N                                                     if BH9U == 1
generate   D_max          = BH9N + 1                                                 if BH9U == 1
replace    D_min          = BH9N*365.25/12                                           if BH9U == 2
replace    D_max          = (BH9N + 1)*365.25/12                                     if BH9U == 2
replace    D_min          = BH9N*365.25                                              if BH9U == 3
replace    D_max          = (BH9N + 1)*365.25                                        if BH9U == 3

replace    D_min          = BH9C*365.25/12                                           if BH9U == 9
replace    D_max          = (BH9C + 1)*365.25/12                                     if BH9U == 9
rename     BH3 sex
rename     BH2 multiple
generate   mother         = 1
sort       WM1 WM2 WM3 B_min
bysort     WM1 WM2 WM3:  generate   bidx      = _n
keep       WM1 WM2 WM3 B_* D_* sex multiple mother bidx
merge m:1  WM1 WM2 WM3 using `mics', nogenerate
sort       caseid B_min
bysort     caseid:  generate   k              = _n
bysort     caseid:  generate   K              = _N 
drop       WM1 WM2 WM3
order      cluster caseid 
label drop _all
sort       cluster caseid k
save      `mics', replace

use       "`pATh'/MWIR7AFL.DTA", clear
keep       caseid v001 v002 v003 v005 v006 v007 v016 v009 v010 v012 v023 v024 v025 v169a vcal_1 v008 v018
save      `temp', replace

use       "`pATh'/MWBR7AFL.DTA", clear
keep       caseid v001 v002 v003 v006 v007 v016 v008 v009 v010 v012 v023 v024 v025 bidx b0 b1 b2 b3 b4 b6 b7 b17
generate   mother         = 1
merge m:1  caseid v001 v002 v003 v006 v007 v016 v008 v009 v010 v012 v023 v024 v025 using `temp', nogenerate
sort       caseid bidx
bysort     caseid: generate   k = _n
bysort     caseid: generate   K = _N

generate   cluster        = v001
generate   HH             = v002
generate   respondent     = v003

sort       cluster HH caseid k
generate   woman          = 1                                               if k     == 1
bysort     cluster: egen       Women          = sum(woman)
bysort     cluster: replace    woman          = sum(woman)
generate   w              = Women                                           if k     == 1 & woman == 1
generate   iNDeX          = sum(w)
replace    iNDeX          = iNDeX - Women
drop       w

generate   DOB            = mdy(v009,1,v010)
generate   interview      = mdy(v006,v016,v007)
generate   Birth          = mdy(b1,b17,b2)
format     %tdDD/NN/CCYY interview Birth DOB
rename     b4 sex
rename     v012 age
rename     v169a mobile

generate   D_min          = b6 - 100                                                 if b6   != .   & b6 <= 200
generate   D_max          = D_min + 1                                                if b6   != .   & b6 <= 200
replace    D_min          = 0                                                        if b6   == 198 | b6 == 199  
replace    D_max          = 365.25/12                                                if b6   == 198 | b6 == 199	
replace    D_min          = (b6 - 200)*365.25/12                                     if b6   != .   & b6 >= 200 & b6  < 300 
replace    D_max          = (b6 - 200 + 1)*365.25/12                                 if b6   != .   & b6 >= 200 & b6  < 300
replace    D_min          = 0                                                        if b6   == 298 | b6 == 299  
replace    D_max          = 24*365.25/12                                             if b6   == 298 | b6 == 299
replace    D_min          = (b6 - 300)*365.25                                        if b6   != .   & b6 >= 300 
replace    D_max          = (b6 - 300 + 1)*365.25                                    if b6   != .   & b6 >= 300
replace    D_min          = 0                                                        if b6   == 398 | b6 == 399  
replace    D_max          = max(year(interview) - year(Birth),0)*365.25              if b6   == 398 | b6 == 399

generate   W              = v005/1000000
generate   CAL            = "C" + vcal_1                                             if k    == 1 | k  == .
generate   row            = v018                                                     if k    == 1 | k  == .
keep       cluster HH respondent age sex interview Birth D_* mobile W bidx caseid DOB k K woman Women iNDeX mother CAL row
order      cluster HH respondent age sex interview Birth D_* mobile W bidx caseid DOB k K woman Women iNDeX mother CAL row
generate   reproductive = 1
save      `dhs', replace

contract   caseid cluster HH respondent age W mobile reproductive
keep       caseid cluster HH respondent age W mobile reproductive
save      `temp2', replace


use       "`pATh'/MWPR7AFL.DTA", clear
keep       hhid hvidx hv000 hv001 hv002 hv003 hv004 hv005 hv006 hv007 hv009 hv012 hv013 hv016 hv023 hv024 hv025 hv101 hv102 hv103 hv104 hv105 hv106 hv107 hv108 hv109 hv140 hv201 hv206 hv215 hv243a
generate   cluster        = hv001
generate   HH             = hv002
generate   respondent     = hvidx
merge 1:1  cluster HH respondent using `temp2', nogenerate
replace    hv105          = age    if age    != .
replace    hv243a         = mobile if mobile != .
drop       age mobile

sort       cluster HH hhid hvidx
bysort     hhid:   generate   k = _n
bysort     hhid:   generate   K = _N
generate   h              = 1                                               if hvidx == 1
bysort     cluster: egen       H              = sum(h)
bysort     cluster: replace    h              = sum(h)
generate   s              = H                                               if hvidx == 1 & h    == 1
generate   iNDeX          = sum(s)
replace    iNDeX          = iNDeX - H
drop       s

generate   interview      = mdy(hv006,hv016,hv007)
format     %tdDD/NN/CCYY interview
generate   mobile         = hv243a
generate   sex            = hv104
generate   Region         = hv024

tabulate   hv109
generate   Education      = hv109
recode     Education   (8 = 1)
recode     Education   (0 = 1) (3 = 2) (4 = 3) (5 = 3)

generate   ageG           = 0     if hv105  < 18
replace    ageG           = 1     if hv105  < 50 & ageG  == .
replace    ageG           = 2     if hv105  < 65 & ageG  == .
generate   age            = hv105
replace    ageG           = 3     if ageG  == .
generate   UR             = hv025
generate   Electricity    = hv206
generate   Roofing        = 0
recode     Roofing     (0 = 1)    if hv215 == 31 | hv215 == 33 | hv215 == 34 | hv215 == 35 | hv215 == 36 /* good material excluding wood*/
generate   Water          = 0
recode     Water       (0 = 1)    if hv201 == 11 | hv201 == 12 | hv201 == 13 | hv201 == 14 | hv201 == 21 | hv201 == 31 | hv201 == 41 | hv201 == 51 | hv201 == 62 | hv201 == 71
generate   jure           = hv102
bysort     hhid: egen   HH_size = sum(jure)
generate   Wh             = hv005/HH_size/1000000
sum        Wh
replace    Wh             = Wh/r(mean)
generate   household      = 1     if HH_si  < 5
replace    household      = 2     if HH_si  < 9  & house == .
recode     household   (. = 3)
egen       GO             = cut(age), at(18,30,40,50,65) icodes
replace    GO             = GO + 1
recode     Electricity Roofing Water (0 = 1) (1 = 2)
keep       cluster HH respondent hhid age* sex UR Region household Education Electricity Roofing Water mobile W* interview caseid jure reproductive GO HH_size k K h H iNDeX
order      cluster HH respondent hhid age* sex UR Region household Education Electricity Roofing Water mobile W* interview caseid jure reproductive GO HH_size k K h H iNDeX
sort       cluster HH hhid k
save      `temp2', replace

drop if    ageG          == 0
drop if    ageG          == 3
keep if    jure          == 1
generate   Total          = 1
tabulate   Total [aw = Wh], matcell(Total)
local      lISt           = "GO sex UR Region household Education Electricity"

local      variable       = ""
local      totals         = "_cons = 10000"
foreach var of local lISt {
	local            variable       = "`variable'" + " i.`var'"
	tabulate        `var' [aw = Wh], matcell(`var')
	forvalues i = 1(1)`= rowsof(`var')' {
		local            number         = `var'[`i',1]/Total[1,1]*10000
		local            totals         = "`totals'" + " `i'.`var' = `number'"
		}	
	}

keep if    reproductive  == 1
tabulate   Total [aw = W], matcell(Total)
local      rep_lISt       = "GO UR Region household Education Electricity"

local      rep_variable   = ""
local      rep_totals     = "_cons = 10000"
foreach var of local rep_lISt {
	local            rep_variable   = "`rep_variable'" + " i.`var'"
	tabulate        `var' [aw = W], matcell(`var')
	forvalues i = 1(1)`= rowsof(`var')' {
		local            number         = `var'[`i',1]/Total[1,1]*10000
		local            rep_totals     = "`rep_totals'" + " `i'.`var' = `number'"
		}	
	}

use       `temp2', clear
keep if    mobile        == 1
drop if    ageG          == 0
drop if    ageG          == 3
keep if    jure          == 1
local      alpha          = 0.10		
		
svycal     rake `rep_variable' if reproductive == 1 [pw = W], force generate(temp) totals(`rep_totals')
xtile      Q_temp         = temp, nq(200)
replace    Q_temp         = min(max(Q_temp/2,100*`alpha'/2),100*(1 - `alpha'/2))
bysort     Q_temp: egen max = max(temp)
bysort     Q_temp: egen min = min(temp)
egen       LB             = min(max)
egen       UB             = max(min)
*replace    temp          = min(max(temp,LB),UB) /**/
*replace    temp          = min(temp,UB) /**/
sum        temp
generate   WR             = temp/r(mean) if reproductive == 1
drop       Q_temp LB UB max min temp
	
svycal     rake `variable' [pw = Wh], force generate(temp) totals(`totals')
xtile      Q_temp         = temp, nq(200)
replace    Q_temp         = min(max(Q_temp/2,100*`alpha'/2),100*(1 - `alpha'/2))
bysort     Q_temp: egen max = max(temp)
bysort     Q_temp: egen min = min(temp)
egen       LB             = min(max)               
egen       UB             = max(min)               
*replace    temp          = min(max(temp,LB),UB) /**/
*replace    temp          = min(temp,UB) /**/      
sum        temp
generate   WT             = temp/r(mean)
drop       Q_temp LB UB max min temp
keep       cluster HH respondent WR WT
merge 1:1  cluster HH respondent using `temp2', nogenerate
order      cluster HH respondent age* GO sex jure UR Region household Education Electricity Roofing Water mobile W WR WT interview 
sort       cluster HH hhid k
export     delimited using "`pATh'/DHSmalawidst.csv", replace
keep       cluster HH respondent GO jure UR Region household Education Electricity Roofing Water W WR WT HH_size
save      `temp2', replace

use       `dhs', clear
merge m:1  cluster HH respondent using `temp2', nogenerate keep(master match)
label drop _all
sort       cluster HH caseid k
export     delimited using "`pATh'/DHSmalawi.csv", replace


use       `mics', clear
keep if    k             == 1
keep if    mobile        == 1
svycal     rake `rep_variable' [pw = W], force generate(temp) totals(`rep_totals')
xtile      Q_temp         = temp, nq(200)
replace    Q_temp         = min(max(Q_temp/2,100*`alpha'/2),100*(1 - `alpha'/2))
bysort     Q_temp: egen max = max(temp)
bysort     Q_temp: egen min = min(temp)
egen       LB             = min(max)
egen       UB             = max(min)
*replace    temp          = min(max(temp,LB),UB) /**/
*replace    temp          = min(temp,UB) /**/
sum        temp
generate   WR             = temp/r(mean)
drop       Q_temp LB UB max min temp
keep       caseid WR
save      `temp2', replace

use       `mics', clear
merge m:1  caseid using `temp2', nogenerate
sort       cluster caseid k
export     delimited using "`pATh'/MICSmalawi.csv", replace

use       `Histories', clear
contract   caseid sample
generate   Hist_sample    = 1
keep       caseid sample Hist_sample
save      `temp2', replace

use       `sample', clear
merge 1:1  caseid using `temp2', nogenerate
generate   missing        = "Electricity "            if Electricity == .
replace    missing        = "HH_size "     + missing  if household   == .
replace    missing        = "Education "   + missing  if Education   == .

egen       GO             = cut(age), at(18,30,40,50,65) icodes
replace    GO             = GO + 1
generate   interview      = dofc(starttime)
format     %tdDD/NN/CCYY interview

generate   temp           = mdy(month(dofc(starttime)),day(dofc(starttime)),year(dofc(starttime)))
format     %tdDD/NN/CCYY temp
generate   group          = 1  if temp <= mdy(5,25,2022)
recode     group       (. = 2) if temp <= mdy(9,13,2022)
recode     group       (. = 3) if temp <= mdy(3,20,2023)
recode     group       (. = 4)
*replace    group          = 1                                       /*Activate not to post-stratify by waves*/ 
drop       temp
save      `temp', replace

foreach var of local lISt {
	drop if    `var' == .
	tabulate   `var'
	}

foreach var of local rep_lISt {
	tabulate   `var' if reproductive == 1
	}
	
generate   WR             = .
generate   WT             = .

forvalues i = 1(1)4 { /*post-stratification by block and sample (all participants vs women in reproductive ages)*/
	svycal     rake `rep_variable' if group == `i' & reproductive == 1, force generate(temp) totals(`rep_totals')
	xtile      Q_temp   = temp, nq(200)
	replace    Q_temp   = min(max(Q_temp/2,100*`alpha'/2),100*(1 - `alpha'/2))  if group == `i' & reproductive == 1
	bysort     Q_temp: egen max = max(temp)
	bysort     Q_temp: egen min = min(temp)
	egen       LB       = min(max)                 if group == `i' & reproductive == 1
	egen       UB       = max(min)                 if group == `i' & reproductive == 1
	*replace    temp    = min(max(temp,LB),UB) /**/
	*replace    temp    = min(temp,UB) /**/
	sum        temp
	replace    WR       = temp/r(mean)             if group == `i' & reproductive == 1
	drop       Q_temp LB UB max min temp
	
	svycal     rake `variable' if group == `i', force generate(temp) totals(`totals')
	xtile      Q_temp   = temp, nq(200)
	replace    Q_temp   = min(max(Q_temp/2,100*`alpha'/2),100*(1 - `alpha'/2))  if group == `i'
	bysort     Q_temp: egen max = max(temp)
	bysort     Q_temp: egen min = min(temp)
	egen       LB       = min(max)                 if group == `i'
	egen       UB       = max(min)                 if group == `i'
	*replace    temp    = min(max(temp,LB),UB) /**/
	*replace    temp    = min(temp,UB) /**/
	sum        temp
	replace    WT       = temp/r(mean)             if group == `i'
	drop       Q_temp LB UB max min temp
	}

merge 1:1  caseid using `temp', nogenerate
drop       SubmissionDate starttime endtime
destring   *_duration, replace
recode     *_duration (0  = .)
export     delimited using "`pATh'/RaMMPSdst.csv", replace
keep       caseid WR WT Region UR household Education* Electricity Roofing Water group interview kinship age ageG sex GO *duration missing
save      `temp', replace
keep       caseid WT WR sex age
save     "`pATh'/weights.dta", replace

use       `Histories', clear
merge m:1  caseid using `temp', nogenerate keep(master match) force

local      N              = `=_N' + 1
set obs   `N'
foreach var of varlist _all {
	if substr("`: type `var''",1,3) == "str" {
		quiet replace   `var' = "0A0" in `N'
		}
	else {
		quiet replace   `var' = 0    in `N'
		}
	}
replace    B_max           = min(B_max,interview)    if B_max             != .
replace    B_max           = .                       if (B_min + B_max)/2 >= interview
replace    B_min           = .                       if B_max             == .	
sort       caseid k K
export     delimited using "`pATh'/RaMMPS.csv", replace
drop if    caseid         == "0A0"
save     "`pATh'/Pregnancy Histories.dta", replace

use       `HHDeaths', clear
bysort     caseid: generate   k = _n
bysort     caseid: generate   K = _N
merge m:1  caseid using `temp', nogenerate keep(master match)

drop       WR
local      N              = `=_N' + 1
set obs   `N'
foreach var of varlist _all {
	if substr("`: type `var''",1,3) == "str" {
		quiet replace   `var' = "0A0" in `N'
		}
	else {
		quiet replace   `var' = 0    in `N'
		}
	}
order      caseid k K
sort       caseid k K
export     delimited using "`pATh'/RaMMPSHH.csv", replace




local      pATh           = "/Users/lshjr3/Documents/RaMMPS/under-five mortality"
tempfile   temp
use      "/Users/lshjr3/Desktop/RaMMPS/MW_AllCallAttempts.dta", clear
generate   temp           = "A" + string(caseid)
drop       caseid
rename     temp caseid
destring   ivr, replace
recode     ivr         (. = 0)
bysort     caseid: egen   IVR = max(ivr)

generate   outcome        = ""
generate   date           = ""
generate   attempts       = .
forvalues i = 1(1)6 { /*to delete duplicated records of multiple attempts and partial interviews*/
	replace    attempts  = `i'                if attempt_`i'_status != "" & outcome   == ""
	replace    date      = attempt_`i'_time   if attempt_`i'_status != "" & outcome   == ""
	replace    outcome   = attempt_`i'_status if attempt_`i'_status != "" & outcome   == ""
	}

generate   temp_date      = clock(date,"YMD hms")
format     %tcMon_dd,_CCYY_hh:MM:SS_AM temp_date
drop       date
rename     temp_date date	

keep       caseid IVR outcome attempts date SubmissionDate starttime endtime C1_1
replace    date           = endtime if date == .

generate   calloutcome    = .
replace    calloutcome    =  1 if outcome     == "Completed interview" & C1_1 == 1
replace    calloutcome    =  2 if outcome     == "Partially complete (no callback)"
replace    calloutcome    =  2 if outcome     == "Incomplete (callback)"
replace    calloutcome    =  3 if outcome     == "Refusal"
replace    calloutcome    =  4 if outcome     == "Answered, but not by the respondent"
replace    calloutcome    =  4 if outcome     == "Busy tone"
replace    calloutcome    =  4 if outcome     == "No answer"
replace    calloutcome    =  4 if outcome     == "Number not accessible (number exists, but (temporarily) not reachable)"
replace    calloutcome    =  5 if outcome     == "Number not in use (not registered on the network)"
replace    calloutcome    =  6 if outcome     == "Already interviewed (contacted on another SIM)"
replace    calloutcome    =  6 if outcome     == "Ineligible(no refer/ defer)"
replace    calloutcome    =  6 if outcome     == "Respondent unavailable (Quota full cases)"
replace    calloutcome    =  6 if outcome     == "Quota full (available /Defer)"
replace    calloutcome    =  6 if outcome     == "Quota full (ineligible)"
replace    calloutcome    =  7 if outcome     == "Incomprehensible (Language issue)"
replace    calloutcome    =  7 if outcome     == "Incomprehensible (Technical issue)"
replace    calloutcome    =  7 if outcome     == "Reassigned (Need female interviewer)"
replace    calloutcome    =  7 if outcome     == "Reassigned to other enumerator (General)"
replace    calloutcome    =  7 if outcome     == "Deferral"
replace    calloutcome    =  7 if outcome     == "Referral"
replace    calloutcome    =  7 if outcome     == "Other (specify)"
replace    calloutcome    =  7 if outcome     == ""
replace    calloutcome    =  3 if calloutcome == .                     & C1_1 != 1
tabulate   calloutcome, m
drop if    calloutcome   == .

generate   temp           = 1  if calloutcome == 1
recode     temp        (. = 0)
gsort      caseid date - calloutcome
bysort     caseid: generate   CC         = sum(temp)
replace    CC             = CC - temp + 1
drop       temp
gsort      caseid CC date - calloutcome
bysort     caseid CC: generate   k = _n
bysort     caseid CC: generate   K = _N
label      define calloutcome 1 "Complete CATI" 2 "Partial CATI" 3 "Refusal" 4 "Number Not Accessible (i.e., no answer, busy tone, not the respondent)" 5 "Number Not in Use (i.e., not registered on the network)" 6 "Ineligible (e.g., quota full, already interviewed using another SIM)" 7 "Other (e.g., technical difficulties, reassigned)"

generate   totalcases     = 1           if k == K
generate   catioutcome    = calloutcome if k == K
generate   catiOUTcome    = catioutcome
label      values catioutcome calloutcome calloutcome
numlabel   _all, add

destring   caseid, replace ignore("A")
save      `temp', replace

preserve
generate   complete       = 1 if catioutcome == 1
duplicates tag caseid         if complete    == 1, gen(dup)
generate   refusal        = 1 if catioutcome == 3
generate   eligible       = 1 if catioutcome == 1 | catioutcome ==  2 | catioutcome ==  3 | catioutcome == 4 | catioutcome == 7


collapse  (sum) complete refusal eligible (count) calls = caseid (sum) totalcases
display    in yellow  "calls: " calls
display    in yellow  "cases:" totalcases
display    in yellow  "response rate: " complete/eligible*100
display    in yellow  "refusal rate: " refusal/eligible*100
display    in yellow  "calls per CATI: " calls/complete
display    in yellow  "% of numbers with CATI: " complete/(totalcases)*100
restore

generate   temp           = "A" + string(caseid)
drop       caseid
rename     temp caseid

local      pATh           = "/Users/lshjr3/Documents/RaMMPS/under-five mortality"
merge m:1  caseid using "`pATh'/weights.dta", keepusing(WT WR sex age)

replace    caseid         = "B" + substr(caseid,2,.) if _merge == 1
replace    caseid         = caseid + "-" + string(CC)
drop       _merge CC

order      caseid
sort       caseid date
format     %-tcDD/NN/CCYY_HH:MM:SS SubmissionDate starttime endtime date
save      `temp', replace

keep if    K             == k
keep       caseid starttime 
generate   temp           = mdy(month(dofc(starttime)),day(dofc(starttime)),year(dofc(starttime)))
format     %tdDD/NN/CCYY temp
generate   group          = 1  if temp <= mdy(5,25,2022)
recode     group       (. = 2) if temp <= mdy(9,13,2022)
recode     group       (. = 3) if temp <= mdy(3,20,2023)
recode     group       (. = 4)
keep       caseid group
merge 1:m  caseid using `temp', nogenerate
sort       caseid k K
export     delimited using "`pATh'/RaMMPScalls.csv", replace

