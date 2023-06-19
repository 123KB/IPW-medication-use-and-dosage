
* Estimating the comparative effectiveness of dynamic treatment strategies for medication use and dosage: application of marginal structural models to emulate a hypothetical target trial using observational data

* Sample code for calculating IPW using different methods

* METHODS 
* A Logistic regression for zero dose and standard linear regression for log dose
* B Logistic regression for zero dose and heteroscedastic linear regression for log dose
* C Logistic regression for zero dose, heteroscedastic linear regression for log dose and multinomial regression for patients who recently received very low or high doses
* D Ordinal logistic regression


****************************************************************************************************************************
****************************************************************************************************************************
* METHOD A Logistic regression for zero dose and standard linear regression for log dose
****************************************************************************************************************************
****************************************************************************************************************************

****************************************************************************************************************************
* LOGISTIC REGRESSION MODELS
****************************************************************************************************************************

* (1)

****************************************************************************************************************************
* Model for those on treatment, coming off (elig_model_on_treat==1)
****************************************************************************************************************************
 
logistic zero_dose i.baseline_age_grp etc timespl* if elig_model_on_treat==1

* Obtain probability of zero dose (i.e. coming off EPO)

predict prob_zero_dose_thoseon, pr 

* Note we could have created a variable for non-zero dose where the predicted probability would be for not receiving darbepoetin 


* (2)

***************************************************************************************************************************
* Model for those off treatment, staying off (elig_model_off_treat==1)
***************************************************************************************************************************

logistic zero_dose i.baseline_age_grp etc timespl* if elig_model_off_treat==1

* Obtain probability of zero dose (i.e. staying off EPO)
predict prob_zero_dose_thoseoff, pr 


****************************************************************************************************************************
* Overall probabiltiy of zero dose

gen prob_zero_dose = prob_zero_dose_thoseon if elig_model_on_treat==1
replace prob_zero_dose = prob_zero_dose_thoseoff if elig_model_off_treat==1


****************************************************************************************************************************


* (3)

***************************************************************************************************************************
* Normal linear regression regression for non-zer doses
***************************************************************************************************************************

regress log_darb i.baseline_age_grp etc timespl* 

* Test for evidence of heteroskedasticity
estat hettest

***********************************************************************
* Obtain the predicted value
predict xb 
predict resid, residuals
sum resid, det
di e(rmse)

*************************************************************
* Calculate strategy 1 probability in range (note: min and max of accceptable range for each strategy is calculated from the dose change decisions protocol and acceptable dose changes table)
gen strat1_lstd = (log_strat1_min - xb ) / e(rmse)
gen strat1_ustd = (log_strat1_max - xb) / e(rmse)

gen strat1_lnorm = normal(strat1_lstd)
gen strat1_unorm = normal(strat1_ustd)

gen strat1_probint = strat1_unorm - strat1_lnorm


*************************************************************
* Calculate strategy 2 probability in range
gen strat2_lstd = (log_strat2_min - xb) / e(rmse) 
gen strat2_ustd = (log_strat2_max - xb) / e(rmse)

gen strat2_lnorm = normal(strat2_lstd)
gen strat2_unorm = normal(strat2_ustd)

gen strat2_probint = strat2_unorm - strat2_lnorm


****************************************************************************************************************************
* Overall probability of adhering for each strategy
* Combine the probability of non-zero dose (from the logistic model) AND the probability of being in the range from the linear regression model

* Strategy 1
* If the strategy says zero dose, pr(adhere) = pr(zero dose)
gen strat1_prob_adhere = prob_zero_dose if strat1_min==0 & strat1_max==0
* If the strategy says 0 or range, pr(adhere) = pr(zero) + [pr(non-zero) x pr(in range)]
replace strat1_prob_adhere = prob_zero_dose + ((1-prob_zero_dose)*strat1_probint) if strat1_min==0 & strat1_max>0 & strat1_max<. & strat1_prob_adhere==. 
* If the strategy says non-zero only, pr(adhere) = [pr(non-zero) x pr(in range)]
replace strat1_prob_adhere = ((1-prob_zero_dose)*strat1_probint) if strat1_min>0 & strat1_min<. & strat1_max>0 & strat1_max<. & strat1_prob_adhere==. 


* Strategy 2
* If the strategy says zero dose, pr(adhere) = pr(zero dose)
gen strat2_prob_adhere = prob_zero_dose if strat2_min==0 & strat2_max==0
* If the strategy says 0 or range, pr(adhere) = pr(zero) + [pr(non-zero) x pr(in range)]
replace strat2_prob_adhere = prob_zero_dose + ((1-prob_zero_dose)*strat2_probint) if strat2_min==0 & strat2_max>0 & strat2_max<. & strat2_prob_adhere==. 
* If the strategy says non-zero only, pr(adhere) = [pr(non-zero) x pr(in range)]
replace strat2_prob_adhere = ((1-prob_zero_dose)*strat2_probint) if strat2_min>0 & strat2_min<. & strat2_max>0 & strat2_max<. & strat2_prob_adhere==. 



****************************************************************************************************************************
* Calculate weights

* Single for each time point, useful for later when want to investigate large weights
gen strat1_weight_single = 1 / strat1_prob_adhere
gen strat2_weight_single = 1 / strat2_prob_adhere

* Cumulative version
* Estimated probability of their complete history up to each month
gen strat1_weight = 1 / strat1_prob_adhere
gen strat2_weight = 1 / strat2_prob_adhere
sort id day
* Day 0, by definition adhere = 1, so weight = 1
by id: replace strat1_weight = 1 if _n ==1
by id: replace strat2_weight = 1 if _n ==1
by id: replace strat1_weight=strat1_weight*strat1_weight[_n-1] if _n!=1
by id: replace strat2_weight=strat2_weight*strat2_weight[_n-1] if _n!=1








****************************************************************************************************************************
****************************************************************************************************************************
* METHOD B Logistic regression for zero dose and heteroscedastic linear regression for log dose
****************************************************************************************************************************
****************************************************************************************************************************

* Logistic regression modesl as per Method A AND...


****************************************************************************************************************************
* Heteroskedastic linear regression 
****************************************************************************************************************************

* Note command doesn't let you choose baseline categories e.g. b3. so need to use i.
hetregress log_darb i(.baseline_age_grp etc) timespl* het(i.baseline_age_grp etc)
predict xb_het 
predict sigma , sigma


*************************************************************

* Calculate strategy 1 probability in range
gen strat1_lstd = (log_strat1_min - xb_het ) / sigma
gen strat1_ustd = (log_strat1_max - xb_het) / sigma

gen strat1_lnorm = normal(strat1_lstd)
gen strat1_unorm = normal(strat1_ustd)

gen strat1_probint = strat1_unorm - strat1_lnorm

*************************************************************

* Calculate strategy 2 probability in range, combine probabilites and caluculate weights as per Method A





****************************************************************************************************************************
****************************************************************************************************************************
* METHOD C Logistic regression for zero dose, heteroscedastic linear regression for log dose and multinomial regression for patients who recently received very low or high doses
****************************************************************************************************************************
****************************************************************************************************************************

* Logistic regression modesl as per Method A AND Heteroskedastic linear regression as per Method B AND

* For example, from a darbepoetin dose of 2.5 mcg/week

* 8 mutually-exclusive groups G for an individual's dose change, for example: (0) go off darbepoetin (i.e. move to zero dose), 
(1) unacceptable decrease in dose for both strategies (i.e. the patient's darbepoetin dose was lowered when the protocol for both strategies said the dose should be constant or increased), 
(2) acceptable decrease for low Hb strategy only (for when the protocol would only say to lower the darbepoetin dose for the low Hb strategy), 
(3) acceptable decrease both strategies, (4) keep the darbepoetin dose constant, (5) acceptable increase for both strategies, (6) acceptable increase for high Hb strategy only, and (7) unacceptable increase for both strategies. 

* Small N in model. Use hb in g/dl and spines with just one knot so the model can cope
mlogit from_darb_2_5 hb_dl_spl1 lag_hb_dl_spl1 timespl* , base(4) rrr

* Predict probabilites, for where categories existed for this variable (from a dose of 2.5)
tab from_darb_2_5
predict p2_5_0 p2_5_4 p2_5_5 p2_5_6 p2_5_7 

*************************************************************
* Calculate probabilities
* Combine the probability of going off darb, staying the same, acceptable increase etc, as appropriate

* Strategy 1

* If the strategy says zero dose, pr(adhere) = pr(zero dose)
gen s1pa = p2_5_0 if (strat1_min==0 & strat1_max ==0) & from_darb_2_5 <. 

* If the strategy says stay the same, pr(adhere) = pr(stay same)
replace s1pa = p2_5_4 if (strat1_min ==2.5 & strat1_max ==2.5) & from_darb_2_5 <. & s1pa==.

* If the strategy says 0 or range above 2.5, pr(adhere) = pr(zero dose) + pr(stay same) + pr(acceptable increase both)
replace s1pa = p2_5_0 + p2_5_4 + p2_5_5 if (strat1_min ==0 & strat1_max >2.5) & from_darb_2_5 <. & s1pa==.

* If the strategy says increase, pr(adhere) = pr(accept inc both)
replace s1pa = p2_5_5 if (strat1_min >2.5 & strat1_max >2.5) & from_darb_2_5 <. & s1pa==.

* If the strategy says stay same or increase, pr(adhere) = pr(stay same) + pr(acceptable increase both)
replace s1pa = p2_5_4 + p2_5_5 if (strat1_min >=2.5 & strat1_max >2.5) & from_darb_2_5 <. & s1pa==.


* Strategy 2

* If the strategy says zero dose, pr(adhere) = pr(zero dose)
gen s2pa = p2_5_0 if (strat2_min ==0 & strat2_max ==0) & from_darb_2_5 <.

* If the strategy says stay the same, pr(adhere) = pr(stay same)
replace s2pa = p2_5_4 if (strat2_min ==2.5 & strat2_max ==2.5) & from_darb_2_5 <. & s2pa==.

* If the strategy says 0 or range above 2.5, pr(adhere) = pr(zero dose) + pr(stay same) + pr(acceptable increase both) + pr(accept inc 2 only)
replace s2pa = p2_5_0 + p2_5_4 + p2_5_5 + p2_5_6 if (strat2_min ==0 & strat2_max >2.5) & from_darb_2_5 <. & s2pa==.

* If the strategy says increase, pr(adhere) = pr(accept inc both) + pr(accept inc 2 only)
replace s2pa = p2_5_5 + p2_5_6 if (strat2_min >2.5 & strat2_max >2.5) & from_darb_2_5 <. & s2pa==.

* If the strategy says stay same or increase, pr(adhere) = pr(stay same) + pr(acceptable increase both) + pr(accept inc 2 only)
replace s2pa = p2_5_4 + p2_5_5 + p2_5_6 if (strat2_min >=2.5 & strat2_max >2.5) & from_darb_2_5 <. & s2pa==.


*************************************************************

* Repeat for other extreme dose levels (2.5, 5, 120 and 150 mcg/week)

* Calculate overall probabilites as per above, but replace with probabilites from mlogit models for the extreme doses
* Caluculate weights as above



****************************************************************************************************************************
****************************************************************************************************************************
* METHOD D Ordinal logistic regression
****************************************************************************************************************************
****************************************************************************************************************************

* Finally, in Method D, we transformed dose into an ordinal variable V, with 17 levels (coded 0-16) to represent the dosing ladder: 
* 0, 2.5, 5, 10, 15, 20, 25, 30, 40, 50, 60, 70, 80, 100, 120, 150, 180 mcg/week. 
* Doses between ladder rungs were coded as the higher rung e.g., a dose of 35 was coded as 40 mcg/week.


version 14.2: quietly ologit dose_ladder baseline_age_grp etc timespl* 
brant

ologit dose_ladder i.baseline_age_grp etc timespl* , or

predict p* 

rename p1 p0
rename p2 p1
rename p3 p2
rename p4 p3
rename p5 p4
rename p6 p5
rename p7 p6
rename p8 p7
rename p9 p8
rename p10 p9
rename p11 p10
rename p12 p11
rename p13 p12
rename p14 p13
rename p15 p14
rename p16 p15
rename p17 p16


****************************************************************************************************************************
* Overall probability of adhering for each strategy

* Strategy 1
gen strat1_prob_adhere = 0
replace strat1_prob_adhere = strat1_prob_adhere + p0 if strat1_min <=0 & strat1_max >=0
replace strat1_prob_adhere = strat1_prob_adhere + p1 if strat1_min <=2.5 & strat1_max >=2.5
replace strat1_prob_adhere = strat1_prob_adhere + p2 if strat1_min <=5 & strat1_max >=5
replace strat1_prob_adhere = strat1_prob_adhere + p3 if strat1_min <=10 & strat1_max >=10
replace strat1_prob_adhere = strat1_prob_adhere + p4 if strat1_min <=15 & strat1_max >=15
replace strat1_prob_adhere = strat1_prob_adhere + p5 if strat1_min <=20 & strat1_max >=20
replace strat1_prob_adhere = strat1_prob_adhere + p6 if strat1_min <=25 & strat1_max >=25
replace strat1_prob_adhere = strat1_prob_adhere + p7 if strat1_min <=30 & strat1_max >=30
replace strat1_prob_adhere = strat1_prob_adhere + p8 if strat1_min <=40 & strat1_max >=40
replace strat1_prob_adhere = strat1_prob_adhere + p9 if strat1_min <=50 & strat1_max >=50
replace strat1_prob_adhere = strat1_prob_adhere + p10 if strat1_min <=60 & strat1_max >=60
replace strat1_prob_adhere = strat1_prob_adhere + p11 if strat1_min <=70 & strat1_max >=70
replace strat1_prob_adhere = strat1_prob_adhere + p12 if strat1_min <=80 & strat1_max >=80
replace strat1_prob_adhere = strat1_prob_adhere + p13 if strat1_min <=100 & strat1_max >=100
replace strat1_prob_adhere = strat1_prob_adhere + p14 if strat1_min <=120 & strat1_max >=120
replace strat1_prob_adhere = strat1_prob_adhere + p15 if strat1_min <=150 & strat1_max >=150


*************************************************************

* Repeat for strategy 2
* Calculate weights as above













 
