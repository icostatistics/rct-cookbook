
/****************
Setup the randomisation with strata and blocks
*****************/
ralloc block bsize trt, nsubj(50)  sav(rct)  ntreat(2) ratio(1) osize(2) init(4) seed(1913) strata(4)

/****************
Clean up and add repeated time points
*****************/

label drop _all
gen pid = _n
replace trt=trt-1
label def trt 0 "Placebo" 1 "Active"
label val trt trt
rename StratID site
drop block bsize Seq*
order pid site trt
sort pid
expand 4
sort pid
by pid: gen time=_n-1


/********************
Simulate the normaly distributed continuous endpoints
********************/

/**************
Generate the random intercept and residual
***************/
gen incept0 = rnormal(0,1) if time == 0
by pid: egen incept = max(incept0)
drop incept0

gen resid = rnormal(0,1)


/******************
Set the mean by time and treatment
*******************/


generate mean = cond(trt==0 & time==0, 4, ///
                cond(trt==0 & time==1, 4.5, ///
                cond(trt==0 & time==2, 5.2, ///
                cond(trt==0 & time==3, 6, ///
                cond(trt==1 & time==0, 4, ///
                cond(trt==1 & time==1, 4.4, ///
                cond(trt==1 & time==2, 5, ///
                cond(trt==1 & time==3, 5.5,. ///
                ))))))))

/***************
Set the mean by site and baseline covariate
****************/
gen mean_site =  site/2 - 1.25

gen covar = round(rnormal(6,2),0.1)

/****************
Generate the outcome
*****************/

gen contout = 10 + incept + mean + mean_site + covar + resid
replace contout = round(contout,0.1)
gen baseline0 = contout if time == 0
by pid: egen contbl = max(baseline0)
drop baseline0


/****************
Simulate the categorical outcome by latent variable
****************/

gen catout = -10 + incept + mean + mean_site + covar + rlogistic() > 0
label def  catout 0 "Negative" 1 "Positive"
label val catout catout
/****************
Remove all observations with a positive categorical outcome at baseline to simulate a selection criteria
*****************/
by pid: egen catoutbl = max(cond(time==0,catout,.))
drop if catoutbl==1

/**************
Simulate a time to event outcome
*****************/
gen timemean =  2 + incept + mean + mean_site + covar if time == 0
gen timeout0 = rweibull(1.5,timemean)
gen censor0 = rnormal(12,1) if time == 0
gen timecens0 = timeout < censor0 if time == 0
replace timeout = censor0 if timecens == 0
by pid: egen timeout = max(timeout0)
by pid: egen timecens = max(timecens0)
label def cens 0 "Censored" 1 "Event"
label val timecens cens

/**********
Clean up and name variables
*************/
drop incept resid mean mean_site timemean - timecens0 catoutbl
order pid - covar contbl

label var pid "Patient identifier"
label var site "Site"
label var trt "Treatment"
label var time "Time point"
label var covar "Continuous baseline covariate"
label var contbl "Baseline value continuous outcome"
label var contout "Continuous outcome"
label var catout "Categorical outcome"
label var timeout "Time to event outcome"
label var timecens "Censoring/event identifier"

save rct, replace
