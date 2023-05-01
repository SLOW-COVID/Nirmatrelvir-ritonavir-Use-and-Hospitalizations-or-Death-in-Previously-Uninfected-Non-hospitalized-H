

/*PS match*/

ods graphics on;
%macro match(class,vars,var1);
proc psmatch data=xx.all region=allobs;
   class group1 &class.;
   psmodel group1(Treated='NMV/r')= &vars.;
match method=greedy(k=1) stat=ps  caliper=0.2 exact=();
assess lps  var=(&var1.)/ weight=none;
output out(obs=match)=xx.OutEx4_ps matchid=_MatchID;
run;
%mend;

%match(Race Sex DistrictNameFY17 vacci,
age5 Race Sex bmi DistrictNameFY17 vacci DiabetesAny2yrs HTN2yrs CVD2yrs CKD2yrs ChrLungDisease2yrs Cancer2yrs,
age5 Sex bmi DiabetesAny2yrs HTN2yrs CVD2yrs CKD2yrs ChrLungDisease2yrs Cancer2yrs);
ods graphics off;

/*IPTW match*/

ods graphics on;
%macro IPTW(class,vars,var1);
proc psmatch data= xx.all region=allobs;
   class group1 &class.;
   psmodel group1(Treated='NMV/r')= &vars.;
assess lps var=(&var1.)/ varinfo nlargestwgt=6 plots=(barchart boxplot(display=(lps)) wgtcloud) weight=ATEWGT;
output out(obs=all)=xx.OutEx4_IPTW atewgt=_ATEWgt_;
run;
ods graphics off;
%mend;

%IPTW(Race Sex DistrictNameFY17 vacci,
age5 Race Sex bmi DistrictNameFY17 vacci DiabetesAny2yrs HTN2yrs CVD2yrs CKD2yrs ChrLungDisease2yrs Cancer2yrs,
age5 Sex bmi DiabetesAny2yrs HTN2yrs CVD2yrs CKD2yrs ChrLungDisease2yrs Cancer2yrs);
ods graphics off;

/* Table 1*/


%macro continue(name);
data table2;set table01;keep ScrSSN &name. group1;where &name. ne .;run;
proc means data=table2 MAXDEC=1 n median p25 p75;var &name.;class group1;run;
data table21;set table11;keep ScrSSN &name. group1;where &name. ne .;run;
proc means data=table21 MAXDEC=1 n median p25 p75;var &name.;class group1;run;
%mend continue;

%macro nocontinue(var);
data table2;set table01;keep ScrSSN &var. group1;run;
proc freq  data=table2;table &var.*group1 / norow nopercent;run;
data table12;set table11;keep ScrSSN &var. group1;run;
proc freq  data=table12;table &var.*group1 / norow nopercent;run;
%mend nocontinue;

%stddiff( inds =table01, 
				groupvar =group2, 
				numvars = age bmi CCI2yrs, 
				charvars =Gender Races bmi30 DiabetesAny2yrs HTN2yrs CVD2yrs CKD2yrs ChrLungDisease2yrs Cancer2yrs vacci Rurality,  
				wtvar = ,
   				stdfmt = 8.4,
				outds = stddiff_result ); 

/* Figure 2*/

%macro t2(where);
data table;set table2;where &where.;run;
proc freq data=table;
title1 "&where.";
table group1*incident/nocol norow nopercent riskdiff (CL=(WALD));
run;
proc genmod data=table;
       class group1 (ref="Control");
       model incident(event='1') = group1 / dist=binomial link=identity;
       lsmeans group1 / diff cl;
       run;
%mend t2;


/*Figure 3*/

%style1(cc1 = blue, cc2 = red, cc3 = black, marker1 = circlefilled, 
 marker2 = trianglefilled, marker3 = square); 
quit;

ods listing gpath="" image_dpi=500;
ods graphics on /LINEPATTERNOBSMAX=23200 imagename='Figure2_KM1_PSM_new';
%kmplot(time = t, survival = surv, censored = censored, 
 tatrisk = t, natrisk = atrisk,  atevent = atevent, group = group, linethick = 1, 
 title1 = %str(), 
 title2 = %str(), 
 ylab = %str(Percentage of free from hospitalize and death), 
 legendloc = inside, legendpos = %str(HALIGN = right VALIGN = top), 
 ystart = 1, yend = 0, yinc = -0.05, 
 xlab = %str(Days since diagnosis), 
 xstart = 0, xend = 30, xinc = (3 5 10 15 20 25 30)); 


 proc sgrender data = ple template = kmplot; 
run; 
ods graphics off;
