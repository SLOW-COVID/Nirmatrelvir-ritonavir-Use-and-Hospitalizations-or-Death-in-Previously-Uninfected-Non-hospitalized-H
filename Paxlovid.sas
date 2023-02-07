
/******1:1 match*******/
ods graphics on;
%macro match(vars);
proc psmatch data=xx.newall12 region=allobs;
   class group &vars.;
   psmodel group(Treated='study')= &vars.;
   match method=greedy(k=1) stat=ps  caliper=0.2 exact=(&vars.);
   assess lps/ stddev=pooled(allobs=no) stdbinvar=no  weight=none;
   output out(obs=match)=xx.OutEx4_0 matchid=_MatchID;
run;
%mend;

%match(age5 Race Sex bmis CCI2yrs DistrictNumberFY17 vaccin);
ods graphics off;


/******table 1*******/
proc freq data=table2;
table group*incident/nocol norow nopercent;
run;
%macro t2(where);
data table;set table2;where &where.;run;
proc freq data=table;
title1 "&where.";
table group*incident/nocol norow nopercent;
run;
proc genmod data=table;
 class group (ref='control');
 model incident(event='1') = group / dist=binomial link=identity;
 lsmeans group / diff cl;
run;
%mend t2;
/******table 2*******/
title;
proc freq data=table2;
table group*incident/nocol norow nopercent;
run;
proc genmod data=table2;
 class id group (ref='control');
 model incident(event='1') = group / dist=binomial link=identity;
 repeated subject=id;
 lsmeans group / diff cl;
run;
%macro t2(where);
data table;set table2;where &where.;run;
proc freq data=table;
title1 "&where.";
table group*incident/nocol norow nopercent;
run;
proc genmod data=table;
 class id group (ref='control');
 model incident(event='1') = group / dist=binomial link=identity;
 repeated subject=id;
 lsmeans group / diff cl;
run;
%mend t2;


%t2(age60 eq 0);
%t2(age60 eq 1);
%t2(races eq 0);
%t2(races eq 1);
%t2(races eq 2);
%t2(Gender eq 0);
%t2(Gender eq 1);
%t2(bmi30 eq 1);
%t2(DiabetesAny2yrs eq 1);
%t2(CVD2yrs eq 1);
%t2(CKD2yrs eq 1);
%t2(COPD2yrs eq 1);
%t2(Cancer2yrs eq 1);
%t2(vaccin0 eq 1);
%t2(vaccin0 eq 2);
%t2(vaccin0 eq 3);
%t2(vaccin1 eq 1);
%t2(Symptoms eq 1);
%t2(NoRecordOfSymptoms30d eq 1);

data forest;
  input Indent Subgroup $3-50 study $52-59 control $61-70 Value $72-94 Mean  Low  High;
  datalines;
0 Total........................................... 96/4886. 1372/63350. -0.2(-0.61 to 0.2),,,, -0.2 -0.61 0.2
0 Age............................................. ........ ........... ,,,,,,,,,,,,,,,,,,,,,, . . .     
2 >60 years....................................... 77/3127. 1097/27608. -1.51(-2.1 to -0.92),, -1.51 -2.1 -0.92
2 <=60 years...................................... 19/1759. 275/35742.. 0.31(-0.18 to 0.8),,,, 0.31 -0.18 0.8
0 Race............................................ ........ ........... ,,,,,,,,,,,,,,,,,,,,,, . . .     
2 White........................................... 61/3223. 972/34518.. -0.92(-1.43 to -0.42), -0.92 -1.43 -0.42
2 Black........................................... 19/1210. 289/12605.. -0.72(-1.47 to 0.03),, -0.72 -1.47 0.03
2 Others/unknown.................................. 16/453.. 111/16227.. 2.85(1.14 to 4.55),,,, 2.85 1.14 4.55
0 Sex............................................. ........ ........... ,,,,,,,,,,,,,,,,,,,,,, . . .     
2 Male............................................ 86/4141. 1273/46258. -0.68(-1.13 to -0.22), -0.68 -1.13 -0.22
2 Female.......................................... 10/745.. 99/17092... 0.76(-0.07 to 1.6),,,, 0.76 -0.07 1.6
0 Risk factors.................................... ........ ........... ,,,,,,,,,,,,,,,,,,,,,, . . .     
2 Obese(BMI>30)................................... 43/2224. 503/22570.. -0.3(-0.9 to 0.31),,,, -0.3 -0.9 0.31
2 Diabetes........................................ 48/1820. 706/15980.. -1.78(-2.58 to -0.98), -1.78 -2.58 -0.98
2 Cardiovascular disease.......................... 63/1941. 901/18224.. -1.7(-2.55 to -0.85),, -1.7 -2.55 -0.85
2 Chronic kidney disease.......................... 18/575.. 482/6775... -3.98(-5.53 to -2.43), -3.98 -5.53 -2.43
2 Chronic lung disease(COPD)...................... 32/855.. 476/8508... -1.85(-3.21 to -0.49), -1.85 -3.21 -0.49
2 Cancer diagnosis................................ 38/1027. 450/8550... -1.56(-2.81 to -0.32), -1.56 -2.81 -0.32
0 Vaccination status(no prior infection).......... ........ ........... ,,,,,,,,,,,,,,,,,,,,,, . . .     
2 Unvaccinated or primary series incomplete....... 17/892.. 431/27067.. 0.31(-0.6 to 1.22),,,, 0.31 -0.6 1.22
2 Primary series complete, no booster............. 22/1080. 389/15663.. -0.45(-1.32 to 0.43),, -0.45 -1.32 0.43
2 Primary series complete, booster................ 57/2914. 552/20620.. -0.72(-1.27 to -0.17), -0.72 -1.27 -0.17
2 Primary +/- booster, last dose >3 months ago.... 72/3764. 722/26712.. -0.79(-1.27 to -0.31), -0.79 -1.27 -0.31
0 Symptomatic at baseline......................... 88/3930. 936/35058.. -0.43(-0.92 to 0.06),, -0.43 -0.92 0.06
0 Asymptomatic at baseline........................ 8/956... 436/28292.. -0.7(-1.3 to -0.11),,, -0.7 -1.3 -0.11
;
run;


data forest;
  input Indent Subgroup $3-50 study $52-59 control $61-68 Value $70-92 Mean  Low  High;
  datalines;
0 Total........................................... 73/4148. 109/4148 -0.87(-1.49 to -0.25),, -0.87 -1.49 -0.25
0 Age............................................. ........ ........ ,,,,,,,,,,,,,,,,,,,,,,, . . .     
2 >60 years....................................... 60/2656. 88/2657. -1.05(-1.93 to -0.18),, -1.05 -1.93 -0.18
2 <=60 years...................................... 13/1492. 21/1491. -0.54(-1.3 to 0.23),,,, -0.54 -1.3 0.23
0 Race............................................ ........ ........ ,,,,,,,,,,,,,,,,,,,,,,, . . .     
2 White........................................... 52/2926. 84/2926. -1.09(-1.86 to -0.33),, -1.09 -1.86 -0.33
2 Black........................................... 15/994.. 21/994.. -0.6(-1.75 to 0.55),,,, -0.6 -1.75 0.55
2 Others/unknown.................................. 6/228... 4/228... 0.88(-1.84 to 3.59),,,, 0.88 -1.84 3.59
0 Sex............................................. ........ ........ ,,,,,,,,,,,,,,,,,,,,,,, . . .     
2 Male............................................ 70/3613. 105/3613 -0.97(-1.67 to -0.27),, -0.97 -1.67 -0.27
2 Female.......................................... 3/535... 4/535... -0.19(-1.16 to 0.78),,, -0.19 -1.16 0.78
0 Risk factors.................................... ........ ........ ,,,,,,,,,,,,,,,,,,,,,,, . . .     
2 Obese(BMI>30)................................... 33/1971. 43/1971. -0.51(-1.35 to 0.34),,, -0.51 -1.35 0.34
2 Diabetes........................................ 34/1485. 51/1375. -1.42(-2.65 to -0.19),, -1.42 -2.65 -0.19
2 Cardiovascular disease.......................... 45/1603. 73/1674. -1.55(-2.82 to -0.29),, -1.55 -2.82 -0.29
2 Chronic kidney disease.......................... 13/431.. 34/534.. -3.35(-6 to -0.7),,,,,, -3.35 -6 -0.7
2 Chronic lung disease(COPD)...................... 23/708.. 34/749.. -1.29(-3.26 to 0.67),,, -1.29 -3.26 0.67
2 Cancer diagnosis................................ 29/827.. 36/737.. -1.38(-3.31 to 0.55),,, -1.38 -3.31 0.55
0 Vaccination status(no prior infection).......... ........ ........ ,,,,,,,,,,,,,,,,,,,,,,, . . .     
2 Unvaccinated or primary series incomplete....... 13/745.. 27/745.. -1.88(-3.54 to -0.22),, -1.88 -3.54 -0.22
2 Primary series complete, no booster............. 16/904.. 26/904.. -1.11(-2.51 to 0.3),,,, -1.11 -2.51 0.3
2 Primary series complete, booster................ 44/2499. 56/2499. -0.48(-1.23 to 0.27),,, -0.48 -1.23 0.27
2 Primary +/- booster, last dose >3 months ago.... 55/3208. 64/2464. -0.88(-1.64 to -0.13),, -0.88 -1.64 -0.13
0 Symptomatic at baseline......................... 66/3337. 65/2587. -0.53(-1.3 to 0.23),,,, -0.53 -1.3 0.23
0 Asymptomatic at baseline........................ 7/811... 44/1561. -1.96(-3 to -0.92),,,,, -1.96 -3 -0.92
;
run;

data forest2;
  set forest;
  subgroup=translate(subgroup, ' ', '.');
  study=translate(study, ' ', '.');
  control=translate(control, ' ', '.');
  Value=translate(Value, ' ', ',');
  indent=ifn(indent eq 2, 1, 0);
  if _N_ in(2,3,4,9,10,11,19,20,21,22,23) then ref=subgroup;
  run;

proc template;
  define statgraph Forest;
  dynamic _color _thk _lab;
    begingraph;
	discreteattrmap name='text';
        value '0' / textattrs=(weight=bold SIZE=_lab);
        value other / textattrs=(SIZE=_lab);
      enddiscreteattrmap;
      discreteattrvar attrvar=type var=indent attrmap='text';

      layout lattice / columns=5 columnweights=(0.35 0.10 0.10 0.36 0.19);

      sidebar / align=top;
        layout lattice / rows=2 columns=3 columnweights=(0.35 0.2 0.55);
          entry textattrs=(size=8) halign=left "Subgroup";
          entry textattrs=(size=8) halign=center "NMV/r         Control";
          entry halign=center textattrs=(size=8) "Absolute risk difference(95% CI)" ;
          entry " "; 
          entry halign=center textattrs=(size=6) "no. of events / no. of participants";
          entry halign=center textattrs=(size=6) "percentage points";
        endlayout;
      endsidebar;

      layout overlay / walldisplay=none xaxisopts=(display=none) 
          yaxisopts=(reverse=true display=none 
                     tickvalueattrs=(weight=bold));
        referenceline y=ref / lineattrs=(thickness=_thk color=_color);
        axistable y=subgroup value=subgroup / indentweight=indent textgroup=type valueattrs=(SIZE=_lab) display=(values);
       endlayout;

       layout overlay / xaxisopts=(display=none) 
            yaxisopts=(reverse=true display=none) walldisplay=none;
         referenceline y=ref / lineattrs=(thickness=_thk color=_color);
         axistable y=subgroup value=study/valueattrs=(SIZE=_lab) display=(values);
       endlayout;

       layout overlay / xaxisopts=(display=none) 
            yaxisopts=(reverse=true display=none) walldisplay=none;
         referenceline y=ref / lineattrs=(thickness=_thk color=_color);
         axistable y=subgroup value=control/valueattrs=(SIZE=_lab) display=(values);
       endlayout;

       layout overlay / xaxisopts=(label='<---NMV/r Better----  ----Control Better--->'  labelattrs=(SIZE=_lab)
           linearopts=(tickvaluepriority=true tickvaluelist=(-10 -5 0 5 10 )) tickvalueattrs=(SIZE=_lab))
           /*linearopts=(tickvaluepriority=true tickvaluelist=(-15 -10 -5 0 5 10 15 20)) tickvalueattrs=(SIZE=_lab))*/
           yaxisopts=(reverse=true display=none) walldisplay=none;
         referenceline y=ref / lineattrs=(thickness=_thk color=_color);
         scatterplot y=subgroup x=mean / xerrorlower=low xerrorupper=high markerattrs=(symbol=squarefilled);
         referenceline x=0;
       endlayout;

       layout overlay / xaxisopts=(display=none) 
            yaxisopts=(reverse=true display=none) walldisplay=none;
         referenceline y=ref / lineattrs=(thickness=_thk color=_color);
         axistable y=subgroup value=Value /display=(values) valueattrs=(SIZE=_lab);
       endlayout;

     endlayout;
   endgraph;
  end;
run;
