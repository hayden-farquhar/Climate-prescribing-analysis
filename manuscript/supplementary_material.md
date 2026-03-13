---
title: |
  **Supplementary Material**

  Ambient Temperature and Population-Level Medication Dispensing in Australia: A Distributed Lag Non-Linear Analysis with Equity Stratification
author:
  - "Hayden Farquhar MBBS MPHTM"
date: "March 2026"
geometry: margin=2.5cm
fontsize: 10pt
linestretch: 1.3
numbersections: false
header-includes:
  - \usepackage{graphicx}
  - \usepackage{booktabs}
  - \usepackage{longtable}
  - \usepackage{float}
  - \usepackage{multirow}
  - \usepackage{caption}
  - \captionsetup[figure]{labelfont=bf}
  - \captionsetup[table]{labelfont=bf}
  - \usepackage{pdflscape}
  - \usepackage{hyperref}
  - \hypersetup{colorlinks=true, linkcolor=blue, citecolor=blue, urlcolor=blue}
  - \renewcommand{\figurename}{Supplementary Figure}
  - \renewcommand{\tablename}{Supplementary Table}
---

\newpage

# Supplementary Tables

## Supplementary Table S1. Medication group definitions

\begin{table}[H]
\centering
\small
\begin{tabular}{p{3.5cm}p{3cm}p{6.5cm}}
\toprule
\textbf{Analysis group} & \textbf{AIHW category label} & \textbf{Likely constituent medications} \\
\midrule
\multicolumn{3}{l}{\textit{Cardiovascular (ATC class C)}} \\
\quad Antithrombotic & Antithrombotic agents & Aspirin (low-dose), clopidogrel, warfarin, rivaroxaban, apixaban \\
\quad BP lowering & Blood pressure lowering & ACE inhibitors, ARBs, calcium channel blockers, beta-blockers, diuretics \\
\quad Lipid modifying & Lipid modifying agents & Statins (atorvastatin, rosuvastatin), ezetimibe, fibrates \\
\quad Other CVD & Other cardiovascular & Digoxin, amiodarone, nitrates, other cardiac agents \\
\midrule
\multicolumn{3}{l}{\textit{Mental Health (ATC classes N05/N06)}} \\
\quad Antidepressants & Antidepressants & SSRIs, SNRIs, TCAs (e.g., sertraline, venlafaxine, amitriptyline) \\
\quad Anxiolytics & Anxiolytics & Benzodiazepines (diazepam, oxazepam), buspirone \\
\quad Other MH & Other mental health & Antipsychotics, mood stabilisers (e.g., quetiapine, lithium) \\
\midrule
\multicolumn{3}{l}{\textit{Respiratory (ATC class R03 + H02\textsuperscript{a})}} \\
\quad Relievers & Respiratory -- relievers & Short-acting beta-2 agonists: salbutamol, terbutaline \\
\quad Preventers & Respiratory -- preventers & Inhaled corticosteroids (ICS): fluticasone, budesonide; ICS/LABA combinations: fluticasone/salmeterol, budesonide/formoterol \\
\quad COPD/Other & Other respiratory (COPD-specific treatments, other asthma) & Long-acting muscarinic antagonists: tiotropium, glycopyrronium; LAMA/LABA combinations: umeclidinium/vilanterol; leukotriene receptor antagonists: montelukast \\
\quad Oral corticosteroids\textsuperscript{a} & Oral corticosteroids & Prednis(ol)one (oral formulations) \\
\bottomrule
\multicolumn{3}{p{13cm}}{\footnotesize Medication group definitions are as provided by the Australian Institute of Health and Welfare (AIHW) in the PHE-338 ``Geography and time-specific health data for environmental analysis'' dataset. The AIHW does not publish the exact ATC codes underlying each therapeutic group; constituent medications listed here are the most commonly dispensed PBS items within each group based on published PBS statistics.} \\
\multicolumn{3}{l}{\footnotesize \textsuperscript{a}Oral corticosteroids are ATC class H02AB, not R03, but are included in the respiratory analysis given their clinical use in managing acute respiratory exacerbations.} \\
\end{tabular}
\end{table}

## Supplementary Table S2. Detailed climate descriptive statistics

\begin{table}[H]
\centering
\small
\begin{tabular}{lr}
\toprule
\textbf{Variable} & \textbf{Mean (SD)} \\
\midrule
Weekly mean Tmax ($^\circ$C) & 22.78 (6.46) \\
Weekly max Tmax ($^\circ$C) & 26.31 (7.00) \\
Weekly mean Tmin ($^\circ$C) & 12.67 (5.64) \\
Weekly total precipitation (mm) & 14.60 (23.13) \\
Days above 90th percentile per week & 0.72 (1.41) \\
Days above 95th percentile per week & 0.37 (0.96) \\
Days above 99th percentile per week & 0.08 (0.39) \\
Heatwave days per week & 0.17 (0.79) \\
\bottomrule
\end{tabular}
\end{table}

## Supplementary Table S3. Equity-stratified cumulative rate ratios for all medication groups

\begin{table}[H]
\centering
\small
\begin{tabular}{llcccc}
\toprule
\textbf{Group} & \textbf{Quintile} & \textbf{n SA4s} & \textbf{RR at p95 (95\% CI)} & \textbf{RR at p99 (95\% CI)} & \textbf{Dispersion} \\
\midrule
\multicolumn{6}{l}{\textit{Cardiovascular (interaction p = 0.66)}} \\
& Q1 (most disadvantaged) & 17 & 1.018 (1.001--1.037) & 1.044 (1.018--1.071) & 68.65 \\
& Q2 & 14 & 1.012 (0.997--1.027) & 1.025 (1.002--1.047) & 65.99 \\
& Q3 & 14 & 1.004 (0.988--1.020) & 1.002 (0.981--1.024) & 91.93 \\
& Q4 & 16 & 1.000 (0.986--1.015) & 0.989 (0.967--1.011) & 90.83 \\
& Q5 (least disadvantaged) & 8 & 0.997 (0.977--1.017) & 1.000 (0.973--1.028) & 92.29 \\
\midrule
\multicolumn{6}{l}{\textit{Mental Health (interaction p = 0.48)}} \\
& Q1 & 17 & 1.012 (0.996--1.028) & 1.030 (1.006--1.054) & 18.96 \\
& Q2 & 14 & 0.991 (0.977--1.005) & 1.006 (0.985--1.026) & 20.85 \\
& Q3 & 14 & 1.002 (0.988--1.016) & 0.998 (0.979--1.017) & 25.41 \\
& Q4 & 16 & 1.006 (0.993--1.019) & 0.995 (0.975--1.015) & 27.85 \\
& Q5 & 8 & 0.996 (0.975--1.016) & 1.006 (0.979--1.034) & 34.27 \\
\midrule
\multicolumn{6}{l}{\textit{Respiratory (interaction p = 6.4 $\times$ 10$^{-11}$)}} \\
& Q1 (most disadvantaged) & 17 & 0.920 (0.892--0.949) & 0.940 (0.899--0.984) & 27.37 \\
& Q2 & 14 & 0.898 (0.873--0.923) & 0.892 (0.856--0.930) & 29.70 \\
& Q3 & 14 & 0.892 (0.865--0.921) & 0.873 (0.836--0.911) & 42.13 \\
& Q4 & 16 & 0.900 (0.872--0.928) & 0.844 (0.804--0.885) & 44.27 \\
& Q5 (least disadvantaged) & 8 & 0.884 (0.847--0.922) & 0.837 (0.789--0.887) & 42.71 \\
\bottomrule
\end{tabular}
\end{table}

## Supplementary Table S4. Remoteness-stratified cumulative rate ratios

\begin{table}[H]
\centering
\small
\begin{tabular}{llccc}
\toprule
\textbf{Group} & \textbf{Remoteness} & \textbf{n SA4s} & \textbf{RR at p95 (95\% CI)} & \textbf{RR at p99 (95\% CI)} \\
\midrule
Cardiovascular & Major Cities & 11 & 1.007 (0.990--1.024) & 0.993 (0.972--1.015) \\
& Inner Regional & 28 & 1.015 (1.004--1.026) & 1.006 (0.992--1.021) \\
& Outer Regional/Remote & 31 & 1.009 (0.999--1.020) & 1.028 (1.011--1.046) \\
\midrule
Mental Health & Major Cities & 11 & 1.012 (0.997--1.026) & 1.010 (0.991--1.029) \\
& Inner Regional & 28 & 1.017 (1.007--1.027) & 1.011 (0.998--1.025) \\
& Outer Regional/Remote & 31 & 1.003 (0.993--1.013) & 1.022 (1.006--1.039) \\
\midrule
Respiratory & Major Cities & 11 & 0.978 (0.942--1.014) & 0.895 (0.853--0.939) \\
& Inner Regional & 28 & 0.907 (0.888--0.926) & 0.848 (0.824--0.873) \\
& Outer Regional/Remote & 31 & 0.907 (0.891--0.924) & 0.910 (0.883--0.938) \\
\bottomrule
\multicolumn{5}{l}{\footnotesize Interaction F-test: CVD p = 0.0002; MH p = 0.045; Respiratory p = 1.4 $\times$ 10$^{-14}$.} \\
\end{tabular}
\end{table}

## Supplementary Table S5. Attributable burden of heat and cold on respiratory dispensing

\begin{table}[H]
\centering
\small
\begin{tabular}{lrrr}
\toprule
\textbf{Outcome / Stratum} & \textbf{Heat-attributable (AF \%)} & \textbf{Extreme heat (AF \%)} & \textbf{Cold-attributable (AF \%)} \\
\midrule
\textit{Overall} & & & \\
\quad Respiratory (total) & $-$1,669,137 ($-$1.73\%) & $-$479,809 ($-$0.50\%) & +692,927 (+0.72\%) \\
\quad Relievers & $-$765,322 ($-$2.70\%) & $-$210,844 ($-$0.74\%) & +318,737 (+1.13\%) \\
\quad Preventers & $-$855,022 ($-$1.90\%) & $-$248,841 ($-$0.55\%) & +324,816 (+0.72\%) \\
\midrule
\textit{By IRSD quintile (respiratory total)} & & & \\
\quad Q1 (most disadvantaged) & $-$367,543 ($-$1.84\%) & $-$46,938 ($-$0.24\%) & +162,039 (+0.81\%) \\
\quad Q2 & $-$444,888 ($-$2.63\%) & $-$135,207 ($-$0.80\%) & +395,366 (+2.34\%) \\
\quad Q3 & $-$483,669 ($-$2.24\%) & $-$171,983 ($-$0.80\%) & +161,464 (+0.75\%) \\
\quad Q4 & $-$253,132 ($-$0.97\%) & $-$129,791 ($-$0.50\%) & +115,993 (+0.44\%) \\
\quad Q5 (least disadvantaged) & $-$280,467 ($-$2.51\%) & $-$129,988 ($-$1.16\%) & +178,875 (+1.60\%) \\
\bottomrule
\multicolumn{4}{l}{\footnotesize Heat-attributable: excess/deficit dispensings attributable to temperatures above the median.} \\
\multicolumn{4}{l}{\footnotesize Extreme heat: attributable to temperatures above the 95th percentile.} \\
\multicolumn{4}{l}{\footnotesize Cold-attributable: excess dispensings attributable to temperatures below the median.} \\
\end{tabular}
\end{table}

## Supplementary Table S6. Cold and heat associations across the full temperature distribution

\begin{table}[H]
\centering
\small
\begin{tabular}{lccccc}
\toprule
\textbf{Group} & \textbf{RR at p01} & \textbf{RR at p05} & \textbf{RR at p10} & \textbf{RR at p90} & \textbf{RR at p95} \\
\midrule
Cardiovascular & 0.986 (0.978--0.994) & 0.990 (0.984--0.996) & 0.993 (0.988--0.998) & 1.001 (0.995--1.007) & 1.004 (0.997--1.011) \\
Mental Health & 1.000 (0.992--1.007) & 0.999 (0.994--1.004) & 0.999 (0.994--1.004) & 1.001 (0.995--1.007) & 1.003 (0.997--1.010) \\
Respiratory & 0.991 (0.975--1.006) & 1.007 (0.995--1.018) & 1.016 (1.005--1.026) & 0.925 (0.914--0.936) & 0.910 (0.898--0.923) \\
Relievers & 1.002 (0.981--1.024) & 1.016 (1.001--1.032) & 1.025 (1.011--1.038) & 0.885 (0.872--0.900) & 0.868 (0.853--0.884) \\
Preventers & 0.986 (0.970--1.002) & 1.004 (0.992--1.016) & 1.015 (1.004--1.026) & 0.915 (0.904--0.927) & 0.898 (0.885--0.911) \\
\bottomrule
\multicolumn{6}{l}{\footnotesize All RRs are cumulative (0--12 week lag), centred at the median temperature (22.8$^\circ$C).} \\
\end{tabular}
\end{table}

## Supplementary Table S7. Black Summer bushfire difference-in-differences (full results)

\begin{table}[H]
\centering
\small
\begin{tabular}{lcccc}
\toprule
\textbf{Medication} & \textbf{Medium smoke RR (95\% CI)} & \textbf{p} & \textbf{High smoke RR (95\% CI)} & \textbf{p} \\
\midrule
Total respiratory & 1.014 (0.997--1.032) & 0.10 & 1.062 (1.046--1.078) & 2.4 $\times$ 10$^{-15}$ \\
Relievers & 1.026 (1.002--1.050) & 0.03 & 1.093 (1.073--1.114) & 4.2 $\times$ 10$^{-20}$ \\
Preventers & 1.008 (0.990--1.025) & 0.40 & 1.071 (1.055--1.087) & 3.0 $\times$ 10$^{-18}$ \\
COPD/Other & 1.017 (1.002--1.032) & 0.03 & 1.013 (1.000--1.026) & 0.047 \\
\bottomrule
\multicolumn{5}{l}{\footnotesize Models include SA4 fixed effects, temperature, precipitation, and seasonal controls.} \\
\multicolumn{5}{l}{\footnotesize N = 5,772 SA4-week observations across 78 SA4 regions.} \\
\end{tabular}
\end{table}

## Supplementary Table S8. Newey-West HAC standard error comparison

\begin{table}[H]
\centering
\small
\begin{tabular}{lccccc}
\toprule
\textbf{Group} & \textbf{Standard SE} & \textbf{HAC SE} & \textbf{SE ratio} & \textbf{Sig (standard)} & \textbf{Sig (HAC)} \\
\midrule
Cardiovascular & 0.00355 & 0.00556 & 1.57 & No & No \\
Mental Health & 0.00327 & 0.00510 & 1.56 & No & No \\
Respiratory & 0.00688 & 0.01177 & 1.71 & Yes & Yes \\
Resp. Relievers & 0.00905 & 0.01577 & 1.74 & Yes & Yes \\
Resp. Preventers & 0.00732 & 0.01297 & 1.77 & Yes & Yes \\
\bottomrule
\multicolumn{6}{l}{\footnotesize Newey-West bandwidth = 12 weeks. SE ratio > 1 indicates model-based SEs underestimate uncertainty.} \\
\end{tabular}
\end{table}

## Supplementary Table S9. Lag constraint sensitivity analysis

\begin{table}[H]
\centering
\small
\begin{tabular}{llccc}
\toprule
\textbf{Group} & \textbf{Lag specification} & \textbf{RR at p95 (95\% CI)} & \textbf{Dispersion} & \textbf{QAIC} \\
\midrule
\multirow{5}{*}{Cardiovascular} & Primary (ns, 0--12w) & 1.004 (0.997--1.011) & 81.09 & 34,234 \\
& Unconstrained indicator & 1.005 (0.998--1.012) & 80.42 & 34,277 \\
& Polynomial (degree 3) & 1.004 (0.997--1.011) & 81.07 & 34,230 \\
& Shorter (0--8w) & 0.998 (0.992--1.004) & 81.18 & 34,232 \\
& Longer (0--16w) & 0.986 (0.978--0.994) & 80.75 & 34,229 \\
\midrule
\multirow{5}{*}{Respiratory} & Primary (ns, 0--12w) & 0.910 (0.898--0.923) & 37.05 & 31,556 \\
& Unconstrained indicator & 0.911 (0.899--0.924) & 36.65 & 31,692 \\
& Polynomial (degree 3) & 0.911 (0.899--0.923) & 37.06 & 31,550 \\
& Shorter (0--8w) & 0.894 (0.883--0.904) & 37.07 & 31,574 \\
& Longer (0--16w) & 0.902 (0.889--0.916) & 36.90 & 31,560 \\
\bottomrule
\end{tabular}
\end{table}

## Supplementary Table S10. Fourier harmonic sensitivity analysis

\begin{table}[H]
\centering
\small
\begin{tabular}{llccc}
\toprule
\textbf{Group} & \textbf{Fourier pairs} & \textbf{RR at p95 (95\% CI)} & \textbf{Dispersion} & \textbf{QAIC} \\
\midrule
Cardiovascular & 2 pairs & 1.004 (0.997--1.011) & 81.09 & 34,234 \\
& 3 pairs & 1.005 (0.999--1.012) & 75.12 & 31,796 \\
Mental Health & 2 pairs & 1.003 (0.997--1.010) & 25.09 & 34,569 \\
& 3 pairs & 1.005 (0.999--1.011) & 23.71 & 32,695 \\
Respiratory & 2 pairs & 0.910 (0.898--0.923) & 37.05 & 31,556 \\
& 3 pairs & 0.912 (0.900--0.923) & 32.74 & 28,027 \\
Relievers & 2 pairs & 0.868 (0.853--0.884) & 18.33 & 31,096 \\
& 3 pairs & 0.870 (0.856--0.885) & 16.24 & 27,730 \\
Preventers & 2 pairs & 0.898 (0.885--0.911) & 19.55 & 30,840 \\
& 3 pairs & 0.899 (0.887--0.911) & 17.52 & 27,798 \\
\bottomrule
\multicolumn{5}{l}{\footnotesize 3 Fourier pairs preferred by QAIC in all groups, but RR estimates change $<$ 0.5\%.} \\
\end{tabular}
\end{table}

## Supplementary Table S11. Influence diagnostics: sensitivity to high-leverage observations

\begin{table}[H]
\centering
\small
\begin{tabular}{lcccc}
\toprule
\textbf{Group} & \textbf{RR (full)} & \textbf{RR (trimmed)} & \textbf{N removed} & \textbf{\% removed} \\
\midrule
Cardiovascular & 1.004 (0.997--1.011) & 0.991 (0.984--0.998) & 1,852 & 5.4\% \\
Mental Health & 1.003 (0.997--1.010) & 0.998 (0.991--1.004) & 2,133 & 6.2\% \\
Respiratory & 0.910 (0.898--0.923) & 0.922 (0.910--0.935) & 1,400 & 4.1\% \\
\bottomrule
\multicolumn{5}{l}{\footnotesize Trimmed models exclude observations with Cook's distance in the top 5\%.} \\
\end{tabular}
\end{table}

## Supplementary Table S12. Spatial autocorrelation in model residuals

\begin{table}[H]
\centering
\small
\begin{tabular}{lcccc}
\toprule
\textbf{Group} & \textbf{Moran's I (overall)} & \textbf{p-value} & \textbf{Weekly mean I} & \textbf{\% weeks significant} \\
\midrule
Cardiovascular & 0.275 & 0.0006 & 0.304 & 88\% \\
Mental Health & 0.288 & 0.0003 & 0.293 & 88\% \\
Respiratory & 0.354 & $<$ 0.0001 & 0.339 & 90\% \\
\bottomrule
\multicolumn{5}{l}{\footnotesize Moran's I computed on deviance residuals using SA4 contiguity weights (queen criterion).} \\
\multicolumn{5}{l}{\footnotesize Weekly significance assessed at $\alpha$ = 0.05 across 50 randomly sampled weeks.} \\
\end{tabular}
\end{table}

## Supplementary Table S13. E-values for primary findings and bushfire analysis

\begin{table}[H]
\centering
\small
\begin{tabular}{llcccc}
\toprule
\textbf{Source} & \textbf{Outcome} & \textbf{RR} & \textbf{E-value (point)} & \textbf{E-value (CI)} \\
\midrule
Primary DLNM & Cardiovascular & 1.004 & 1.07 & 1.00 \\
& Mental Health & 1.003 & 1.06 & 1.00 \\
& Respiratory & 0.910 & 1.43 & 1.39 \\
& Relievers & 0.868 & 1.57 & 1.52 \\
& Preventers & 0.898 & 1.47 & 1.43 \\
\midrule
Bushfire DiD & All Respiratory & 1.062 & 1.32 & 1.27 \\
& Relievers & 1.093 & 1.41 & 1.35 \\
& Preventers & 1.071 & 1.35 & 1.29 \\
& COPD/Other & 1.013 & 1.13 & 1.01 \\
\bottomrule
\multicolumn{5}{l}{\footnotesize E-value: the minimum RR an unmeasured confounder would need with both exposure} \\
\multicolumn{5}{l}{\footnotesize and outcome to explain away the observed association.} \\
\end{tabular}
\end{table}

## Supplementary Table S14. Temporal adaptation analysis (early vs late study period)

\begin{table}[H]
\centering
\small
\begin{tabular}{llcc}
\toprule
\textbf{Group} & \textbf{Period} & \textbf{RR at p95 (95\% CI)} & \textbf{Dispersion} \\
\midrule
Cardiovascular & 2013--2017 & 0.999 (0.991--1.007) & 53.23 \\
& 2018--2023 & 0.982 (0.971--0.992) & 90.35 \\
Mental Health & 2013--2017 & 1.004 (0.997--1.012) & 17.86 \\
& 2018--2023 & 0.985 (0.976--0.994) & 26.37 \\
Respiratory & 2013--2017 & 0.909 (0.897--0.922) & 21.10 \\
& 2018--2023 & 0.885 (0.865--0.905) & 47.55 \\
\bottomrule
\multicolumn{4}{l}{\footnotesize Interaction F-test: CVD p = 9.5 $\times$ 10$^{-13}$; MH p = 0.005; Respiratory p = 0.0001.} \\
\end{tabular}
\end{table}

## Supplementary Table S15. Model diagnostics summary

\begin{table}[H]
\centering
\small
\begin{tabular}{lcccccc}
\toprule
\textbf{Group} & \textbf{N obs} & \textbf{Dispersion} & \textbf{Mean resid} & \textbf{SD resid} & \textbf{Skewness} & \textbf{Mean ACF lag-1} \\
\midrule
Cardiovascular & 34,568 & 81.09 & $-$0.08 & 8.93 & 0.82 & 0.34 \\
Mental Health & 34,568 & 25.09 & $-$0.04 & 4.99 & 0.25 & 0.29 \\
Respiratory & 34,568 & 37.05 & $-$0.09 & 5.79 & 3.16 & 0.55 \\
\bottomrule
\multicolumn{7}{l}{\footnotesize Deviance residuals from quasi-Poisson DLNM models. ACF = autocorrelation function.} \\
\end{tabular}
\end{table}

## Supplementary Table S16. Leave-one-out jurisdiction analysis (respiratory dispensing)

\begin{table}[H]
\centering
\small
\begin{tabular}{lccc}
\toprule
\textbf{Jurisdiction excluded} & \textbf{n SA4s remaining} & \textbf{RR at p95 (95\% CI)} & \textbf{RR at p99 (95\% CI)} \\
\midrule
ACT & 69 & 0.912 (0.900--0.925) & 0.895 (0.873--0.917) \\
NSW & 51 & 0.896 (0.881--0.912) & 0.892 (0.865--0.920) \\
NT & 68 & 0.917 (0.905--0.929) & 0.894 (0.874--0.915) \\
Queensland & 54 & 0.921 (0.906--0.936) & 0.886 (0.860--0.912) \\
SA & 64 & 0.909 (0.896--0.922) & 0.892 (0.869--0.915) \\
Tasmania & 66 & 0.911 (0.898--0.923) & 0.892 (0.870--0.915) \\
Victoria & 57 & 0.917 (0.904--0.930) & 0.907 (0.884--0.931) \\
WA & 61 & 0.903 (0.890--0.916) & 0.886 (0.867--0.907) \\
\bottomrule
\multicolumn{4}{l}{\footnotesize All estimates remain statistically significant regardless of which jurisdiction is excluded.} \\
\end{tabular}
\end{table}

\newpage

# Supplementary Figures

\begin{figure}[H]
\centering
\includegraphics[width=\textwidth]{figures/figS1_lag_response.png}
\caption{Lag-response curves for temperature--dispensing associations. Panels show the lag-specific rate ratio at the 95th percentile of temperature for each medication group over 0--12 weeks. Shaded areas represent 95\% confidence intervals.}
\label{fig:s1_lag}
\end{figure}

\begin{figure}[H]
\centering
\includegraphics[width=\textwidth]{figures/figS2_remoteness.png}
\caption{Remoteness-stratified exposure-response curves. Cumulative rate ratios for temperature and medication dispensing stratified by remoteness classification (Major Cities, Inner Regional, Outer Regional/Remote).}
\label{fig:s2_remoteness}
\end{figure}

\begin{figure}[H]
\centering
\includegraphics[width=\textwidth]{figures/figS3_negative_control.png}
\caption{Negative control analysis. Exposure-response curves for antidepressants and anxiolytics (negative controls) compared with respiratory medications (positive control). The null associations for negative controls support the specificity of the temperature--respiratory dispensing relationship.}
\label{fig:s3_negcontrol}
\end{figure}

\begin{figure}[H]
\centering
\includegraphics[width=\textwidth]{figures/figS4_temporal.png}
\caption{Temporal adaptation analysis. Comparison of exposure-response curves between the early (2013--2017) and late (2018--2023) study periods. The respiratory temperature--dispensing association appears to strengthen over time.}
\label{fig:s4_temporal}
\end{figure}

\begin{figure}[H]
\centering
\includegraphics[width=\textwidth]{figures/newey_west_forest.png}
\caption{Newey-West HAC standard error comparison. Forest plot comparing cumulative rate ratios at the 95th percentile with model-based standard errors versus Newey-West heteroskedasticity and autocorrelation-consistent (HAC) standard errors (bandwidth = 12 weeks). HAC confidence intervals are wider but all respiratory findings remain significant.}
\label{fig:s5_nw}
\end{figure}

\begin{figure}[H]
\centering
\includegraphics[width=\textwidth]{figures/fourier_sensitivity_forest.png}
\caption{Fourier harmonic sensitivity analysis. Forest plot comparing cumulative rate ratios at the 95th percentile using 2 versus 3 Fourier harmonic pairs for seasonal adjustment. Point estimates are virtually identical across specifications.}
\label{fig:s6_fourier}
\end{figure}

\begin{figure}[H]
\centering
\includegraphics[width=\textwidth]{figures/loo_jurisdiction_panel.png}
\caption{Leave-one-out jurisdiction analysis. Cumulative rate ratios at the 95th percentile for each medication group when sequentially excluding each Australian state/territory. No single jurisdiction drives the primary findings.}
\label{fig:s7_loo}
\end{figure}

\begin{figure}[H]
\centering
\includegraphics[width=\textwidth]{figures/lag_constraint_comparison.png}
\caption{Lag constraint sensitivity analysis. Cumulative rate ratios at the 95th percentile across five lag specifications: primary natural cubic spline (0--12 weeks), unconstrained indicator lags, polynomial degree 3, shorter window (0--8 weeks), and longer window (0--16 weeks). All specifications produce similar estimates for respiratory medications.}
\label{fig:s8_lag}
\end{figure}

\begin{figure}[H]
\centering
\includegraphics[width=\textwidth]{figures/crossbasis_sensitivity_panel.png}
\caption{Cross-basis sensitivity analysis. Cumulative rate ratios across 9 alternative cross-basis specifications varying the number and placement of knots for both the temperature and lag dimensions. Respiratory estimates are robust across all specifications.}
\label{fig:s9_cb}
\end{figure}

\begin{figure}[H]
\centering
\includegraphics[width=\textwidth]{figures/qaic_comparison.png}
\caption{QAIC model comparison. Delta-QAIC values (relative to the best model) for 9 model specifications across medication groups. The full model (with all covariates) consistently achieves the lowest QAIC.}
\label{fig:s10_qaic}
\end{figure}

\begin{figure}[H]
\centering
\includegraphics[width=\textwidth]{figures/mmt_curves_panel.png}
\caption{Minimum morbidity temperature (MMT) analysis. Full exposure-response curves with MMT identified for each medication group. For respiratory medications, the MMT is at the upper extreme of the temperature distribution, confirming a monotonically decreasing relationship.}
\label{fig:s11_mmt}
\end{figure}

\begin{figure}[H]
\centering
\includegraphics[width=\textwidth]{figures/diagnostics_residuals_resp_total.png}
\caption{Model diagnostic plots for respiratory medications. Deviance residuals versus fitted values (left) and temporal residual pattern (right) from the primary quasi-Poisson DLNM model.}
\label{fig:s12_diag}
\end{figure}

\begin{figure}[H]
\centering
\includegraphics[width=\textwidth]{figures/diagnostics_acf_resp_total.png}
\caption{Autocorrelation function (ACF) of deviance residuals for respiratory medications. The significant lag-1 autocorrelation (mean 0.55) motivates the Newey-West HAC standard error analysis.}
\label{fig:s13_acf}
\end{figure}

\begin{figure}[H]
\centering
\includegraphics[width=\textwidth]{figures/descriptive_temperature_distribution.png}
\caption{Distribution of weekly mean maximum temperature across SA4-week observations. Vertical dashed lines indicate the 5th, 50th, and 95th percentiles of the distribution.}
\label{fig:s14_tempdist}
\end{figure}

\begin{figure}[H]
\centering
\includegraphics[width=\textwidth]{figures/descriptive_time_series.png}
\caption{Time series of weekly dispensing counts by medication group across the study period (2013--2022). The seasonal pattern in respiratory dispensing (winter peaks, summer troughs) is clearly visible.}
\label{fig:s15_timeseries}
\end{figure}
