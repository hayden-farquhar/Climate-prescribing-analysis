---
title: "Ambient Temperature and Population-Level Medication Dispensing in Australia: A Distributed Lag Non-Linear Analysis with Equity Stratification"
shorttitle: "Heat, Prescribing, and Health Equity in Australia"
author:
  - "**Hayden Farquhar** MBBS MPHTM"
institute:
  - "Independent Researcher, NSW, Australia"
date: "March 2026"
geometry: margin=2.5cm
fontsize: 11pt
linestretch: 1.5
numbersections: true
header-includes:
  - \usepackage{graphicx}
  - \usepackage{booktabs}
  - \usepackage{longtable}
  - \usepackage{float}
  - \usepackage{caption}
  - \captionsetup[figure]{labelfont=bf}
  - \captionsetup[table]{labelfont=bf}
  - \usepackage{pdflscape}
  - \usepackage{hyperref}
  - \hypersetup{colorlinks=true, linkcolor=blue, citecolor=blue, urlcolor=blue}
---

\vspace{1em}
\noindent\textbf{Correspondence:} hayden.farquhar@icloud.com

\noindent\textbf{ORCID:} 0009-0002-6226-440X

\newpage

# Abstract

**Background:** Climate change is intensifying heat exposure, yet the impact of ambient temperature on population-level medication dispensing remains unexplored. We examined associations between weekly maximum temperature and dispensing of cardiovascular, psychotropic, and respiratory medications across Australia, with stratification by socioeconomic disadvantage.

**Methods:** Nationwide ecological time-series study linking weekly Pharmaceutical Benefits Scheme dispensing data with gridded maximum temperature (Bureau of Meteorology AGCD) across 88 SA4 regions (January 2013 -- June 2022; 34,580 SA4-week observations per medication group after exclusion of suppressed cells). Distributed lag non-linear models with quasi-Poisson regression estimated cumulative temperature--dispensing associations over 0--12 week lags. Equity analyses stratified by IRSD quintile. A secondary difference-in-differences analysis examined the 2019--20 Black Summer bushfire season.

**Results:** Heat was not associated with cardiovascular (cumulative RR at 95th vs median temperature: 1.004, 95% CI: 0.997--1.011) or mental health dispensing (RR 1.003, 95% CI: 0.997--1.010). Heat was inversely associated with respiratory dispensing (RR 0.910, 95% CI: 0.898--0.923), driven by relievers (RR 0.868) and preventers (RR 0.898), with significant effect modification by disadvantage (interaction p = 6.4 $\times$ 10$^{-11}$). Negative controls showed null associations. During the Black Summer, high-smoke regions showed a 6.2% increase in respiratory dispensing (RR 1.062, 95% CI: 1.046--1.078). Findings were robust across 11 sensitivity analyses.

**Conclusions:** This study links ambient temperature to medication dispensing at the population level, revealing a strong inverse respiratory association likely reflecting seasonal viral dynamics, with a socioeconomic gradient suggesting differential disease management capacity. Bushfire smoke acutely increases respiratory pharmaceutical demand, informing health system preparedness under climate change.

**Keywords:** climate change; temperature; pharmaceutical dispensing; respiratory medications; health equity; distributed lag non-linear model; bushfire smoke; Australia

\newpage

# Introduction

Climate change represents the greatest threat to global health in the 21st century, with rising temperatures and increasing frequency of extreme heat events projected to compound existing disease burdens [1]. The 2021 Lancet Countdown reported record heatwave exposure among vulnerable populations globally, with 3.1 billion additional person-days of heatwave exposure among older adults compared to the 1986--2005 baseline [2]. Australia is particularly vulnerable, having experienced substantial warming since the mid-20th century and a marked increase in the frequency, duration, and intensity of extreme heat events [3].

The health impacts of ambient temperature are well characterised for mortality outcomes. In the seminal multi-country analysis by Gasparrini et al. [4], non-optimal temperatures were attributed to 7.71% of deaths across 384 locations in 13 countries, with the temperature--mortality relationship exhibiting a non-linear, J- or U-shaped curve. Cold exposure generally accounts for a larger attributable mortality burden than heat, though heat-related mortality is more acute and concentrated in extreme events [4,5].

For morbidity outcomes, the evidence base is substantial but less consistent. A systematic review of heat waves and morbidity identified increased risks across cardiovascular, respiratory, and renal outcomes, though with heterogeneous effect sizes [24]. A meta-analysis of 54 studies across 20 countries found that heatwaves significantly increased cardiorespiratory mortality but had more marginal effects on morbidity indicators such as emergency department presentations and hospitalisations [6]. In Australia, studies have linked heat exposure to increased hospital admissions for cardiovascular and renal conditions [7,8], emergency department presentations [9], and ambulance call-outs [10], though effect sizes are generally modest and vary by health outcome, age group, and region.

No study has investigated whether ambient temperature systematically alters population-level medication dispensing patterns. Prior work has focused exclusively on how specific medications increase individual heat vulnerability---for example, antipsychotics impairing thermoregulation [11] and heat-sensitising drugs elevating hospitalisation risk [12]---but the reverse direction remains unexplored: how does temperature affect the volume and pattern of medication use across populations? This matters because dispensing data reflects the full healthcare response---from new diagnoses to worsening symptoms to changed prescribing---not just the acute presentations captured by emergency department or hospital data.

Recent studies have documented socioeconomic gradients in heat-related health impacts across diverse settings. In Western Australia, dose-response relationships between heatwave intensity and health service usage were steeper in disadvantaged and remote populations [13,14]. In the Paris region, women aged 65 years and over in the most deprived municipalities had heat-related mortality risk ratios up to 24% higher than those in less deprived areas [15], and a multi-city study in Latin America found that lower socioeconomic status amplified temperature-related mortality [16]. Given the large socioeconomic and geographic disparities in healthcare access across Australia, the equity dimension of climate-health impacts is particularly salient. As extreme heat events become more frequent and intense [25], understanding these disparities becomes increasingly urgent.

The 2019--20 Black Summer bushfire season in Australia was unprecedented in scale, burning approximately 24 million hectares and producing smoke plumes that blanketed major population centres for weeks [17]. Bushfire smoke contains fine particulate matter (PM$_{2.5}$) and other toxic compounds known to cause respiratory inflammation and exacerbate pre-existing cardiorespiratory conditions [18,19]. While studies have documented increased emergency department presentations and hospital admissions during bushfire smoke events [10,20], the impact on pharmaceutical dispensing---which may capture less severe but more widespread effects on respiratory health---has not been investigated.

The distributed lag non-linear model (DLNM) framework, developed by Gasparrini et al. [21], provides a flexible methodology for simultaneously modelling the non-linear shape of the exposure-response relationship and the lagged structure of delayed effects. Originally applied to temperature--mortality associations, DLNMs have been widely adopted in environmental epidemiology to capture the complex temporal dynamics of climate-health relationships [22].

In this study, we present the first population-level analysis of ambient temperature and medication dispensing, using weekly PBS data across 88 SA4 regions in Australia from 2013 to 2022. We estimate cumulative associations between maximum temperature and dispensing of cardiovascular, mental health, and respiratory medications using DLNMs, with stratification by socioeconomic disadvantage to identify equity-relevant effect modification. In a secondary analysis, we use a difference-in-differences design to estimate the impact of the 2019--20 Black Summer bushfire smoke on respiratory medication dispensing.

# Methods

## Study design and setting

We conducted a nationwide ecological time-series study using a panel of 88 Australian Statistical Area Level 4 (SA4) regions observed weekly from 7 January 2013 to 20 June 2022 (494 weeks), yielding 608,608 potential SA4-week observations. SA4 regions are large geographic areas defined by the Australian Bureau of Statistics (ABS) that represent functional labour markets or distinct rural/regional areas, with populations generally between 100,000 and 500,000 people [26]. The study period was determined by the availability of geographically disaggregated PBS dispensing data from the AIHW.

## Medication dispensing data

Weekly counts of PBS-subsidised prescriptions were obtained from the Australian Institute of Health and Welfare (AIHW) PHE-338 dataset [27] for the following Anatomical Therapeutic Chemical (ATC) classification groups: cardiovascular (C), mental health (N05/N06), and respiratory (R03). Respiratory medications were further classified into relievers, preventers, COPD/other, and oral corticosteroids based on the AIHW therapeutic groupings (Supplementary Table S1). Over the study period, dispensing totals were approximately 992 million cardiovascular prescriptions, 356 million mental health prescriptions, and 123 million respiratory prescriptions (Table 1). To protect confidentiality, the AIHW suppresses cells with fewer than 5 dispensings; after excluding these suppressed SA4-weeks, 34,580 observations per medication group were available for analysis (approximately 79% of the maximum possible 88 $\times$ 494 = 43,472 SA4-weeks). Suppressed cells were concentrated in smaller, more remote SA4 regions with low dispensing volumes. Antidepressants (N06A) and anxiolytics (N05B) were designated *a priori* as negative control outcomes [33], as there is no plausible biological mechanism linking ambient temperature to the underlying conditions for which these medications are prescribed.

## Climate data

Daily maximum temperature (Tmax) and precipitation were obtained from the Bureau of Meteorology's Australian Gridded Climate Data (AGCD) product [28], which provides interpolated daily observations on an approximately 5 km $\times$ 5 km grid derived from station data using Barnes successive-correction analysis. For each SA4 region, area-weighted spatial means were computed using the regionmask Python package for polygon masking, then aggregated to weekly means (temperature) or totals (precipitation).

## Socioeconomic and geographic data

Area-level socioeconomic disadvantage was measured using the Australian Bureau of Statistics Index of Relative Socio-economic Disadvantage (IRSD), derived from the 2016 Census [29]. SA4 regions were classified into quintiles from Q1 (most disadvantaged) to Q5 (least disadvantaged). Remoteness was classified according to the ABS Remoteness Structure (Supplementary Table S4).

## Air quality data

For the bushfire analysis, PM$_{2.5}$ data were obtained from the OpenAQ platform and spatially interpolated to SA4 centroids using inverse distance weighting. Bushfire smoke exposure was classified by comparing mean PM$_{2.5}$ during the Black Summer period (October 2019 to February 2020) with pre-bushfire baselines: SA4 regions with PM$_{2.5}$ ratios $\geq$ 2.0 were classified as "high" smoke exposure, 1.5--2.0 as "medium," and < 1.5 as "low."

## Statistical analysis

### Primary analysis: Distributed lag non-linear models

For each medication group, we fitted quasi-Poisson generalised linear models with a DLNM cross-basis specification [30]. The cross-basis used natural cubic splines for the exposure dimension (knots at the 10th, 50th, and 90th percentiles of temperature, following standard recommendations [30]) and natural cubic splines for the lag dimension (0--12 weeks, with 3 log-spaced internal knots). Models adjusted for long-term trends (natural cubic spline with degrees of freedom set to twice the number of study years), seasonality (2 Fourier harmonic pairs), weekly total precipitation, and SA4 fixed effects. Cumulative rate ratios (RR) and 95% confidence intervals were computed at the 90th, 95th, and 99th percentiles relative to the median temperature (22.8$^\circ$C).

### Equity stratification

Models were re-fitted separately within each IRSD quintile. Effect modification was formally tested using an interaction model with cross-basis $\times$ IRSD quintile product terms, assessed via an F-test.

### Bushfire difference-in-differences

The impact of 2019--20 Black Summer bushfire smoke on respiratory dispensing was estimated using a difference-in-differences framework comparing high-smoke and medium-smoke SA4 regions with low-smoke regions during the Black Summer period, adjusting for temperature, precipitation, and SA4 fixed effects. The parallel trends assumption was assessed by comparing respiratory dispensing trajectories across smoke exposure groups in the 12 months preceding the Black Summer (October 2018 -- September 2019); trends were visually parallel and a formal interaction test for differential pre-period trends was not significant (p = 0.31).

### Sensitivity and validation analyses

We conducted 11 pre-specified sensitivity analyses: (1) Benjamini-Hochberg and Bonferroni multiple testing corrections; (2) negative control outcomes; (3) Newey-West HAC standard errors (bandwidth = 12 weeks); (4) five alternative lag specifications; (5) Fourier harmonic sensitivity; (6) Cook's distance influence diagnostics; (7) cross-basis sensitivity; (8) Moran's I spatial autocorrelation tests; (9) E-values for unmeasured confounding; (10) leave-one-out jurisdiction analysis; and (11) QAIC model comparison. Full details of all sensitivity analyses are provided in the Supplementary Material (Supplementary Tables S6--S13, Supplementary Figures S1--S8).

Attributable burden was computed following the backward perspective framework of Gasparrini and Leone [31]. All analyses were conducted in R (version 4.3+) using the dlnm, splines, sandwich, and lmtest packages.

# Results

## Study population and descriptive statistics

The analysis encompassed 88 SA4 regions across all eight Australian states and territories, observed over 494 weeks (January 2013 to June 2022). After exclusion of AIHW-suppressed cells (counts < 5), 34,580 observations per medication group were available for the primary analysis (79.5% of maximum possible; Table 1). Excluded SA4-weeks were concentrated in remote and sparsely populated regions.

Weekly mean maximum temperature across all SA4-weeks was 22.8$^\circ$C (SD 6.5), with the 5th, 50th, and 95th percentiles at 12.3$^\circ$C, 22.8$^\circ$C, and 33.7$^\circ$C, respectively. Mean weekly precipitation was 14.6 mm (SD 23.1). Medication group definitions are provided in Supplementary Table S1, with additional descriptive statistics in Supplementary Table S2.

\begin{table}[H]
\centering
\caption{Study characteristics and descriptive statistics}
\label{tab:descriptive}
\begin{tabular}{ll}
\toprule
\textbf{Characteristic} & \textbf{Value} \\
\midrule
Study period & 7 January 2013 -- 20 June 2022 \\
SA4 regions & 88 \\
Total weeks & 494 \\
SA4-week observations & 608,608 \\
Observations per model (after exclusions) & 34,580 \\
Jurisdictions & 8 (NSW, VIC, QLD, SA, WA, TAS, NT, ACT) \\
\midrule
\textbf{Weekly dispensing counts, mean (SD)} & \\
\quad Cardiovascular (total) & 22,821 (11,338) \\
\quad Mental health (total) & 8,198 (4,140) \\
\quad Respiratory (total) & 2,839 (1,496) \\
\quad\quad Relievers & 818 (504) \\
\quad\quad Preventers & 1,352 (715) \\
\quad\quad COPD/Other & 487 (263) \\
\quad\quad Oral corticosteroids & 182 (112) \\
\midrule
\textbf{Climate variables, mean (SD)} & \\
\quad Weekly mean Tmax ($^\circ$C) & 22.8 (6.5) \\
\quad Weekly total precipitation (mm) & 14.6 (23.1) \\
\quad Heatwave days per week & 0.17 (0.79) \\
\quad Days above 95th percentile per week & 0.37 (0.96) \\
\bottomrule
\end{tabular}
\end{table}

## Primary temperature--dispensing associations

Heat exposure was not significantly associated with cardiovascular medication dispensing (cumulative RR at 95th percentile: 1.004, 95% CI: 0.997--1.011; p = 0.25) or mental health dispensing (RR 1.003, 95% CI: 0.997--1.010; p = 0.32) (Table 2, Figure 1). Both showed flat exposure-response curves across the temperature distribution.

Heat was strongly and inversely associated with respiratory medication dispensing (RR at 95th percentile: 0.910, 95% CI: 0.898--0.923; p = 1.5 $\times$ 10$^{-42}$), representing a 9.0% reduction. This was dose-dependent: RR 0.925 (95% CI: 0.914--0.936) at the 90th percentile and RR 0.893 (95% CI: 0.872--0.915) at the 99th percentile (Figure 1). Reliever dispensing showed the largest reduction (RR 0.868, 95% CI: 0.853--0.884; 13.2% decrease), while preventer dispensing decreased by 10.3% (RR 0.898, 95% CI: 0.885--0.911).

At the cold end, respiratory reliever dispensing at the 5th percentile was modestly elevated (RR 1.016, 95% CI: 1.001--1.032), consistent with winter increases in respiratory illness. The minimum morbidity temperature (MMT)---the temperature at which dispensing is lowest---for respiratory medications was at the upper extreme of the distribution (42.7$^\circ$C, 100th percentile), confirming a monotonically decreasing relationship across the observed temperature range.

\begin{table}[H]
\centering
\caption{Cumulative rate ratios for temperature--dispensing associations (primary DLNM models)}
\label{tab:primary_rr}
\small
\begin{tabular}{lcccc}
\toprule
\textbf{Medication group} & \textbf{RR at p90 (95\% CI)} & \textbf{RR at p95 (95\% CI)} & \textbf{RR at p99 (95\% CI)} & \textbf{p-value (p95)} \\
\midrule
Cardiovascular & 1.001 (0.995--1.007) & 1.004 (0.997--1.011) & 1.013 (1.000--1.025) & 0.25 \\
Mental health & 1.001 (0.995--1.007) & 1.003 (0.997--1.010) & 1.010 (0.998--1.021) & 0.32 \\
\textbf{Respiratory} & \textbf{0.925 (0.914--0.936)} & \textbf{0.910 (0.898--0.923)} & \textbf{0.893 (0.872--0.915)} & \textbf{1.5 $\times$ 10$^{-42}$} \\
\quad Relievers & 0.885 (0.872--0.900) & 0.868 (0.853--0.884) & 0.856 (0.830--0.883) & 2.5 $\times$ 10$^{-55}$ \\
\quad Preventers & 0.915 (0.904--0.927) & 0.898 (0.885--0.911) & 0.874 (0.852--0.897) & 1.3 $\times$ 10$^{-49}$ \\
Antidepressants\textsuperscript{a} & 0.999 (0.993--1.004) & 0.999 (0.993--1.005) & 1.001 (0.990--1.012) & 0.82 \\
Anxiolytics\textsuperscript{a} & 1.008 (0.999--1.017) & 1.007 (0.997--1.018) & 1.004 (0.984--1.023) & 0.18 \\
\bottomrule
\multicolumn{5}{l}{\footnotesize Rate ratios are cumulative over 0--12 week lags, centred at the median temperature (22.8$^\circ$C).} \\
\multicolumn{5}{l}{\footnotesize Bold indicates p < 0.05 after Bonferroni correction. \textsuperscript{a}Negative control outcomes.} \\
\end{tabular}
\end{table}

\begin{figure}[H]
\centering
\includegraphics[width=\textwidth]{figures/fig1_exposure_response.png}
\caption{Cumulative exposure-response curves for temperature and medication dispensing. Panels show the cumulative rate ratio (0--12 week lag) for cardiovascular, mental health, total respiratory, respiratory relievers, and respiratory preventers medications as a function of weekly mean maximum temperature, centred at the median (22.8$^\circ$C). Shaded areas represent 95\% confidence intervals. Dashed vertical lines indicate the 5th, 50th, and 95th percentiles of the temperature distribution.}
\label{fig:er_curves}
\end{figure}

## Equity stratification

The inverse temperature--respiratory dispensing association showed significant effect modification by socioeconomic disadvantage (interaction F-test: F = 11.3, p = 6.4 $\times$ 10$^{-11}$) (Figure 2). The magnitude of the heat-associated reduction varied across IRSD quintiles: Q1 (most disadvantaged) RR 0.920 (95% CI: 0.892--0.949); Q3 RR 0.892 (95% CI: 0.865--0.921); Q5 (least disadvantaged) RR 0.884 (95% CI: 0.847--0.922). Full quintile-specific results are in Supplementary Table S3. No significant equity interaction was observed for cardiovascular (p = 0.66) or mental health (p = 0.48) medications.

Among respiratory medications, heat exposure above the median was associated with an estimated 1.67 million fewer dispensed prescriptions over the study period (attributable fraction --1.73%). Attributable burden by quintile is reported in Supplementary Table S5.

\begin{figure}[H]
\centering
\includegraphics[width=0.85\textwidth]{figures/fig2_equity_forest.png}
\caption{Equity stratification of the temperature--respiratory dispensing association. Forest plot showing cumulative rate ratios at the 95th percentile of temperature for total respiratory medication dispensing, stratified by IRSD quintile (Q1 = most disadvantaged, Q5 = least disadvantaged). Error bars represent 95\% confidence intervals.}
\label{fig:equity}
\end{figure}

## Black Summer bushfire analysis

The DiD analysis revealed significant increases in respiratory medication dispensing in high-smoke regions during the Black Summer (Figure 3): total respiratory RR 1.062 (95% CI: 1.046--1.078; p = 2.4 $\times$ 10$^{-15}$), relievers RR 1.093 (95% CI: 1.073--1.114; p = 4.2 $\times$ 10$^{-20}$), preventers RR 1.071 (95% CI: 1.055--1.087; p = 3.0 $\times$ 10$^{-18}$), and COPD/Other RR 1.013 (95% CI: 1.000--1.026; p = 0.047). Medium-smoke regions showed smaller, mostly non-significant effects (Supplementary Table S7).

\begin{figure}[H]
\centering
\includegraphics[width=0.85\textwidth]{figures/fig3_bushfire.png}
\caption{Black Summer bushfire impact on respiratory medication dispensing. Difference-in-differences rate ratios comparing high-smoke and medium-smoke SA4 regions with low-smoke regions during the 2019--20 Black Summer period (October 2019 -- February 2020). Error bars represent 95\% confidence intervals.}
\label{fig:bushfire}
\end{figure}

## Sensitivity and validation analyses

The respiratory finding was robust across all 11 sensitivity analyses (Table 3, Supplementary Tables S6--S13, Supplementary Figures S5--S8). All three primary respiratory findings survived both Benjamini-Hochberg and Bonferroni corrections (adjusted p < 10$^{-40}$). Negative control outcomes (antidepressants, anxiolytics) showed null associations at all temperature percentiles (Supplementary Figure S3).

Newey-West HAC standard errors were 1.57--1.77 times larger than model-based SEs. All respiratory findings remained significant with HAC-adjusted CIs (total respiratory HAC 95% CI: 0.890--0.932; Supplementary Table S8). RR estimates were stable across five lag specifications (range 0.894--0.911; Supplementary Table S9) and two Fourier pair settings (Supplementary Table S10). After excluding high-influence observations (Cook's D top 5%), the respiratory RR shifted modestly from 0.910 to 0.922 (Supplementary Table S11). Moran's I tests indicated moderate spatial autocorrelation in residuals (I = 0.28--0.35; Supplementary Table S12). E-values for the respiratory finding were 1.43 (point estimate) and 1.39 (CI bound; Supplementary Table S13), indicating moderate robustness to unmeasured confounding --- an unmeasured confounder associated with both temperature and respiratory dispensing at RR $\geq$ 1.43 could explain the observed association, though it would need to be unrelated to the covariates already adjusted for. Leave-one-out jurisdiction analysis showed respiratory RR ranged from 0.896 to 0.917 with no single state driving the result (Supplementary Figure S7). When the study period was split into early (2013--2017) and late (2018--2022) halves, the inverse temperature--respiratory association was stronger in the later period (RR 0.885 vs 0.909; interaction p = 0.0001; Supplementary Table S14), potentially reflecting changes in respiratory disease epidemiology or prescribing practices over time.

\begin{table}[H]
\centering
\caption{Sensitivity analysis summary for respiratory medication dispensing (RR at 95th percentile)}
\label{tab:sensitivity}
\small
\begin{tabular}{llr}
\toprule
\textbf{Analysis} & \textbf{RR (95\% CI)} & \textbf{Notes} \\
\midrule
Primary model & 0.910 (0.898--0.923) & NS lag, 0--12 weeks \\
Newey-West HAC SE & 0.910 (0.890--0.932) & SE ratio 1.71 \\
Unconstrained indicator lags & 0.911 (0.899--0.924) & QAIC +135 \\
Polynomial DL (degree 3) & 0.911 (0.899--0.923) & QAIC $-$7 \\
Shorter lag (0--8 weeks) & 0.894 (0.883--0.904) & QAIC +18 \\
Longer lag (0--16 weeks) & 0.902 (0.889--0.916) & QAIC +3 \\
3 Fourier pairs & 0.912 (0.900--0.923) & QAIC $-$3,530 \\
Excluding high-influence obs & 0.922 (0.910--0.935) & 4.1\% removed \\
Bonferroni-corrected p & --- & p = 2.8 $\times$ 10$^{-41}$ \\
E-value (point / CI bound) & 1.43 / 1.39 & \\
\bottomrule
\end{tabular}
\end{table}

\begin{figure}[H]
\centering
\includegraphics[width=\textwidth]{figures/fig4_sensitivity_forest.png}
\caption{Sensitivity analysis forest plot for respiratory medication dispensing. Cumulative rate ratios at the 95th percentile of temperature across primary and sensitivity model specifications. Error bars represent 95\% confidence intervals. The dashed vertical line indicates RR = 1 (null).}
\label{fig:sensitivity}
\end{figure}

# Discussion

## Principal findings

This study provides the first population-level evidence linking ambient temperature to pharmaceutical dispensing patterns. Using nearly a decade of nationwide PBS data across 88 SA4 regions, we found that heat exposure was strongly associated with reduced respiratory medication dispensing but not with cardiovascular or mental health prescribing. The inverse temperature--respiratory association was dose-dependent, consistent across lag specifications, robust to multiple testing corrections and HAC-adjusted standard errors, and specific to respiratory medications (with null negative control outcomes). Socioeconomic disadvantage significantly modified the respiratory association, with the most disadvantaged areas showing the smallest heat-related reductions. In a secondary analysis, bushfire smoke exposure during the 2019--20 Black Summer was associated with a 6.2% increase in respiratory dispensing in high-smoke regions, with the largest effects on reliever medications.

## Novelty and positioning in the literature

To our knowledge, this is the first study to systematically examine temperature--prescribing associations at the population level. Previous studies have investigated how specific medications alter individual heat vulnerability---for example, antipsychotics increasing heat-related mortality risk [11], or heat-sensitising medications elevating hospitalisation risk during heat events [12]. A pharmacovigilance analysis identified heat stroke and exhaustion cases associated with various drug classes [23]. However, all prior work has focused on the medication $\rightarrow$ heat vulnerability pathway, whereas our study examines the reverse direction: how temperature affects population-level medication dispensing as a marker of disease burden.

## Interpretation of findings

### Respiratory medications and temperature

The strong inverse association (9% reduction at the 95th percentile) likely reflects the well-established seasonality of respiratory viral infections rather than a direct heat effect on airway pathology. Respiratory syncytial virus, influenza, and other common viral pathogens show reduced transmission in warmer conditions [32], and temperature is one of the key environmental predictors of epidemic timing across both temperate and tropical climates. The observation that reliever dispensing (reflecting acute symptoms) showed the largest reduction (13.2%) while preventer dispensing (more maintenance-oriented) showed a smaller reduction (10.3%) is consistent with reduced acute exacerbation frequency. The RR estimate was stable across 2 and 3 Fourier harmonic pairs (0.910 vs 0.912), suggesting the association is not an artefact of inadequate seasonal adjustment, though we cannot fully exclude residual seasonal confounding (see Limitations). A supply-side mechanism --- whereby seasonal patterns in GP consultations, medication reviews, or repeat prescription timing contribute to the observed variation --- is also plausible and cannot be distinguished from demand-side changes using aggregate dispensing data alone.

### Null cardiovascular and mental health findings

Cardiovascular medications are primarily prescribed for chronic conditions (hypertension, hyperlipidaemia) where dispensing reflects long-term management rather than acute events. The weekly temporal resolution may be too coarse to capture acute cardiovascular events. For mental health medications, the null finding is expected given the chronic nature of prescribing, and their successful use as negative controls supports this interpretation.

### Equity gradient in respiratory dispensing

The significant interaction (p = 6.4 $\times$ 10$^{-11}$) revealed a gradient where the most disadvantaged areas (Q1) showed the smallest heat-related dispensing reductions and the least disadvantaged (Q5) showed the largest. The most likely explanation is that higher baseline chronic respiratory disease prevalence in disadvantaged areas makes dispensing less temperature-responsive --- if a larger proportion of prescriptions reflects persistent disease management rather than seasonal exacerbations, the seasonal temperature signal is attenuated. Healthcare access barriers may also reduce symptom-responsive medication adjustment in disadvantaged populations, while environmental co-exposures (indoor air quality, occupational exposures, tobacco smoke) could maintain respiratory burden regardless of ambient temperature. These findings suggest that climate warming may widen existing inequities in respiratory disease burden.

### Temporal strengthening of the respiratory association

The inverse temperature--respiratory dispensing association was stronger in the later study period (2018--2022 RR 0.885 vs 2013--2017 RR 0.909; interaction p = 0.0001). Several factors may contribute. COVID-19 public health measures (lockdowns, mask mandates, social distancing) from early 2020 substantially reduced respiratory viral transmission, which could amplify the apparent temperature sensitivity of dispensing during the later period. Changes in prescribing practices, including the shift toward telehealth consultations during 2020--2022, may also have altered seasonal dispensing patterns. Whether this temporal trend reflects genuine climate adaptation or period-specific confounding cannot be determined from these data.

### Black Summer bushfire impacts

The 6.2% increase in respiratory dispensing (9.3% for relievers) in high-smoke regions provides population-level evidence complementing prior studies of emergency presentations [10,17,18]. The dose-response gradient from low to medium to high smoke zones supports a causal interpretation. The dominant reliever effect is consistent with acute bronchospasm from particulate matter inhalation.

## Strengths and limitations

The principal strengths of this study are its novel use of population-level dispensing data as a climate-health endpoint and the scale of the dataset (88 SA4 regions, 494 weeks, approximately 1.47 billion dispensed prescriptions). The DLNM framework provides flexible characterisation of non-linear, delayed exposure-response relationships, and the comprehensive sensitivity programme---including HAC standard errors, negative controls, E-values, and leave-one-out jurisdiction analysis---provides confidence in the robustness of the primary findings. The pre-specified equity stratification and the secondary bushfire natural experiment add policy-relevant dimensions.

Several limitations warrant consideration. As an ecological study, we cannot make individual-level causal inferences; the associations we observe operate at the SA4-region level and may not reflect individual-level exposure-response relationships. The weekly temporal resolution, while appropriate for DLNM modelling, may obscure sub-weekly dynamics such as acute heat events lasting 2--3 days. Residual spatial autocorrelation in model residuals (Moran's I 0.28--0.35, significant in approximately 90% of weeks sampled) is only partially addressed by SA4 fixed effects and Newey-West temporal HAC standard errors; spatially correlated errors inflate the effective sample size and may produce confidence intervals that are narrower than warranted. Future work should consider spatially explicit error structures such as spatial random effects or generalised estimating equations with spatial correlation.

The respiratory finding, while robust to seasonal adjustment with Fourier harmonics, may partly capture residual seasonality not fully absorbed by 2--3 harmonic pairs; the substantially improved QAIC with 3 pairs ($\Delta$ = $-$3,530) indicates that seasonal variation is complex and incompletely modelled. We note, however, that additional Fourier pairs would absorb variation at frequencies that overlap with the temperature signal, making it difficult to fully separate temperature effects from seasonality in any observational design. PBS data exclude private (non-subsidised) prescriptions, which may introduce differential capture across socioeconomic groups, and we cannot distinguish new prescriptions from repeats. Sparse PM$_{2.5}$ monitoring in rural Australia limits the precision of bushfire exposure classification. Attributable burden estimates assume a causal relationship and should be interpreted as indicative rather than definitive.

## Implications

If the inverse temperature--respiratory dispensing association reflects genuine temperature sensitivity of viral transmission, rising mean temperatures could compress seasonal dispensing peaks, requiring adjustment to PBS stockpiling protocols that are currently calibrated to historical demand cycles. However, this projection depends on whether inter-annual warming replicates the within-year seasonal pattern, which remains uncertain. The equity gradient indicates that climate adaptation strategies for respiratory health should prioritise socioeconomically disadvantaged communities, where baseline disease burden is higher and temperature-responsive dispensing reductions are smaller. The Black Summer bushfire analysis demonstrates that smoke events generate acute pharmaceutical demand surges --- particularly for reliever medications --- that should be factored into emergency pharmaceutical supply planning. More broadly, population-level dispensing data from pharmaceutical benefits schemes represents an underutilised and readily available resource for environmental health surveillance and early warning systems.

# Conclusions

Ambient temperature is strongly and inversely associated with population-level respiratory medication dispensing, likely mediated by the temperature dependence of viral respiratory pathogen transmission. This association is modified by socioeconomic disadvantage, with the smallest seasonal variation in the most disadvantaged areas where chronic respiratory disease burden is highest. Bushfire smoke produces acute, dose-dependent increases in respiratory dispensing. These findings establish medication dispensing data as a novel and readily accessible endpoint for climate-health research. Replication in other countries with pharmaceutical benefits schemes and investigation of whether long-term warming trends alter dispensing patterns are priorities for future work.

\newpage

# Declarations

**Ethics approval:** This study used publicly available, aggregated administrative data with no individual-level information. Ethical approval was not required under Australian National Health and Medical Research Council guidelines for research involving collections of non-identifiable data.

**Consent to participate:** Not applicable (aggregate data with no individual participants).

**Data availability:** PBS dispensing data are available from the AIHW PHE-338 dataset (aihw.gov.au). AGCD climate data are available from the Bureau of Meteorology (bom.gov.au). PM$_{2.5}$ data are available from OpenAQ (openaq.org). SEIFA indices and SA4 boundaries are available from the Australian Bureau of Statistics (abs.gov.au). Analysis code is available at [repository URL to be inserted].

**Author contributions:** HF conceived the study, acquired and processed the data, conducted all analyses, and wrote the manuscript.

**Funding:** This research received no specific funding.

**Conflicts of interest:** The author declares no conflicts of interest.

**AI disclosure:** AI tools (Claude, Anthropic) were used to assist with manuscript drafting, code development, and literature search. The author takes full responsibility for the content and has verified all analyses, interpretations, and factual claims.

\newpage

# References

\small

1. Watts N, Amann M, Arnell N, et al. The 2019 report of The Lancet Countdown on health and climate change: ensuring that the health of a child born today is not defined by a changing climate. *Lancet.* 2019;394(10211):1836-1878.

2. Romanello M, McGushin A, Di Napoli C, et al. The 2021 report of the Lancet Countdown on health and climate change: code red for a healthy future. *Lancet.* 2021;398(10311):1619-1662.

3. Amoatey P, Xu Z, Odebeatu CC, Singh N, Osborne NJ, Phung D. Impact of extreme heat on health in Australia: a scoping review. *BMC Public Health.* 2025;25:522.

4. Gasparrini A, Guo Y, Hashizume M, et al. Mortality risk attributable to high and low ambient temperature: a multicountry observational study. *Lancet.* 2015;386(9991):369-375.

5. Tobias A, Hashizume M, Honda Y, et al. Geographical variations of the minimum mortality temperature at a global scale: a multicountry study. *Environ Epidemiol.* 2021;5(5):e169.

6. Cheng J, Xu Z, Bambrick H, et al. Cardiorespiratory effects of heatwaves: a systematic review and meta-analysis of global epidemiological evidence. *Environ Res.* 2019;177:108610.

7. Wilson LA, Morgan GG, Hanigan IC, et al. The impact of heat on mortality and morbidity in the Greater Metropolitan Sydney Region: a case crossover analysis. *Environ Health.* 2013;12:98.

8. Khalaj B, Lloyd G, Sheppeard V, Dear K. The health impacts of heat waves in five regions of New South Wales, Australia: a case-only analysis. *Int Arch Occup Environ Health.* 2010;83(7):833-842.

9. Goldie J, Alexander L, Lewis SC, Sherwood SC. Comparative evaluation of human heat stress indices on selected hospital admissions in Sydney, Australia. *Aust N Z J Public Health.* 2017;41(4):381-387.

10. Mason H, King J, Peden AE, Franklin RC. Systematic review of the impact of heatwaves on health service demand in Australia. *BMC Health Serv Res.* 2022;22:960.

11. Chen SX, Lee MJ, McVea DA, et al. Antipsychotics and other risk factors for mortality among people with schizophrenia during an extreme heat event: a population-based case-control study. *Sci Rep.* 2025;15:34505.

12. Layton JB, Li W, Yuan J, Gilman JP, Horton DB, Setoguchi S. Heatwaves, medications, and heat-related hospitalization in older Medicare beneficiaries with chronic conditions. *PLoS One.* 2020;15(12):e0243665.

13. Xiao J, Spicer T, Jian L, et al. Variation in population vulnerability to heat wave in Western Australia. *Front Public Health.* 2017;5:64.

14. Patel D, Jian L, Xiao J, et al. Joint effect of heatwaves and air quality on emergency department attendances for vulnerable population in Perth, Western Australia, 2006-2015. *Environ Res.* 2019;174:80-87.

15. Pascal M, Goria S, Forceville G, et al. Analyzing effect modifiers of the temperature-mortality relationship in the Paris Region to identify social and environmental levers for more effective adaptation to heat. *Health Place.* 2024;89:103325.

16. Bakhtsiyarava M, Schinasi LH, Sanchez BN, et al. Modification of temperature-related human mortality by area-level socioeconomic and demographic characteristics in Latin American cities. *Soc Sci Med.* 2022;317:115526.

17. Walter CM, Schneider-Futschik EK, Knibbs LD, Irving LB. Health impacts of bushfire smoke exposure in Australia. *Respirology.* 2020;25(5):495-501.

18. Rodney RM, Swaminathan A, Calear AL, et al. Physical and mental health effects of bushfire and smoke in the Australian Capital Territory 2019-20. *Front Public Health.* 2021;9:682402.

19. Zhang Y, Ye T, Huang W, et al. Health impacts of wildfire smoke on children and adolescents: a systematic review and meta-analysis. *Curr Environ Health Rep.* 2024;11(1):46-60.

20. Heaney E, Hunter L, Clulow A, Bowles D, Vardoulakis S. Efficacy of communication techniques and health outcomes of bushfire smoke exposure: a scoping review. *Int J Environ Res Public Health.* 2021;18(20):10889.

21. Gasparrini A, Armstrong B, Kenward MG. Distributed lag non-linear models. *Stat Med.* 2010;29(21):2224-2234.

22. He C, Breitner S, Zhang S, et al. Nocturnal heat exposure and stroke risk. *Eur Heart J.* 2024;45(24):2158-2166.

23. Gamboa L, Sáez Lafuente A, Morera-Herreras T, Garcia M, Aguirre C, Lertxundi U. Analysis of heat stroke and heat exhaustion cases in EudraVigilance pharmacovigilance database. *Eur J Clin Pharmacol.* 2023;79(5):679-685.

24. Li M, Gu S, Bi P, Yang J, Liu Q. Heat waves and morbidity: current knowledge and further direction --- a comprehensive literature review. *Int J Environ Res Public Health.* 2015;12(5):5256-5283.

25. Zhou S, Wu Y, Liu Y, et al. Lethal heat and humidity events across the globe. *Annu Rev Environ Resour.* 2025;50(1):247-272.

26. Australian Bureau of Statistics. Australian Statistical Geography Standard (ASGS): volume 1 --- main structure and greater capital city statistical areas, July 2021. ABS cat. no. 1270.0.55.001. Canberra: ABS; 2021.

27. Australian Institute of Health and Welfare. Geography and time-specific health data for environmental analysis (PHE-338). Canberra: AIHW; 2023. Available from: https://www.aihw.gov.au/reports/burden-of-disease/health-data-environmental-analysis

28. Jones DA, Wang W, Fawcett R. High-quality spatial climate data-sets for Australia. *Aust Meteorol Oceanogr J.* 2009;58(4):233-248.

29. Australian Bureau of Statistics. Socio-Economic Indexes for Areas (SEIFA), Technical Paper, 2016. ABS cat. no. 2033.0.55.001. Canberra: ABS; 2018.

30. Gasparrini A. Modeling exposure-lag-response associations with distributed lag non-linear models. *Stat Med.* 2014;33(5):881-899.

31. Gasparrini A, Leone M. Attributable risk from distributed lag models. *BMC Med Res Methodol.* 2014;14:55.

32. Tamerius JD, Shaman J, Alonso WJ, et al. Environmental predictors of seasonal influenza epidemics across temperate and tropical climates. *PLoS Pathog.* 2013;9(3):e1003194.

33. Lipsitch M, Tchetgen Tchetgen E, Cohen T. Negative controls: a tool for detecting confounding and bias in observational studies. *Epidemiology.* 2010;21(3):383-388.
