# Sources used in the research

**Every entry below was used in the research; to be read and verified by the author before citation.**

This is a record, not a bibliography. It lists each published method, model, dataset and paper that the code
implements, reproduces, or uses for a parameter or a verification target, with the identifier exactly as
recorded in the project's source files (`final_pipeline.py` docstrings, `experiments/*/README.md`,
`experiments/phase8_cardiac/MODELS.md`, `experiments/phase9_waveform/METHODS.md`,
`experiments/phase10_final/METHODS.md`, and module docstrings). Nothing was added that the code does not use.
Several primary papers were **not accessible** when the code was written; the project then worked from an open
restatement, as recorded in the "access" column. The author should read every primary source before citing it.
Full bibliographic details (journal, pages) were not always recorded; they are given as recorded.

## 1. Datasets (PhysioNet; licence: Open Data Commons Attribution License v1.0, per each database page)

| dataset | identifier | used in | citation PhysioNet asks for (from the database page, 2026-10-03) |
|---|---|---|---|
| MIT-BIH Arrhythmia Database 1.0.0 (mitdb) | doi:10.13026/C2F305 | Phase 7 (all), Phase 9 exploratory, Phase 10 development | Moody GB, Mark RG. The impact of the MIT-BIH Arrhythmia Database. IEEE Eng Med Biol 20(3):45-50 (2001), PMID 11446209 |
| MIT-BIH Normal Sinus Rhythm Database 1.0.0 (nsrdb) | doi:10.13026/C2NK5R | Phase 9 Part F, Phase 10 development | standard PhysioNet citation |
| BIDMC Congestive Heart Failure Database 1.0.0 (chfdb) | doi:10.13026/C29G60 | Phase 9 Part F, Phase 10 development | Baim DS et al. Survival of patients with severe congestive heart failure treated with oral milrinone. J Am Coll Cardiol 7(3):661-670 (1986) |
| Normal Sinus Rhythm RR Interval Database 1.0.0 (nsr2db) | doi:10.13026/C2S881 | Phase 10 confirmatory | standard PhysioNet citation |
| Congestive Heart Failure RR Interval Database 1.0.0 (chf2db) | doi:10.13026/C2F598 | Phase 10 confirmatory | standard PhysioNet citation |
| MIT-BIH Noise Stress Test Database 1.0.0 (nstdb) | doi:10.13026/C2HS3T | Phase 9 recording noise | Moody GB, Muldrow WE, Mark RG. A noise stress test for arrhythmia detectors. Computers in Cardiology 11:381-384 (1984) |
| QT Database 1.0.0 (qtdb) | doi:10.13026/C24K53 | Phase 9 delineator check | Laguna P, Mark RG, Goldberger AL, Moody GB. A database for evaluation of algorithms for measurement of QT and other waveform intervals in the ECG. Computers in Cardiology 24:673-676 (1997) |
| PhysioNet (standard citation for every database) | | all of the above | Goldberger AL et al. PhysioBank, PhysioToolkit, and PhysioNet. Circulation 101(23):e215-e220 (2000) |

## 2. ECG processing and beat detection

| source (as recorded) | used for | where recorded | access |
|---|---|---|---|
| Pan J, Tompkins WJ (1985), QRS detection algorithm | R-peak detector (`fp.detect_r_peaks`) | `final_pipeline.py` line ~395; `docs/PHASE1_PAN_TOMPKINS.md` | (Phase 1) |
| ANSI/AAMI EC57 (150 ms beat-match window) | R-peak quality and label transfer (Phase 7) | `experiments/phase7_mitbih/data.py`; Phase 7 report 2.1 | standard |
| McSharry PE, Clifford GD, Tarassenko L, Smith LA. A dynamical model for generating synthetic electrocardiogram signals. IEEE TBME 50(3):289-294 (2003); ECGSYN 1.0.0 on PhysioNet (`ecgsyn.c`) | Phase 5 RR spectrum (N1 generator); Phase 9 ECGSYN waveform generator | `experiments/phase5_rr/README.md` [McSharry 2003]; Phase 9 METHODS 1 | open (authors' HTML on PhysioNet) |
| LITFL "Premature Ventricular Complex" and "Premature Atrial Complex" (web pages) | ectopic beat morphology criteria (Phase 9) | Phase 9 METHODS 2 | open |
| Moody GB, Muldrow WE, Mark RG (1984); WFDB `nst` documentation (physionet.org/physiotools/wag/nst-1.htm) | SNR definition for added noise (Phase 9) | Phase 9 METHODS | open |
| Martínez JP et al. IEEE TBME 51:570 (2004), wavelet delineator, via NeuroKit2 0.2.13; Di Marco LY, Chiari L (2011) PMC3076264, Table 6 | delineator check on QTDB (Phase 9; **failed, dropped**) | Phase 9 METHODS 7 | primary blocked; restatement open |
| Laguna P et al. (1994), ecgpuwave | considered for delineation (Phase 9; **not run**: no Fortran compiler) | Phase 9 METHODS 7 | |

## 3. Nonlinear time-series methods

| source (as recorded) | used for | where recorded | access |
|---|---|---|---|
| So P, Ott E, Schiff SJ, Kaplan DT, Sauer T, Grebogi C. Phys Rev Lett 76:4705 (1996); So P, Ott E, Sauer T, Gluckman BJ, Grebogi C, Schiff SJ. Phys Rev E 55:5398 (1997) | UPO transform, surrogate significance, stability (pipeline, Phases 2E-4) | `final_pipeline.py` header and section 4 | (Phase 2) |
| Rosenstein MT, Collins JJ, De Luca CJ (1993) | largest Lyapunov exponent (pipeline; Phase 3 C1) | `final_pipeline.py` header; Phase 3 estimators | (Phase 2) |
| Cao L (1997), E1/E2 embedding dimension | embedding dimension (pipeline) | `final_pipeline.py` | (Phase 2) |
| Fraser AM, Swinney HL (1986), mutual information delay | TDMI delay (pipeline; Phase 9 surrogates) | `experiments/phase9_waveform/surrogates_wave.py` | |
| Theiler exclusion window | nearest-neighbour exclusion (pipeline) | `final_pipeline.py` header | |
| Kantz H (1994)-style divergence | Phase 3 candidate C2 | `experiments/phase3_lle/estimators.py` | |
| Gottwald GA, Melbourne I (2009), 0-1 test, modified correlation method | Phase 3 candidate C4 | `experiments/phase3_lle/estimators.py` | |
| Schreiber T, Schmitz A (1996), iterative AAFT surrogates | IAAFT surrogates (Phase 3 C1 = `lle_chaos_test`; K3) | `final_pipeline.py`; Phase 3 estimators | |
| AAFT surrogates (amplitude-adjusted Fourier transform) | production LLE test (pipeline baseline) | `final_pipeline.py` (`fp.run_surrogate_analysis`) | |
| Schreiber T, Schmitz A. Physica D 142:346 (2000), arXiv:chao-dyn/9909037, Sect. 4.6 | multivariate IAAFT (Phase 9 C2) | Phase 9 METHODS 7 | open (arXiv) |
| Sugihara G, May RM. Nature 344:734 (1990), restated in PMC9760897 | nonlinear prediction and forecast-error growth (Phase 8 C2/C3; K3 growth gate) | `experiments/phase9_waveform/candidates9.py`; Phase 9 METHODS 8b | primary not accessible; restatement open |
| Peltola MA (2012), PMC3358711 | ectopy masking/editing rule (K3; Phase 9) | Phase 9 METHODS 8b, `candidates9.py` | open |
| Gao J et al. Front Physiol (2011), PMC3264951 (SDLE; first defined PRE 74:066204, 2006) | D1 SDLE (Phase 9) | Phase 9 METHODS D1 | open; PRE blocked |
| Aurell E et al. (1997) arXiv:chao-dyn/9606014; Boffetta G et al. Phys Rep (2002) arXiv:nlin/0101029, Eq. 3.37; Cencini M et al. (2000) arXiv:nlin/0002018 | D2 FSLE; D3 (ε, τ)-entropy (Phase 9; K5 FSLE gate) | Phase 9 METHODS D2-D3 | open (arXiv) |
| Gaspard P, Wang XJ. Phys Rep 235:291 (1993) | D3 ε-entropy definition (Phase 9) | `experiments/phase9_waveform/measures.py` | blocked |
| Bandt C, Pompe B (2002); Rosso OA et al. (2007); ordpy (arXiv:2102.06786) | D4 permutation entropy / complexity-entropy plane (Phase 9) | Phase 9 METHODS D4 | primaries blocked; ordpy open |
| Marwan N. IJBC 21:1003 (2011), arXiv:1007.2215 | D5 RQA determinism (Phase 9) | Phase 9 METHODS D5 | open (arXiv) |
| Grassberger P et al. Chaos 3:127 (1993); Hegger R, Kantz H, Schreiber T. Chaos 9:413 (1999), arXiv:chao-dyn/9810005; TISEAN 3.0.1 (`ghkss.c`, `fsle`, `d2`) | D6 GHKSS noise reduction; reference programs for D2-D3 (Phase 9) | Phase 9 METHODS D6 | 1993 blocked; 1999 open |
| Small M, Yu D, Harrison RG. PRL 87:188101 (2001), via Luo X, Nakamura T, Small M, arXiv:nlin/0404054; TimeseriesSurrogates.jl `noiseradius` | pseudo-periodic surrogates (Phase 9 B2; **not verified, unusable**) | Phase 9 METHODS 5 | primary blocked |
| Theiler J. Phys Lett A 196:335 (1995), via Luo et al. | cycle-shuffle surrogates (Phase 9 B2; **not verified, unusable**) | Phase 9 METHODS 5 | primary blocked |
| Thiel M, Romano MC, Kurths J, Rolfs M, Kliegl R. Europhys Lett 75:535 (2006) and 2008 restatement (arXiv) | twin surrogates (Phase 9 B2; verified, uninformative for chaos) | Phase 9 METHODS 5 | open restatement |

## 4. Noise titration (Phases 8 and 10)

| source (as recorded) | used for | where recorded | access |
|---|---|---|---|
| Poon CS, Barahona M. Titration of chaos with added noise. PNAS 98:7107-7112 (2001), PMC34630 | titration algorithm, noise limit, verification targets (Figs. 1-5) | Phase 10 METHODS 2.1 [PB01]; Phase 8 DIAGNOSIS D4 | full text via PMC page |
| Wu GQ et al. Chaotic signatures of heart rate variability and its power spectrum in health, aging and heart failure. PLoS ONE 4:e4323 (2009) | VWK model (Eq. 1), criterion (Eq. 2), study design, Wu-style preprocessing, 12-min segments | Phase 10 METHODS 2.1 [Wu09]; copy in `experiments/phase10_final/sources/` (CC-BY) | open |
| Wysocki M et al. Respir Physiol Neurobiol 153:54 (2006), arXiv:nlin/0606032 | model selection, F-test at 1 %, routine defaults κ = 6, d = 3 | Phase 10 METHODS 2.1 [Wy06]; Phase 8 DIAGNOSIS D4 | open (arXiv) |
| Poon CS, Li C, Wu GQ. A unified theory of chaos ... arXiv:1004.1427 | restatement of Eqs. 1-2 | Phase 10 METHODS 2.1 [P10] | open |
| Barahona M, Poon CS. Detection of nonlinear dynamics in short, noisy time series. Nature 381:215 (1996) | primary definition of the indicator (taken from the restatements) | Phase 10 METHODS 2.1 [BP96] | **not accessible** (abstract only) |
| Poon CS, Merrill CK. Decrease of cardiac chaos in congestive heart failure. Nature 389:492 (1997) | the CHF claim tested in P2 | Phase 10 METHODS 2.1 [PM97] | **not accessible** (abstract only) |
| Korenberg's recursive orthogonal estimation (Volterra-Wiener-Korenberg model fitting) | least-squares fitting of the VWK terms | `experiments/phase10_final/titration10.py` | named in the restatements; no Korenberg paper recorded |

## 5. Synthetic systems and cardiac models

| source (as recorded) | used for | where recorded | access |
|---|---|---|---|
| Hénon M (1976). A two-dimensional mapping with a strange attractor | Hénon map (Phases 2E-10) | `experiments/phase5_rr/README.md` | |
| Logistic map (r = 3.5, 3.58-4) | null/periodic/chaotic test system | Phase 2E and later configs | textbook |
| Lorenz EN (1963). Deterministic nonperiodic flow. J Atmos Sci | G1 Lorenz maxima; titration verification | `experiments/phase5_rr/README.md` | |
| Rössler OE (1976). An equation for continuous chaos. Phys Lett A | G2 Rössler flow; Phase 9 surrogate checks | `experiments/phase5_rr/README.md` | |
| Mackey MC, Glass L (1977); Glass L, Mackey MC, "Mackey-Glass equation", Scholarpedia 5(3):6908; parameter set of Farmer (1982); exponents from arXiv:1810.01016 (Table 1) and arXiv:1509.06057 | G3; Phase 8 DEV family; titration verification | Phase 5 README; Phase 8 MODELS 3 | Science 1977 blocked; Scholarpedia open |
| Tong H, Lim KS (1980). Threshold autoregression, limit cycles and cyclical data | S1 SETAR null | `experiments/phase5_rr/README.md` | |
| Diagne K, Bury TM, Pettebone ME, ..., Glass L et al. Phase resetting in human stem cell derived cardiomyocytes explains complex cardiac arrhythmias. PLoS Comput Biol (2026), doi:10.1371/journal.pcbi.1013935, PMC12900431 (S1/S2 Text, S1-S3 Tables) | phase-resetting map (Phase 8 TEST family; Phase 10 spike-ins) | Phase 8 MODELS 1; supplements in `experiments/phase8_cardiac/sources/` | open |
| Sun J, Amellal F, Glass L, Billette J. J Theor Biol 173:79 (1995), equations via arXiv:math/0609106 Eq. (51) | AV-nodal conduction model (Phase 8 non-chaotic family) | Phase 8 MODELS 2 | primary blocked |
| da Silva Lima G, Savi MA, Bessa WM. Adaptive control of cardiac rhythms. Sci Rep 14:23446 (2024), PMC11458860 | three coupled modified van der Pol oscillators (Phase 8; Phase 10 spike-ins) | Phase 8 MODELS 4 | open |
| Gois & Savi (2009) model via PMC9938421 (Heliyon 2023) | candidate coupled-oscillator model (Phase 8; **dropped**: parameter table inconsistent) | Phase 8 MODELS 4-5 | open |
| Gall WM et al. arXiv:2202.12406 (KTz map; Kinouchi-Tragtenberg) | KTz action-potential morphology family (Phase 9 TEST) | Phase 9 METHODS 3 | open |
| Qu Z et al. Phys Rep (2014), PMC4175480 | basis for mapping APD alternans to the T wave (Phase 9) | Phase 9 METHODS 3 | open |
| Tran DX et al. PRL 102:258103 (2009), PMC2726623; Sato D et al. PNAS (2009) PMC2651322; Biophys J (2010) PMC2913181 | modified Luo-Rudy I EAD family (Phase 9; **dropped**: chaos not reproduced) | Phase 9 METHODS 3.2 | Tran open; Sato full texts unavailable |
| Arnold VI (1965); Glass L, Perez R. PRL 49:1782 (1982); Guevara MR, Glass L. J Math Biol 14:1 (1982) | sine circle map (Phase 9 development-only controls) | Phase 9 METHODS 8; `devmaps.py` | |

## 6. Physiology norms and clinical definitions

| source (as recorded) | used for | where recorded |
|---|---|---|
| Task Force of the European Society of Cardiology and the North American Society of Pacing and Electrophysiology (1996). Heart rate variability: standards of measurement ... | HRV bands and LF/HF norms (realism checks, Phases 5, 8, 9) | `experiments/phase5_rr/README.md`, `experiments/phase6_robust/README.md` |
| Nunan D, Sandercock GRH, Brodie DA (2010). A quantitative systematic review of normal values for short-term heart rate variability | SDNN / RMSSD targets (realism checks) | `experiments/phase5_rr/README.md` |
| Al-Khatib SM et al. 2017 AHA/ACC/HRS guideline (ventricular arrhythmias) | definition of ventricular tachycardia runs (Phase 6 E4) | `experiments/phase6_robust/README.md` |
| Tarvainen MP, Ranta-aho PO, Karjalainen PA (2002). An advanced detrending method ... (smoothness priors) | detrending candidate D3 and `rr_detrend = "smoothness_priors"` | `experiments/phase6_robust/README.md`; `final_pipeline.py` |

## 7. Statistics

| source (as recorded) | used for | where recorded |
|---|---|---|
| Wilson score interval | rate intervals (all phases) | phase analysis modules |
| Clopper-Pearson exact interval | real-rate bound U (Phase 10 P1); Phase 9 Part F | `experiments/phase10_final/analysis10.py` |
| Firth D (1993), bias-reduced (Jeffreys-penalised) logistic regression | P1 detection curves (descriptive) | Phase 10 METHODS 4; `analysis10.py` |
| Liang KY, Zeger SL (1986), GEE, as implemented in statsmodels 0.15.0 | P3 dose-response models | Phase 10 METHODS 5 |
| DeLong ER et al. AUC comparison, fast algorithm of Sun X, Xu W (2014) | secondary AUC intervals (Phase 7) | `experiments/phase7_mitbih/analysis.py` |
| Cohen's kappa; exact McNemar test | Q4 agreement (Phase 10) | `analysis10.py` |

## 8. Software (versions in `requirements-lock.txt`)

numpy, scipy, pandas, scikit-learn, statsmodels, matplotlib, wfdb (WFDB Python package, PhysioNet), numba,
NeuroKit2 0.2.13 (delineator check only), ordpy 1.2.3 (verification only), pytest; TISEAN 3.0.1 and
TimeseriesSurrogates.jl were used as reference implementations for verification (Phase 9), not as dependencies.
