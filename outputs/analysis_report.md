# Analysis report

## Part 2 - Relative frequency summary table

`outputs/part2_cell_frequencies.csv` - 52,500 rows (10,500 samples x 5 populations). Columns: sample, total_count, population, count, percentage.

First rows:

```
     sample  total_count population  count  percentage
sample00000        93214     b_cell  10908   11.702105
sample00000        93214 cd8_t_cell  24440   26.219237
sample00000        93214 cd4_t_cell  20491   21.982749
sample00000        93214    nk_cell  13864   14.873302
sample00000        93214   monocyte  23511   25.222606
sample00001       100824     b_cell   6777    6.721614
sample00001       100824 cd8_t_cell  19407   19.248393
sample00001       100824 cd4_t_cell  33459   33.185551
sample00001       100824    nk_cell  18170   18.021503
sample00001       100824   monocyte  23011   22.822939
```

## Part 3 - Responders vs non-responders

Cohort: melanoma patients treated with miraclib, PBMC samples only, response recorded - **1,968 samples from 656 subjects**.

Primary test: two-sided Mann-Whitney U per population; p-values corrected across the five populations with Benjamini-Hochberg (FDR 5%). Effect size is the rank-biserial correlation.

```
population_label  n_responder  n_non_responder  median_responder  median_non_responder  median_difference  p_mannwhitney  q_value_bh  rank_biserial_r  significant
     CD4+ T cell          993              975           30.2206               29.6576             0.5630         0.0133      0.0667           0.0644        False
          B cell          993              975            9.4323                9.7884            -0.3561         0.0557      0.1393          -0.0498        False
         NK cell          993              975           14.5093               14.7991            -0.2899         0.1211      0.2018          -0.0404        False
        Monocyte          993              975           19.6097               19.9423            -0.3326         0.1632      0.2039          -0.0363        False
     CD8+ T cell          993              975           24.7280               24.6031             0.1250         0.6391      0.6391          -0.0122        False
```

**Conclusion.** No immune population showed a statistically significant difference in relative frequency between responders and non-responders after Benjamini-Hochberg correction (all q >= 0.05).

**Sensitivity analysis.** Each subject contributes up to three samples (days 0/7/14), so the pooled test treats correlated observations as independent. Repeating the analysis on subject-level means (656 independent observations) reproduces the same set of significant populations, so the finding is not an artefact of repeated measures.

```
population_label  n_responder  n_non_responder  median_responder  median_non_responder  median_difference  p_mannwhitney  q_value_bh  rank_biserial_r  significant
     CD4+ T cell          331              325           30.2098               29.8225             0.3873         0.0124      0.0621           0.1128        False
         NK cell          331              325           14.7397               14.9598            -0.2201         0.1267      0.3169          -0.0689        False
          B cell          331              325            9.6714                9.8446            -0.1732         0.3458      0.4322          -0.0425        False
        Monocyte          331              325           19.7945               20.2768            -0.4823         0.2645      0.4322          -0.0504        False
     CD8+ T cell          331              325           24.8969               25.0097            -0.1128         0.6221      0.6221          -0.0222        False
```

Figures: `outputs/figures/part3_boxplot_responders_vs_nonresponders.png`, `part3_boxplot_small_multiples.png`, `part3_boxplot_subject_level.png`.

## Part 4 - Baseline subset analysis

Melanoma PBMC samples at `time_from_treatment_start = 0` from miraclib-treated patients: **656 samples from 656 subjects** (`outputs/part4_baseline_samples.csv`).

### 4.2a - Samples per project

```
breakdown category  samples  subjects
  project     prj1      384       384
  project     prj3      272       272
```

### 4.2b - Subjects by response

```
breakdown category  samples  subjects
 response       no      325       325
 response      yes      331       331
```

### 4.2c - Subjects by sex

```
breakdown category  samples  subjects
      sex        F      312       312
      sex        M      344       344
```

## Closing question

> Considering melanoma males of all sample and treatment types, what is the average number of B cells for responders at time = 0?

**10206.15** B cells, averaged over 485 samples from 485 subjects (condition = melanoma, sex = M, response = yes, time_from_treatment_start = 0; all sample types and treatments included).
