# T-Zone Enrichment in a Stochastic Model of Cell Death

This repository investigates when cell deaths occur more often in a predefined T-shaped region than its share of the observation area would predict. The simulations use an event-driven Gillespie algorithm with thinning.

Four factors vary across the experiments.

1. The contrast in activation intensity between the T-zone and the surrounding region.
2. The spatial extent of the active regions.
3. The intensity of cell death.
4. The radius and duration of ERK-mediated protection after a death.

## 1. Burn-in and experimental design

Each completed replicate starts from an empty system and uses an independent random seed. The first **500 observed deaths** serve as burn-in to reduce the influence of the initial state. The **next 1,000 observed deaths** enter the spatial analysis. Each completed replicate therefore contains 1,500 accepted death events.

## 2. Testing for enrichment in the T-zone

Let $T$ denote the T-zone and $W$ the observation window. Of the $N$ deaths included in the analysis, let $K$ be the number that fall in $T$. Write $\pi_T$ for the probability that a death falls in $T$.

The one-sided exact binomial test uses the area fraction $p_0$ as its reference probability. Under the binomial null model, each death independently falls in the T-zone with probability

$$
p_0 = \frac{\lvert T\rvert}{\lvert W\rvert},
\qquad K \sim \mathrm{Binomial}(N,p_0).
$$

The hypotheses are

$$
H_0:\ \pi_T = p_0
\qquad \text{versus} \qquad
H_1:\ \pi_T > p_0.
$$

For an observed count $K$, the one-sided $p$-value is $\Pr(X\ge K)$ for $X\sim\mathrm{Binomial}(N,p_0)$. The analyses use a nominal significance level of $\alpha=0.05$.

## 3. T-zone density ratio

The density ratio compares deaths per unit area inside and outside the T-zone.

$$
R_T =
\frac{K/\lvert T\rvert}
     {(N-K)/(\lvert W\rvert-\lvert T\rvert)}.
$$

- $R_T>1$ means deaths are more concentrated inside the T-zone.
- $R_T=1$ means the two regions have the same death density.
- $R_T<1$ means deaths are less concentrated inside the T-zone.

## 4. Simulation experiments

### Step 1. Activation contrast and localization

`step1code.py` explores how T-zone enrichment responds to four activation parameters.

$$
\lambda_{a,T},\qquad \lambda_{a,c},\qquad \beta_{a,R},\qquad \beta_{a,T}.
$$

The sweep varies both the absolute activation intensities and their ratio $\lambda_{a,T}/\lambda_{a,c}$, as well as the radii and lifetimes of active regions. Each parameter setting is simulated in multiple independent replicates.

The reference regime showed no significant enrichment. Increasing the activation contrast alone was insufficient in the settings tested. Enrichment emerged when activation was sufficiently localized in space and time.

| File | Purpose |
| --- | --- |
| `step1code.py` | Runs the Gillespie simulations and sweeps the activation and T-zone parameters |
| `step1result.csv` | Simulation results for individual replicates |
| `step1resultsummary.csv` | Descriptive statistics by parameter setting |
| `step1pooled_binomial_analysis.py` | Pools independent replicates and runs the one-sided exact binomial test |
| `step1pooled_binomial_summary.csv` | Pooled T-zone enrichment results |

### Step 2. Factorial sweep of death and ERK parameters

`step2code.py` fixes activation at an intermediate level of localization that showed enrichment in Step 1.

$$
\lambda_{a,1}=5.0,\quad \lambda_{a,T}=0.5,\quad \lambda_{a,c}=0.005,\quad
\beta_{a,R}=2.5,\quad \beta_{a,T}=3.0.
$$

The experiment crosses three death intensities with three ERK radius rates and three ERK lifetime rates.

$$
\lambda_d\in\{0.5,1,2\},\qquad
\beta_{d,R}\in\{1,2,4\},\qquad
\beta_{d,T}\in\{0.4,0.8,1.6\}.
$$

This gives $3\times3\times3=27$ settings, with five independent replicates per setting. T-zone enrichment persisted across all 27 settings. Varying the death and ERK parameters mainly changed the time scale of observed deaths and, to a lesser extent, the magnitude of enrichment.

| File | Purpose |
| --- | --- |
| `step2code.py` | Runs the full factorial sweep of death and ERK feedback parameters |
| `step2result.csv` | Simulation results for individual replicates |
| `step2resultsummary.csv` | Summary statistics by death and ERK parameter setting |
| `step2binomialtestcode.py` | Runs the one-sided exact binomial test on pooled counts |
| `step2pooled_binomial_summary.csv` | Pooled binomial test results |

### Step 3. Isolating the effect of activation localization

`step3code.py` compares three activation regimes while keeping the activation intensities fixed.

$$
\lambda_{a,1}=5.0,\qquad \lambda_{a,T}=0.5,\qquad \lambda_{a,c}=0.005.
$$

| Activation regime | $\beta_{a,R}$ | $\beta_{a,T}$ | Interpretation |
| --- | ---: | ---: | --- |
| `non_local_activation` | 1.0 | 1.2 | Large, long-lived active disks |
| `medium_activation` | 2.5 | 3.0 | Smaller, shorter-lived disks than in the nonlocal regime |
| `local_activation` | 5.0 | 3.0 | The smallest disks, with the same mean lifetime as in the intermediate regime |

For each activation regime, the experiment crosses five death intensities, $\lambda_d\in\{0.25,0.5,1,2,4\}$, with three levels of ERK protection (strong, reference, and weak). This gives $5\times3=15$ settings per regime. Each setting pools **10 independent replicates**, for a total of **10,000 analyzed deaths** after burn-in.

| File | Purpose |
| --- | --- |
| `step3code.py` | Runs the experiment across activation, death, and ERK regimes |
| `step3result.csv` | Pooled binomial test results for each parameter setting |
| `step3resultsummary.csv` | Summary of enrichment across the 15 settings for each activation regime |

## 5. Main finding

The T-zone occupies $p_0=4/9\approx0.444$ of the observation window. With activation intensities held fixed, the mean fraction of deaths inside it ranged from 0.424 under nonlocal activation to 0.907 under local activation.

| Activation regime | Settings with significant enrichment | Mean fraction of deaths in the T-zone | Mean density ratio |
| --- | ---: | ---: | ---: |
| Local | 15 / 15 | 0.907 | 12.21 |
| Intermediate | 15 / 15 | 0.622 | 2.07 |
| Nonlocal | 0 / 15 | 0.424 | 0.92 |

Enrichment persisted across every death and ERK setting tested under local and intermediate activation. The nonlocal regime showed no significant enrichment in any setting. Across these experiments, the spatial extent of activation was the main determinant of the T-shaped pattern. Death intensity and ERK feedback primarily affected the time scale and the strength of the observed signal.

## 6. Model parameters

| Parameter | Meaning |
| --- | --- |
| `lambda_a_1` | Activation intensity inside active regions, also used as the upper bound for thinning |
| `lambda_a_T` | Baseline activation intensity inside the T-zone |
| `lambda_a_c` | Baseline activation intensity outside the T-zone |
| `lambda_d` | Intensity of candidate death events |
| `beta_a_R` | Exponential rate for the radius of an active disk |
| `beta_a_T` | Exponential rate for the lifetime of an active region |
| `beta_d_R` | Exponential rate for the radius of an ERK protection zone |
| `beta_d_T` | Exponential rate for the lifetime of an ERK protection zone |

Under the exponential model, the mean radii and lifetimes are

$$
\mathbb E[R^a]=\frac{1}{\beta_{a,R}},\qquad
\mathbb E[T^a]=\frac{1}{\beta_{a,T}},\qquad
\mathbb E[R^d]=\frac{1}{\beta_{d,R}},\qquad
\mathbb E[T^d]=\frac{1}{\beta_{d,T}}.
$$

A higher radius rate gives a smaller mean radius. A higher lifetime rate gives a shorter mean lifetime.

## 7. Software requirements

The code requires Python 3. Install the dependencies with

```bash
pip install numpy pandas scipy matplotlib
```
