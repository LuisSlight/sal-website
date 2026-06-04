# Chapter 2.5: A Complete Manual Derivation of One-Factor Analysis

> **From \(Q_i = \lambda_iF + \varepsilon_i\) to \(\Sigma = \Lambda\Lambda' + \Psi\)**

At this point, we have established the measurement model (Section 2.1), derived how a latent factor leaves its fingerprint in the covariance structure (Sections 2.2-2.3), and expressed the full result in matrix form (Section 2.4). The purpose of this section is to bring everything together in one place.

We will walk through a complete worked example from start to finish: simulating data with a *known* underlying structure, discarding the hidden mechanism, and then *recovering* it by hand using only the observable correlations. At the end, we verify the result against Stata's built-in estimator.

The goal is not to introduce a new technique. It is to make the whole derivation visible as one coherent object.

After completing this section, you should be able to:

- explain how a latent factor generates observable correlations;
- derive the covariance equation \(\operatorname{Cov}(Q_i,Q_j)=\lambda_i\lambda_j\);
- derive the variance decomposition \(\operatorname{Var}(Q_i)=\lambda_i^2+\psi_i\);
- explain why \(\Sigma=\Lambda\Lambda' + \Psi\) is a compact representation of many covariance equations;
- reproduce a one-factor solution manually; and
- verify the derivation using Stata.

---

## Step 0. The Setup: A Survey of Statistical Anxiety

Suppose we are designing a short questionnaire to measure **statistical anxiety**, an unobservable psychological disposition. We write three items:

| Item | Wording |
|---|---|
| \(q_1\) | "Statistics formulas make me nervous." |
| \(q_2\) | "I worry that I cannot learn statistics." |
| \(q_3\) | "I feel anxious when doing statistics problems." |

We assume that all three items are caused by a single latent factor \(F\), interpreted as anxiety. The measurement model is

$$
Q_i = \lambda_iF + \varepsilon_i,\qquad i=1,2,3 .
$$

The true parameter values, known to us as data generators but unknown to the analyst, are

$$
\lambda_1=0.83,\qquad \lambda_2=0.87,\qquad \lambda_3=0.78,
$$

and

$$
\psi_1=0.31,\qquad \psi_2=0.24,\qquad \psi_3=0.39 .
$$

Here \(\lambda_i\) is the loading for item \(i\), and \(\psi_i=\operatorname{Var}(\varepsilon_i)\) is its unique variance.

## Step 1. Simulate Data from the Known Model

We generate \(N=100\) observations. First, we draw the latent factor

$$
F \sim \mathcal{N}(0,1),
$$

which satisfies the normalization assumption \(\operatorname{Var}(F)=1\). Then we draw independent errors

$$
\varepsilon_i \sim \mathcal{N}(0,\psi_i),
$$

ensuring \(\operatorname{Cov}(F,\varepsilon_i)=0\) and \(\operatorname{Cov}(\varepsilon_i,\varepsilon_j)=0\) for \(i\ne j\).

```stata
clear all
set seed 2026
set obs 100

* Latent factor: Var(F) = 1
gen F = rnormal(0, 1)

* Unique errors: independent of F and of each other
gen epsilon1 = rnormal(0, sqrt(0.31))
gen epsilon2 = rnormal(0, sqrt(0.24))
gen epsilon3 = rnormal(0, sqrt(0.39))

* Observed variables: q_i = lambda_i * F + epsilon_i
gen q1 = 0.83 * F + epsilon1
gen q2 = 0.87 * F + epsilon2
gen q3 = 0.78 * F + epsilon3
```

We can verify the three assumptions directly in the simulated data:

```stata
* Assumption 1: Cov(F, epsilon_i) is approximately 0
correlate F epsilon1 epsilon2 epsilon3, covariance

* Assumption 2: Cov(epsilon_i, epsilon_j) is approximately 0 for i != j
correlate epsilon1 epsilon2 epsilon3, covariance

* Assumption 3: Var(F) is approximately 1
summarize F
```

Now we discard \(F\) and all \(\varepsilon_i\). From this point on, we pretend we are the analyst: we see only \(q_1,q_2,q_3\).

```stata
drop F epsilon1 epsilon2 epsilon3
```

## Step 2. Compute the Sample Covariance Matrix \(S\)

We standardize the three items so that all variances equal 1. Once the variables are standardized, the covariance matrix is the correlation matrix. This is standard practice in exploratory factor analysis when items are measured on comparable but not identical scales.

```stata
foreach var of varlist q1 q2 q3 {
    summarize `var'
    gen z_`var' = (`var' - r(mean)) / r(sd)
}

correlate z_q1 z_q2 z_q3
matrix S = r(C)
```

The resulting sample correlation matrix is

$$
S=
\begin{pmatrix}
1     & 0.677 & 0.485\\
0.677 & 1     & 0.608\\
0.485 & 0.608 & 1
\end{pmatrix}.
$$

This is now **all we have to work with**. The task of factor analysis is to find \(\Lambda\) and \(\Psi\) such that

$$
\Lambda\Lambda' + \Psi \approx S .
$$

The latent factor \(F\) has been discarded. We must recover its structure from the correlations alone.

## Step 3. Estimate the Loadings \(\lambda\) from Off-Diagonal Elements

Recall from Section 2.2 that, under the three assumptions,

$$
\operatorname{Cov}(q_i,q_j)=\lambda_i\lambda_j,\qquad i\ne j .
$$

Because a one-factor model with three standardized variables contains exactly three unknown loadings and three off-diagonal correlations, the system is exactly identified. This allows us to solve the loadings algebraically, without iterative optimization.

The three equations are

$$
r_{12}=\lambda_1\lambda_2,\qquad
r_{13}=\lambda_1\lambda_3,\qquad
r_{23}=\lambda_2\lambda_3 .
$$

Divide the first equation by the third:

$$
\frac{r_{12}}{r_{23}}
=
\frac{\lambda_1\lambda_2}{\lambda_2\lambda_3}
=
\frac{\lambda_1}{\lambda_3}.
$$

Similarly,

$$
\frac{r_{12}}{r_{13}}=\frac{\lambda_2}{\lambda_3}.
$$

Substitute this relation back into \(r_{23}=\lambda_2\lambda_3\):

$$
r_{23}
=
\left(\frac{r_{12}}{r_{13}}\lambda_3\right)\lambda_3
=
\frac{r_{12}}{r_{13}}\lambda_3^2 .
$$

Solving for the third loading gives

$$
\lambda_3
=
\sqrt{\frac{r_{23}r_{13}}{r_{12}}}.
$$

The remaining two loadings follow immediately:

$$
\lambda_1=\frac{r_{13}}{\lambda_3},
\qquad
\lambda_2=\frac{r_{23}}{\lambda_3}.
$$

Equivalently, using the ratios above,

$$
\lambda_1=\left(\frac{r_{12}}{r_{23}}\right)\lambda_3,
\qquad
\lambda_2=\left(\frac{r_{12}}{r_{13}}\right)\lambda_3 .
$$

In Stata:

```stata
local r12 = S[1,2]    // 0.6766
local r13 = S[1,3]    // 0.4848
local r23 = S[2,3]    // 0.6080

local ratio_1 = `r12' / `r23'    // lambda_1 / lambda_3
local ratio_2 = `r12' / `r13'    // lambda_2 / lambda_3

local lambda3 = sqrt(`r23' / `ratio_2')
local lambda1 = `ratio_1' * `lambda3'
local lambda2 = `ratio_2' * `lambda3'
```

The estimated loadings are:

| Variable | Manual \(\hat{\lambda}\) | True \(\lambda\) |
|---|---:|---:|
| \(q_1\) | 0.7345 | 0.83 |
| \(q_2\) | 0.9212 | 0.87 |
| \(q_3\) | 0.6601 | 0.78 |

The estimates deviate somewhat from the true values. This is expected: with \(N=100\) and a purely algebraic method, sampling variability matters. But the ordering and relative magnitudes are preserved correctly.

## Step 4. Estimate the Uniqueness \(\psi\) from Diagonal Elements

Recall from Section 2.3 that

$$
\operatorname{Var}(q_i)=\lambda_i^2+\psi_i .
$$

Since the variables are standardized, \(\operatorname{Var}(q_i)=1\) for every item. Therefore

$$
\psi_i=1-\lambda_i^2 .
$$

The quantity \(\lambda_i^2\) is called the **communality** of item \(i\): the proportion of the item's variance explained by the latent factor.

```stata
local psi1 = 1 - `lambda1'^2
local psi2 = 1 - `lambda2'^2
local psi3 = 1 - `lambda3'^2
```

| Variable | \(\hat{\lambda}\) | Communality \(\hat{\lambda}^2\) | Uniqueness \(\hat{\psi}\) |
|---|---:|---:|---:|
| \(q_1\) | 0.7345 | 0.5394 | 0.4606 |
| \(q_2\) | 0.9212 | 0.8485 | 0.1515 |
| \(q_3\) | 0.6601 | 0.4357 | 0.5643 |

**Interpretation.** Item \(q_2\), "I worry that I cannot learn statistics," has the highest communality, approximately 0.85. Nearly all of its variance is accounted for by the latent anxiety factor, making it the purest indicator. Item \(q_3\) has the most unique variance, which means it appears to capture something that the common factor alone does not fully explain.

## Step 5. Reconstruct \(\hat{\Sigma}=\Lambda\Lambda' + \Psi\)

We now verify that our estimated parameters reproduce the original correlation matrix. Construct the loading vector and uniqueness matrix:

$$
\hat{\Lambda}
=
\begin{pmatrix}
0.7345\\
0.9212\\
0.6601
\end{pmatrix},
\qquad
\hat{\Psi}
=
\begin{pmatrix}
0.4606 & 0      & 0\\
0      & 0.1515 & 0\\
0      & 0      & 0.5643
\end{pmatrix}.
$$

In Stata:

```stata
matrix LAMBDA = [`lambda1' \ `lambda2' \ `lambda3']
matrix PSI    = diag((`psi1', `psi2', `psi3'))

matrix SIGMA_HAT = LAMBDA * LAMBDA' + PSI
matrix RESIDUAL  = S - SIGMA_HAT

matrix list RESIDUAL
```

The residual matrix is

$$
S-\hat{\Sigma}
\approx
\begin{pmatrix}
0 & -4.4\times10^{-16} & -1.7\times10^{-16}\\
\cdot & 0 & -1.1\times10^{-16}\\
\cdot & \cdot & 0
\end{pmatrix}.
$$

The off-diagonal residuals are numerically zero; the small values are only floating-point rounding. This confirms that the algebraic solution exactly reproduces \(S\). In real data, non-zero residuals indicate sampling variability, model misspecification, or too few factors.

## Step 6. Verify Against Stata's Built-In Estimator

Stata's `factor` command uses iterative principal factors (IPF), a more sophisticated algorithm than our algebraic solution. We compare:

```stata
factor z_q1 z_q2 z_q3, ipf factors(1)
```

Stata's output is:

```text
Factor loadings (pattern matrix) and unique variances

    Variable |  Factor1 |   Uniqueness
    ---------+----------+--------------
        z_q1 |   0.7346 |      0.4604
        z_q2 |   0.9210 |      0.1518
        z_q3 |   0.6601 |      0.5642
```

Comparing Stata to the manual estimates:

| Parameter | Manual | Stata IPF |
|---|---:|---:|
| \(\hat{\lambda}_1\) | 0.7345 | 0.7346 |
| \(\hat{\lambda}_2\) | 0.9212 | 0.9210 |
| \(\hat{\lambda}_3\) | 0.6601 | 0.6601 |
| \(\hat{\psi}_1\) | 0.4606 | 0.4604 |
| \(\hat{\psi}_2\) | 0.1515 | 0.1518 |
| \(\hat{\psi}_3\) | 0.5643 | 0.5642 |

The results are virtually identical. The small differences arise because Stata's IPF algorithm iterates to a refined solution, while our algebraic approach solves the exactly identified system directly. The underlying logic is the same: find the loading matrix \(\Lambda\) and uniqueness matrix \(\Psi\) that reproduce the observed covariance structure.

## Step 7. Computing Factor Scores

Once loadings are estimated, we can compute a **factor score** for each respondent: a weighted composite of responses that serves as the best reconstruction of the respondent's underlying anxiety level.

```stata
predict factor_score
label variable factor_score "Statistical Anxiety (Factor Score)"

summarize factor_score
```

The scoring coefficients printed by Stata indicate the weight assigned to each item:

| Item | Scoring weight |
|---|---:|
| \(z_{q_1}\) | 0.187 |
| \(z_{q_2}\) | 0.711 |
| \(z_{q_3}\) | 0.137 |

Notice that \(q_2\), the item with the highest loading, receives the greatest weight. This is not arbitrary. Items that are more strongly connected to the latent factor carry more information about it, and the regression-based scoring method formalizes this intuition.

Factor scores are standardized: \(\bar{F}\approx0\) and \(\operatorname{SD}(F)\approx1\). A score of \(+1\) indicates a respondent one standard deviation above average in statistical anxiety; \(-1\) indicates one standard deviation below average.

## Summary: The Complete Logic Chain

The derivation in this section follows a single continuous thread:

$$
\underbrace{Q_i=\lambda_iF+\varepsilon_i}_{\text{measurement model}}
\xrightarrow{\text{assumptions}}
\underbrace{\operatorname{Cov}(q_i,q_j)=\lambda_i\lambda_j}_{\text{off-diagonal fingerprint}}
\xrightarrow{\text{solve}}
\underbrace{\hat{\lambda}_i}_{\text{estimated loadings}} .
$$

The diagonal elements provide the second part:

$$
\underbrace{\operatorname{Var}(q_i)=\lambda_i^2+\psi_i}_{\text{variance decomposition}}
\xrightarrow{\text{solve}}
\underbrace{\hat{\psi}_i=1-\hat{\lambda}_i^2}_{\text{estimated uniqueness}}
\xrightarrow{\text{combine}}
\underbrace{\hat{\Sigma}=\hat{\Lambda}\hat{\Lambda}' + \hat{\Psi}\approx S}_{\text{model reproduction}} .
$$

The matrix equation

$$
\Sigma=\Lambda\Lambda' + \Psi
$$

is not a new model. It is simply a compact notation for the covariance equations derived from the measurement model. Every element in that matrix corresponds to a statement about how observable variables co-vary because they share a common cause.

> **Core insight.** The latent factor \(F\) is unobservable, but it is not undetectable. Its presence leaves a systematic fingerprint in the off-diagonal elements of the correlation matrix. Factor analysis is the method for reading that fingerprint: recovering the hidden structure that generated the correlations we observe.
