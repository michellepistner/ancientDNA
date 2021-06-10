# Analysis of Covariation

Welcome! This repo has the code for the analysis of covariation between temporal and spatial effects in ancient DNA. With ancient DNA, much of the metadata is confounded which each other as certain samples are typically only collected from a specific location at a specific time period. In addition, these covariates can be further confounded with other information, as burial location/time might be related to other factors such as class, diet, etc.

We wanted to assess what percentage of signal that is explainable by temporal factors is also *potentially* explainable by geographic location. To do this, we used a multinomial-logistic normal model to model the effect of time on our data set. Letting $Y_{ij}$ denote the observed abundance for taxa $i = 1,...,D$ in sample $j = 1,...N$ and $X$ denote a $Q \times N$ matrix of covariate information, the general form for these models can be written as:


\begin{equation}
  \begin{split}
  Y_j &\sim Multinomial (\pi_j)\\
  \pi_j &= \phi^{-1}(\eta_j)\\
  \eta_j &\sim N(\Lambda X_j, \Sigma)
  \end{split}
\end{equation}

with priors

\begin{equation}
  \begin{split}
  \Lambda &\sim MN_{(D-1) \times P} (\Theta, \Sigma, \Gamma)\\
  \Sigma &\sim IW(\Xi, v),
  \end{split}
\end{equation}

where $\eta_j = ALR(\pi_j)$ (\textit{e.g.}, $\eta_j$ denotes the additive log-ratio coordinates of the proportions which exist on the simplex) and $MN_{(D-1) \times P}$ represents the Matrix-Normal distribution of correct dimension. These methods are efficiently implemented in the *fido* R package (Silverman, 2019) of which we used  version 0.1.13 for our analyses.

**References**

Silverman, J. D., Roche, K., Holmes, Z. C., David, L. A., & Mukherjee, S. (2019). Bayesian multinomial logistic normal models through marginally latent matrix-T processes. *arXiv preprint arXiv:1903.11695*.
