# Merwe Scaled Sigma Points (MSSP)
The MSSP repository provides MATLAB code for calculating the Van der Merwe scaled sigma points and weights described by a mean $x \in \mathbb{R}^{d\times 1}$, covariance $P \in \mathbb{R}^{d\times d}$, spread parameter $\alpha$, prior knowledge parameter $\beta$, and secondary scaling $\kappa$, as defined in [1]. The sigma points { $\Sigma_i$ } and weights { $\theta_i^m, \theta_i^c$ } are calculated as follows:<br>

$$
\begin{align}
   \Sigma\_0 &= x \\
   \Sigma\_i &= x - {(\sqrt{(\lambda + d)P})}\_i; \quad i=1,\dots,d\\
   \Sigma\_{i+d} &= x + {(\sqrt{(\lambda + d)P})}\_i; \quad i=1,\dots,d\\
   \theta_0^m &= \frac{\lambda}{d+\lambda}, \quad \theta_0^c = \theta_0^m + (1-\alpha^2 + \beta); \quad i=0\\
   \theta_i^m &= \theta_i^c = \frac{1}{2(d+\lambda)}; \quad i=1,2,\dots,2d,
\end{align}
$$

where $\lambda = \alpha^2(d + \kappa) - d$, and $0 \leq \alpha \leq 1$, $0 \leq \beta$ and $\kappa$ are tuning parameters. These $2d+1$ can be used to exactly construct the mean $x$ and covariance $P$ in the following way: <br>

$$
\begin{align}
    x &= \sum_{i=0}^{2d} \theta_i^m \Sigma_i \\
    P &= \sum_{i=0}^{2d} \theta_i^c(\Sigma_i - x)(\Sigma_i - x)^T
\end{align}
$$

Please direct any questions to blhanson@ucsd.edu. <br><br>

## References
[1] R. Van der Merwe "P-Point Kalman Filters for Probabilitic Inference in Dynamic State-Space Models" (Doctoral dissertation)

![mssp_test](https://github.com/user-attachments/assets/74121144-be06-41b5-8d18-e676bba692f3)
