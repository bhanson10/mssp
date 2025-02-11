# Merwe Scaled Sigma Points (MSSP)
The MSSP repository provides MATLAB code for calculating the returns the Van der Merwe scaled sigma points and weights described by a mean $x \in \mathbb{R}^{d\times 1}$, covariance $P \in \mathbb{R}^{d\times d}$, spread parameter $\alpha$, prior knowledge parameter $\beta$, and secondary scaling $\kappa$, as defined in [1]. The sigma points $\{\Sigma}$ <br>

$$
\begin{gather}
    p(\mathbf{x}|\boldsymbol{\mu},\Sigma,\beta) =  A(\beta,d)
    \exp(-[B(\beta,d)(\mathbf{x}-\boldsymbol{\mu})^T \Sigma^{-1}(\mathbf{x}-\boldsymbol{\mu})]^{\beta}),
    \\ 
    \text{where}\quad A(\beta,d)= \Big(\frac{B(\beta,d)}{\pi}\Big)^{\frac{d}{2}}\cdot\frac{\Gamma(\frac{d}{2})\beta}{\Gamma(\frac{d}{2\beta})|\Sigma|^{\frac{1}{2}}}\quad\text{and}\quad B(\beta,d)= \frac{\Gamma(\frac{d+2}{2\beta})}{d\Gamma(\frac{d}{2\beta})}, 
\end{gather}
$$

When $\beta = 1$, *mvggd.m* becomes *mvnpdf.m*. Please direct any questions to blhanson@ucsd.edu. <br><br>

## References
[1] R. Van der Merwe "P-Point Kalman Filters for Probabilitic Inference in Dynamic State-Space Models" (Doctoral dissertation)

![mssp_test](https://github.com/user-attachments/assets/74121144-be06-41b5-8d18-e676bba692f3)
