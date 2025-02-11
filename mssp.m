function [Sigma, Theta] = mssp(x, P, alpha, beta, kappa, n)
%MSSP Merwe Scaled Sigma Points
%   [SIGMA, THETA] = MSSP() returns the Van der Merwe scaled sigma points SIGMA 
%   and weights THETA for zero mean, identity covariance, spread alpha=1e-3, 
%   prior knowledge beta=2, secondary scaling kappa=0, and number of pts n=3. 
%
%   [SIGMA, THETA] = MSSP(X) returns the Van der Merwe scaled sigma points 
%   SIGMA and weights THETA for mean X, identity covariance, spread alpha=1e-3, 
%   prior knowledge beta=2, and secondary scaling kappa=0. X is a D-by-1 
%   vector, and number of pts n=2*D+1. 
%
%   [SIGMA, THETA] = MSSP(X,P) returns the Van der Merwe scaled sigma 
%   points SIGMA and weights THETA for mean X, covariance P, spread alpha=1e-3, 
%   prior knowledge beta=2, and secondary scaling kappa=0. P is a
%   D-by-D matrix, and number of pts n=2*D+1. 
%
%   [SIGMA, THETA] = MSSP(X,P,ALPHA) returns the Van der Merwe scaled 
%   sigma points SIGMA and weights THETA for mean X, covariance P, spread
%   ALPHA, prior knowledge beta=2, secondary scaling kappa=0, and number of 
%   pts n=2*D+1. 
%
%   [SIGMA, THETA] = MSSP(X,P,ALPHA,BETA) returns the Van der Merwe scaled 
%   sigma points SIGMA and weights THETA for mean X, covariance P, spread
%   ALPHA, prior knowledge BETA, secondary scaling kappa=0, and number of 
%   pts n=2*D+1. 
%
%   [SIGMA, THETA] = MSSP(X,P,ALPHA,BETA,KAPPA) returns the Van der Merwe 
%   scaled sigma points SIGMA and weights THETA for mean X, covariance P, 
%   spread ALPHA, prior knowledge BETA, secondary scaling KAPPA, and number of 
%   pts n=2*D+1. 
%
%   [SIGMA, THETA] = MSSP(X,P,ALPHA,BETA,KAPPA,N) returns the Van der Merwe 
%   scaled sigma points SIGMA and weights THETA for mean X, covariance P, 
%   spread ALPHA, prior knowledge BETA, secondary scaling KAPPA, and number of 
%   pts N=2*N+1. 
%
%   SIGMA is an N-by-D matrix and THETA is a N-by-2 matrix where THETA(:,1) 
%   are the mean weights and THETA(:,2) are the covariance weights.
%
%   Example:
%
%     x = [-1; 0.5]; P = [2.5 5; 5 11]; 
%     figure; hold on; legend; 
%     npts=100; 
%     tt = linspace(0, 2*pi, npts);
%     x1 = cos(tt); x2=sin(tt);
%     ap = [x1(:) x2(:)]';
%     for i = 1:3
%         [v,d]=eig(P); 
%         d = i * sqrt(d); 
%         bp = (v*d*ap) + repmat(x, 1, size(ap,2)); 
%         plot(bp(1,:), bp(2,:), "k-", 'LineWidth',2, "HandleVisibility", "off");  % [1,2,3] std. dev. ellipse
%     end
%     [Sigma, Theta] = mssp(x, P, 0.8, 2, 0); % sigma points
%     scatter(Sigma(:,1), Sigma(:,2), 100, "k", "filled", "Marker", "d", "DisplayName", "Merwe Scaled Sigma Points");
%
%   References:
%      R. Van der Merwe "P-Point Kalman Filters for Probabilitic 
%      Inference in Dynamic State-Space Models" (Doctoral dissertation)

if nargin < 1
    x    = 0; 
    P = 1; 
    alpha = 1e-3;
    beta  = 2;
    kappa = 0;
    n     = 3; 
else
    [d, ~] = size(x); 
end

if nargin<2
    P = eye(d); 
    alpha = 1e-3;
    beta  = 2;
    kappa = 0;
    n     = 2*d+1; 
elseif nargin<3
    [d1_P, d2_P] = size(P);
    if (d1_P ~= d)||(d2_P ~= d)
        error("BadPDimension");
    end
    alpha = 1e-3;
    beta  = 2;
    kappa = 0;
    n     = 2*d+1; 
elseif nargin<4
    [d1_P, d2_P] = size(P);
    if (d1_P ~= d)||(d2_P ~= d)
        error("BadPDimension");
    end
    if ~isscalar(alpha)
        error("BadAlphaDimension");
    else
        if (alpha < 0)||(alpha > 1)
            error("BadAlpha");
        end
    end
    beta  = 2;
    kappa = 0;
    n     = 2*d+1; 
elseif nargin<5
    [d1_P, d2_P] = size(P);
    if (d1_P ~= d)||(d2_P ~= d)
        error("BadPDimension");
    end
    if ~isscalar(alpha)
        error("BadAlphaDimension");
    else
        if (alpha < 0)||(alpha > 1)
            error("BadAlpha");
        end
    end
    if ~isscalar(beta)
        error("BadBetaDimension");
    else
        if (beta < 0)
            error("BadBeta");
        end
    end
    kappa = 0;
    n     = 2*d+1; 
elseif nargin<6
    [d1_P, d2_P] = size(P);
    if (d1_P ~= d)||(d2_P ~= d)
        error("BadPDimension");
    end
    if ~isscalar(alpha)
        error("BadAlphaDimension");
    else
        if (alpha < 0)||(alpha > 1)
            error("BadAlpha");
        end
    end
    if ~isscalar(beta)
        error("BadBetaDimension");
    else
        if (beta < 0)
            error("BadBeta");
        end
    end
    if ~isscalar(kappa)
        error("BadKappaDimension");
    else
        if (kappa < 0)
            error("BadKappa");
        end
    end
    n     = 2*d+1; 
elseif nargin<7
    [d1_P, d2_P] = size(P);
    if (d1_P ~= d)||(d2_P ~= d)
        error("BadPDimension");
    end
    if ~isscalar(alpha)
        error("BadAlphaDimension");
    else
        if (alpha < 0)||(alpha > 1)
            error("BadAlpha");
        end
    end
    if ~isscalar(beta)
        error("BadBetaDimension");
    else
        if (beta < 0)
            error("BadBeta");
        end
    end
    if ~isscalar(kappa)
        error("BadKappaDimension");
    else
        if (kappa < 0)
            error("BadKappa");
        end
    end
    if ~isscalar(n)
        error("BadNDimension");
    else
        if (n < 0)
            error("BadN");
        else
            n = 2*n + 1;
        end
    end
else
    error("BadNumberOfInputs");
end

lambda = alpha^2 * (d + kappa) - d; 
c = .5 / (d + lambda); 
Theta = c .* ones(n, 2);
Theta(1,1) = lambda / (d + lambda); 
Theta(1,2) = lambda / (d + lambda) + (1 - alpha^2 + beta); 

U = chol((lambda + d) .* P);
Sigma = NaN(n,d);
Sigma(1,:) = x';
for i = 1:d
    Sigma(i+1,:)   = x' + U(i,:); 
    Sigma(d+i+1,:) = x' - U(i,:); 
end

end