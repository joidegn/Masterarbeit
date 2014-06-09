# Here is how Dynamic Factor models with principal components work:

X are rotated such that they are uncorrelated to estimate an equation of the form
X = Lambda\*F
using PCA. The resulting rotation matrix is the matrix that rotates X into the space spanned by F. In other words the inverse of the rotation matrix is what would give us back X if multiplied with the factors.
Lambda = (rotation matrix)^-1  The rotation matrix consists of the Eigenvectors of X'X (i.e. it is orthogonal) or XX' scaled by T or N respectively (say Bai and Ng 2002)

Summary:
Singular Value decomposition of X yields:
X = U\*S\*Vt
Vt are the Eigenvectors: V (untransposed) gives the rotation matrix (also principal component directions and hopefully factor loadings)
U\*Diag(S) = X\*V (because V^-1 = V') which gives the rotated X or the Scores (principal components)

## Dynamic PCA
Factors might change over time --> make barplots of the eigenvectors -> eyeball if there are important variables (==scores/factors)

Find out how VAR models figure in that: Giannone, Reichlin and Sala (2002) or Forni, Lippi and Reichlin (2003) 


How are dynamic Factor models different from PCR?

$F'F = I_N$ I think due to the normalization performed by the julia function pca (see also BAi Ng, 2002 p.198)

Why is BIC inconsistent as a criterion to estimate the number of factors?

An alternative to principal components to estimate the factors is e.g. to write the rewrite model in a state-space setting and use a Kalman-Filter to maximize a gaussian likelihood. (Bai 2003)


If one is interested in the original factors then the static factors are inappropriate. For forecasting we can simply include more static factors (more static factors account for lags of the original factors). (Breitung, Eickmeier 2011 p72)
