# Optimized-Bayesian-Nonlocal-means-with-block(OBNLM)
Optimized bayesian nonlocal-means algorithm for denoising ultrasound image

# Description
- BayesianNLM.m : the OBNLM algorithm
- getPearsonDistance.m: compute the pearson distance based bayesian framework
- ImgNormalize.m: pre-processing the image (histogram stretching) 
- testBayesianNLM.m: test our algorithm using the input (noisyImage.png) and get the output (despeckledImage.png)

# Reference Paper
Coupé, Pierrick, et al. "Nonlocal means-based speckle filtering for ultrasound images." IEEE transactions on image processing 18.10 (2009): 2221-2229.

# Basic Principle
The blockwise Nonlocal Means algorithm is completed. The basic principle is shown below:
![basic principle](https://github.com/Xingorno/Figures/blob/master/NLM.png?raw=true)

**Object Function**: <img src="http://latex.codecogs.com/svg.latex?NL(u)(B_j) = \sum_{i\in\Delta_j}w(B_i,B_j)u(Bi)" border="0"/> 
with <img src="http://latex.codecogs.com/svg.latex? w(B_i,B_j)=\frac{1}{Z^j}e^{-\frac{dp(u(B_i),u(B_j))}{h^2}}" border="0"/>

<img src="http://latex.codecogs.com/svg.latex? B" border="0"/>: the blockwise of image

<img src="http://latex.codecogs.com/svg.latex? dp(u(B_i),u(B_j))" border="0"/>:the pearson distance between blockwise



