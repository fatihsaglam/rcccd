# rcccd
Class Cover Catch Digraph Classification (CCCD). In CCCD, a ball set, $S$ is determined for each class. Center of these balls are dominant points and radius, $r$ determines the cover area of the ball. Two methods are implemented to determine ball sets: pure and proper cover (PCCCD) and random walk (RWCCCD). In PCCCD, $S$ covers all samples in the class it tries to cover, i.e. target class. This set also never covers any non-target class(es). This methood is vulnerable to the noise and tends to overfit. In RCCCD, dominant points and radii are determined by moving afar from the points step by step, just like walking. Each step is made to the next nearest neighbor and according to the class of neighbor, a score is calculated. There are two methods to calculate scores. One is equal weights for each class, the other is Manukyan and Ceyhan's (2016) step weights to deal with imbalancced classes. Maximum scored sample is the dominant point. This process is repeated until conditions are met. Prediction is made similar to k nearest neighbors but in radii unit. Nearest dominant sample from each class using radii unit is detected. Prediction is made as the closer one. More detail can be found in Marchette and Socolinsky (2003) and Manukyan and Ceyhan (2016).

# R installation
devtools::install_github("https://github.com/fatihsaglam/rcccd")

# References
Marchette, C. E. P. D. J., & Socolinsky, J. G. D. D. A. (2003). Classiﬁcation Using Class Cover Catch Digraphs. Journal of Classiﬁcation, 20, 3-23.

Manukyan, A., & Ceyhan, E. (2016). Classification of imbalanced data with a geometric digraph family. The Journal of Machine Learning Research, 17(1), 6504-6543.
