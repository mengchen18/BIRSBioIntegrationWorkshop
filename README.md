# BIRSBioIntegrationWorkshop

This repository stores the code to solve the questions listed in:
https://github.com/BIRSBiointegration/Hackathon/tree/master/sc-targeted-proteomics


*	What integrative data analysis question have you addressed with the selected data and why?

"How should we approach integrating partially overlapping proteomic data collected on different patients with similar phenotypes?"

I believe the multiple-block PCA, the methods I have implemented, can be extended to solve these questions. I am curious to see 
how well its performance is compared to some other methods. 

*	What are the advantages and performance of your approach?

I have mainly worked on two methods: simple linear regression (as a baseline/control) and 
multi-block PCA (MBPCA; including multiple co-inertia, multiple canonical correspondence 
analysis as special cases). In theory, MBPCA should outperform simple linear regression 
because it finds the correlated pattern across multiple datasets, preventing the potential 
problem of overfitting to one dataset. But I find this is not always the case (also see next 
point). I apply the linear regression and MBPCA to two scenarios. (1) multiple datasets
representing the different gene expression profiles over the same set of samples. The 
question I asked is if one dataset has some missed samples, could we use 
the information from other datasets to predict the gene expression profile 
of the missed samples. (2) is closer to the question of we are asking here, could 
we predict protein expression of single cells based on partially overlapping proteomic 
data. In fact, the only difference between the two questions depends on how the 
matrices are transposed. However, I find that MBPCA only outperforms linear regression 
in scenario 1 (figure below), whereas the performance of the two methods is comparable in scenario 2. 
The main reason is discussed in the next point.


<img src="Fig/benchmark_NCI60_train40_test20.png" alt="benchmark_NCI60_train40_test20" width="400" />

Other remarks 
1) two other categories of methods, multiple-factorial analysis and partial least square, 
could also be applied to solve the question, but I haven’t got the time to include all 
them in the evaluation, and they share the same problems as MBPCA which I will discuss 
in the next point. 
2) MBPCA is also 10-100 times faster than linear regression methods. 


*	What were the specific challenges you have encountered ?

The MBPCA was originally developed as a method to integrate multiple omics datasets where 
different variables (features) measured on the same set of samples (similar to scenarios 1 above). 
In practice, these data are often normalized in a way that the variable distributions are 
comparable between samples (column). As a result, in MBPCA, we can directly center/scale each 
variable so that the divergence among samples is highlighted. However, scenario 2 treats 
proteins as samples and individual cells as variables. However, the overall abundances 
of different proteins are different in orders of magnitude, leading to highly correlated 
variables after scaling/centering. The accuracy of prediction is negatively influenced by 
the collinearity of scale variables. 


*	How are you going to address those challenges?

I have tried to scale the protein expression and it improves the prediction accuracy. 
I believe there are more spaces to improve by applying different preprocessing or different 
metrics to transform the data. A starting point of the next step would try to transform the 
values to marginal probabilities as in correspondence analysis. Of course, I would also happy 
to discuss this with others if I would have a chance to join the workshop.

