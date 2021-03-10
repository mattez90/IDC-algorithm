Implementation of the ID and IDC algorithm that allow to identify causal quantities in complex graph and also in presence of hidden confounding
In this programs I implement the IDC and ID algorithm. This algorithm was proposed by Pearl and Shpitser (reference papers are 
https://www.aaai.org/Papers/AAAI/2006/AAAI06-191.pdf and https://ftp.cs.ucla.edu/pub/stat_ser/r329-uai.pdf) and allow to identify a causal effect in case we are dealing with: observational datas;  complex graph (since generalize the well known frontdoor and backdoor criterion and do() calculus); hidden confounding. Moreover,
it is possible specify a conditional set z of variable and measures the causal effect on more than one variable putting in place an ideal intervention on more than 
one variable. As far as I know at the momento there is no python implementation for this algorithm.  At the moment the code return the probability structure. I'm working to produce a function that "write in LateX form" the identify expressoion. Another interesting development would be harmonize this output with DoWhy framework in order to have a possibility to estimate the identify causal quantities


