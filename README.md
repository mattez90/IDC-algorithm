In this programs I implement the IDC and ID algorithms. This algorithms were proposed by Pearl and Shpitser (reference papers are 
https://www.aaai.org/Papers/AAAI/2006/AAAI06-191.pdf and https://ftp.cs.ucla.edu/pub/stat_ser/r329-uai.pdf) and allow to identify a causal effect in case we are dealing with: observational datas;  complex graph (since generalize the well known frontdoor and backdoor criterion and do() calculus); hidden confounding. Moreover,
it is possible specify a conditional set z of variable and measures the causal effect on more than one variable putting in place an ideal intervention on more than 
one variable. As far as I know at there are no python implementation for this algorithm (there is an R package, causaleffect, that implement this method of identification and from which I got inspired for the output structure)  At the moment the code return the probability structure. I'm working to produce a function that "write in LateX form" the identified expression. Another interesting development would be harmonize this output with DoWhy framework in order to have a possibility to estimate the identify causal quantities


