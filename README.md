# AUMIC: MATLAB Toolbox for the Cybernetic Modeling and More
Hyun-Seob Song<sup>1</sup> and Doraiswami Ramkrishna<sup>2</sup> 

<sup>1</sup> Biological Sciences Division, Pacific Northwest National Laboratory, Richland, WA 99352  
<sup>2</sup> Davidson School of Chemical Engineering, Purdue University, West Lafayette, IN 47907

We introduce a new tool for dynamic metabolic modeling – AUMIC (AUtomated tool for Metabolic modeling Integrated with the Cybernetic regulatory mechanisms)! 

AUMIC is a MATLAB toolbox that enables automatic construction of multiple types of cybernetic models (Ramkrishna and Song, 2012;Ramkrishna and Song, 2018), including: 
* Lumped cybernetic model (LCM) (e.g., Kompala et al., 1986)
* Hybrid cybernetic model (HCM) (Kim et al., 2008; Song et al., 2009)
* Lumped hybrid cybernetic model (L-HCM) (Song and Ramkrishna, 2010; 2011)

The AUMIC can also implement other dynamic modeling approaches for comparative studies, e.g., lumped kinetic model (e.g., Bijkerk and Hall, 1977), and macroscopic bioreaction model (Provost and Bastin, 2004). 

In AUMIC, all of the foregoing approaches are subsumed as quasi-steady-state models which are formulated using a unifying framework through the following step-by-step procedures, i.e., 
*	Network decomposition into elementary modes (EMs) using efmtool (Terzer and Stelling, 2008) or metatool (von Kamp and Schuster, 2006) 
*	EM classification and processing
*	Design of kinetics for uptake fluxes (and kinetics for enzyme synthesis)

Through user-defined functions, one can flexibly account for any functional forms of kinetic equations (without being limited to Monod-type equations), and various types of reactor configurations (including batch, fed-batch, and continuous reactors). Moreover, the AUMIC provides the parameter identification routine so that one can develop their own models using specific experimental data. Finally, the results are quickly analyzed using the post-processing module.

The upcoming version of AUMIC will include dynamic predictions of cellular response to gene knock-outs, a key feature for the rational design of metabolic pathways in microorganisms. We expect to provide updates in the future to make AUMIC as effective and versatile as possible.

Please contact Dr. Song (hyunsbsong@gmail.com) or Dr. Ramkrishna (ramkrish@ecn.purdue.edu) for any questions regarding the use of AUMIC. 

## References:  
- Bijkerk, A.H.E., and Hall, R.J. (1977). Mechanistic Model of Aerobic Growth of Saccharomyces cerevisiae. Biotechnology and Bioengineering 19, 267-296. doi: DOI 10.1002/bit.260190209.  
- Kim, J.I., Varner, J.D., and Ramkrishna, D. (2008). A Hybrid Model of Anaerobic E. coli GJT001: Combination of Elementary Flux Modes and Cybernetic Variables. Biotechnology Progress 24, 993-1006. doi: 10.1002/btpr.73.  
- Kompala, D.S., Ramkrishna, D., Jansen, N.B., and Tsao, G.T. (1986). Investigation of Bacterial-Growth on Mixed Substrates - Experimental Evaluation of Cybernetic Models. Biotechnology and Bioengineering 28, 1044-1055. doi: DOI 10.1002/bit.260280715.  
- Provost, A., and Bastin, G. (2004). Dynamic metabolic modelling under the balanced growth condition. Journal of Process Control 14, 717-728. doi: 10.1016/j.jprocont.2003.12.004.  
- Ramkrishna, D., and Song, H.-S. (2018). Cybernetic Modeling for Bioreaction Engineering. Cambridge University Press.  
- Ramkrishna, D., and Song, H.S. (2012). Dynamic models of metabolism: Review of the cybernetic approach. Aiche Journal 58, 986-997. doi: 10.1002/aic.13734.  
- Song, H.S., Morgan, J.A., and Ramkrishna, D. (2009). Systematic Development of Hybrid Cybernetic Models: Application to Recombinant Yeast Co-Consuming Glucose and Xylose. Biotechnology and Bioengineering 103, 984-1002. doi: DOI 10.1002/bit.22332.  
- Song, H.S., and Ramkrishna, D. (2010). Prediction of Metabolic Function From Limited Data: Lumped Hybrid Cybernetic Modeling (L-HCM). Biotechnology and Bioengineering 106, 271-284. doi: 10.1002/bit.22692.  
- Song, H.S., and Ramkrishna, D. (2011). Cybernetic Models Based on Lumped Elementary Modes Accurately Predict Strain-Specific Metabolic Function. Biotechnology and Bioengineering 108, 127-140. doi: 10.1002/bit.22922.  
- Terzer, M., and Stelling, J. (2008). Large-scale computation of elementary flux modes with bit pattern trees. Bioinformatics 24, 2229-2235. doi: 10.1093/bioinformatics/btn401.  
- Von Kamp, A., and Schuster, S. (2006). Metatool 5.0: fast and flexible elementary modes analysis. Bioinformatics 22, 1930-1931. doi: 10.1093/bioinformatics/btl267.  
