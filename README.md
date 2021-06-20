# MergeRatePBH 

Code in phyton to plot and obtain the primordial black hole merge rate considering different models.

These codes are based on the work developed in Abram Pérez Herrero's final thesis. The data of the article 2010.14533 is freely available on https://dcc.ligo.org/public/0171/P2000434/002/Fig-3-m1-ppd.tar.gz. and it is necessary to obtain this data in order for the codes to work correctly. Once downloaded, it is recommended to download the data with name 'o1o2o3_mass_c_iid_mag_two_comp_iid_tilt_powerlaw_redshift_mass_data' and insert in a new folder namely DataLVO3. 

Once this action is done, the following codes generate a txt with the data of the constraints on the abundance of primordial black holes: 

 ConstraintsCC.py  # To obtain the constraints for a Critical Collapse mass function\
 ConstraintsLN.py  # To obtain the constraints for a Log-Normal mass function\
 ConstraintsPL.py  # To obtain the constraints for a Power-Law mass function\
 ConstraintsMonochromatic.py # To obtain the constraints for a monochromatic mass function\
 Monochromaticlowmass.py  # To obtain the constraints for a monochromatic mass function and lower masses\
 
Also, all the mass function described in the work are defined in MassFunction.py

The Figures are obtained with the following codes:

**Figure 3.2**\
FigureTotalMargeRateMonochromatic.py  ![TMergeRateMonochromatic](https://github.com/AbramPerezHerrero/MergeRatePBH/blob/850d26bf4e370f766f7fe2db536c4f2a2698a4c5/Plots/MergerateMonocrhomatic.pdf)\

**Figure 3.3**\
FiguredRdmdm2D.py  ![DifferentialMergeRate](https://github.com/AbramPerezHerrero/MergeRatePBH/blob/850d26bf4e370f766f7fe2db536c4f2a2698a4c5/Plots/dRdmdmPlot2D.pdf)\
\
**Figure 3.4**\
FiguratotalmergeExtendedMassFunction.py  ![TotalMergeRateExtended](https://github.com/AbramPerezHerrero/MergeRatePBH/blob/7205f45581f6835d1986202d884e20ee98221a00/Plots/TotalFigure.pdf)\
\
**Figure 4.2**\
FigureMassFunctionSelected.py  ![MassFunctionSelected](https://github.com/AbramPerezHerrero/MergeRatePBH/blob/7205f45581f6835d1986202d884e20ee98221a00/Plots/MassFunctionPlotselected.pdf)\
\
**Figure 4.3**\
FiguratotalmergeExtendedMassFunction.py ![TotalMergeClustering](https://github.com/AbramPerezHerrero/MergeRatePBH/blob/7205f45581f6835d1986202d884e20ee98221a00/Plots/TotalFigure+clustering.pdf)\
\
**Figure 4.4**\
FigureConstraintsMonochromatic.py  ![Monochromatic](https://github.com/AbramPerezHerrero/MergeRatePBH/blob/7205f45581f6835d1986202d884e20ee98221a00/Plots/constraintsMonochromatic.pdf)\
\
**Figure 4.5**\
FigureConstraintsMonochromatic.py  ![MonochromaticTotal](https://github.com/AbramPerezHerrero/MergeRatePBH/blob/7205f45581f6835d1986202d884e20ee98221a00/Plots/constraintsMonochromaticTotal.pdf)\
\
**Figure 4.6**\
FigureConstraintsCC.py ![CC](https://github.com/AbramPerezHerrero/MergeRatePBH/blob/1867d2af4edae471a63beb479bc6b445b1ee59a4/Plots/constraintsCC.pdf)\
FigureConstraintsLN.py ![LN](https://github.com/AbramPerezHerrero/MergeRatePBH/blob/1867d2af4edae471a63beb479bc6b445b1ee59a4/Plots/constraintsLN.pdf)\
FigureConstraintsPL.py ![PL](https://github.com/AbramPerezHerrero/MergeRatePBH/blob/1867d2af4edae471a63beb479bc6b445b1ee59a4/Plots/constraintsPL.pdf)\
\
**Figure 4.7**\
FigureConstraintsCC+Clustering.py ![CC+clustering](https://github.com/AbramPerezHerrero/MergeRatePBH/blob/1867d2af4edae471a63beb479bc6b445b1ee59a4/Plots/constraintsCC+Clustering.pdf)\
FigureConstraintsLN+Clustering.py  ![LN+clustering](https://github.com/AbramPerezHerrero/MergeRatePBH/blob/fa0a2b4606c7df0fcf1de74bd539a58e74f493e9/Plots/constrainsLNCl.pdf)\
FigureConstraintsPL+Cluestering.py  ![PL+clustering](https://github.com/AbramPerezHerrero/MergeRatePBH/blob/fa0a2b4606c7df0fcf1de74bd539a58e74f493e9/Plots/constrainsPLCL.pdf)\
\
**Figure 4.8**\
ComparisonMonochromatic.py  ![MonoComparison](https://github.com/AbramPerezHerrero/MergeRatePBH/blob/fa0a2b4606c7df0fcf1de74bd539a58e74f493e9/Plots/constrainsPLCL.pdf)\
ComparisonLN.py  ![LNComparison](https://github.com/AbramPerezHerrero/MergeRatePBH/blob/fa0a2b4606c7df0fcf1de74bd539a58e74f493e9/Plots/constrainsPLCL.pdf)\
## References
We have also used the data obtained from the article https://arxiv.org/abs/2010.14533 , https://arxiv.org/abs/2012.02786 and https://arxiv.org/abs/2102.00868 . Some of this data is obtained by https://apps.automeris.io/wpd/.
