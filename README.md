# MergeRatePBH

Code in phyton to plot and obtain the primordial black hole merge rate considering different models.

These codes are based on the work developed in Abram Pérez Herrero's final thesis. We have also used the data obtained from the article https://arxiv.org/abs/2010.14533 , https://arxiv.org/abs/2012.02786 and https://arxiv.org/abs/2102.00868 . The data of the article 2010.14533 is freely available on https://dcc.ligo.org/public/0171/P2000434/002/Fig-3-m1-ppd.tar.gz. and it is necessary to obtain this data in order for the codes to work correctly. Once downloaded, it is recommended to download the data with name 'o1o2o3_mass_c_iid_mag_two_comp_iid_tilt_powerlaw_redshift_mass_data' and insert in a new folder namely DataLVO3. 

Once this action is done, the following codes generate a txt with the data of the constraints on the abundance of primordial black holes: 

 ConstraintsCC.py  # To obtain the constraints for a Critical Collapse mass function\
 ConstraintsLN.py  # To obtain the constraints for a Log-Normal mass function\
 ConstraintsPL.py  # To obtain the constraints for a Power-Law mass function\
 ConstraintsMonochromatic.py # To obtain the constraints for a monochromatic mass function\
 Monochromaticlowmass.py  # To obtain the constraints for a monochromatic mass function and lower masses\
 
Also, all the mass function described in the work are defined in MassFunction.py

The Figures are obtained with the following codes:

**Figure 3.2**\
FigureTotalMargeRateMonochromatic.py\
\
**Figure 3.3**\
FiguredRdmdm2D.py\
\
**Figure 3.4**\
FiguratotalmergeExtendedMassFunction.py\
\
**Figure 4.2**\
FigureMassFunctionSelected.py\
\
**Figure 4.4**\
FigureConstraintsMonochromatic.py\
\
**Figure 4.5**\
FigureConstraintsCC.py\
FigureConstraintsLN.py\
FigureConstraintsPL.py\
\
**Figure 4.6**\
FigureConstraintsCC+Clusterin.py\
FigureConstraintsLN+Clustering.py\
FigureConstraintsPL+Cluestering.py\
\
**Figure 4.7**\
ComparisonMonochromatic.py;\
ComparisonLN.py
