# UMAT for an aging viscoelastic material						   

Description: UMAT file with semi-analytical aging viscoelastic implementation with constant stress between time increments (Zienkiewicz et al., 1968) with an aging factor (Bazant, 1977).The algorithm treats the volumetric part of the constitutive tensor as elastic,and the deviatoric part of the constitutive tensor as a solid standard.

Layout of the standard solid:



                                           TAU_K/TAU_G  
                              G0/K0      |-----||-------|
                           ---/\/\/\-----|              |-----(...)
                                         |-----/\/\/\---| 
                                               G1/K1	