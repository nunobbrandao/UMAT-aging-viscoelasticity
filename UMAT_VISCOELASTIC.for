                                                                     
C***********************************************************************
C Name: UMAT_VISCOELASTIC_EXPLICIT 											   
C Data:    13/09/2016												   
c Author: NUNO BANDARRINHA BRANDÁO								   
c***********************************************************************
c Descrição:  Arquivo com a implementacao viscoelastica semi-analitica
c             explicita com variacao de tensao constante entre os 
c			  incrementos de tempo(Zienkiewicz et al., 1968) com fator
c             de aging dado pela teoria da solidificacao (Bazant, 1977).
c             O algoritmo trata separadamente a parcela volumetrica e 
c             desviatoria como solido padrao. 
c 
C Esquematizaçao do solido padrao:
c                                           TAU_K/TAU_G  
c                              G0/K0      |-----||-------|
c                           ---/\/\/\-----|              |-----(...)
c										  |-----/\/\/\---|           
C												G1/K1
cc Subrotinas do usuário Abaqus:  UMAT
c
C  Formulação do elemento de interface 2D quadrático com grau de       
C   liberdade de poro, com comportamento do tipo traction-separation,  
C   modelo constitutivo de Mohr-Coulomb com cut-off zero. Inicialização!
C   das tensões in-situ.									           
      SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     1 RPL,DDSDDT,DRPLDE,DRPLDT,
     2 STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,
     3 NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,
     4 CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,JSTEP,KINC)
C
C	DDSDDE - Constitutive Jacobian
C	NDI - Number of direct stress components at this point.
C	NSHR - Number of engineering shear stress components at this point.
C	NTENS - Size of the stress or strain component array (NDI + NSHR).
C	DSTRAN - Relative-displacements as “strains” (STRAN and DSTRAN)
C
      INCLUDE 'ABA_PARAM.INC'
C - 
      CHARACTER*80 CMNAME
      DIMENSION STRESS(NTENS),STATEV(NSTATV),
     1 DDSDDE(NTENS,NTENS),DDSDDT(NTENS),DRPLDE(NTENS),
     2 STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1),
     3 PROPS(NPROPS),COORDS(3),DROT(3,3),DFGRD0(3,3),DFGRD1(3,3),
     4 JSTEP(4),DEVSTRESS(NTENS),DECRSTRAN(NTENS)
C
      INTEGER KELDIM
      REAL*8 RTG,AGE,I1,S1,v,E0,Gv,Kv,K0,G0
      REAL*8, DIMENSION (1) :: E, RT
      REAL*8, ALLOCATABLE, DIMENSION(:) :: DEPSILON, D_VOL_DEPSILON,
     1 D_DEV_DEPSILON, VOL_STRESS, DEV_STRESS,D_STATEV,G,K,GAMMA
C
      IF (NDI.NE.3) THEN
      WRITE (6, *) 'ERROR: THIS SUBROUTINE CANNOT BE USED IN '// 
     1 'PLANE-STRESS CONDITION!'
      CALL XIT
      ENDIF
C  

C --------------------------------------------------------------------
C  PARAMETERS
C --------------------------------------------------------------------
	v = 0.15 
	E0 = 6.852
	E=(/8.04/)
      RT=(/1/)
C	
	G0 = E0/(2*(1+v))
	K0 = E0/(3*(1-2*v))
      KELDIM = SIZE(E)
      ALLOCATE(DEPSILON(NTENS))
      ALLOCATE(D_VOL_DEPSILON(NTENS))
      ALLOCATE(D_DEV_DEPSILON(NTENS))
      ALLOCATE(VOL_STRESS(NTENS)) 
      ALLOCATE(DEV_STRESS(NTENS)) 
      ALLOCATE(D_STATEV(2*NTENS)) 
      ALLOCATE(G(KELDIM))
      ALLOCATE(K(KELDIM))
      ALLOCATE(GAMMA(KELDIM))
	DO J= 1,KELDIM
	  G(J)=E(J)/(2*(1+v))
	  K(J)=E(J)/(3*(1-2*v))
	  GAMMA(J)=E0/E(J)
	END DO                                                          
	AGE = 1.d0
c	0.45959789239450355*EXP(-0.0028069922021957444*(TIME(2)+DTIME/2))
c	1 + 0.7252619855721698*EXP(-0.30491164042326163*(TIME(2)+DTIME/2))
C --------------------------------------------------------------------
C	TRANSFORM THE ENGINEERING DSTRAN INTO NORMAL ONES
C -------------------------------------------------------------------- 
      DO K1 =  1, NDI
        DEPSILON(K1) = DSTRAN(K1) 
      END DO
	DO K1 = NDI + 1, NTENS
        DEPSILON(K1) = DSTRAN(K1)/2 
      END DO
C --------------------------------------------------------------------
C INCREMENT OF VOLUMETRIC STRAIN TENSOR
C --------------------------------------------------------------------
      S1=0.d0
      DO K1 = 1,NDI
         S1 = S1 + DEPSILON(K1) 
      END DO
      DO K1 = 1,NDI
         D_VOL_DEPSILON(K1) = S1/3
      END DO
      DO K1 = NDI + 1, NTENS
         D_VOL_DEPSILON(K1) = 0
      END DO
C --------------------------------------------------------------------
C INCREMENT OF DEVIATORIC STRAIN TENSOR
C --------------------------------------------------------------------	  
      DO K1 = 1,NDI
         D_DEV_DEPSILON(K1) = DEPSILON(K1) - S1/3
      END DO
      DO K1 = NDI + 1, NTENS
         D_DEV_DEPSILON(K1) = DEPSILON(K1)
      END DO
C --------------------------------------------------------------------
C VOLUMETRIC STRESS TENSOR
C --------------------------------------------------------------------
      I1=0.d0
      DO K1 = 1,NDI
         I1 = I1 + STRESS(K1) 
      END DO
      DO K1 = 1,NDI
         VOL_STRESS(K1) = I1/3
      END DO
      DO K1 = NDI + 1, NTENS
         VOL_STRESS(K1) = 0 
      END DO
C --------------------------------------------------------------------
C DEVIATORIC STRESS TENSOR
C --------------------------------------------------------------------
      DO K1 = 1,NDI
         DEV_STRESS(K1) = STRESS(K1) - I1/3
      END DO
      DO K1 = NDI + 1, NTENS
         DEV_STRESS(K1) = STRESS(K1) 
      END DO
C --------------------------------------------------------------------
C INCREMENT OF VISCOELASTIC STRAIN - NEWTON-RAPHSON  
C -------------------------------------------------------------------- 
      D_STATEV=0.0D0
       WRITE(6,*),'-------------------------'
       WRITE(6,*),'INCREMENT=' ,KINC
	 CALL NEWTON_V(Kv,D_STATEV,NTENS,GAMMA,3*K0,RT,DTIME, 
     1 D_VOL_DEPSILON, STATEV,VOL_STRESS,KELDIM)
	 CALL NEWTON_D(Gv,D_STATEV,NTENS,GAMMA,2*G0,RT,DTIME,
     1 D_DEV_DEPSILON,STATEV,DEV_STRESS,KELDIM)
       WRITE(6,*),'D_STATEV = ' ,D_STATEV
C --------------------------------------------------------------------
C ACTUALAZING STATEV
C --------------------------------------------------------------------
      DO K1 = 1, NTENS
        STATEV(K1) = D_STATEV(K1) + STATEV(K1) 
      END DO
      DO K1 = NTENS+1, 2*NTENS
        STATEV(K1) = D_STATEV(K1) + STATEV(K1) 
      END DO
C --------------------------------------------------------------------
C ACTUALAZING DEVIATORIC AND VOLUMETRIC STRESSES
C --------------------------------------------------------------------
      DO K1 = 1, NTENS
         VOL_STRESS(K1)=VOL_STRESS(K1)+
     1 3*K0/AGE*(D_VOL_DEPSILON(K1)-D_STATEV(K1)*AGE) 
      END DO 
	DO K1 = NTENS+1, 2*NTENS
         DEV_STRESS(K1-NTENS)=DEV_STRESS(K1-NTENS)+
     1 2*G0/AGE*(D_DEV_DEPSILON(K1-NTENS)-D_STATEV(K1)*AGE)
      END DO 
C --------------------------------------------------------------------
C ACTUALAZING TOTAL STRESSES
C --------------------------------------------------------------------
      DO K1 = 1, NTENS
         STRESS(K1) = VOL_STRESS(K1) + DEV_STRESS(K1)
      END DO   
C--------------------------------------------------------------------
C ISOTROPIC ELASTIC CONSTITUTIVE MATRIX [C]
C --------------------------------------------------------------------
      DO K1 = 1, NDI
        DO K2 = 1, NDI
            DDSDDE(K2,K1) = K0 - (2/3)*G0
C			DDSDDE(K2,K1) = Kv/3 - (2/3)*Gv/2
        END DO
        DDSDDE(K1,K1) = K0 + (4/3)*G0
C		DDSDDE(K1,K1) = Kv/3 + (4/3)*Gv/2
      END DO
      DO K1 = NDI + 1, NTENS
        DDSDDE(K1, K1) = G0
c        DDSDDE(K1, K1) = Gv/2
      END DO	  
C --------------------------------------------------------------------
      IF (NOEL == 1) THEN
       IF (NPT == 1) THEN
c        IF (KINC == 1) THEN
c        END IF 
       END IF
      END IF
      RETURN
      END
!
      SUBROUTINE NEWTON_V(Cve,D_STATEV,NTENS,GAMMA,E0,RT,DTIME, 
     1 DEPSILON, STATEV,STRESS,KELDIM)
!
      INCLUDE 'ABA_PARAM.INC'
!
      REAL(8) :: ITER_max, TOL2
      PARAMETER (ITER_max=20, TOL2=1e-10)
!   
      INTEGER, INTENT(IN)::NTENS,KELDIM
      REAL(8), INTENT(IN)::E0,DTIME, DEPSILON(NTENS),STATEV(2*NTENS),
     1 STRESS(NTENS),GAMMA(KELDIM),RT(KELDIM)
      REAL(8), INTENT(INOUT)::D_STATEV(2*NTENS)
      REAL(8), INTENT(OUT)::Cve
!   
      REAL(8) ::RES, REST, Tn1(KELDIM), Tn0(KELDIM), R21(NTENS),
     1 R22(NTENS),R23(NTENS),DSTRESS_TRY(NTENS),delta_sigma(NTENS),
     2 delta_h(NTENS),R1(NTENS), R2(NTENS),SUMCv
      integer::k                                    
!    
      k = 1                                     
	REST = 1.0D0
	SUMCv = 0.0D0
	R21=0.0D0
	R22=0.0D0
	R23=0.0D0  
      DO I = 1,KELDIM                            		
	  Tn1(I) = (1-RT(I)/DTIME*(1-EXP(-DTIME/RT(I))))
	  Tn0(I) = 1-EXP(-DTIME/RT(I))-Tn1(I)
	  SUMCv = SUMCv + Tn1(I)*GAMMA(I)
	  DO J = 1,NTENS 
	      R21(J) = R21(J)+(EXP(-DTIME/RT(I))-1)*STATEV(J)           
		  R22(J) = R22(J)+ GAMMA(I)*Tn0(I)*STRESS(J)/E0 
      END DO
	END DO   
      WRITE(6,*) 'DTIME = ', DTIME
      WRITE(6,*) 'RT = ', RT
      WRITE(6,*) 'STATEV = ', STATEV 
      WRITE(6,*) 'STRESS = ', STRESS
	DO WHILE (k<ITER_max .AND. REST > TOL2)
		RES=0.0D0   
		DO J = 1,NTENS 
			DSTRESS_TRY(J)=E0*(DEPSILON(J)-D_STATEV(J))    
			R23(J)=(STRESS(J)+DSTRESS_TRY(J))*SUMCv/E0
			R1(J)=DSTRESS_TRY(J)/E0-(DEPSILON(J)-D_STATEV(J))
			R2(J)=R21(J)+R22(J)+R23(J)-D_STATEV(J)                
			Cve = E0/(1+SUMCv) 
			delta_sigma(J) = -(Cve)*(R1(J)+R2(J))       
			delta_h(J)= R2(J)+ SUMCv*delta_sigma(J)/E0
			DSTRESS_TRY(J)=DSTRESS_TRY(J)+delta_sigma(J)  
			D_STATEV(J)=D_STATEV(J)+delta_h(J)             
			RES = RES + R1(J)**2 + R2(J)**2                        
		END	DO											
	  REST=sqrt(RES)                                      
	  k = k + 1          
	END DO                
      RETURN
      END 
!---------------------------------------------------------      
!---------------------------------------------------------
      SUBROUTINE NEWTON_D(Cve,D_STATEV,NTENS,GAMMA,E0,RT,KELDIM,DTIME, 
     1 DEPSILON, STATEV,STRESS)
!
      INCLUDE 'ABA_PARAM.INC'
!
      REAL(8) :: ITER_max, TOL2
      PARAMETER (ITER_max=20, TOL2=1e-10)
!     
      INTEGER, INTENT(IN)::NTENS,KELDIM
      REAL(8), INTENT(IN)::E0,DTIME, DEPSILON(NTENS), STATEV(2*NTENS),
     1 STRESS(NTENS),GAMMA(KELDIM),RT(KELDIM)
      REAL(8), INTENT(INOUT)::D_STATEV(2*NTENS)
      REAL(8), INTENT(OUT)::Cve
!   
      REAL(8) ::RES, REST, Tn1(KELDIM), Tn0(KELDIM), R21(NTENS),
     1 R22(NTENS),R23(NTENS),DSTRESS_TRY(NTENS),delta_sigma(NTENS),
     2 delta_h(NTENS),R1(NTENS), R2(NTENS),SUMCv
      integer::k                                   
!
      k = 1                                     
	REST = 1.D0
	SUMCv = 0.D0
	DO J = 1,NTENS
	  R21(J)=0.D0 
	  R22(J)=0.D0
	  R23(J)=0.D0
	END DO     
      DO I = 1,KELDIM                            		
	  Tn1(I) = (1-RT(I)/DTIME*(1-EXP(-DTIME/RT(I))))
	  Tn0(I) = 1-EXP(-DTIME/RT(I))-Tn1(I)
	  SUMCv = SUMCv + Tn1(I)*GAMMA(I)
	  DO J = 1,NTENS 
	      R21(J) = R21(J)+(EXP(-DTIME/RT(I))-1)*STATEV(J+NTENS)           
		    R22(J) = R22(J)+ GAMMA(I)*Tn0(I)*STRESS(J)/E0 
        END DO
	END DO
	DO WHILE (k<ITER_max .AND. REST > TOL2)
		RES=0.D0   
		DO J = 1,NTENS 
			DSTRESS_TRY(J)=(K0)*(DEPSILON(J)-D_STATEV(J+NTENS))    
			R23(J)= (STRESS(J)+DSTRESS_TRY(J))*SUMCv/E0  
			R1(J)=DSTRESS_TRY(J)/E0-(DEPSILON(J)-D_STATEV(J+NTENS))
			R2(J)=R21(J)+R22(J)+R23(J)-D_STATEV(J+NTENS)                
			Cve= E0/(1+SUMCv)                                
			delta_sigma(J)=-(Cve)*(R1(J)+R2(J))       
			delta_h(J)=R2(J)+ SUMCv*delta_sigma(J)/E0
			DSTRESS_TRY(J)=DSTRESS_TRY(J)+delta_sigma(J)  
			D_STATEV(J+NTENS)=D_STATEV(J+NTENS)+delta_h(J)             
			RES = RES + R1(J)**2 + R2(J)**2                      
		END	DO											
	  REST=sqrt(RES)                                      
	  k = k + 1                          
	END DO
      RETURN
      END 