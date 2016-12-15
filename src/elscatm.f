C      INCLUDE 'elsepa.f'

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
C    +++++++++++++++++++++++++++++
C    +++    PROGRAM ELSCATM    +++
C    +++++++++++++++++++++++++++++
C
C
C                              F. Salvat, A. Jablonski and C.J. Powell
C                              November 30, 2003
C
C
C  ELastic SCATtering of electrons and positrons by Molecules
C  --      ----                                     -
C
C     This program calculates differential cross sections (DCS) and spin
C  polarization functions for elastic scattering of electrons and
C  positrons by molecules using a single-scattering independent-atom
C  approximation.
C
C  ****  The input data file.
C
C     Data are read from an input file (unit 5); the first line is the
C  name of the target molecule (up to 35 characters, used only for the
C  user's information). Each of the subsequent lines consists of a set
C  of numerical values (in free format) followed by a short text
C  describing the data. This text is a reminder for the user and is not
C  read by the program. The following example corresponds to the water
C  molecule.
C
C
C ----+----1----+----2----+----3----+----4----+----5----+----6----+----7
C Water (vapour)                        name of the molecule
C 3                                     number of atoms in a molecule
C 8     0.0E+0     0.0E+0     0.0E+0    oxygen, coord. (cm)
C 1  9.5720E-9     0.0E+0     0.0E+0    hydrogen, coord. (cm)
C 1 -2.3999E-9  9.2663E-9     0.0E+0    hydrogen, coord. (cm)
C 1                                     MEXCH
C 2  1.457E-24                          MCPOL, mol. polarizab. (cm**3)
C 1  2.00  6.20                         MABS, Aabs, Delta (eV)
C -1                                    electron (-1) / positron (+1)
C 100                                   kinetic energy (eV)
C 200                                   optionally more energies ...
C ----+----1----+----2----+----3----+----4----+----5----+----6----+----7
C
C     The code stops when it finds an inconsistent input datum. The
C  conflicting quantity appears in the last line written on the screen.
C
C     The computed information is delivered in various files with the
C  extension '.dat'. The output file 'dcs_xpyyyezz.dat' contains the
C  calculated DCS for the energy x.yyyEzz (E format, in eV) in a format
C  ready for visualization with a plotting program.
C
C     This program uses the subroutine package 'elsepa.f', which has
C  been inserted into this source file by using an INCLUDE statement
C  (see above).
C
C
      IMPLICIT DOUBLE PRECISION (A-B,D-H,O-Z), COMPLEX*16 (C),
     1   INTEGER*4 (I-N)
      CHARACTER*35 MNAME
      CHARACTER*12 BUFFER,OFILE
      PARAMETER (PI=3.1415926535897932D0)
C
      PARAMETER (A0B=5.291772083D-9)  ! Bohr radius (cm)
      PARAMETER (HREV=27.2113834D0)  ! Hartree energy (eV)
      PARAMETER (SL=137.03599976D0)  ! Speed of light (1/alpha)
      PARAMETER (A0B2=A0B*A0B)
C
      PARAMETER (NATM=64)
      DIMENSION IZA(NATM),XA(NATM),YA(NATM),ZA(NATM)
C  ****  Results from the partial wave calculation (molecules).
      PARAMETER (NGT=650)
      COMMON/DCSMOL/ECSM,TCS1M,TCS2M,ECSI,TCS1I,TCS2I,THM(NGT),XTM(NGT),
     1              DCSC(NGT),DCSI(NGT),SPOLC(NGT),ERRM(NGT),NTABM
C
C  ************  Input data.
C
      READ(5,'(A35)') MNAME
      WRITE(6,'(2X,A35)') MNAME
C  ****  Number of atoms in a molecule.
      READ(5,*) NA
      WRITE(6,*) '  ',NA
      IF(NA.GT.NATM) STOP 'Too many atoms.'
C  ****  Atomic numbers and position coordinates (cm).
      DO I=1,NA
        READ(5,*) IZA(I),XA(I),YA(I),ZA(I)
        WRITE(6,*) IZA(I),XA(I),YA(I),ZA(I)
      ENDDO
C  ****  Electron exchange potential.
      READ(5,*) MEXCH
      IF(MEXCH.LT.1.OR.MEXCH.GT.3) MEXCH=0
      WRITE(6,*) MEXCH
C  ****  Correlation-polarization potential, molecular polarizability.
      READ(5,*) MCPOL,AMPOL
      IF(MCPOL.LT.1.OR.MCPOL.GT.2) MCPOL=0
      IF(MCPOL.EQ.0) AMPOL=1.0D-35
      WRITE(6,*) MCPOL,AMPOL
      IF(AMPOL.LT.-1.0D-35) STOP
C  ****  Absorption potential, first excitation energy.
      READ(5,*) MABS,VABSA,FEXCE
      IF(MABS.NE.1) MABS=0
      IF(VABSA.LT.-1.0D-35) VABSA=0.0D0
      IF(FEXCE.LT.-1.0D-35) FEXCE=0.0D0
      WRITE(6,*) MABS,VABSA,FEXCE
C  ****  Projectile.
      READ(5,*) IELEC
      IF(IELEC.NE.+1) IELEC=-1
      WRITE(6,*) '  ',IELEC
C
C  ************  Partial wave analysis.
C
      OPEN(25,FILE='tcstable.dat')
      WRITE(25,1001)
 1001 FORMAT(1X,'#  Total cross section table (last run of',
     1  ' ''elscatm'').', /1X,'#')
      WRITE(25,'(1X,''#  Molecule: '',A35,
     2  /'' #'',/'' #  Structure:'' )') MNAME
      DO I=1,NA
        WRITE(25,1002) IZA(I),XA(I),YA(I),ZA(I)
 1002   FORMAT(1X,'#  IZ =',I3,',    X,Y,Z = ',1P,3E12.5,' cm')
      ENDDO
      IF(IELEC.EQ.-1) THEN
        WRITE(25,1003)
 1003   FORMAT(1X,'#  Projectile: electron')
      ELSE
        WRITE(25,1004)
 1004   FORMAT(1X,'#  Projectile: positron')
      ENDIF
      WRITE(25,1005)
 1005 FORMAT(1X,'#',/1X,'#   Energy',9X,'ECS',9X,'TCS1',9X,
     1  'TCS2',/1X,'#',4X,'(eV)',8X,'(cm**2)',6X,'(cm**2)',
     2  6X,'(cm**2)',/1X,'#',51('-'),
     3  2X,'<------  incoherent addition  ------>')
C
C
    1 CONTINUE
      READ(5,*,END=9999,ERR=9999) EV
      WRITE(6,*) '   '
      WRITE(6,*) 'E (eV)=',EV
      IF(EV.LT.10.0D0) STOP 'The kinetic energy is too small'
C  ****  DCS output file.
      WRITE(BUFFER,'(1P,E12.5)') EV
      OFILE=BUFFER(2:2)//'p'//BUFFER(4:6)//'e'//BUFFER(11:12)
      OPEN(8,FILE='dcs_'//OFILE(1:8)//'.dat')
C
      CALL ELSEPM(MNAME,IELEC,EV,IZA,XA,YA,ZA,NA,MEXCH,MCPOL,
     1            AMPOL,MABS,VABSA,FEXCE,8)
C
      WRITE(25,1006) EV,ECSM,TCS1M,TCS2M,ECSI,TCS1I,TCS2I
 1006 FORMAT(1X,1P,7E13.5)
      CLOSE(8)
      GO TO 1
C
 9999 CONTINUE
      CLOSE(25)
      END
C  *********************************************************************
C                       SUBROUTINE ELSEPM
C  *********************************************************************
      SUBROUTINE ELSEPM(MNAME,IELEC,EV,IZA,XA,YA,ZA,NA,MEXCH,MCPOL,
     1                  AMPOL,MABS,VABSA,FEXCE,IW)
C
C     This subroutine computes elastic scattering of electrons and
C  positrons by molecules. It uses a single-scattering independent-atom
C  approximation, with scattering amplitudes calculated by using the
C  subroutine package 'elsepa'.
C
C  Input arguments:
C    MNAME ..... name of the molecule (character*35)
C    IELEC ..... electron-positron flag;
C                =-1 for electrons,
C                =+1 for positrons.
C    EV ........ projectile's kinetic energy (in eV).
C    IZA(I) .... atomic number of the I-th atom in the molecule.
C    XA(I), YA(I), ZA(I) ... coordinates (cm) of the I-th atom.
C    NA ........ number of atoms in the molecule (.le.64)
C    MEXCH ..... exchange correction for electrons.
C                  0 --> no exchange correction,
C                  1 --> Furness-McCarthy (FM),
C                  2 --> Thomas-Fermi (TF),
C                  3 --> Riley-Truhlar (RT).
C    MCPOL ..... correlation-polarization correction.
C                  0 --> no correlation-polarization correction,
C                  1 --> Buckingham potential (B),
C                  2 --> Local density approximation (LDA).
C      AMPOL ... molecular induced polarizability (in cm**3).
C    MABS ...... absorption correction (imaginary potential).
C                  0 --> no absorption correction,
C                  1 --> LDA.
C      VABSA ... strength of the absorption potential.
C      FEXCE ... first excitation energy (eV) of the molecule.
C    IW ........ output unit (to be defined in the main program).
C
C  Output (through the common block /DCSMOL/):
C     ECSM ....... total cross section (cm**2).
C     TCS1M ...... 1st transport cross section (cm**2).
C     TCS2M ...... 2nd transport cross section (cm**2).
C     ECSI ....... ICA total cross section (cm**2).
C     TCS1I ...... ICA 1st transport cross section (cm**2).
C     TCS2I ...... ICA 2nd transport cross section (cm**2).
C     THM(I) ...... scattering angles (in deg)
C     XTM(I) ...... values of (1-COS(TH(I)))/2.
C     DCSC(I) .... differential cross section per unit solid
C                  angle at TH(I) (in cm**2/sr).
C     DCSI(I) .... ICA differential cross section per unit solid
C                  angle at TH(I) (in cm**2/sr).
C     SPOLC(I) ... Sherman spin-polarization function at TH(I).
C     ERRM(I) .... relative uncertainty of the computed DCS
C                  values. Estimated from the convergence of the
C                  series.
C     NTABM ...... number of angles in the table.
C
      IMPLICIT DOUBLE PRECISION (A-B,D-H,O-Z), COMPLEX*16 (C),
     1   INTEGER*4 (I-N)
      CHARACTER*35 MNAME
      CHARACTER*12 BUFFER
      PARAMETER (PI=3.1415926535897932D0)
C
      PARAMETER (A0B=5.291772083D-9)  ! Bohr radius (cm)
      PARAMETER (HREV=27.2113834D0)  ! Hartree energy (eV)
      PARAMETER (SL=137.03599976D0)  ! Speed of light (1/alpha)
      PARAMETER (A0B2=A0B*A0B)
C  ****  Molecular structure.
      PARAMETER (NATM=64)
      DIMENSION IZA(NATM),XA(NATM),YA(NATM),ZA(NATM),D(NATM,NATM)
C  ****  Results from the partial wave calculation (atoms).
      PARAMETER (NGT=650)
      COMMON/DCSTAB/ECS,TCS1,TCS2,TH(NGT),XT(NGT),DCST(NGT),SPOL(NGT),
     1              ERROR(NGT),NTAB
      DIMENSION CFA(NATM,NGT),CGA(NATM,NGT),DCSA(NATM,NGT),
     1          ECSA(NATM),TCS1A(NATM),TCS2A(NATM)
C  ****  Results from the partial wave calculation (molecules).
      COMMON/DCSMOL/ECSM,TCS1M,TCS2M,ECSI,TCS1I,TCS2I,THM(NGT),XTM(NGT),
     1              DCSC(NGT),DCSI(NGT),SPOLC(NGT),ERRM(NGT),NTABM
C  ****  Phase shifts.
      PARAMETER (NDM=25000)
      COMMON/PHASES/DP(NDM),DM(NDM),NPH,ISUMP
C
C  ****  Atomic polarizabilities of free atoms (in cm**3), from
C    Thomas M. Miller, 'Atomic and molecular polarizabilities' in
C    CRC Handbook of Chemistry and Physics, Editor-in-chief David
C    R. Linde. 79th ed., 1998-1999, pp. 10-160 to 10-174.
      DIMENSION ATPOL(103)
      DATA ATPOL/   0.666D-24, 0.205D-24, 24.30D-24, 5.600D-24,
     A   3.030D-24, 1.760D-24, 1.100D-24, 0.802D-24, 0.557D-24,
     1   3.956D-25, 24.08D-24, 10.06D-24, 6.800D-24, 5.380D-24,
     A   3.630D-24, 2.900D-24, 2.180D-24, 1.641D-24, 43.40D-24,
     2   22.80D-24, 17.80D-24, 14.60D-24, 12.40D-24, 11.60D-24,
     A   9.400D-24, 8.400D-24, 7.500D-24, 6.800D-24, 6.100D-24,
     3   7.100D-24, 8.120D-24, 6.070D-24, 4.310D-24, 3.770D-24,
     A   3.050D-24, 2.484D-24, 47.30D-24, 27.60D-24, 22.70D-24,
     4   17.90D-24, 15.70D-24, 12.80D-24, 11.40D-24, 9.600D-24,
     A   8.600D-24, 4.800D-24, 7.200D-24, 7.200D-24, 10.20D-24,
     5   7.700D-24, 6.600D-24, 5.500D-24, 5.350D-24, 4.044D-24,
     A   59.60D-24, 39.70D-24, 31.10D-24, 29.60D-24, 28.20D-24,
     6   31.40D-24, 30.10D-24, 28.80D-24, 27.70D-24, 23.50D-24,
     A   25.50D-24, 24.50D-24, 23.60D-24, 22.70D-24, 21.80D-24,
     7   21.00D-24, 21.90D-24, 16.20D-24, 13.10D-24, 11.10D-24,
     A   9.700D-24, 8.500D-24, 7.600D-24, 6.500D-24, 5.800D-24,
     8   5.100D-24, 7.600D-24, 6.800D-24, 7.400D-24, 6.800D-24,
     A   6.000D-24, 5.300D-24, 48.70D-24, 38.30D-24, 32.10D-24,
     9   32.10D-24, 25.40D-24, 24.90D-24, 24.80D-24, 24.50D-24,
     A   23.30D-24, 23.00D-24, 22.70D-24, 20.50D-24, 19.70D-24,
     1   23.80D-24, 18.20D-24, 17.50D-24, 20.00D-24/
C
C  ****  Ionization energies of neutral atoms (in eV).
C        NIST Physical Reference Data.
C        http://sed.nist.gov/PhysRefData/IonEnergy/tblNew.html
C  For astatine (Z=85), the value given below was calculated with
C  the DHFXA code.
      DIMENSION EIONZ(103)
      DATA EIONZ/  13.5984D0, 24.5874D0, 5.39170D0, 9.32270D0,
     A  8.29800D0, 11.2603D0, 14.5341D0, 13.6181D0, 17.4228D0,
     1  21.5646D0, 5.13910D0, 7.64620D0, 5.98580D0, 8.15170D0,
     A  10.4867D0, 10.3600D0, 12.9676D0, 15.7596D0, 4.34070D0,
     2  6.11320D0, 6.56150D0, 6.82810D0, 6.74620D0, 6.76650D0,
     A  7.43400D0, 7.90240D0, 7.88100D0, 7.63980D0, 7.72640D0,
     3  9.39420D0, 5.99930D0, 7.89940D0, 9.78860D0, 9.75240D0,
     A  11.8138D0, 13.9996D0, 4.17710D0, 5.69490D0, 6.21710D0,
     4  6.63390D0, 6.75890D0, 7.09240D0, 7.28000D0, 7.36050D0,
     A  7.45890D0, 8.33690D0, 7.57620D0, 8.99380D0, 5.78640D0,
     5  7.34390D0, 8.60840D0, 9.00960D0, 10.4513D0, 12.1298D0,
     A  3.89390D0, 5.21170D0, 5.57690D0, 5.53870D0, 5.47300D0,
     6  5.52500D0, 5.58200D0, 5.64360D0, 5.67040D0, 6.15010D0,
     A  5.86380D0, 5.93890D0, 6.02150D0, 6.10770D0, 6.18430D0,
     7  6.25420D0, 5.42590D0, 6.82510D0, 7.54960D0, 7.86400D0,
     A  7.83350D0, 8.43820D0, 8.96700D0, 8.95870D0, 9.22550D0,
     8  10.4375D0, 6.10820D0, 7.41670D0, 7.28560D0, 8.41700D0,
     A  9.50000D0, 10.7485D0, 4.07270D0, 5.27840D0, 5.17000D0,
     9  6.30670D0, 5.89000D0, 6.19410D0, 6.26570D0, 6.02620D0,
     A  5.97380D0, 5.99150D0, 6.19790D0, 6.28170D0, 6.42000D0,
     1  6.50000D0, 6.58000D0, 6.65000D0, 4.90000D0/
C
C  ****  Coherent average parameters.
C
      TPOL=0.0D0
      DO I=1,NA
        TPOL=TPOL+ATPOL(IZA(I))
      ENDDO
      FPOL=AMPOL/TPOL
C
      DO I=1,NA
        DO J=1,NA
          D(I,J)=SQRT((XA(I)-XA(J))**2+(YA(I)-YA(J))**2
     1          +(ZA(I)-ZA(J))**2)
        ENDDO
      ENDDO
      IF(FEXCE.LT.-1.0D-35) FEXCE=0.0D0
C
C  ****  Scattering potential model parameters.
C
      IF(IELEC.EQ.-1) THEN
        MEXCH0=MEXCH
      ELSE
        MEXCH0=0
      ENDIF
      MNUCL=3           ! Fermi nuclear charge distributions
      MELEC=4           ! DF atomic electron densities
      MUFFIN=0          ! free atoms
      RMUF=200.0D0      !
      IHEF=0            ! no high-energy factorization
C
      E=EV/HREV
      RK=DSQRT(E*(E+2.0D0*SL*SL))/SL  ! Wavenumber, atomic units
C
      DO 2 I=1,NA
        IZ=IZA(I)
        IF(I.GT.1) THEN
          DO J=1,I-1
            IF(IZ.EQ.IZA(J)) THEN
              ECSA(I)=ECSA(J)
              TCS1A(I)=TCS1A(J)
              TCS2A(I)=TCS2A(J)
              DO K=1,NTAB
                DCSA(I,K)=DCSA(J,K)
                CFA(I,K)=CFA(J,K)
                CGA(I,K)=CGA(J,K)
              ENDDO
              GO TO 2
            ENDIF
          ENDDO
        ENDIF
C
        NELEC=IZ
C  ****  Correlation-polarization potential.
        IF(EV.GT.10.0D3.OR.MCPOL.EQ.0) THEN
          MCPOL0=0
          VPOLA=0.0D0
          VPOLB=0.0D0
        ELSE
          MCPOL0=MCPOL
          VPOLA=FPOL*ATPOL(IZ)
          VPOLB=0.25D0*SQRT(MAX(EV-50.0D0,0.0D0))
        ENDIF
C  ****  Absorption potential,
        IF(EV.GT.1.0D6.OR.MABS.EQ.0) THEN
          MABS0=0
          AABSA=0.0D0
          VABSD=0.0D0
        ELSE
          MABS0=1
          IF(VABSA.GT.-1.0D-15) THEN
            AABSA=VABSA
          ELSE
            AABSA=2.0D0
          ENDIF
          IF(IELEC.EQ.-1) THEN
            VABSD=FEXCE
            IF(VABSD.LT.-1.0D-35) VABSD=0.0D0
          ELSE
            VABSD=MAX(0.0D0,EIONZ(IZ)-6.8D0)
          ENDIF
        ENDIF
C
        IW1=33
        FIW=1000+IZ
        WRITE(BUFFER,'(1P,E12.5)') FIW
        OPEN(IW1,FILE='dcs_'//BUFFER(4:6)//'.dat')
        CALL ELSEPA(IELEC,EV,IZ,NELEC,MNUCL,MELEC,MUFFIN,RMUF,
     1    MEXCH0,MCPOL0,VPOLA,VPOLB,MABS0,AABSA,VABSD,IHEF,IW1)
        CLOSE(IW1)
        ECSA(I)=ECS
        TCS1A(I)=TCS1
        TCS2A(I)=TCS2
        DO K=1,NTAB
          THRAD=TH(K)*PI/180.0D0
          CALL DPWA(THRAD,CCF,CCG,DCS,SPL,ERRF,ERRG)
          CFA(I,K)=CCF
          CGA(I,K)=CCG
          DCSA(I,K)=DCS
        ENDDO
    2 CONTINUE
C
C  ****  Molecular DCS and integrated cross sections.
C
      CI=DCMPLX(0.0D0,1.0D0)
      DO K=1,NTAB
        DCSC(K)=0.0D0
        DCSI(K)=0.0D0
        CSPL1=0.0D0
        CSPL2=0.0D0
C  ----  Momentum transfer in atomic units.
        Q=RK*2.0D0*DSIN(0.5D0*TH(K)*PI/180.0D0)
        DO I=1,NA
          DO J=1,NA
            ARG=Q*D(I,J)/A0B
            IF(ARG.GT.1.0D-9) THEN
              FACT=SIN(ARG)/ARG
            ELSE
              FACT=1.0D0
            ENDIF
            DCSC(K)=DCSC(K)+FACT*(CFA(I,K)*DCONJG(CFA(J,K))
     1             +CGA(I,K)*DCONJG(CGA(J,K)))
            CSPL1=CSPL1+FACT*CFA(I,K)*DCONJG(CGA(J,K))
            CSPL2=CSPL2+FACT*DCONJG(CFA(I,K))*CGA(J,K)
          ENDDO
          DCSI(K)=DCSI(K)+DCSA(I,K)
        ENDDO
        CSPL1=CSPL1*CI
        CSPL2=CSPL2*CI
        TST=CDABS(CSPL1-CSPL2)/(CDABS(CSPL1)+1.0D-30)
        IF(TST.GT.1.0D-3) THEN
          SPOLC(K)=(CSPL1-CSPL2)/DCSC(K)
        ELSE
          SPOLC(K)=0.0D0
        ENDIF
      ENDDO
C  ****  Total cross sections (coherent sum).
      ECS0=4.0D0*PI*RMOM(XT,DCSC,NTAB,0)
      ECS1=4.0D0*PI*RMOM(XT,DCSC,NTAB,1)
      ECS2=4.0D0*PI*RMOM(XT,DCSC,NTAB,2)
      ECSM=ECS0
      TCS1M=2.0D0*ECS1
      TCS2M=6.0D0*(ECS1-ECS2)
C  ****  Total cross sections (incoherent sum).
      ECSI=0.0D0
      TCS1I=0.0D0
      TCS2I=0.0D0
      DO I=1,NA
        ECSI=ECSI+ECSA(I)
        TCS1I=TCS1I+TCS1A(I)
        TCS2I=TCS2I+TCS2A(I)
      ENDDO
C
      WRITE(IW,2000)
 2000 FORMAT(1X,'#',/1X,'# Subroutine ELSEPM. Elastic scattering of',
     1  'electrons and positrons by',/1X,'#',7X,
     2  'molecules (single-scattering independent-atom approximation)')
      WRITE(IW,'('' #'',/'' # Molecule: '',A35,
     2  /'' #'',/'' # Structure:'' )') MNAME
      DO I=1,NA
        WRITE(IW,2001) IZA(I),XA(I),YA(I),ZA(I)
 2001   FORMAT(1X,'# IZ =',I3,',    X,Y,Z = ',1P,3E13.5,' cm')
      ENDDO
C
      WRITE(IW,'(1X,''#'')')
      IF(IELEC.EQ.-1) THEN
        WRITE(IW,2002)
 2002   FORMAT(1X,'# Projectile: electron')
      ELSE
        WRITE(IW,2003)
 2003   FORMAT(1X,'# Projectile: positron')
      ENDIF
      WRITE(IW,2004) EV
 2004 FORMAT(1X,'# Kinetic energy =',1P,E14.7,' eV')
C
C
      IF(MEXCH0.EQ.1) THEN
        WRITE(IW,2005)
 2005   FORMAT(1X,'#',/1X,'# Furness-McCarthy exchange',
     1    ' potential')
      ELSE IF(MEXCH0.EQ.2) THEN
        WRITE(IW,2006)
 2006   FORMAT(1X,'#',/1X,'# Thomas-Fermi exchange potential')
      ELSE IF(MEXCH0.EQ.3) THEN
        WRITE(IW,2007)
 2007   FORMAT(1X,'#',/1X,'# Riley-Truhlar exchange potential')
      ENDIF
C
      IF(MCPOL0.EQ.1) THEN
        WRITE(IW,2008) AMPOL
 2008   FORMAT(1X,'#',/1X,'# Correlation-polarization potential (Buc',
     1    'kingham):',/1X,'#    molecular polarizability =',1P,E12.5,
     2    ' cm**3')
      ELSE IF(MCPOL0.EQ.2) THEN
C  ****  LDA correlation-polarization potential.
        WRITE(IW,2009) AMPOL
 2009   FORMAT(1X,'#',/1X,'# Correlation-polarization potential (LDA',
     1    '):',/1X,'#    molecular polarizability =',1P,E12.5,' cm**3')
      ENDIF
C
      IF(MABS0.NE.0) THEN
        WRITE(IW,2010) VABSD,VABSA
 2010   FORMAT(1X,'#',/1X,'# LDA absorption potential: Delta =',
     1         1P,E12.5,' eV',/1X,'#',28X,'Aabs =',E12.5)
      ENDIF
C
      WRITE(IW,2011)
 2011 FORMAT(1X,'#',/1X,'# Integrated cross sections (cm**2):',/1X,
     1  '#',21X,'coherent',6X,'incoherent')
      WRITE(IW,2012) ECSM,ECSI
 2012 FORMAT(1X,'#         total cs:',1P,E14.7,1X,E14.7)
      WRITE(IW,2013) TCS1M,TCS1I
 2013 FORMAT(1X,'# 1st transport cs:',1P,E14.7,1X,E14.7)
      WRITE(IW,2014) TCS2M,TCS2I
 2014 FORMAT(1X,'# 2nd transport cs:',1P,E14.7,1X,E14.7)
C
      WRITE(IW,'(1X,''#'',/1X,''# Differential cross section:'',
     1  21X,''MU=(1-COS(THETA))/2'')')
      WRITE(IW,'(1X,''#'',/1X,''#   THETA'',11X,''MU'',9X,
     1  ''DCS coh'',6X,''DCS incoh'',6X,''Sherman'',/1X,
     2  ''#   (deg)'',21X,''(cm**2/sr)'',4X,''(cm**2/sr)'',
     3  5X,''function'',/1X,''#'',68(''-''))')
C
      NTABM=NTAB
      DO I=1,NTABM
        THM(I)=TH(I)
        XTM(I)=XT(I)
      ENDDO
      DO K=1,NTAB
        WRITE(IW,'(1P,5E14.6)')
     1    THM(K),XTM(K),DCSC(K),DCSI(K),SPOLC(K)
      ENDDO
C
      RETURN
      END
