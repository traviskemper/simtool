!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!     Version 5.0 06/29/2012 T. W. Kemper                      !
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!     Calculate dipole momment distribution for an index
!
      SUBROUTINE ATYPNOPLSAA(STN)
      USE structure
      USE const
      USE specify 
      USE potential 
!
      IMPLICIT none
!
      INTEGER :: STN,I
!
      DO I=1,NA(STN)
        IF( TYP(I,STN).EQ."F") TYPN(I,STN) =  1  !  "Fluoride -CH2-F (UA)" 9 18.998 1 
        IF( TYP(I,STN).EQ."C2") TYPN(I,STN) =  2  !  "Fluoride -CH2-F (UA)" 6 14.027 2 
        IF( TYP(I,STN).EQ."C") TYPN(I,STN) =  3  !  "Acetic Acid -COOH (UA)" 6 12.011 3
        IF( TYP(I,STN).EQ."O") TYPN(I,STN) =  4  !  "Acetic Acid >C=O (UA)" 8 15.999 1
        IF( TYP(I,STN).EQ."OH") TYPN(I,STN) =  5  !  "Acetic Acid -OH (UA)" 8 15.999 2
        IF( TYP(I,STN).EQ."C3") TYPN(I,STN) =  6  !  "Acetic Acid CH3- (UA)" 6 15.035 1
        IF( TYP(I,STN).EQ."HO") TYPN(I,STN) =  7  !  "Acetic Acid -OH (UA)" 1 1.008 1
        IF( TYP(I,STN).EQ."C4") TYPN(I,STN) =  8  !  "Methane CH4 (UA)" 6 16.043 0 
        IF( TYP(I,STN).EQ."C3") TYPN(I,STN) =  6  !  "Ethane CH3- (UA)" 6 15.035 1 
        IF( TYP(I,STN).EQ."C3") TYPN(I,STN) =  6  !  "N-Alkane CH3- (UA)" 6 15.035 1 
        IF( TYP(I,STN).EQ."C3") TYPN(I,STN) =  6  !  "Isobutane CH3- (UA)" 6 15.035 1 
        IF( TYP(I,STN).EQ."C3") TYPN(I,STN) =  6  !  "Neopentane CH3- (UA)" 6 15.035 1 
        IF( TYP(I,STN).EQ."C2") TYPN(I,STN) =  2  !  "Alkanes -CH2- (UA)" 6 14.027 2 
        IF( TYP(I,STN).EQ."C9") TYPN(I,STN) =  9  !  "1-Alkene CH2= (UA)" 6 14.027 1 
        IF( TYP(I,STN).EQ."CH") TYPN(I,STN) =  10  !  "Isobutane CH (UA)" 6 13.019 3 
        IF( TYP(I,STN).EQ."C8") TYPN(I,STN) =  11  !  "2-Alkene -CH= (UA)" 6 13.019 2 
        IF( TYP(I,STN).EQ."CD") TYPN(I,STN) =  12  !  "Aromatic CH (UA)" 6 13.019 2 
        IF( TYP(I,STN).EQ."CT") TYPN(I,STN) =  13  !  "Neopentane C (UA)" 6 12.011 4 
        IF( TYP(I,STN).EQ."C7") TYPN(I,STN) =  14  !  "Isobutene >C= (UA)" 6 12.011 3 
        IF( TYP(I,STN).EQ."OH") TYPN(I,STN) =  5  !  "Alcohol OH (UA)" 8 15.999 2 
        IF( TYP(I,STN).EQ."HO") TYPN(I,STN) =  7  !  "Alcohol OH (UA)" 1 1.008 1 
        IF( TYP(I,STN).EQ."C3") TYPN(I,STN) =  6  !  "Methanol CH3- (UA)" 6 15.035 1 
        IF( TYP(I,STN).EQ."C2") TYPN(I,STN) =  2  !  "Ethanol -CH2OH (UA)" 6 14.027 2 
        IF( TYP(I,STN).EQ."SH") TYPN(I,STN) =  15  !  "Hydrogen Sulfide H2S" 16 32.060 2 
        IF( TYP(I,STN).EQ."SH") TYPN(I,STN) =  15  !  "Alkyl Sulfide RSH (UA)" 16 32.060 2
        IF( TYP(I,STN).EQ."S") TYPN(I,STN) =  16  !  "Thioether RSR (UA)" 16 32.060 2 
        IF( TYP(I,STN).EQ."S") TYPN(I,STN) =  16  !  "Disulfide RSSR (UA)" 16 32.060 2 
        IF( TYP(I,STN).EQ."HS") TYPN(I,STN) =  17  !  "Hydrogen Sulfide H2S" 1 1.008 1 
        IF( TYP(I,STN).EQ."HS") TYPN(I,STN) =  17  !  "Alkyl Sulfide RSH (UA)" 1 1.008 1
        IF( TYP(I,STN).EQ."C3") TYPN(I,STN) =  6  !  "Methyl Sulfide CH3 (UA)" 6 15.035 1
        IF( TYP(I,STN).EQ."C2") TYPN(I,STN) =  2  !  "Alkyl Sulfide CH2 (UA)" 6 14.027 2
        IF( TYP(I,STN).EQ."C3") TYPN(I,STN) =  6  !  "Thioether CH3 (UA)" 6 15.035 1 
        IF( TYP(I,STN).EQ."C2") TYPN(I,STN) =  2  !  "Thioether CH2 (UA)" 6 14.027 2 
        IF( TYP(I,STN).EQ."C3") TYPN(I,STN) =  6  !  "Disulfide CH3 (UA)" 6 15.035 1 
        IF( TYP(I,STN).EQ."C2") TYPN(I,STN) =  2  !  "Disulfide CH2 (UA)" 6 14.027 2 
        IF( TYP(I,STN).EQ."NZ") TYPN(I,STN) =  18  !  "Acetonitrile -CN (UA)" 7 14.007 1 
        IF( TYP(I,STN).EQ."CZ") TYPN(I,STN) =  19  !  "Acetonitrile -CN (UA)" 6 12.011 2 
        IF( TYP(I,STN).EQ."C3") TYPN(I,STN) =  6  !  "Acetonitrile CH3 (UA)" 6 15.035 1 
        IF( TYP(I,STN).EQ."OW") TYPN(I,STN) =  20  !  "TIP5P Water O" 8 15.999 4 
        IF( TYP(I,STN).EQ."HW") TYPN(I,STN) =  21  !  "TIP5P Water H" 1 1.008 1 
        IF( TYP(I,STN).EQ."LP") TYPN(I,STN) =  22  !  "TIP5P Water M" 99 0.000 1 
        IF( TYP(I,STN).EQ."DM") TYPN(I,STN) =  23  !  "Dummy Atom" 99 0.000 0  
        IF( TYP(I,STN).EQ."He") TYPN(I,STN) =  24  !  "Helium Atom" 2 4.003 0  
        IF( TYP(I,STN).EQ."Ne") TYPN(I,STN) =  25  !  "Neon Atom" 10 20.179 0  
        IF( TYP(I,STN).EQ."Ar") TYPN(I,STN) =  26  !  "Argon Atom" 18 39.948 0  
        IF( TYP(I,STN).EQ."Kr") TYPN(I,STN) =  27  !  "Krypton Atom" 36 83.800 0  
        IF( TYP(I,STN).EQ."Xe") TYPN(I,STN) =  28  !  "Xenon Atom" 54 131.300 0  
        IF( TYP(I,STN).EQ."CH") TYPN(I,STN) =  10  !  "Isopropanol >CHOH (UA)" 6 13.019 3 
        IF( TYP(I,STN).EQ."CT") TYPN(I,STN) =  13  !  "t-Butanol COH (UA)" 6 12.011 4 
        IF( TYP(I,STN).EQ."OS") TYPN(I,STN) =  29  !  "Ether ROR (UA)" 8 15.999 2 
        IF( TYP(I,STN).EQ."C3") TYPN(I,STN) =  6  !  "Ether CH3-OR (UA)" 6 15.035 1 
        IF( TYP(I,STN).EQ."C2") TYPN(I,STN) =  2  !  "Ether -CH2-OR (UA)" 6 14.027 2 
        IF( TYP(I,STN).EQ."OW") TYPN(I,STN) =  20  !  "TIP3P Water O" 8 15.999 2 
        IF( TYP(I,STN).EQ."HW") TYPN(I,STN) =  21  !  "TIP3P Water H" 1 1.008 1 
        IF( TYP(I,STN).EQ."OW") TYPN(I,STN) =  20  !  "TIP4P Water O" 8 15.999 3 
        IF( TYP(I,STN).EQ."HW") TYPN(I,STN) =  21  !  "TIP4P Water H" 1 1.008 1 
        IF( TYP(I,STN).EQ."LP") TYPN(I,STN) =  22  !  "TIP4P Water M" 99 0.000 1 
        IF( TYP(I,STN).EQ."OW") TYPN(I,STN) =  20  !  "TIP3F Water O" 8 15.999 2 
        IF( TYP(I,STN).EQ."HW") TYPN(I,STN) =  21  !  "TIP3F Water H" 1 1.008 1 
        IF( TYP(I,STN).EQ."C2") TYPN(I,STN) =  2  !  "Methylene Chloride (UA)" 6 14.027 2 
        IF( TYP(I,STN).EQ."Cl") TYPN(I,STN) =  30  !  "Methylene Chloride (UA)" 17 35.453 1 
        IF( TYP(I,STN).EQ."CH") TYPN(I,STN) =  10  !  "Chloroform CHCl3 (UA)" 6 12.011 3 
        IF( TYP(I,STN).EQ."Cl") TYPN(I,STN) =  30  !  "Chloroform CHCl3 (UA)" 17 35.453 1 
        IF( TYP(I,STN).EQ."CT") TYPN(I,STN) =  13  !  "Carbon Tetrachloride" 6 12.011 4  
        IF( TYP(I,STN).EQ."Cl") TYPN(I,STN) =  30  !  "Carbon Tetrachloride" 17 35.453 1  
        IF( TYP(I,STN).EQ."SZ") TYPN(I,STN) =  31  !  "DMSO >S=O (UA)" 16 32.060 3 
        IF( TYP(I,STN).EQ."OY") TYPN(I,STN) =  32  !  "DMSO >S=O (UA)" 8 15.999 1 
        IF( TYP(I,STN).EQ."C3") TYPN(I,STN) =  6  !  "DMSO CH3- (UA)" 6 15.035 1 
        IF( TYP(I,STN).EQ."NT") TYPN(I,STN) =  33  !  "Ammonia NH3" 7 14.007 3  
        IF( TYP(I,STN).EQ."H") TYPN(I,STN) =  34  !  "Ammonia NH3" 1 1.008 1  
        IF( TYP(I,STN).EQ."O") TYPN(I,STN) =  4  !  "DMF C=O (UA)" 8 15.999 1 
        IF( TYP(I,STN).EQ."N") TYPN(I,STN) =  35  !  "DMF CON< (UA)" 7 14.007 3 
        IF( TYP(I,STN).EQ."C") TYPN(I,STN) =  3  !  "DMF C=O (UA)" 6 12.011 3 
        IF( TYP(I,STN).EQ."C3") TYPN(I,STN) =  6  !  "DMF CH3- (UA)" 6 15.035 1 
        IF( TYP(I,STN).EQ."OW") TYPN(I,STN) =  20  !  "SPC Water O" 8 15.999 2 
        IF( TYP(I,STN).EQ."HW") TYPN(I,STN) =  21  !  "SPC Water H" 1 1.008 1 
        IF( TYP(I,STN).EQ."CT") TYPN(I,STN) =  13  !  "Alkane CH3-" 6 12.011 4  
        IF( TYP(I,STN).EQ."CT") TYPN(I,STN) =  13  !  "Alkane -CH2-" 6 12.011 4  
        IF( TYP(I,STN).EQ."CT") TYPN(I,STN) =  13  !  "Alkane >CH-" 6 12.011 4  
        IF( TYP(I,STN).EQ."CT") TYPN(I,STN) =  13  !  "Methane CH4" 6 12.011 4  
        IF( TYP(I,STN).EQ."CT") TYPN(I,STN) =  13  !  "Alkane >C<" 6 12.011 4  
        IF( TYP(I,STN).EQ."HC") TYPN(I,STN) =  36  !  "Alkane H-C" 1 1.008 1  
        IF( TYP(I,STN).EQ."CM") TYPN(I,STN) =  37  !  "Alkene R2-C=" 6 12.011 3  
        IF( TYP(I,STN).EQ."CM") TYPN(I,STN) =  37  !  "Alkene RH-C=" 6 12.011 3  
        IF( TYP(I,STN).EQ."CM") TYPN(I,STN) =  37  !  "Alkene H2-C=" 6 12.011 3  
        IF( TYP(I,STN).EQ."HC") TYPN(I,STN) =  36  !  "Alkene H-C=" 1 1.008 1  
        IF( TYP(I,STN).EQ."CA") TYPN(I,STN) =  38  !  "Aromatic C" 6 12.011 3  
        IF( TYP(I,STN).EQ."HA") TYPN(I,STN) =  39  !  "Aromatic H-C" 1 1.008 1  
        IF( TYP(I,STN).EQ."CA") TYPN(I,STN) =  38  !  "Naphthalene Fusion C" 6 12.011 3 
        IF( TYP(I,STN).EQ."CT") TYPN(I,STN) =  13  !  "Ethyl Benzene CH3-" 6 12.011 4 
        IF( TYP(I,STN).EQ."CT") TYPN(I,STN) =  13  !  "Ethyl Benzene -CH2-" 6 12.011 4 
        IF( TYP(I,STN).EQ."C=") TYPN(I,STN) =  40  !  "Diene =CH-CH=" 6 12.011 3  
        IF( TYP(I,STN).EQ."OH") TYPN(I,STN) =  5  !  "Alcohol -OH" 8 15.999 2  
        IF( TYP(I,STN).EQ."HO") TYPN(I,STN) =  7  !  "Alcohol -OH" 1 1.008 1  
        IF( TYP(I,STN).EQ."HC") TYPN(I,STN) =  36  !  "Methanol CH3-" 1 1.008 1  
        IF( TYP(I,STN).EQ."CT") TYPN(I,STN) =  13  !  "Alcohol CH3OH & RCH2OH" 6 12.011 4
        IF( TYP(I,STN).EQ."CT") TYPN(I,STN) =  13  !  "Alcohol R2CHOH" 6 12.011 4  
        IF( TYP(I,STN).EQ."CT") TYPN(I,STN) =  13  !  "Alcohol R3COH" 6 12.011 4  
        IF( TYP(I,STN).EQ."CT") TYPN(I,STN) =  13  !  "Trifluoroethanol -CH2-" 6 12.011 4  
        IF( TYP(I,STN).EQ."CT") TYPN(I,STN) =  13  !  "Trifluoroethanol CF3-" 6 12.011 4  
        IF( TYP(I,STN).EQ."OH") TYPN(I,STN) =  5  !  "Trifluoroethanol -OH" 8 15.999 2  
        IF( TYP(I,STN).EQ."HO") TYPN(I,STN) =  7  !  "Trifluoroethanol -OH" 1 1.008 1  
        IF( TYP(I,STN).EQ."F") TYPN(I,STN) =  1  !  "Trifluoroethanol F" 9 18.998 1  
        IF( TYP(I,STN).EQ."HC") TYPN(I,STN) =  36  !  "Trifluoroethanol -CH2-" 1 1.008 1  
        IF( TYP(I,STN).EQ."CA") TYPN(I,STN) =  38  !  "Phenol C-OH" 6 12.011 3  
        IF( TYP(I,STN).EQ."OH") TYPN(I,STN) =  5  !  "Phenol -OH" 8 15.999 2  
        IF( TYP(I,STN).EQ."HO") TYPN(I,STN) =  7  !  "Phenol -OH" 1 1.008 1  
        IF( TYP(I,STN).EQ."OH") TYPN(I,STN) =  5  !  "Diol -OH" 8 15.999 2  
        IF( TYP(I,STN).EQ."HO") TYPN(I,STN) =  7  !  "Diol -OH" 1 1.008 1  
        IF( TYP(I,STN).EQ."OH") TYPN(I,STN) =  5  !  "Triol -OH" 8 15.999 2  
        IF( TYP(I,STN).EQ."HO") TYPN(I,STN) =  7  !  "Triol -OH" 1 1.008 1  
        IF( TYP(I,STN).EQ."CT") TYPN(I,STN) =  13  !  "Diol & Triol -CH2OH" 6 12.011 4
        IF( TYP(I,STN).EQ."CT") TYPN(I,STN) =  13  !  "Diol & Triol -CHROH" 6 12.011 4
        IF( TYP(I,STN).EQ."CT") TYPN(I,STN) =  13  !  "Diol & Triol -CR2OH" 6 12.011 4
        IF( TYP(I,STN).EQ."HC") TYPN(I,STN) =  36  !  "Diol & Triol H-COH" 1 1.008 1
        IF( TYP(I,STN).EQ."OS") TYPN(I,STN) =  29  !  "Diphenyl Ether" 8 15.999 2  
        IF( TYP(I,STN).EQ."C=") TYPN(I,STN) =  40  !  "Diene =CR-CR=" 6 12.011 3  
        IF( TYP(I,STN).EQ."OS") TYPN(I,STN) =  29  !  "Anisole -OCH3" 8 15.999 1  
        IF( TYP(I,STN).EQ."OS") TYPN(I,STN) =  29  !  "Dialkyl Ether -O-" 8 15.999 2 
        IF( TYP(I,STN).EQ."CT") TYPN(I,STN) =  13  !  "Methyl Ether CH3OR" 6 12.011 4 
        IF( TYP(I,STN).EQ."CT") TYPN(I,STN) =  13  !  "Ethyl Ether -CH2OR" 6 12.011 4 
        IF( TYP(I,STN).EQ."CT") TYPN(I,STN) =  13  !  "Isopropyl Ether >CHOR" 6 12.011 4 
        IF( TYP(I,STN).EQ."CT") TYPN(I,STN) =  13  !  "t-Butyl Ether COR" 6 12.011 4 
        IF( TYP(I,STN).EQ."HC") TYPN(I,STN) =  36  !  "Alkyl Ether H-COR" 1 1.008 1 
        IF( TYP(I,STN).EQ."OS") TYPN(I,STN) =  29  !  "Acetal RO-CR2OX" 8 15.999 4  
        IF( TYP(I,STN).EQ."OH") TYPN(I,STN) =  5  !  "Hemiacetal -OH" 8 15.999 2  
        IF( TYP(I,STN).EQ."HO") TYPN(I,STN) =  7  !  "Hemiacetal -OH" 1 1.008 1  
        IF( TYP(I,STN).EQ."CO") TYPN(I,STN) =  41  !  "Acetal RO-CH2-OR" 6 12.011 4  
        IF( TYP(I,STN).EQ."HC") TYPN(I,STN) =  36  !  "Acetal RO-CH2-OR" 1 1.008 1  
        IF( TYP(I,STN).EQ."CO") TYPN(I,STN) =  41  !  "Hemiacetal RO-CH2-OH" 6 12.011 4  
        IF( TYP(I,STN).EQ."HC") TYPN(I,STN) =  36  !  "Hemiacetal RO-CH2-OH" 1 1.008 1  
        IF( TYP(I,STN).EQ."CO") TYPN(I,STN) =  41  !  "Acetal RO-CHR-OR" 6 12.011 4  
        IF( TYP(I,STN).EQ."HC") TYPN(I,STN) =  36  !  "Acetal RO-CHR-OR" 1 1.008 1  
        IF( TYP(I,STN).EQ."CO") TYPN(I,STN) =  41  !  "Hemiacetal RO-CHR-OH" 6 12.011 4  
        IF( TYP(I,STN).EQ."HC") TYPN(I,STN) =  36  !  "Hemiacetal RO-CHR-OH" 1 1.008 1  
        IF( TYP(I,STN).EQ."CO") TYPN(I,STN) =  41  !  "Acetal RO-CR2-OR" 6 12.011 4  
        IF( TYP(I,STN).EQ."CO") TYPN(I,STN) =  41  !  "Hemiacetal RO-CR2-OH" 6 12.011 4  
        IF( TYP(I,STN).EQ."CA") TYPN(I,STN) =  38  !  "Anisole C-OCH3" 6 12.011 4  
        IF( TYP(I,STN).EQ."SH") TYPN(I,STN) =  15  !  "Thiol -SH" 16 32.060 2  
        IF( TYP(I,STN).EQ."SH") TYPN(I,STN) =  15  !  "Hydrogen Sulfide H2S" 16 32.060 2 
        IF( TYP(I,STN).EQ."S") TYPN(I,STN) =  16  !  "Sulfide -S-" 16 32.060 2  
        IF( TYP(I,STN).EQ."S") TYPN(I,STN) =  16  !  "Disulfide -S-S-" 16 32.060 2  
        IF( TYP(I,STN).EQ."HS") TYPN(I,STN) =  17  !  "Thiol -SH" 1 1.008 1  
        IF( TYP(I,STN).EQ."HS") TYPN(I,STN) =  17  !  "Hydrogen Sulfide H2S" 1 1.008 1 
        IF( TYP(I,STN).EQ."CT") TYPN(I,STN) =  13  !  "Thiol -CH2-SH" 6 12.011 4  
        IF( TYP(I,STN).EQ."CT") TYPN(I,STN) =  13  !  "Thiol >CH-SH" 6 12.011 4  
        IF( TYP(I,STN).EQ."CT") TYPN(I,STN) =  13  !  "Thiol C-SH" 6 12.011 4  
        IF( TYP(I,STN).EQ."CT") TYPN(I,STN) =  13  !  "Methyl Sulfide CH3-SR" 6 12.011 4 
        IF( TYP(I,STN).EQ."CT") TYPN(I,STN) =  13  !  "Sulfide RCH2-SR" 6 12.011 4  
        IF( TYP(I,STN).EQ."CT") TYPN(I,STN) =  13  !  "Sulfide R2CH-SR" 6 12.011 4  
        IF( TYP(I,STN).EQ."CT") TYPN(I,STN) =  13  !  "Sulfide R3C-SR" 6 12.011 4  
        IF( TYP(I,STN).EQ."CT") TYPN(I,STN) =  13  !  "Disulfide CH3-S-SR" 6 12.011 4  
        IF( TYP(I,STN).EQ."CT") TYPN(I,STN) =  13  !  "Disulfide RCH2-S-SR" 6 12.011 4  
        IF( TYP(I,STN).EQ."CT") TYPN(I,STN) =  13  !  "Disulfide R2CH-S-SR" 6 12.011 4  
        IF( TYP(I,STN).EQ."CT") TYPN(I,STN) =  13  !  "Disulfide R3C-S-SR" 6 12.011 4  
        IF( TYP(I,STN).EQ."CT") TYPN(I,STN) =  13  !  "Methanethiol CH3-SH" 6 12.011 4  
        IF( TYP(I,STN).EQ."CT") TYPN(I,STN) =  13  !  "Benzyl Alcohol -CH2OH" 6 12.011 4 
        IF( TYP(I,STN).EQ."CT") TYPN(I,STN) =  13  !  "Benzyl Alcohol -CHROH" 6 12.011 4 
        IF( TYP(I,STN).EQ."CT") TYPN(I,STN) =  13  !  "Benzyl Alcohol -CR2OH" 6 12.011 4 
        IF( TYP(I,STN).EQ."CA") TYPN(I,STN) =  38  !  "Benzyl Alcohol/Nitrile" 6 12.011 3  
        IF( TYP(I,STN).EQ."S") TYPN(I,STN) =  16  !  "Thioanisole -SCH3" 16 32.060 2  
        IF( TYP(I,STN).EQ."CT") TYPN(I,STN) =  13  !  "RCH2-NH2 & GLY CA" 6 12.011 4
        IF( TYP(I,STN).EQ."CT") TYPN(I,STN) =  13  !  "RCHR-NH2 & ALA CA" 6 12.011 4
        IF( TYP(I,STN).EQ."CT") TYPN(I,STN) =  13  !  "R3C-NH2 & AIB CA" 6 12.011 4
        IF( TYP(I,STN).EQ."Cl") TYPN(I,STN) =  30  !  "Chloroalkene Cl-CH=" 17 35.453 1  
        IF( TYP(I,STN).EQ."CM") TYPN(I,STN) =  37  !  "Chloroalkene Cl-CH=" 6 12.011 3  
        IF( TYP(I,STN).EQ."CA") TYPN(I,STN) =  38  !  "Thioanisole C-SCH3" 6 12.011 3  
        IF( TYP(I,STN).EQ."CT") TYPN(I,STN) =  13  !  "Amide -NH-CHR2" 6 12.011 4  
        IF( TYP(I,STN).EQ."CT") TYPN(I,STN) =  13  !  "Amide -NH-CR3" 6 12.011 4  
        IF( TYP(I,STN).EQ."C") TYPN(I,STN) =  3  !  "Benzophenone C=O" 6 12.011 3  
        IF( TYP(I,STN).EQ."C") TYPN(I,STN) =  3  !  "Benzaldehyde C=O" 6 12.011 3  
        IF( TYP(I,STN).EQ."C") TYPN(I,STN) =  3  !  "Acetophenone C=O" 6 12.011 3  
        IF( TYP(I,STN).EQ."C") TYPN(I,STN) =  3  !  "Benzamide C=O" 6 12.011 3  
        IF( TYP(I,STN).EQ."C") TYPN(I,STN) =  3  !  "Amide C=O" 6 12.011 3  
        IF( TYP(I,STN).EQ."O") TYPN(I,STN) =  4  !  "Amide C=O" 8 15.999 1  
        IF( TYP(I,STN).EQ."N") TYPN(I,STN) =  35  !  "Amide -CO-NH2" 7 14.007 3  
        IF( TYP(I,STN).EQ."N") TYPN(I,STN) =  35  !  "Amide -CO-NHR" 7 14.007 3  
        IF( TYP(I,STN).EQ."N") TYPN(I,STN) =  35  !  "Amide -CO-NR2" 7 14.007 3  
        IF( TYP(I,STN).EQ."H") TYPN(I,STN) =  34  !  "Amide -CO-NH2" 1 1.008 1  
        IF( TYP(I,STN).EQ."H") TYPN(I,STN) =  34  !  "Amide -CO-NHR" 1 1.008 1  
        IF( TYP(I,STN).EQ."CT") TYPN(I,STN) =  13  !  "Amide -NH-CH3" 6 12.011 4  
        IF( TYP(I,STN).EQ."CT") TYPN(I,STN) =  13  !  "Amide -NR-CH3" 6 12.011 4  
        IF( TYP(I,STN).EQ."CT") TYPN(I,STN) =  13  !  "Amide -NH-CH2R" 6 12.011 4  
        IF( TYP(I,STN).EQ."CT") TYPN(I,STN) =  13  !  "Amide -NR-CH2R & PRO CD" 6 12.011
        IF( TYP(I,STN).EQ."CT") TYPN(I,STN) =  13  !  "Amide -NR-CHR2 & PRO CA" 6 12.011
        IF( TYP(I,STN).EQ."C") TYPN(I,STN) =  3  !  "Urea C=O" 6 12.011 3  
        IF( TYP(I,STN).EQ."O") TYPN(I,STN) =  4  !  "Urea C=O" 8 15.999 1  
        IF( TYP(I,STN).EQ."N") TYPN(I,STN) =  35  !  "Urea -NH2" 7 14.007 3  
        IF( TYP(I,STN).EQ."H") TYPN(I,STN) =  34  !  "Urea -NH2" 1 1.008 1  
        IF( TYP(I,STN).EQ."N") TYPN(I,STN) =  35  !  "Imide -NH-" 7 14.007 3  
        IF( TYP(I,STN).EQ."C") TYPN(I,STN) =  3  !  "Imide C=O" 6 12.011 3  
        IF( TYP(I,STN).EQ."O") TYPN(I,STN) =  4  !  "Imide C=O" 8 15.999 1  
        IF( TYP(I,STN).EQ."H") TYPN(I,STN) =  34  !  "Imide -NH-" 1 1.008 1  
        IF( TYP(I,STN).EQ."HC") TYPN(I,STN) =  36  !  "Formimide H-C=O" 1 1.008 1  
        IF( TYP(I,STN).EQ."CT") TYPN(I,STN) =  13  !  "Imide CH3-CONHCO-" 6 12.011 4  
        IF( TYP(I,STN).EQ."CT") TYPN(I,STN) =  13  !  "Imide -CH2-CONHCO-" 6 12.011 4  
        IF( TYP(I,STN).EQ."CT") TYPN(I,STN) =  13  !  "Imide >CH-CONHCO-" 6 12.011 4  
        IF( TYP(I,STN).EQ."CT") TYPN(I,STN) =  13  !  "Imide C-CONHCO-" 6 12.011 4  
        IF( TYP(I,STN).EQ."CA") TYPN(I,STN) =  38  !  "Benzonitrile C-CN" 6 12.011 3  
        IF( TYP(I,STN).EQ."CZ") TYPN(I,STN) =  19  !  "Benzonitrile -CN" 6 12.011 2  
        IF( TYP(I,STN).EQ."NZ") TYPN(I,STN) =  18  !  "Benzonitrile -CN" 7 14.007 1  
        IF( TYP(I,STN).EQ."CA") TYPN(I,STN) =  38  !  "Chlorobenzene C-Cl" 6 12.011 3  
        IF( TYP(I,STN).EQ."Cl") TYPN(I,STN) =  30  !  "Chlorobenzene C-Cl" 17 35.453 1  
        IF( TYP(I,STN).EQ."N") TYPN(I,STN) =  35  !  "N-Phenylacetamide N" 7 14.007 3  
        IF( TYP(I,STN).EQ."CA") TYPN(I,STN) =  38  !  "N-Phenylacetamide N-CA" 6 12.011 3  
        IF( TYP(I,STN).EQ."C") TYPN(I,STN) =  3  !  "Carboxylic Acid -COOH" 6 12.011 3 
        IF( TYP(I,STN).EQ."OH") TYPN(I,STN) =  5  !  "Carboxylic Acid C=O" 8 15.999 1 
        IF( TYP(I,STN).EQ."O") TYPN(I,STN) =  4  !  "Carboxylic Acid -OH" 8 15.999 2 
        IF( TYP(I,STN).EQ."HO") TYPN(I,STN) =  7  !  "Carboxylic Acid -COOH" 1 1.008 1 
        IF( TYP(I,STN).EQ."C") TYPN(I,STN) =  3  !  "Carboxylate COO-" 6 12.011 3  
        IF( TYP(I,STN).EQ."O2") TYPN(I,STN) =  42  !  "Carboxylate COO-" 8 15.999 1  
        IF( TYP(I,STN).EQ."CT") TYPN(I,STN) =  13  !  "Carboxylate CH3-COO-" 6 12.011 4  
        IF( TYP(I,STN).EQ."CT") TYPN(I,STN) =  13  !  "Carboxylate RCH2-COO-" 6 12.011 4  
        IF( TYP(I,STN).EQ."CT") TYPN(I,STN) =  13  !  "Carboxylate R2CH-COO-" 6 12.011 4  
        IF( TYP(I,STN).EQ."CT") TYPN(I,STN) =  13  !  "Carboxylate R3C-COO-" 6 12.011 4  
        IF( TYP(I,STN).EQ."C") TYPN(I,STN) =  3  !  "Aldehyde/Acyl Halide C=O" 6 12.011 3 
        IF( TYP(I,STN).EQ."O") TYPN(I,STN) =  4  !  "Aldyhyde/Acyl Halide C=O" 8 15.999 1 
        IF( TYP(I,STN).EQ."HC") TYPN(I,STN) =  36  !  "Aldehyde/Formamide HCO-" 1 1.008 1  
        IF( TYP(I,STN).EQ."C") TYPN(I,STN) =  3  !  "Ketone C=O" 6 12.011 3  
        IF( TYP(I,STN).EQ."O") TYPN(I,STN) =  4  !  "Ketone C=O" 8 15.999 1  
        IF( TYP(I,STN).EQ."HC") TYPN(I,STN) =  36  !  "Acyl H-C-COX" 1 1.008 1  
        IF( TYP(I,STN).EQ."CT") TYPN(I,STN) =  13  !  "C-Terminal ALA CA" 6 12.011 4 
        IF( TYP(I,STN).EQ."CT") TYPN(I,STN) =  13  !  "C-Terminal GLY CA" 6 12.011 4 
        IF( TYP(I,STN).EQ."CT") TYPN(I,STN) =  13  !  "C-Terminal AIB CA" 6 12.011 4 
        IF( TYP(I,STN).EQ."CT") TYPN(I,STN) =  13  !  "C-Terminal PRO CA" 6 12.011 4 
        IF( TYP(I,STN).EQ."N3") TYPN(I,STN) =  43  !  "Ammonium NH4+" 7 14.007 4  
        IF( TYP(I,STN).EQ."N3") TYPN(I,STN) =  43  !  "Ammonium RNH3+" 7 14.007 4  
        IF( TYP(I,STN).EQ."N3") TYPN(I,STN) =  43  !  "Ammonium R4N+" 7 14.007 4  
        IF( TYP(I,STN).EQ."H3") TYPN(I,STN) =  44  !  "Ammonium NH4+" 1 1.008 1  
        IF( TYP(I,STN).EQ."H3") TYPN(I,STN) =  44  !  "Ammonium RNH3+" 1 1.008 1  
        IF( TYP(I,STN).EQ."CT") TYPN(I,STN) =  13  !  "Ammonium CH3-NH3+" 6 12.011 4  
        IF( TYP(I,STN).EQ."CT") TYPN(I,STN) =  13  !  "CH3NH3+/N-Term GLY CA" 6 12.011 4 
        IF( TYP(I,STN).EQ."CT") TYPN(I,STN) =  13  !  "RCH2NH3+/N-Term ALA CA" 6 12.011 4 
        IF( TYP(I,STN).EQ."CT") TYPN(I,STN) =  13  !  "R3C-NH3+/N-Term AIB CA" 6 12.011 4 
        IF( TYP(I,STN).EQ."CT") TYPN(I,STN) =  13  !  "N-Terminal PRO CA" 6 12.011 4 
        IF( TYP(I,STN).EQ."CT") TYPN(I,STN) =  13  !  "N-Terminal PRO CD" 6 12.011 4 
        IF( TYP(I,STN).EQ."CT") TYPN(I,STN) =  13  !  "Ammonium CH3-NH2R+" 6 12.011 4  
        IF( TYP(I,STN).EQ."CT") TYPN(I,STN) =  13  !  "GLY Zwitterion CA" 6 12.011 4 
        IF( TYP(I,STN).EQ."CT") TYPN(I,STN) =  13  !  "ALA Zwitterion CA" 6 12.011 4 
        IF( TYP(I,STN).EQ."N2") TYPN(I,STN) =  45  !  "Guanidinium -NH2" 7 14.007 3  
        IF( TYP(I,STN).EQ."H3") TYPN(I,STN) =  44  !  "Guanidinium -NH2" 1 1.008 1  
        IF( TYP(I,STN).EQ."CA") TYPN(I,STN) =  38  !  "Guanidinium C+" 6 12.011 3  
        IF( TYP(I,STN).EQ."N2") TYPN(I,STN) =  45  !  "Guanidinium -NHR" 7 14.007 3  
        IF( TYP(I,STN).EQ."H3") TYPN(I,STN) =  44  !  "Guanidinium -NHR" 1 1.008 1  
        IF( TYP(I,STN).EQ."CT") TYPN(I,STN) =  13  !  "Me Guanidinium CH3-" 6 12.011 4 
        IF( TYP(I,STN).EQ."CT") TYPN(I,STN) =  13  !  "Et Guanidinium CH3-" 6 12.011 4 
        IF( TYP(I,STN).EQ."CT") TYPN(I,STN) =  13  !  "Et Guan -CH2- & ARG CD" 6
        IF( TYP(I,STN).EQ."CT") TYPN(I,STN) =  13  !  "Pr Guan -CH2- & ARG CG" 6
        IF( TYP(I,STN).EQ."N3") TYPN(I,STN) =  43  !  "Ammonium R2NH2+" 7 14.007 4  
        IF( TYP(I,STN).EQ."H3") TYPN(I,STN) =  44  !  "Ammonium R2NH2+" 1 1.008 1  
        IF( TYP(I,STN).EQ."NC") TYPN(I,STN) =  46  !  "Diaminopyridine N1" 7 14.007 2  
        IF( TYP(I,STN).EQ."CA") TYPN(I,STN) =  38  !  "Diaminopyridine C2" 6 12.011 3  
        IF( TYP(I,STN).EQ."N2") TYPN(I,STN) =  45  !  "Diaminopyridine -NH2" 7 14.007 3  
        IF( TYP(I,STN).EQ."H") TYPN(I,STN) =  34  !  "Diaminopyridine -NH2" 1 1.008 3  
        IF( TYP(I,STN).EQ."CA") TYPN(I,STN) =  38  !  "Diaminopyridine C3" 6 12.011 3  
        IF( TYP(I,STN).EQ."HA") TYPN(I,STN) =  39  !  "Diaminopyridine H3" 1 1.008 1  
        IF( TYP(I,STN).EQ."CA") TYPN(I,STN) =  38  !  "Diaminopyridine C4" 6 12.011 3  
        IF( TYP(I,STN).EQ."HA") TYPN(I,STN) =  39  !  "Diaminopyridine H4" 1 1.008 1  
        IF( TYP(I,STN).EQ."NA") TYPN(I,STN) =  47  !  "Uracil & Thymine N1" 7 14.007 3
        IF( TYP(I,STN).EQ."C") TYPN(I,STN) =  3  !  "Uracil & Thymine C2" 6 12.011 3
        IF( TYP(I,STN).EQ."NA") TYPN(I,STN) =  47  !  "Uracil & Thymine N3" 7 14.007 3
        IF( TYP(I,STN).EQ."C") TYPN(I,STN) =  3  !  "Uracil & Thymine C4" 6 12.011 3
        IF( TYP(I,STN).EQ."CM") TYPN(I,STN) =  37  !  "Uracil & Thymine C5" 6 12.011 3
        IF( TYP(I,STN).EQ."CM") TYPN(I,STN) =  37  !  "Uracil & Thymine C6" 6 12.011 3
        IF( TYP(I,STN).EQ."H") TYPN(I,STN) =  34  !  "Uracil & Thymine HN1" 1 1.008 1
        IF( TYP(I,STN).EQ."O") TYPN(I,STN) =  4  !  "Uracil & Thymine O2" 8 15.999 1
        IF( TYP(I,STN).EQ."H") TYPN(I,STN) =  34  !  "Uracil & Thymine HN3" 1 1.008 1
        IF( TYP(I,STN).EQ."O") TYPN(I,STN) =  4  !  "Uracil & Thymine O4" 8 15.999 1
        IF( TYP(I,STN).EQ."HC") TYPN(I,STN) =  36  !  "Uracil & Thymine HC5" 1 1.008 1
        IF( TYP(I,STN).EQ."HC") TYPN(I,STN) =  36  !  "Uracil & Thymine HC6" 1 1.008 1
        IF( TYP(I,STN).EQ."CT") TYPN(I,STN) =  13  !  "Thymine CH3-" 6 12.011 4  
        IF( TYP(I,STN).EQ."HC") TYPN(I,STN) =  36  !  "Thymine CH3-" 1 1.008 1  
        IF( TYP(I,STN).EQ."NA") TYPN(I,STN) =  47  !  "Cytosine N1" 7 14.007 3  
        IF( TYP(I,STN).EQ."C") TYPN(I,STN) =  3  !  "Cytosine C2" 6 12.011 3  
        IF( TYP(I,STN).EQ."NC") TYPN(I,STN) =  46  !  "Cytosine N3" 7 14.007 2  
        IF( TYP(I,STN).EQ."CA") TYPN(I,STN) =  38  !  "Cytosine C4" 6 12.011 3  
        IF( TYP(I,STN).EQ."CM") TYPN(I,STN) =  37  !  "Cytosine C5" 6 12.011 3  
        IF( TYP(I,STN).EQ."CM") TYPN(I,STN) =  37  !  "Cytosine C6" 6 12.011 3  
        IF( TYP(I,STN).EQ."H") TYPN(I,STN) =  34  !  "Cytosine HN1" 1 1.008 1  
        IF( TYP(I,STN).EQ."O") TYPN(I,STN) =  4  !  "Cytosine O2" 8 15.999 1  
        IF( TYP(I,STN).EQ."N2") TYPN(I,STN) =  45  !  "Cytosine NH2-" 7 14.007 3  
        IF( TYP(I,STN).EQ."H") TYPN(I,STN) =  34  !  "Cytosine NH2- (N3)" 1 1.008 1 
        IF( TYP(I,STN).EQ."H") TYPN(I,STN) =  34  !  "Cytosine NH2- (C5)" 1 1.008 1 
        IF( TYP(I,STN).EQ."HC") TYPN(I,STN) =  36  !  "Cytosine HC5" 1 1.008 1  
        IF( TYP(I,STN).EQ."H4") TYPN(I,STN) =  48  !  "Cytosine HC6" 1 1.008 1  
        IF( TYP(I,STN).EQ."NC") TYPN(I,STN) =  46  !  "Adenine N1" 7 14.007 2  
        IF( TYP(I,STN).EQ."CQ") TYPN(I,STN) =  49  !  "Adenine C2" 6 12.011 3  
        IF( TYP(I,STN).EQ."NC") TYPN(I,STN) =  46  !  "Adenine N3" 7 14.007 2  
        IF( TYP(I,STN).EQ."CB") TYPN(I,STN) =  50  !  "Adenine C4" 6 12.011 3  
        IF( TYP(I,STN).EQ."CB") TYPN(I,STN) =  50  !  "Adenine C5" 6 12.011 3  
        IF( TYP(I,STN).EQ."CA") TYPN(I,STN) =  38  !  "Adenine C6" 6 12.011 3  
        IF( TYP(I,STN).EQ."NB") TYPN(I,STN) =  51  !  "Adenine & Guanine N7" 7 14.007 2
        IF( TYP(I,STN).EQ."CK") TYPN(I,STN) =  52  !  "Adenine & Guanine C8" 6 12.011 3
        IF( TYP(I,STN).EQ."NA") TYPN(I,STN) =  47  !  "Adenine & Guanine N9" 7 14.007 2
        IF( TYP(I,STN).EQ."H5") TYPN(I,STN) =  53  !  "Adenine HC2" 1 1.008 1  
        IF( TYP(I,STN).EQ."N2") TYPN(I,STN) =  45  !  "Adenine NH2-" 7 14.007 3  
        IF( TYP(I,STN).EQ."H") TYPN(I,STN) =  34  !  "Adenine NH2- (N1)" 1 1.008 1 
        IF( TYP(I,STN).EQ."H") TYPN(I,STN) =  34  !  "Adenine NH2- (C5)" 1 1.008 1 
        IF( TYP(I,STN).EQ."H5") TYPN(I,STN) =  53  !  "Adenine & Guanine HC8" 1 1.008 1
        IF( TYP(I,STN).EQ."H") TYPN(I,STN) =  34  !  "Adenine & Guanine HN9" 1 1.008 1
        IF( TYP(I,STN).EQ."NA") TYPN(I,STN) =  47  !  "Guanine N1" 7 14.007 3  
        IF( TYP(I,STN).EQ."CA") TYPN(I,STN) =  38  !  "Guanine C2" 6 12.011 3  
        IF( TYP(I,STN).EQ."NC") TYPN(I,STN) =  46  !  "Guanine N3" 7 14.007 2  
        IF( TYP(I,STN).EQ."CB") TYPN(I,STN) =  50  !  "Guanine C4" 6 12.011 3  
        IF( TYP(I,STN).EQ."CB") TYPN(I,STN) =  50  !  "Guanine C5" 6 12.011 3  
        IF( TYP(I,STN).EQ."C") TYPN(I,STN) =  3  !  "Guanine C6" 6 12.011 3  
        IF( TYP(I,STN).EQ."H") TYPN(I,STN) =  34  !  "Guanine HN1" 1 1.008 1  
        IF( TYP(I,STN).EQ."N2") TYPN(I,STN) =  45  !  "Guanine NH2-" 7 14.007 3  
        IF( TYP(I,STN).EQ."H") TYPN(I,STN) =  34  !  "Guanine NH2-" 1 1.008 1  
        IF( TYP(I,STN).EQ."O") TYPN(I,STN) =  4  !  "Guanine O6" 8 15.999 1  
        IF( TYP(I,STN).EQ."CT") TYPN(I,STN) =  13  !  "9-Me A & 9-Me-G CH3-" 6 12.011
        IF( TYP(I,STN).EQ."HC") TYPN(I,STN) =  36  !  "9-Me-A & 9-Me-G CH3-" 1 1.008 1
        IF( TYP(I,STN).EQ."CT") TYPN(I,STN) =  13  !  "1-Me-U & 1-Me-T CH3-" 6 12.011 4
        IF( TYP(I,STN).EQ."HC") TYPN(I,STN) =  36  !  "1-Me-U & 1-Me-T CH3-" 1 1.008 1
        IF( TYP(I,STN).EQ."CT") TYPN(I,STN) =  13  !  "1-Me-Cytosine CH3-" 6 12.011 4  
        IF( TYP(I,STN).EQ."HC") TYPN(I,STN) =  36  !  "1-Me-Cytosine CH3-" 1 1.008 1  
        IF( TYP(I,STN).EQ."NA") TYPN(I,STN) =  47  !  "CytosineH+ N1" 7 14.007 3  
        IF( TYP(I,STN).EQ."C") TYPN(I,STN) =  3  !  "CytosineH+ C2" 6 12.011 3  
        IF( TYP(I,STN).EQ."NA") TYPN(I,STN) =  47  !  "CytosineH+ N3" 7 14.007 3  
        IF( TYP(I,STN).EQ."CA") TYPN(I,STN) =  38  !  "CytosineH+ C4" 6 12.011 3  
        IF( TYP(I,STN).EQ."CM") TYPN(I,STN) =  37  !  "CytosineH+ C5" 6 12.011 3  
        IF( TYP(I,STN).EQ."CM") TYPN(I,STN) =  37  !  "CytosineH+ C6" 6 12.011 3  
        IF( TYP(I,STN).EQ."H") TYPN(I,STN) =  34  !  "CytosineH+ HN1" 1 1.008 1  
        IF( TYP(I,STN).EQ."O") TYPN(I,STN) =  4  !  "CytosineH+ O2" 8 15.999 1  
        IF( TYP(I,STN).EQ."H") TYPN(I,STN) =  34  !  "CytosineH+ HN3" 1 1.008 1  
        IF( TYP(I,STN).EQ."N2") TYPN(I,STN) =  45  !  "CytosineH+ NH2-" 7 14.007 3  
        IF( TYP(I,STN).EQ."H") TYPN(I,STN) =  34  !  "CytosineH+ NH2- (N3)" 1 1.008 1 
        IF( TYP(I,STN).EQ."H") TYPN(I,STN) =  34  !  "CytosineH+ NH2- (C5)" 1 1.008 1 
        IF( TYP(I,STN).EQ."HA") TYPN(I,STN) =  39  !  "CytosineH+ HC5" 1 1.008 1  
        IF( TYP(I,STN).EQ."H4") TYPN(I,STN) =  48  !  "CytosineH+ HC6" 1 1.008 1  
        IF( TYP(I,STN).EQ."CT") TYPN(I,STN) =  13  !  "1-Me-CytosineH+ CH3-" 6 12.011 4  
        IF( TYP(I,STN).EQ."HC") TYPN(I,STN) =  36  !  "1-Me-CytosineH+ CH3-" 1 1.008 1  
        IF( TYP(I,STN).EQ."P") TYPN(I,STN) =  54  !  "DiMePhosphate P (UA)" 15 30.974 4 
        IF( TYP(I,STN).EQ."O2") TYPN(I,STN) =  42  !  "DiMePhosphate O=P-O (UA)" 8 15.999 1 
        IF( TYP(I,STN).EQ."OS") TYPN(I,STN) =  29  !  "DiMePhosphate CH3-O (UA)" 8 15.999 2 
        IF( TYP(I,STN).EQ."CT") TYPN(I,STN) =  13  !  "DiMePhosphate CH3-O (UA)" 6 15.035 1 
        IF( TYP(I,STN).EQ."CM") TYPN(I,STN) =  37  !  "Trifluorothymine CF3-" 6 12.011 4  
        IF( TYP(I,STN).EQ."Cl") TYPN(I,STN) =  30  !  "Chloroalkene Cl2-C=" 17 35.453 1  
        IF( TYP(I,STN).EQ."CM") TYPN(I,STN) =  37  !  "Chloroalkene Cl2-C=" 6 12.011 3  
        IF( TYP(I,STN).EQ."F") TYPN(I,STN) =  1  !  "Fluoride Ion F-" 9 18.998 0 
        IF( TYP(I,STN).EQ."Cl") TYPN(I,STN) =  30  !  "Chloride Ion Cl-" 17 35.453 0 
        IF( TYP(I,STN).EQ."Br") TYPN(I,STN) =  55  !  "Bromide Ion Br-" 35 79.904 0 
        IF( TYP(I,STN).EQ."I") TYPN(I,STN) =  56  !  "Iodide Ion I-" 53 126.905 0 
        IF( TYP(I,STN).EQ."N4") TYPN(I,STN) =  57  !  "Ammonium Ion NH4+ (UA)" 7 18.039 0
        IF( TYP(I,STN).EQ."Li") TYPN(I,STN) =  58  !  "Lithium Ion Li+" 3 6.941 0 
        IF( TYP(I,STN).EQ."Na") TYPN(I,STN) =  59  !  "Sodium Ion Na+" 11 22.990 0 
        IF( TYP(I,STN).EQ."K") TYPN(I,STN) =  60  !  "Potassium Ion K+" 19 39.098 0 
        IF( TYP(I,STN).EQ."Rb") TYPN(I,STN) =  61  !  "Rubidium Ion Rb+" 37 85.468 0 
        IF( TYP(I,STN).EQ."Cs") TYPN(I,STN) =  62  !  "Cesium Ion Cs+" 55 132.905 0 
        IF( TYP(I,STN).EQ."Mg") TYPN(I,STN) =  63  !  "Magnesium Ion Mg+2" 12 24.305 0 
        IF( TYP(I,STN).EQ."Ca") TYPN(I,STN) =  64  !  "Calcium Ion Ca+2" 20 40.080 0 
        IF( TYP(I,STN).EQ."Sr") TYPN(I,STN) =  65  !  "Strontium Ion Sr+2" 38 87.620 0 
        IF( TYP(I,STN).EQ."Ba") TYPN(I,STN) =  66  !  "Barium Ion Ba+2" 56 137.330 0 
        IF( TYP(I,STN).EQ."C3") TYPN(I,STN) =  6  !  "Methyl Thiolate CH3S-" 6 12.011 4 
        IF( TYP(I,STN).EQ."HC") TYPN(I,STN) =  36  !  "Methyl Thiolate CH3S-" 1 1.008 1 
        IF( TYP(I,STN).EQ."SH") TYPN(I,STN) =  15  !  "Methyl Thiolate CH3S-" 16 32.060 1 
        IF( TYP(I,STN).EQ."C3") TYPN(I,STN) =  6  !  "Methoxide CH3O-" 6 12.011 4  
        IF( TYP(I,STN).EQ."HC") TYPN(I,STN) =  36  !  "Methoxide CH3O-" 1 1.008 1  
        IF( TYP(I,STN).EQ."OH") TYPN(I,STN) =  5  !  "Methoxide CH3O-" 8 15.999 1  
        IF( TYP(I,STN).EQ."CT") TYPN(I,STN) =  13  !  "Nitrile Anion CNCH2-" 6 12.011 3 
        IF( TYP(I,STN).EQ."HC") TYPN(I,STN) =  36  !  "Nitrile Anion CNCH2-" 1 1.008 1 
        IF( TYP(I,STN).EQ."CZ") TYPN(I,STN) =  19  !  "Nitrile Anion CNCH2-" 6 12.011 2 
        IF( TYP(I,STN).EQ."NZ") TYPN(I,STN) =  18  !  "Nitrile Anion CNCH2-" 7 14.007 1 
        IF( TYP(I,STN).EQ."C3") TYPN(I,STN) =  6  !  "Me Amine Anion CH3NH-" 6 12.011 4
        IF( TYP(I,STN).EQ."HC") TYPN(I,STN) =  36  !  "Me Amine Anion CH3NH-" 1 1.008 1
        IF( TYP(I,STN).EQ."N3") TYPN(I,STN) =  43  !  "Me Amine Anion CH3NH-" 7 14.007 2
        IF( TYP(I,STN).EQ."H") TYPN(I,STN) =  34  !  "Methyl Amine Anion" 1 1.008 1 
        IF( TYP(I,STN).EQ."C3") TYPN(I,STN) =  6  !  "Ethyl Anion CH3-CH2-" 6 12.011 4 
        IF( TYP(I,STN).EQ."HC") TYPN(I,STN) =  36  !  "Ethyl Anion CH3-CH2-" 1 1.008 1 
        IF( TYP(I,STN).EQ."CT") TYPN(I,STN) =  13  !  "Ethyl Anion CH3-CH2-" 6 12.011 4 
        IF( TYP(I,STN).EQ."HC") TYPN(I,STN) =  36  !  "Ethyl Anion CH3-CH2-" 1 1.008 1 
        IF( TYP(I,STN).EQ."LP") TYPN(I,STN) =  22  !  "Ethyl Anion CH3-CH2-" 99 0.000 1 
        IF( TYP(I,STN).EQ."OH") TYPN(I,STN) =  5  !  "Hydroxide Ion OH-" 8 15.999 1 
        IF( TYP(I,STN).EQ."HO") TYPN(I,STN) =  7  !  "Hydroxide Ion OH-" 1 1.008 1 
        IF( TYP(I,STN).EQ."U") TYPN(I,STN) =  67  !  "Uranyl Ion UO2+" 92 238.029 2 
        IF( TYP(I,STN).EQ."OU") TYPN(I,STN) =  68  !  "Uranyl Ion UO2+" 8 15.999 1 
        IF( TYP(I,STN).EQ."OS") TYPN(I,STN) =  29  !  "GTP O-(POn)2" 8 15.999 2  
        IF( TYP(I,STN).EQ."P") TYPN(I,STN) =  54  !  "DiMe Phosphate P" 15 30.974 4 
        IF( TYP(I,STN).EQ."O2") TYPN(I,STN) =  42  !  "DiMe Phosphate O=P-O" 8 15.999 1 
        IF( TYP(I,STN).EQ."OS") TYPN(I,STN) =  29  !  "DiMe Phosphate CH3-O" 8 15.999 1 
        IF( TYP(I,STN).EQ."CT") TYPN(I,STN) =  13  !  "DiMe Phosphate CH3-O" 6 12.011 4 
        IF( TYP(I,STN).EQ."HC") TYPN(I,STN) =  36  !  "DiMe Phosphate CH3-O" 1 1.008 1 
        IF( TYP(I,STN).EQ."P") TYPN(I,STN) =  54  !  "Me Phosphate P" 15 30.974 4 
        IF( TYP(I,STN).EQ."O2") TYPN(I,STN) =  42  !  "Me Phosphate O=PO2" 8 15.999 1 
        IF( TYP(I,STN).EQ."OS") TYPN(I,STN) =  29  !  "Me Phosphate CH3-O" 8 15.999 1 
        IF( TYP(I,STN).EQ."CT") TYPN(I,STN) =  13  !  "Me Phosphate CH3-O" 6 12.011 4 
        IF( TYP(I,STN).EQ."HC") TYPN(I,STN) =  36  !  "Me Phosphate CH3-O" 1 1.008 1 
        IF( TYP(I,STN).EQ."P") TYPN(I,STN) =  54  !  "Me MePhosphonate P" 15 30.974 4 
        IF( TYP(I,STN).EQ."O2") TYPN(I,STN) =  42  !  "Me MePhosphonate O=P-O" 8 15.999 1 
        IF( TYP(I,STN).EQ."OS") TYPN(I,STN) =  29  !  "Me MePhosphonate CH3-O" 8 15.999 1 
        IF( TYP(I,STN).EQ."CT") TYPN(I,STN) =  13  !  "Me MePhosphonate CH3-O" 6 12.011 4 
        IF( TYP(I,STN).EQ."HC") TYPN(I,STN) =  36  !  "Me MePhosphonate CH3-O" 1 1.008 1 
        IF( TYP(I,STN).EQ."CT") TYPN(I,STN) =  13  !  "Me MePhosphonate CH3-P" 6 12.011 4 
        IF( TYP(I,STN).EQ."HC") TYPN(I,STN) =  36  !  "Me MePhosphonate CH3-P" 1 1.008 1 
        IF( TYP(I,STN).EQ."CA") TYPN(I,STN) =  38  !  "Bz MePhosphonate Cipso" 6 12.011 3 
        IF( TYP(I,STN).EQ."CT") TYPN(I,STN) =  13  !  "Bz MePhosphonate CH3-O" 6 12.011 4 
        IF( TYP(I,STN).EQ."HC") TYPN(I,STN) =  36  !  "Bz MePhosphonate CH3-O" 1 1.008 1 
        IF( TYP(I,STN).EQ."CA") TYPN(I,STN) =  38  !  "Me BzPhosphonate Cipso" 6 12.011 3 
        IF( TYP(I,STN).EQ."CT") TYPN(I,STN) =  13  !  "Me BzPhosphonate CH3-P" 6 12.011 4 
        IF( TYP(I,STN).EQ."HC") TYPN(I,STN) =  36  !  "Me BzPhosphonate CH3-P" 1 1.008 1 
        IF( TYP(I,STN).EQ."CA") TYPN(I,STN) =  38  !  "Ph Phosphate Cipso" 6 12.011 3 
        IF( TYP(I,STN).EQ."CT") TYPN(I,STN) =  13  !  "Barbiturate C6(R2)" 6 12.011 4  
        IF( TYP(I,STN).EQ."C") TYPN(I,STN) =  3  !  "Ester -COOR" 6 12.011 3  
        IF( TYP(I,STN).EQ."O") TYPN(I,STN) =  4  !  "Ester C=O" 8 15.999 1  
        IF( TYP(I,STN).EQ."OS") TYPN(I,STN) =  29  !  "Ester CO-O-R" 8 15.999 2  
        IF( TYP(I,STN).EQ."CT") TYPN(I,STN) =  13  !  "Methyl Ester -OCH3" 6 12.011 4 
        IF( TYP(I,STN).EQ."HC") TYPN(I,STN) =  36  !  "Methyl Ester -OCH3" 1 1.008 1 
        IF( TYP(I,STN).EQ."C") TYPN(I,STN) =  3  !  "Benzoic Acid -COOH" 6 12.011 3 
        IF( TYP(I,STN).EQ."C") TYPN(I,STN) =  3  !  "Aryl Ester -COOR" 6 12.011 3 
        IF( TYP(I,STN).EQ."CA") TYPN(I,STN) =  38  !  "Phenyl Ester Cipso" 6 12.011 3 
        IF( TYP(I,STN).EQ."OS") TYPN(I,STN) =  29  !  "Phenyl Ester -OPh" 8 15.999 2 
        IF( TYP(I,STN).EQ."SY") TYPN(I,STN) =  69  !  "Sulfonamide -SO2N<" 16 32.060 4  
        IF( TYP(I,STN).EQ."OY") TYPN(I,STN) =  32  !  "Sulfonamide -SO2N<" 8 15.999 1  
        IF( TYP(I,STN).EQ."CT") TYPN(I,STN) =  13  !  "Sulfonamide CH3-S" 6 12.011 4  
        IF( TYP(I,STN).EQ."HC") TYPN(I,STN) =  36  !  "Sulfonamide CH3-S" 1 1.008 1  
        IF( TYP(I,STN).EQ."N") TYPN(I,STN) =  35  !  "Sulfonamide -SO2NH2" 7 14.007 3  
        IF( TYP(I,STN).EQ."H") TYPN(I,STN) =  34  !  "Sulfonamide -SO2NH2" 1 1.008 1  
        IF( TYP(I,STN).EQ."N") TYPN(I,STN) =  35  !  "Sulfonamide -SO2NHR" 7 14.007 3  
        IF( TYP(I,STN).EQ."H") TYPN(I,STN) =  34  !  "Sulfonamide -SO2NHR" 1 1.008 1  
        IF( TYP(I,STN).EQ."CT") TYPN(I,STN) =  13  !  "N-Me Sulfonamide CH3-" 6 12.011 4 
        IF( TYP(I,STN).EQ."HC") TYPN(I,STN) =  36  !  "N-Me Sulfonamide CH3-" 1 1.008 1 
        IF( TYP(I,STN).EQ."CT") TYPN(I,STN) =  13  !  "Sulfonamide N-CH2-R" 6 12.011 4  
        IF( TYP(I,STN).EQ."HC") TYPN(I,STN) =  36  !  "Sulfonamide N-CH2-R" 1 1.008 1  
        IF( TYP(I,STN).EQ."CT") TYPN(I,STN) =  13  !  "N-Et Sulfonamide CH3-" 6 12.011 4 
        IF( TYP(I,STN).EQ."HC") TYPN(I,STN) =  36  !  "N-Et Sulfonamide CH3-" 1 1.008 1 
        IF( TYP(I,STN).EQ."CA") TYPN(I,STN) =  38  !  "Aryl Sulfonamide C-SO2N" 6 12.011 3 
        IF( TYP(I,STN).EQ."CA") TYPN(I,STN) =  38  !  "Aryl Sulfoxide C-S=O" 6 12.011 3 
        IF( TYP(I,STN).EQ."CT") TYPN(I,STN) =  13  !  "Et Ester -OCH2R" 6 12.011 4 
        IF( TYP(I,STN).EQ."CT") TYPN(I,STN) =  13  !  "i-Pr Ester -OCHR2" 6 12.011 4 
        IF( TYP(I,STN).EQ."CT") TYPN(I,STN) =  13  !  "t-Bu Ester -OCR3" 6 12.011 4 
        IF( TYP(I,STN).EQ."SY") TYPN(I,STN) =  69  !  "Sulfone R-SO2-R" 16 32.060 4  
        IF( TYP(I,STN).EQ."OY") TYPN(I,STN) =  32  !  "Sulfone R-SO2-R" 8 15.999 1  
        IF( TYP(I,STN).EQ."SZ") TYPN(I,STN) =  31  !  "Alkyl Aryl Sulfoxide" 16 32.060 3 
        IF( TYP(I,STN).EQ."SZ") TYPN(I,STN) =  31  !  "Dialkyl Sulfoxide" 16 32.060 3  
        IF( TYP(I,STN).EQ."OY") TYPN(I,STN) =  32  !  "Sulfoxide R-SO-R" 8 15.999 1  
        IF( TYP(I,STN).EQ."CT") TYPN(I,STN) =  13  !  "Sulfoxide CH3-SO-R" 6 12.011 4  
        IF( TYP(I,STN).EQ."CT") TYPN(I,STN) =  13  !  "Sulfoxide -CH2-SO-R" 6 12.011 4  
        IF( TYP(I,STN).EQ."C*") TYPN(I,STN) =  70  !  "TRP CG" 6 12.011 3  
        IF( TYP(I,STN).EQ."CB") TYPN(I,STN) =  50  !  "TRP CD" 6 12.011 3  
        IF( TYP(I,STN).EQ."CN") TYPN(I,STN) =  71  !  "TRP CE" 6 12.011 3  
        IF( TYP(I,STN).EQ."NA") TYPN(I,STN) =  47  !  "TRP NE, HID ND & HIE NE"
        IF( TYP(I,STN).EQ."H") TYPN(I,STN) =  34  !  "TRP HNE & HID/HIE HN" 1 1.008
        IF( TYP(I,STN).EQ."CT") TYPN(I,STN) =  13  !  "HIS CB" 6 12.011 4  
        IF( TYP(I,STN).EQ."CR") TYPN(I,STN) =  72  !  "HID & HIE CE1" 6 12.011 3
        IF( TYP(I,STN).EQ."CV") TYPN(I,STN) =  73  !  "HID CD2 & HIE CG" 6 12.011
        IF( TYP(I,STN).EQ."CW") TYPN(I,STN) =  74  !  "HID CG & HIE CD2" 6 12.011
        IF( TYP(I,STN).EQ."CR") TYPN(I,STN) =  72  !  "HIP CE1" 6 12.011 3  
        IF( TYP(I,STN).EQ."CX") TYPN(I,STN) =  75  !  "HIP CG & CD2" 6 12.011 3
        IF( TYP(I,STN).EQ."NB") TYPN(I,STN) =  51  !  "HID NE & HIE ND" 7 14.007
        IF( TYP(I,STN).EQ."NA") TYPN(I,STN) =  47  !  "HIP ND & NE" 7 14.007 3
        IF( TYP(I,STN).EQ."H") TYPN(I,STN) =  34  !  "HIP HND & HNE" 1 1.008 1
        IF( TYP(I,STN).EQ."CW") TYPN(I,STN) =  74  !  "TRP CD1" 6 12.011 3  
        IF( TYP(I,STN).EQ."CT") TYPN(I,STN) =  13  !  "i-Pr Benzene -CHMe2" 6 12.011 4 
        IF( TYP(I,STN).EQ."CT") TYPN(I,STN) =  13  !  "t-Bu Benzene -CMe3" 6 12.011 4 
        IF( TYP(I,STN).EQ."CM") TYPN(I,STN) =  37  !  "Vinyl Ether =CH-OR" 6 12.011 3 
        IF( TYP(I,STN).EQ."CM") TYPN(I,STN) =  37  !  "Vinyl Ether =CR-OR" 6 12.011 3 
        IF( TYP(I,STN).EQ."C!") TYPN(I,STN) =  76  !  "Biphenyl C1" 6 12.011 3  
        IF( TYP(I,STN).EQ."NC") TYPN(I,STN) =  46  !  "Pyridine N" 7 14.007 2  
        IF( TYP(I,STN).EQ."CA") TYPN(I,STN) =  38  !  "Pyridine C1" 6 12.011 3  
        IF( TYP(I,STN).EQ."CA") TYPN(I,STN) =  38  !  "Pyridine C2" 6 12.011 3  
        IF( TYP(I,STN).EQ."CA") TYPN(I,STN) =  38  !  "Pyridine C3" 6 12.011 3  
        IF( TYP(I,STN).EQ."HA") TYPN(I,STN) =  39  !  "Pyridine H1" 1 1.008 1  
        IF( TYP(I,STN).EQ."HA") TYPN(I,STN) =  39  !  "Pyridine H2" 1 1.008 1  
        IF( TYP(I,STN).EQ."HA") TYPN(I,STN) =  39  !  "Pyridine H3" 1 1.008 1  
        IF( TYP(I,STN).EQ."NC") TYPN(I,STN) =  46  !  "Pyrazine N" 7 14.007 2  
        IF( TYP(I,STN).EQ."CA") TYPN(I,STN) =  38  !  "Pyrazine CH" 6 12.011 3  
        IF( TYP(I,STN).EQ."HA") TYPN(I,STN) =  39  !  "Pyrazine CH" 1 1.008 1  
        IF( TYP(I,STN).EQ."NC") TYPN(I,STN) =  46  !  "Pyrimidine N" 7 14.007 2  
        IF( TYP(I,STN).EQ."CQ") TYPN(I,STN) =  49  !  "Pyrimidine C2" 6 12.011 3  
        IF( TYP(I,STN).EQ."CA") TYPN(I,STN) =  38  !  "Pyrimidine C4" 6 12.011 3  
        IF( TYP(I,STN).EQ."CA") TYPN(I,STN) =  38  !  "Pyrimidine C5" 6 12.011 3  
        IF( TYP(I,STN).EQ."HA") TYPN(I,STN) =  39  !  "Pyrimidine HC2" 1 1.008 1  
        IF( TYP(I,STN).EQ."HA") TYPN(I,STN) =  39  !  "Pyrimidine HC4" 1 1.008 1  
        IF( TYP(I,STN).EQ."HA") TYPN(I,STN) =  39  !  "Pyrimidine HC5" 1 1.008 1  
        IF( TYP(I,STN).EQ."NC") TYPN(I,STN) =  46  !  "Pyridazine N" 7 14.007 2  
        IF( TYP(I,STN).EQ."CA") TYPN(I,STN) =  38  !  "Pyridazine C3" 6 12.011 3  
        IF( TYP(I,STN).EQ."CA") TYPN(I,STN) =  38  !  "Pyridazine C4" 6 12.011 3  
        IF( TYP(I,STN).EQ."HA") TYPN(I,STN) =  39  !  "Pyridazine HC3" 1 1.008 1  
        IF( TYP(I,STN).EQ."HA") TYPN(I,STN) =  39  !  "Pyridazine HC4" 1 1.008 1  
        IF( TYP(I,STN).EQ."NA") TYPN(I,STN) =  47  !  "Pyrrole N" 7 14.007 3  
        IF( TYP(I,STN).EQ."CW") TYPN(I,STN) =  74  !  "Pyrrole C2" 6 12.011 3  
        IF( TYP(I,STN).EQ."CS") TYPN(I,STN) =  77  !  "Pyrrole C3" 6 12.011 3  
        IF( TYP(I,STN).EQ."H") TYPN(I,STN) =  34  !  "Pyrrole HN" 1 1.008 1  
        IF( TYP(I,STN).EQ."HA") TYPN(I,STN) =  39  !  "Pyrrole HC2" 1 1.008 1  
        IF( TYP(I,STN).EQ."HA") TYPN(I,STN) =  39  !  "Pyrrole HC3" 1 1.008 1  
        IF( TYP(I,STN).EQ."NA") TYPN(I,STN) =  47  !  "Pyrazole N1" 7 14.007 3  
        IF( TYP(I,STN).EQ."NB") TYPN(I,STN) =  51  !  "Pyrazole N2" 7 14.007 2  
        IF( TYP(I,STN).EQ."CU") TYPN(I,STN) =  78  !  "Pyrazole C3" 6 12.011 3  
        IF( TYP(I,STN).EQ."CS") TYPN(I,STN) =  77  !  "Pyrazole C4" 6 12.011 3  
        IF( TYP(I,STN).EQ."CW") TYPN(I,STN) =  74  !  "Pyrazole C5" 6 12.011 3  
        IF( TYP(I,STN).EQ."H") TYPN(I,STN) =  34  !  "Pyrazole HN1" 1 1.008 1  
        IF( TYP(I,STN).EQ."HA") TYPN(I,STN) =  39  !  "Pyrazole HC3" 1 1.008 1  
        IF( TYP(I,STN).EQ."HA") TYPN(I,STN) =  39  !  "Pyrazole HC4" 1 1.008 1  
        IF( TYP(I,STN).EQ."HA") TYPN(I,STN) =  39  !  "Pyrazole HC5" 1 1.008 1  
        IF( TYP(I,STN).EQ."NA") TYPN(I,STN) =  47  !  "Imidazole N1" 7 14.007 3  
        IF( TYP(I,STN).EQ."CR") TYPN(I,STN) =  72  !  "Imidazole C2" 6 12.011 3  
        IF( TYP(I,STN).EQ."NB") TYPN(I,STN) =  51  !  "Imidazole N3" 7 14.007 2  
        IF( TYP(I,STN).EQ."CV") TYPN(I,STN) =  73  !  "Imidazole C4" 6 12.011 3  
        IF( TYP(I,STN).EQ."CW") TYPN(I,STN) =  74  !  "Imidazole C5" 6 12.011 3  
        IF( TYP(I,STN).EQ."H") TYPN(I,STN) =  34  !  "Imidazole HN1" 1 1.008 1  
        IF( TYP(I,STN).EQ."HA") TYPN(I,STN) =  39  !  "Imidazole HC2" 1 1.008 1  
        IF( TYP(I,STN).EQ."HA") TYPN(I,STN) =  39  !  "Imidazole HC4" 1 1.008 1  
        IF( TYP(I,STN).EQ."HA") TYPN(I,STN) =  39  !  "Imidazole HC5" 1 1.008 1  
        IF( TYP(I,STN).EQ."OS") TYPN(I,STN) =  29  !  "Furan O" 8 15.999 2  
        IF( TYP(I,STN).EQ."CW") TYPN(I,STN) =  74  !  "Furan C2" 6 12.011 3  
        IF( TYP(I,STN).EQ."CS") TYPN(I,STN) =  77  !  "Furan C3" 6 12.011 3  
        IF( TYP(I,STN).EQ."HA") TYPN(I,STN) =  39  !  "Furan HC2" 1 1.008 1  
        IF( TYP(I,STN).EQ."HA") TYPN(I,STN) =  39  !  "Furan HC3" 1 1.008 1  
        IF( TYP(I,STN).EQ."OS") TYPN(I,STN) =  29  !  "Oxazole O" 8 15.999 2  
        IF( TYP(I,STN).EQ."CR") TYPN(I,STN) =  72  !  "Oxazole C2" 6 12.011 3  
        IF( TYP(I,STN).EQ."NB") TYPN(I,STN) =  51  !  "Oxazole N" 7 14.007 2  
        IF( TYP(I,STN).EQ."CV") TYPN(I,STN) =  73  !  "Oxazole C4" 6 12.011 3  
        IF( TYP(I,STN).EQ."CW") TYPN(I,STN) =  74  !  "Oxazole C5" 6 12.011 3  
        IF( TYP(I,STN).EQ."HA") TYPN(I,STN) =  39  !  "Oxazole HC2" 1 1.008 1  
        IF( TYP(I,STN).EQ."HA") TYPN(I,STN) =  39  !  "Oxazole HC4" 1 1.008 1  
        IF( TYP(I,STN).EQ."HA") TYPN(I,STN) =  39  !  "Oxazole HC5" 1 1.008 1  
        IF( TYP(I,STN).EQ."OS") TYPN(I,STN) =  29  !  "Isoxazole O" 8 15.999 2  
        IF( TYP(I,STN).EQ."NB") TYPN(I,STN) =  51  !  "Isoxazole N" 7 14.007 2  
        IF( TYP(I,STN).EQ."CU") TYPN(I,STN) =  78  !  "Isoxazole C3" 6 12.011 3  
        IF( TYP(I,STN).EQ."CS") TYPN(I,STN) =  77  !  "Isoxazole C4" 6 12.011 3  
        IF( TYP(I,STN).EQ."CW") TYPN(I,STN) =  74  !  "Isoxazole C5" 6 12.011 3  
        IF( TYP(I,STN).EQ."HA") TYPN(I,STN) =  39  !  "Isoxazole HC3" 1 1.008 1  
        IF( TYP(I,STN).EQ."HA") TYPN(I,STN) =  39  !  "Isoxazole HC4" 1 1.008 1  
        IF( TYP(I,STN).EQ."HA") TYPN(I,STN) =  39  !  "Isoxazole HC5" 1 1.008 1  
        IF( TYP(I,STN).EQ."NA") TYPN(I,STN) =  47  !  "Indole N1" 7 14.007 3  
        IF( TYP(I,STN).EQ."CW") TYPN(I,STN) =  74  !  "Indole C2" 6 12.011 3  
        IF( TYP(I,STN).EQ."CS") TYPN(I,STN) =  77  !  "Indole C3" 6 12.011 3  
        IF( TYP(I,STN).EQ."CA") TYPN(I,STN) =  38  !  "Indole C4" 6 12.011 3  
        IF( TYP(I,STN).EQ."CA") TYPN(I,STN) =  38  !  "Indole C5" 6 12.011 3  
        IF( TYP(I,STN).EQ."CA") TYPN(I,STN) =  38  !  "Indole C6" 6 12.011 3  
        IF( TYP(I,STN).EQ."CA") TYPN(I,STN) =  38  !  "Indole C7" 6 12.011 3  
        IF( TYP(I,STN).EQ."CN") TYPN(I,STN) =  71  !  "Indole C8" 6 12.011 3  
        IF( TYP(I,STN).EQ."CB") TYPN(I,STN) =  50  !  "Indole C9" 6 12.011 3  
        IF( TYP(I,STN).EQ."H") TYPN(I,STN) =  34  !  "Indole HN1" 1 1.008 1  
        IF( TYP(I,STN).EQ."HA") TYPN(I,STN) =  39  !  "Indole HC2" 1 1.008 1  
        IF( TYP(I,STN).EQ."HA") TYPN(I,STN) =  39  !  "Indole HC3" 1 1.008 1  
        IF( TYP(I,STN).EQ."HA") TYPN(I,STN) =  39  !  "Indole HC4" 1 1.008 1  
        IF( TYP(I,STN).EQ."HA") TYPN(I,STN) =  39  !  "Indole HC5" 1 1.008 1  
        IF( TYP(I,STN).EQ."HA") TYPN(I,STN) =  39  !  "Indole HC6" 1 1.008 1  
        IF( TYP(I,STN).EQ."HA") TYPN(I,STN) =  39  !  "Indole HC7" 1 1.008 1  
        IF( TYP(I,STN).EQ."NC") TYPN(I,STN) =  46  !  "Quinoline N1" 7 14.007 2  
        IF( TYP(I,STN).EQ."CA") TYPN(I,STN) =  38  !  "Quinoline C2" 6 12.011 3  
        IF( TYP(I,STN).EQ."CA") TYPN(I,STN) =  38  !  "Quinoline C3" 6 12.011 3  
        IF( TYP(I,STN).EQ."CA") TYPN(I,STN) =  38  !  "Quinoline C4" 6 12.011 3  
        IF( TYP(I,STN).EQ."CA") TYPN(I,STN) =  38  !  "Quinoline C5" 6 12.011 3  
        IF( TYP(I,STN).EQ."CA") TYPN(I,STN) =  38  !  "Quinoline C6" 6 12.011 3  
        IF( TYP(I,STN).EQ."CA") TYPN(I,STN) =  38  !  "Quinoline C7" 6 12.011 3  
        IF( TYP(I,STN).EQ."CA") TYPN(I,STN) =  38  !  "Quinoline C8" 6 12.011 3  
        IF( TYP(I,STN).EQ."CA") TYPN(I,STN) =  38  !  "Quinoline C9" 6 12.011 3  
        IF( TYP(I,STN).EQ."CA") TYPN(I,STN) =  38  !  "Quinoline C10" 6 12.011 3  
        IF( TYP(I,STN).EQ."HA") TYPN(I,STN) =  39  !  "Quinoline HC2" 1 1.008 1  
        IF( TYP(I,STN).EQ."HA") TYPN(I,STN) =  39  !  "Quinoline HC3" 1 1.008 1  
        IF( TYP(I,STN).EQ."HA") TYPN(I,STN) =  39  !  "Quinoline HC4" 1 1.008 1  
        IF( TYP(I,STN).EQ."HA") TYPN(I,STN) =  39  !  "Quinoline HC5" 1 1.008 1  
        IF( TYP(I,STN).EQ."HA") TYPN(I,STN) =  39  !  "Quinoline HC6" 1 1.008 1  
        IF( TYP(I,STN).EQ."HA") TYPN(I,STN) =  39  !  "Quinoline HC7" 1 1.008 1  
        IF( TYP(I,STN).EQ."HA") TYPN(I,STN) =  39  !  "Quinoline HC8" 1 1.008 1  
        IF( TYP(I,STN).EQ."NC") TYPN(I,STN) =  46  !  "Purine N1" 7 14.007 2  
        IF( TYP(I,STN).EQ."CQ") TYPN(I,STN) =  49  !  "Purine C2" 6 12.011 3  
        IF( TYP(I,STN).EQ."NC") TYPN(I,STN) =  46  !  "Purine N3" 7 14.007 2  
        IF( TYP(I,STN).EQ."CB") TYPN(I,STN) =  50  !  "Purine C4" 6 12.011 3  
        IF( TYP(I,STN).EQ."CB") TYPN(I,STN) =  50  !  "Purine C5" 6 12.011 3  
        IF( TYP(I,STN).EQ."CA") TYPN(I,STN) =  38  !  "Purine C6" 6 12.011 3  
        IF( TYP(I,STN).EQ."NB") TYPN(I,STN) =  51  !  "Purine N7" 7 14.007 2  
        IF( TYP(I,STN).EQ."CK") TYPN(I,STN) =  52  !  "Purine C8" 6 12.011 3  
        IF( TYP(I,STN).EQ."NA") TYPN(I,STN) =  47  !  "Purine N9" 7 14.007 3  
        IF( TYP(I,STN).EQ."HA") TYPN(I,STN) =  39  !  "Purine HC2" 1 1.008 1  
        IF( TYP(I,STN).EQ."HA") TYPN(I,STN) =  39  !  "Purine HC6" 1 1.008 1  
        IF( TYP(I,STN).EQ."HA") TYPN(I,STN) =  39  !  "Purine HC8" 1 1.008 1  
        IF( TYP(I,STN).EQ."H") TYPN(I,STN) =  34  !  "Purine HN9" 1 1.008 1  
        IF( TYP(I,STN).EQ."S") TYPN(I,STN) =  16  !  "Thiazole S" 16 32.060 2  
        IF( TYP(I,STN).EQ."CR") TYPN(I,STN) =  72  !  "Thiazole C2" 6 12.011 3  
        IF( TYP(I,STN).EQ."NB") TYPN(I,STN) =  51  !  "Thiazole N" 7 14.007 2  
        IF( TYP(I,STN).EQ."CV") TYPN(I,STN) =  73  !  "Thiazole C4" 6 12.011 3  
        IF( TYP(I,STN).EQ."CW") TYPN(I,STN) =  74  !  "Thiazole C5" 6 12.011 3  
        IF( TYP(I,STN).EQ."HA") TYPN(I,STN) =  39  !  "Thiazole HC2" 1 1.008 1  
        IF( TYP(I,STN).EQ."HA") TYPN(I,STN) =  39  !  "Thiazole HC4" 1 1.008 1  
        IF( TYP(I,STN).EQ."HA") TYPN(I,STN) =  39  !  "Thiazole HC5" 1 1.008 1  
        IF( TYP(I,STN).EQ."NC") TYPN(I,STN) =  46  !  "1,3,5-Triazine N" 7 14.007 2  
        IF( TYP(I,STN).EQ."CQ") TYPN(I,STN) =  49  !  "1,3,5-Triazine CH" 6 12.011 3  
        IF( TYP(I,STN).EQ."HA") TYPN(I,STN) =  39  !  "1,3,5-Triazine CH" 1 1.008 1  
        IF( TYP(I,STN).EQ."CA") TYPN(I,STN) =  38  !  "Serotonin C5-OH" 6 12.011 3  
        IF( TYP(I,STN).EQ."CT") TYPN(I,STN) =  13  !  "Serotonin CH2 on C3" 6 12.011 4
        IF( TYP(I,STN).EQ."NC") TYPN(I,STN) =  46  !  "1,10-Phenanthroline N" 7 14.007 2  
        IF( TYP(I,STN).EQ."CA") TYPN(I,STN) =  38  !  "1,10-Phenanthroline C2" 6 12.011 3  
        IF( TYP(I,STN).EQ."CA") TYPN(I,STN) =  38  !  "1,10-Phenanthroline C3" 6 12.011 3  
        IF( TYP(I,STN).EQ."CA") TYPN(I,STN) =  38  !  "1,10-Phenanthroline C4" 6 12.011 3  
        IF( TYP(I,STN).EQ."CA") TYPN(I,STN) =  38  !  "1,10-Phenanthroline C12" 6 12.011 3  
        IF( TYP(I,STN).EQ."CA") TYPN(I,STN) =  38  !  "1,10-Phenanthroline C11" 6 12.011 3  
        IF( TYP(I,STN).EQ."CA") TYPN(I,STN) =  38  !  "1,10-Phenanthroline C5" 6 12.011 3  
        IF( TYP(I,STN).EQ."HA") TYPN(I,STN) =  39  !  "1,10-Phenanthroline HC2" 1 1.008 1  
        IF( TYP(I,STN).EQ."HA") TYPN(I,STN) =  39  !  "1,10-Phenanthroline HC3" 1 1.008 1  
        IF( TYP(I,STN).EQ."HA") TYPN(I,STN) =  39  !  "1,10-Phenanthroline HC4" 1 1.008 1  
        IF( TYP(I,STN).EQ."HA") TYPN(I,STN) =  39  !  "1,10-Phenanthroline HC5" 1 1.008 1  
        IF( TYP(I,STN).EQ."NA") TYPN(I,STN) =  47  !  "1-Methylimidazole N1" 7 14.007 3  
        IF( TYP(I,STN).EQ."CR") TYPN(I,STN) =  72  !  "1-Methylimidazole C2" 6 12.011 3  
        IF( TYP(I,STN).EQ."NB") TYPN(I,STN) =  51  !  "1-Methylimidazole N3" 7 14.007 2  
        IF( TYP(I,STN).EQ."CV") TYPN(I,STN) =  73  !  "1-Methylimidazole C4" 6 12.011 3  
        IF( TYP(I,STN).EQ."CW") TYPN(I,STN) =  74  !  "1-Methylimidazole C5" 6 12.011 3  
        IF( TYP(I,STN).EQ."CT") TYPN(I,STN) =  13  !  "1-Methylimidazole CH3-" 6 12.011 4  
        IF( TYP(I,STN).EQ."HA") TYPN(I,STN) =  39  !  "1-Methylimidazole HC2" 1 1.008 1  
        IF( TYP(I,STN).EQ."HA") TYPN(I,STN) =  39  !  "1-Methylimidazole HC4" 1 1.008 1  
        IF( TYP(I,STN).EQ."HA") TYPN(I,STN) =  39  !  "1-Methylimidazole HC5" 1 1.008 1  
        IF( TYP(I,STN).EQ."HC") TYPN(I,STN) =  36  !  "1-Methylimidazole CH3-" 1 1.008 1  
        IF( TYP(I,STN).EQ."CT") TYPN(I,STN) =  13  !  "1-Et Imidazole RCH2-" 6 12.011 4 
        IF( TYP(I,STN).EQ."CT") TYPN(I,STN) =  13  !  "1-iPr Imidazole R2CH-" 6 12.011 4 
        IF( TYP(I,STN).EQ."CT") TYPN(I,STN) =  13  !  "1-MeO-Me-Imidazole CH2" 6 12.011 4  
        IF( TYP(I,STN).EQ."CT") TYPN(I,STN) =  13  !  "2-Me Pyridine CH3" 6 12.011 4 
        IF( TYP(I,STN).EQ."CT") TYPN(I,STN) =  13  !  "2-Et Pyridine CH2" 6 12.011 4 
        IF( TYP(I,STN).EQ."CT") TYPN(I,STN) =  13  !  "3-Me Pyridazine CH3" 6 12.011 4 
        IF( TYP(I,STN).EQ."CT") TYPN(I,STN) =  13  !  "3-Et Pyridazine CH2" 6 12.011 4 
        IF( TYP(I,STN).EQ."CT") TYPN(I,STN) =  13  !  "4-Me Pyrimidine CH3" 6 12.011 4 
        IF( TYP(I,STN).EQ."CT") TYPN(I,STN) =  13  !  "4-Et Pyrimidine CH2" 6 12.011 4 
        IF( TYP(I,STN).EQ."CT") TYPN(I,STN) =  13  !  "2-Me Pyrazine CH3" 6 12.011 4 
        IF( TYP(I,STN).EQ."CT") TYPN(I,STN) =  13  !  "2-Et Pyrazine CH2" 6 12.011 4 
        IF( TYP(I,STN).EQ."CT") TYPN(I,STN) =  13  !  "2-Me Pyrrole CH3" 6 12.011 4 
        IF( TYP(I,STN).EQ."CT") TYPN(I,STN) =  13  !  "2-Et Pyrrole CH2" 6 12.011 4 
        IF( TYP(I,STN).EQ."CT") TYPN(I,STN) =  13  !  "2-Me Furan CH3" 6 12.011 4 
        IF( TYP(I,STN).EQ."CT") TYPN(I,STN) =  13  !  "2-Et Furan CH2" 6 12.011 4 
        IF( TYP(I,STN).EQ."SH") TYPN(I,STN) =  15  !  "6-Mercaptopurine SH" 16 32.060 2  
        IF( TYP(I,STN).EQ."HS") TYPN(I,STN) =  17  !  "6-Mercaptopurine SH" 1 1.008 1  
        IF( TYP(I,STN).EQ."CA") TYPN(I,STN) =  38  !  "6-Mercaptopurine C6" 6 12.011 3  
        IF( TYP(I,STN).EQ."C$") TYPN(I,STN) =  79  !  "Beta-Lactam N-C=O" 6 12.011 3  
        IF( TYP(I,STN).EQ."N$") TYPN(I,STN) =  80  !  "Beta-Lactam N-C=O" 7 14.007 4  
        IF( TYP(I,STN).EQ."CY") TYPN(I,STN) =  81  !  "Penicillin CH-N" 6 12.011 4  
        IF( TYP(I,STN).EQ."CY") TYPN(I,STN) =  81  !  "Penicillin CH-CO" 6 12.011 4  
        IF( TYP(I,STN).EQ."CT") TYPN(I,STN) =  13  !  "3-Me Indole CH3" 6 12.011 4 
        IF( TYP(I,STN).EQ."C!") TYPN(I,STN) =  76  !  "2-Phenyl Pyridine C2" 6 12.011 3 
        IF( TYP(I,STN).EQ."C!") TYPN(I,STN) =  76  !  "2-Phenyl Pyridine C2'" 6 12.011 3 
        IF( TYP(I,STN).EQ."C!") TYPN(I,STN) =  76  !  "2-Phenyl Pyridine C3" 6 12.011 3 
        IF( TYP(I,STN).EQ."C!") TYPN(I,STN) =  76  !  "2-Phenyl Pyridine C3'" 6 12.011 3 
        IF( TYP(I,STN).EQ."C!") TYPN(I,STN) =  76  !  "2-Phenyl Pyridine C4" 6 12.011 3 
        IF( TYP(I,STN).EQ."C!") TYPN(I,STN) =  76  !  "2-Phenyl Pyridine C4'" 6 12.011 3 
        IF( TYP(I,STN).EQ."S") TYPN(I,STN) =  16  !  "Diphenyl Thioether S" 16 32.060 2 
        IF( TYP(I,STN).EQ."Ac") TYPN(I,STN) =  82  !  "Actinium Ion Ac+3" 89 227.000 0 
        IF( TYP(I,STN).EQ."Th") TYPN(I,STN) =  83  !  "Thorium Ion Th+4" 90 232.038 0 
        IF( TYP(I,STN).EQ."Am") TYPN(I,STN) =  84  !  "Americium Ion Am+3" 95 243.000 0 
        IF( TYP(I,STN).EQ."C+") TYPN(I,STN) =  85  !  "t-Butyl Cation C+" 6 12.011 3 
        IF( TYP(I,STN).EQ."CT") TYPN(I,STN) =  13  !  "t-Butyl Cation CH3-" 6 12.011 4 
        IF( TYP(I,STN).EQ."HC") TYPN(I,STN) =  36  !  "t-Butyl Cation CH3-" 1 1.008 1 
        IF( TYP(I,STN).EQ."La") TYPN(I,STN) =  86  !  "Lanthanum Ion La+3" 57 138.906 0 
        IF( TYP(I,STN).EQ."Nd") TYPN(I,STN) =  87  !  "Neodymium Ion Nd+3" 60 144.240 0 
        IF( TYP(I,STN).EQ."Eu") TYPN(I,STN) =  88  !  "Europium Ion Eu+3" 63 151.960 0 
        IF( TYP(I,STN).EQ."Gd") TYPN(I,STN) =  89  !  "Gadolinium Ion Gd+3" 64 157.250 0 
        IF( TYP(I,STN).EQ."Yb") TYPN(I,STN) =  90  !  "Ytterbium Ion Yb+3" 70 173.040 0 
        IF( TYP(I,STN).EQ."CM") TYPN(I,STN) =  37  !  "Cl..CH3..Cl- Sn2 TS" 6 12.011 5 
        IF( TYP(I,STN).EQ."Cl") TYPN(I,STN) =  30  !  "Cl..CH3..Cl- Sn2 TS" 17 35.453 1 
        IF( TYP(I,STN).EQ."HC") TYPN(I,STN) =  36  !  "Cl..CH3..Cl- Sn2 TS" 1 1.008 1 
        IF( TYP(I,STN).EQ."CY") TYPN(I,STN) =  81  !  "Cyclopropane -CH2-" 6 12.011 4  
        IF( TYP(I,STN).EQ."CY") TYPN(I,STN) =  81  !  "Cyclopropane -CHR-" 6 12.011 4  
        IF( TYP(I,STN).EQ."CY") TYPN(I,STN) =  81  !  "Cyclopropane -CR2-" 6 12.011 4  
        IF( TYP(I,STN).EQ."CA") TYPN(I,STN) =  38  !  "Cyclopentadienyl Anion" 6 12.011 3  
        IF( TYP(I,STN).EQ."HA") TYPN(I,STN) =  39  !  "Cyclopentadienyl Anion" 1 1.008 1  
        IF( TYP(I,STN).EQ."CA") TYPN(I,STN) =  38  !  "Cyclopentadienyl Radical" 6 12.011 3  
        IF( TYP(I,STN).EQ."HA") TYPN(I,STN) =  39  !  "Cyclopentadienyl Radical" 1 1.008 1  
        IF( TYP(I,STN).EQ."CA") TYPN(I,STN) =  38  !  "Fluorobenzene CF" 6 12.011 3  
        IF( TYP(I,STN).EQ."F") TYPN(I,STN) =  1  !  "Fluorobenzene CF" 9 18.998 1  
        IF( TYP(I,STN).EQ."CA") TYPN(I,STN) =  38  !  "Hexafluorobenzene CF" 6 12.011 3  
        IF( TYP(I,STN).EQ."F") TYPN(I,STN) =  1  !  "Hexafluorobenzene CF" 9 18.998 1  
        IF( TYP(I,STN).EQ."Br") TYPN(I,STN) =  55  !  "Bromide -CH2-Br (UA)" 35 79.904 1 
        IF( TYP(I,STN).EQ."C2") TYPN(I,STN) =  2  !  "Bromide -CH2-Br (UA)" 6 14.027 0 
        IF( TYP(I,STN).EQ."CA") TYPN(I,STN) =  38  !  "TrifluoroMeBenzene C-CF3" 6 12.011 3  
        IF( TYP(I,STN).EQ."CT") TYPN(I,STN) =  13  !  "TrifluoroMeBenzene CF3-" 6 12.011 4  
        IF( TYP(I,STN).EQ."F") TYPN(I,STN) =  1  !  "TrifluoroMeBenzene CF3-" 9 18.998 1  
        IF( TYP(I,STN).EQ."CA") TYPN(I,STN) =  38  !  "Difluorobenzene CF" 6 12.011 3  
        IF( TYP(I,STN).EQ."F") TYPN(I,STN) =  1  !  "Difluorobenzene CF" 9 18.998 1  
        IF( TYP(I,STN).EQ."CA") TYPN(I,STN) =  38  !  "Bromobenzene CBr" 6 12.011 3  
        IF( TYP(I,STN).EQ."Br") TYPN(I,STN) =  55  !  "Bromobenzene CBr" 35 79.904 1  
        IF( TYP(I,STN).EQ."CA") TYPN(I,STN) =  38  !  "Iodobenzene CI" 6 12.011 3  
        IF( TYP(I,STN).EQ."I") TYPN(I,STN) =  56  !  "Iodobenzene CI" 53 126.905 1  
        IF( TYP(I,STN).EQ."CY") TYPN(I,STN) =  81  !  "cProp/cBut Benzene C-Ar" 6 12.011 4 
        IF( TYP(I,STN).EQ."SH") TYPN(I,STN) =  15  !  "Thiophenol SH" 16 32.060 2  
        IF( TYP(I,STN).EQ."CA") TYPN(I,STN) =  38  !  "Thiophenol C-SH" 6 12.011 3  
        IF( TYP(I,STN).EQ."CA") TYPN(I,STN) =  38  !  "Benzamidine CG" 6 12.011 3  
        IF( TYP(I,STN).EQ."CA") TYPN(I,STN) =  38  !  "Benzamidine CD" 6 12.011 3  
        IF( TYP(I,STN).EQ."CA") TYPN(I,STN) =  38  !  "Benzamidine CE" 6 12.011 3  
        IF( TYP(I,STN).EQ."CA") TYPN(I,STN) =  38  !  "Benzamidine CZ" 6 12.011 3  
        IF( TYP(I,STN).EQ."HA") TYPN(I,STN) =  39  !  "Benzamidine HCD" 1 1.008 1  
        IF( TYP(I,STN).EQ."HA") TYPN(I,STN) =  39  !  "Benzamidine HCE" 1 1.008 1  
        IF( TYP(I,STN).EQ."CA") TYPN(I,STN) =  38  !  "Benzamidine C+" 6 12.011 3  
        IF( TYP(I,STN).EQ."N2") TYPN(I,STN) =  45  !  "Benzamidine -NH2" 7 14.007 3  
        IF( TYP(I,STN).EQ."H") TYPN(I,STN) =  34  !  "Benzamidine H1-N" 1 1.008 1  
        IF( TYP(I,STN).EQ."H") TYPN(I,STN) =  34  !  "Benzamidine H2-N" 1 1.008 1  
        IF( TYP(I,STN).EQ."HA") TYPN(I,STN) =  39  !  "Benzamidine HCG" 1 1.008 1  
        IF( TYP(I,STN).EQ."CT") TYPN(I,STN) =  13  !  "Neutral MeGdn CH3-" 6 12.011 4 
        IF( TYP(I,STN).EQ."CT") TYPN(I,STN) =  13  !  "Neutral ARG CD" 6 12.011 4 
        IF( TYP(I,STN).EQ."NY") TYPN(I,STN) =  91  !  "Neutral ARG NE" 7 14.007 3 
        IF( TYP(I,STN).EQ."NC") TYPN(I,STN) =  46  !  "Neutral ARG N1 (HN=C)" 7 14.007 2
        IF( TYP(I,STN).EQ."NY") TYPN(I,STN) =  91  !  "Neutral ARG N2 (H2N-C)" 7 14.007 3
        IF( TYP(I,STN).EQ."CA") TYPN(I,STN) =  38  !  "Neutral ARG CZ (>C=)" 6 12.011 3
        IF( TYP(I,STN).EQ."NZ") TYPN(I,STN) =  18  !  "Alkyl Nitrile -CN" 7 14.007 1 
        IF( TYP(I,STN).EQ."CZ") TYPN(I,STN) =  19  !  "Alkyl Nitrile -CN" 6 12.011 2 
        IF( TYP(I,STN).EQ."CT") TYPN(I,STN) =  13  !  "Acetonitrile CH3-CN" 6 12.011 4  
        IF( TYP(I,STN).EQ."CT") TYPN(I,STN) =  13  !  "Alkyl Nitrile RCH2-CN" 6 12.011 4 
        IF( TYP(I,STN).EQ."CT") TYPN(I,STN) =  13  !  "Alkyl Nitrile R2CH-CN" 6 12.011 4 
        IF( TYP(I,STN).EQ."CT") TYPN(I,STN) =  13  !  "Alkyl Nitrile R3C-CN" 6 12.011 4 
        IF( TYP(I,STN).EQ."HC") TYPN(I,STN) =  36  !  "Alkyl Nitrile H-C-CN" 1 1.008 1 
        IF( TYP(I,STN).EQ."NO") TYPN(I,STN) =  92  !  "Nitroalkane -NO2" 7 14.007 3  
        IF( TYP(I,STN).EQ."ON") TYPN(I,STN) =  93  !  "Nitroalkane -NO2" 8 15.999 1  
        IF( TYP(I,STN).EQ."CT") TYPN(I,STN) =  13  !  "Nitromethane CH3-NO2" 6 12.011 4  
        IF( TYP(I,STN).EQ."HC") TYPN(I,STN) =  36  !  "Nitroalkane H-C-NO2" 1 1.008 1  
        IF( TYP(I,STN).EQ."CT") TYPN(I,STN) =  13  !  "Nitroalkane RCH2-NO2" 6 12.011 4  
        IF( TYP(I,STN).EQ."CT") TYPN(I,STN) =  13  !  "Nitroalkane R2CH-NO2" 6 12.011 4  
        IF( TYP(I,STN).EQ."CT") TYPN(I,STN) =  13  !  "Nitroalkane R3C-NO2" 6 12.011 4  
        IF( TYP(I,STN).EQ."NO") TYPN(I,STN) =  92  !  "Nitrobenzene -NO2" 7 14.007 3  
        IF( TYP(I,STN).EQ."CA") TYPN(I,STN) =  38  !  "Nitrobenzene C-NO2" 6 12.011 3  
        IF( TYP(I,STN).EQ."CT") TYPN(I,STN) =  13  !  "Benzonitrile -CH2-" 6 12.011 0  
        IF( TYP(I,STN).EQ."NC") TYPN(I,STN) =  46  !  "Neutral Benzamidine N" 7 14.007 2 
        IF( TYP(I,STN).EQ."O") TYPN(I,STN) =  4  !  "Propylene Carbonate C=O" 8 15.999 1 
        IF( TYP(I,STN).EQ."C") TYPN(I,STN) =  3  !  "Propylene Carbonate C=O" 6 12.011 3 
        IF( TYP(I,STN).EQ."OS") TYPN(I,STN) =  29  !  "Propylene Carbonate C-O" 8 15.999 2 
        IF( TYP(I,STN).EQ."CT") TYPN(I,STN) =  13  !  "Propylene Carbonate CH2" 6 12.011 4 
        IF( TYP(I,STN).EQ."CT") TYPN(I,STN) =  13  !  "Propylene Carbonate CH" 6 12.011 4 
        IF( TYP(I,STN).EQ."CT") TYPN(I,STN) =  13  !  "Propylene Carbonate CH3" 6 12.011 4 
        IF( TYP(I,STN).EQ."HC") TYPN(I,STN) =  36  !  "Propylene Carbonate CH2" 1 1.008 1 
        IF( TYP(I,STN).EQ."HC") TYPN(I,STN) =  36  !  "Propylene Carbonate CH" 1 1.008 1 
        IF( TYP(I,STN).EQ."HC") TYPN(I,STN) =  36  !  "Propylene Carbonate CH3" 1 1.008 1 
        IF( TYP(I,STN).EQ."OS") TYPN(I,STN) =  29  !  "GTP O-(POn)2" 8 15.999 2  
        IF( TYP(I,STN).EQ."P+") TYPN(I,STN) =  94  !  "Phosphonium R4P+" 15 30.974 4  
        IF( TYP(I,STN).EQ."CT") TYPN(I,STN) =  13  !  "Phosphonium CH3-PR3+" 6 12.011 4  
        IF( TYP(I,STN).EQ."CT") TYPN(I,STN) =  13  !  "Phosphonium RCH2-PR3+" 6 12.011 4  
        IF( TYP(I,STN).EQ."HC") TYPN(I,STN) =  36  !  "Phosphonium CH3-PR3+" 1 1.008 1  
        IF( TYP(I,STN).EQ."P") TYPN(I,STN) =  54  !  "Hexafluorophosphate Ion" 15 30.974 6  
        IF( TYP(I,STN).EQ."F") TYPN(I,STN) =  1  !  "Hexafluorophosphate Ion" 9 18.998 1  
        IF( TYP(I,STN).EQ."N") TYPN(I,STN) =  35  !  "Nitrate Ion NO3-" 7 14.007 3 
        IF( TYP(I,STN).EQ."O") TYPN(I,STN) =  4  !  "Nitrate Ion NO3-" 8 15.999 1 
        IF( TYP(I,STN).EQ."OW") TYPN(I,STN) =  20  !  "TIP4F Water O" 8 15.999 3 
        IF( TYP(I,STN).EQ."HW") TYPN(I,STN) =  21  !  "TIP4F Water H" 1 1.008 1 
        IF( TYP(I,STN).EQ."LP") TYPN(I,STN) =  22  !  "TIP4F Water M" 99 0.000 1 
        IF( TYP(I,STN).EQ."NT") TYPN(I,STN) =  33  !  "Amine RNH2" 7 14.007 3  
        IF( TYP(I,STN).EQ."NT") TYPN(I,STN) =  33  !  "Amine R2NH" 7 14.007 3  
        IF( TYP(I,STN).EQ."NT") TYPN(I,STN) =  33  !  "Amine R3N" 7 14.007 3  
        IF( TYP(I,STN).EQ."CT") TYPN(I,STN) =  13  !  "Amine CH3-NH2" 6 12.011 4  
        IF( TYP(I,STN).EQ."CT") TYPN(I,STN) =  13  !  "Amine CH3-NHR" 6 12.011 4  
        IF( TYP(I,STN).EQ."CT") TYPN(I,STN) =  13  !  "Amine CH3-NR2" 6 12.011 4  
        IF( TYP(I,STN).EQ."CT") TYPN(I,STN) =  13  !  "Amine RCH2-NH2" 6 12.011 4  
        IF( TYP(I,STN).EQ."CT") TYPN(I,STN) =  13  !  "Amine RCH2-NHR" 6 12.011 4  
        IF( TYP(I,STN).EQ."CT") TYPN(I,STN) =  13  !  "Amine RCH2-NR2" 6 12.011 4  
        IF( TYP(I,STN).EQ."H") TYPN(I,STN) =  34  !  "Amine RNH2" 1 1.008 1  
        IF( TYP(I,STN).EQ."H") TYPN(I,STN) =  34  !  "Amine R2NH" 1 1.008 1  
        IF( TYP(I,STN).EQ."HC") TYPN(I,STN) =  36  !  "Amine H-C-N" 1 1.008 1  
        IF( TYP(I,STN).EQ."CT") TYPN(I,STN) =  13  !  "Amine R2CH-NH2" 6 12.011 4  
        IF( TYP(I,STN).EQ."CT") TYPN(I,STN) =  13  !  "Amine R3C-NH2" 6 12.011 4  
        IF( TYP(I,STN).EQ."CT") TYPN(I,STN) =  13  !  "Amine R2CH-NHR" 6 12.011 4  
        IF( TYP(I,STN).EQ."CT") TYPN(I,STN) =  13  !  "Amine R2CH-NR2" 6 12.011 4  
        IF( TYP(I,STN).EQ."CA") TYPN(I,STN) =  38  !  "Aniline C-NH2" 6 12.011 3  
        IF( TYP(I,STN).EQ."CA") TYPN(I,STN) =  38  !  "N-Me Aniline C-NHR" 6 12.011 3 
        IF( TYP(I,STN).EQ."CA") TYPN(I,STN) =  38  !  "N-DiMe Aniline C-NR2" 6 12.011 3 
        IF( TYP(I,STN).EQ."CT") TYPN(I,STN) =  13  !  "Benzyl Amine -CH2NH2" 6 12.011 4 
        IF( TYP(I,STN).EQ."CT") TYPN(I,STN) =  13  !  "Benzyl Amine -CHRNH2" 6 12.011 4 
        IF( TYP(I,STN).EQ."CT") TYPN(I,STN) =  13  !  "Benzyl Amine -CR2NH2" 6 12.011 4 
        IF( TYP(I,STN).EQ."CT") TYPN(I,STN) =  13  !  "Benzyl Ether -CH2OR" 6 12.011 4 
        IF( TYP(I,STN).EQ."CT") TYPN(I,STN) =  13  !  "Benzyl Sulfide -CH2SH" 6 12.011 4 
        IF( TYP(I,STN).EQ."CT") TYPN(I,STN) =  13  !  "Benzyl Amine -CH2NHR" 6 12.011 4 
        IF( TYP(I,STN).EQ."CZ") TYPN(I,STN) =  19  !  "Alkyne HCC-" 6 12.011 2  
        IF( TYP(I,STN).EQ."HC") TYPN(I,STN) =  36  !  "Alkyne HCC-" 1 1.008 1  
        IF( TYP(I,STN).EQ."CZ") TYPN(I,STN) =  19  !  "Alkyne RCCH R w/ 2/3 H" 6
        IF( TYP(I,STN).EQ."CZ") TYPN(I,STN) =  19  !  "Alkyne RCCH R w/ 1 H" 6
        IF( TYP(I,STN).EQ."CZ") TYPN(I,STN) =  19  !  "Alkyne RCCH R w/ O H/Ph" 6
        IF( TYP(I,STN).EQ."HC") TYPN(I,STN) =  36  !  "Alkyne H-C-CC-" 1 1.008 1  
        IF( TYP(I,STN).EQ."CO") TYPN(I,STN) =  41  !  "A & G Sugar C1'" 6 12.011
        IF( TYP(I,STN).EQ."CO") TYPN(I,STN) =  41  !  "C Sugar C1'" 6 12.011 4 
        IF( TYP(I,STN).EQ."CO") TYPN(I,STN) =  41  !  "U & T Sugar C1'" 6 12.011
        IF( TYP(I,STN).EQ."OH") TYPN(I,STN) =  5  !  "Sugar O5'" 8 15.999 5  
        IF( TYP(I,STN).EQ."HO") TYPN(I,STN) =  7  !  "Sugar H3' (-OH)" 1 1.008 1 
        IF( TYP(I,STN).EQ."N*") TYPN(I,STN) =  95  !  "A & G Nucleoside N9" 7 14.007
        IF( TYP(I,STN).EQ."N*") TYPN(I,STN) =  95  !  "C Nucleoside N1" 7 14.007 3 
        IF( TYP(I,STN).EQ."N*") TYPN(I,STN) =  95  !  "U & T Nucleoside N1" 7 14.007
        IF( TYP(I,STN).EQ."CZ") TYPN(I,STN) =  19  !  "Alkyne RCCR" 6 12.011 2  
        IF( TYP(I,STN).EQ."N3") TYPN(I,STN) =  43  !  "Ammonium R3NH+" 7 14.007 4  
        IF( TYP(I,STN).EQ."H3") TYPN(I,STN) =  44  !  "Ammonium R3NH+" 1 1.008 1  
        IF( TYP(I,STN).EQ."CT") TYPN(I,STN) =  13  !  "Ammonium CH3-NHR2+" 6 12.011 4  
        IF( TYP(I,STN).EQ."CT") TYPN(I,STN) =  13  !  "Ammonium RCH2-NHR2+" 6 12.011 4  
        IF( TYP(I,STN).EQ."CT") TYPN(I,STN) =  13  !  "Ammonium R2CH-NHR2+" 6 12.011 4  
        IF( TYP(I,STN).EQ."CT") TYPN(I,STN) =  13  !  "Ammonium R3C-NHR2+" 6 12.011 4  
        IF( TYP(I,STN).EQ."CW") TYPN(I,STN) =  74  !  "2-Phenyl Furan C2" 6 12.011 3 
        IF( TYP(I,STN).EQ."CS") TYPN(I,STN) =  77  !  "2-Phenyl Furan C3" 6 12.011 3 
        IF( TYP(I,STN).EQ."C!") TYPN(I,STN) =  76  !  "2-Phenyl Furan C2'" 6 12.011 3 
        IF( TYP(I,STN).EQ."C!") TYPN(I,STN) =  76  !  "2-Phenyl Furan C3'" 6 12.011 3 
        IF( TYP(I,STN).EQ."HC") TYPN(I,STN) =  36  !  "GLY Zwitterion HA" 1 1.008 1 
        IF( TYP(I,STN).EQ."CT") TYPN(I,STN) =  13  !  "GLY Zwitterion CA" 6 12.011 4 
        IF( TYP(I,STN).EQ."C") TYPN(I,STN) =  3  !  "GLY Zwitterion C" 6 12.011 3 
        IF( TYP(I,STN).EQ."N3") TYPN(I,STN) =  43  !  "GLY Zwitterion N" 7 14.007 4 
        IF( TYP(I,STN).EQ."O2") TYPN(I,STN) =  42  !  "GLY Zwitterion O" 8 15.999 1 
        IF( TYP(I,STN).EQ."H3") TYPN(I,STN) =  44  !  "GLY Zwitterion HN" 1 1.008 1 
        IF( TYP(I,STN).EQ."F") TYPN(I,STN) =  1  !  "Alkyl Fluoride C-F" 9 18.998 1 
        IF( TYP(I,STN).EQ."CT") TYPN(I,STN) =  13  !  "Alkyl Fluoride RCH2-F" 6 12.011 4 
        IF( TYP(I,STN).EQ."HC") TYPN(I,STN) =  36  !  "Alkyl Fluoride H-C-F" 1 1.008 1 
        IF( TYP(I,STN).EQ."CT") TYPN(I,STN) =  13  !  "Alkyl Fluoride R2CH-F" 6 12.011 4 
        IF( TYP(I,STN).EQ."CT") TYPN(I,STN) =  13  !  "Alkyl Fluoride R3C-F" 6 12.011 4 
        IF( TYP(I,STN).EQ."CT") TYPN(I,STN) =  13  !  "Perfluoroalkane CF3-" 6 12.011 4  
        IF( TYP(I,STN).EQ."CT") TYPN(I,STN) =  13  !  "Perfluoroalkane -CF2-" 6 12.011 4  
        IF( TYP(I,STN).EQ."CT") TYPN(I,STN) =  13  !  "Perfluoroalkane >CF-" 6 12.011 4  
        IF( TYP(I,STN).EQ."CT") TYPN(I,STN) =  13  !  "Tetrafluoromethane CF4" 6 12.011 4  
        IF( TYP(I,STN).EQ."F") TYPN(I,STN) =  1  !  "Perfluoroalkane C-F" 9 18.998 1  
        IF( TYP(I,STN).EQ."CT") TYPN(I,STN) =  13  !  "DifluoroMeBenzene -CHF2" 6 12.011 4  
        IF( TYP(I,STN).EQ."HC") TYPN(I,STN) =  36  !  "DifluoroMeBenzene -CHF2" 1 1.008 1  
        IF( TYP(I,STN).EQ."CT") TYPN(I,STN) =  13  !  "Fluoroacetate FCH2-COO-" 6 12.011 4  
        IF( TYP(I,STN).EQ."CT") TYPN(I,STN) =  13  !  "Chloroacetate ClCH2-COO-" 6 12.011 4  
        IF( TYP(I,STN).EQ."Cl") TYPN(I,STN) =  30  !  "Alkyl Chloride C-Cl" 17 35.453 1 
        IF( TYP(I,STN).EQ."CT") TYPN(I,STN) =  13  !  "Alkyl Chloride RCH2-Cl" 6 12.011 4 
        IF( TYP(I,STN).EQ."HC") TYPN(I,STN) =  36  !  "Alkyl Chloride H-C-Cl" 1 1.008 1 
        IF( TYP(I,STN).EQ."CT") TYPN(I,STN) =  13  !  "Alkyl Chloride R2CH-Cl" 6 12.011 4 
        IF( TYP(I,STN).EQ."CT") TYPN(I,STN) =  13  !  "Alkyl Chloride R3C-Cl" 6 12.011 4 
        IF( TYP(I,STN).EQ."Br") TYPN(I,STN) =  55  !  "Alkyl Bromide C-Br" 35 79.904 1 
        IF( TYP(I,STN).EQ."CT") TYPN(I,STN) =  13  !  "Alkyl Bromide RCH2-Br" 6 12.011 4 
        IF( TYP(I,STN).EQ."HC") TYPN(I,STN) =  36  !  "Alkyl Bromide H-C-Br" 1 1.008 1 
        IF( TYP(I,STN).EQ."CT") TYPN(I,STN) =  13  !  "Alkyl Bromide R2CH-Br" 6 12.011 4 
        IF( TYP(I,STN).EQ."CT") TYPN(I,STN) =  13  !  "Alkyl Bromide R3C-Br" 6 12.011 4 
        IF( TYP(I,STN).EQ."F") TYPN(I,STN) =  1  !  "Acyl Fluoride F-C=O" 9 18.998 1 
        IF( TYP(I,STN).EQ."Cl") TYPN(I,STN) =  30  !  "Acyl Chloride Cl-C=O" 17 35.453 1 
        IF( TYP(I,STN).EQ."Br") TYPN(I,STN) =  55  !  "Acyl Bromide Br-C=O" 35 79.904 1 
        IF( TYP(I,STN).EQ."CA") TYPN(I,STN) =  38  !  "Trifluoroanisole C-OCF3" 6 12.011 3  
        IF( TYP(I,STN).EQ."OS") TYPN(I,STN) =  29  !  "Trifluoroanisole -OCF3" 8 15.999 2  
        IF( TYP(I,STN).EQ."CT") TYPN(I,STN) =  13  !  "Trifluoroanisole -OCF3" 6 12.011 4  
        IF( TYP(I,STN).EQ."F") TYPN(I,STN) =  1  !  "Trifluoroanisole -OCF3" 9 18.998 1  
        IF( TYP(I,STN).EQ."N") TYPN(I,STN) =  35  !  "N-Me,N-PhAcetamide N" 7 14.007 3  
        IF( TYP(I,STN).EQ."CA") TYPN(I,STN) =  38  !  "N-Me,N-PhAcetamide Cipso" 6 12.011 3  
        IF( TYP(I,STN).EQ."CT") TYPN(I,STN) =  13  !  "Benzyl Amine -CH2NR2" 6 12.011 4 
        IF( TYP(I,STN).EQ."C") TYPN(I,STN) =  3  !  "Alkyl Hydroxamic Acid C" 6 12.011 3
        IF( TYP(I,STN).EQ."C") TYPN(I,STN) =  3  !  "Aryl Hydroxamic Acid C" 6 12.011 3
        IF( TYP(I,STN).EQ."O") TYPN(I,STN) =  4  !  "Hydroxamic Acid C=O" 8 15.999 1 
        IF( TYP(I,STN).EQ."N") TYPN(I,STN) =  35  !  "Hydroxamic Acid N" 7 14.007 3 
        IF( TYP(I,STN).EQ."H") TYPN(I,STN) =  34  !  "Hydroxamic Acid HN" 1 1.008 1 
        IF( TYP(I,STN).EQ."OH") TYPN(I,STN) =  5  !  "Hydroxamic Acid OH" 8 15.999 2 
        IF( TYP(I,STN).EQ."HO") TYPN(I,STN) =  7  !  "Hydroxamic Acid OH" 1 1.008 1 
        IF( TYP(I,STN).EQ."CT") TYPN(I,STN) =  13  !  "Benzyl Ether -CHROR" 6 12.011 4 
        IF( TYP(I,STN).EQ."CT") TYPN(I,STN) =  13  !  "Benzyl Ether -CR2OR" 6 12.011 4 
        IF( TYP(I,STN).EQ."C!") TYPN(I,STN) =  76  !  "3-Phenyl Pyrrole C3" 6 12.011 3 
        IF( TYP(I,STN).EQ."C!") TYPN(I,STN) =  76  !  "3-Phenyl Pyrrole C3'" 6 12.011 3 
        IF( TYP(I,STN).EQ."C!") TYPN(I,STN) =  76  !  "4-Phenyl Imidazole C4" 6 12.011 3 
        IF( TYP(I,STN).EQ."C!") TYPN(I,STN) =  76  !  "4-Phenyl Imidazole C4'" 6 12.011 3 
        IF( TYP(I,STN).EQ."CA") TYPN(I,STN) =  38  !  "Diphenylmethane Cipso" 6 12.011 3  
        IF( TYP(I,STN).EQ."Zn") TYPN(I,STN) =  96  !  "Zinc Ion Zn+2" 30 0.000 0 
        IF( TYP(I,STN).EQ."CT") TYPN(I,STN) =  13  !  "Alkyl Iodide RCH2-I" 6 12.011 4 
        IF( TYP(I,STN).EQ."CT") TYPN(I,STN) =  13  !  "Alkyl Iodide R2CH-I" 6 12.011 4 
        IF( TYP(I,STN).EQ."CT") TYPN(I,STN) =  13  !  "Alkyl Iodide R3C-I" 6 12.011 4 
        IF( TYP(I,STN).EQ."I") TYPN(I,STN) =  56  !  "Alkyl Iodide C-I" 53 126.905 1 
        IF( TYP(I,STN).EQ."HC") TYPN(I,STN) =  36  !  "Alkyl Iodide H-C-I" 1 1.008 1 
        IF( TYP(I,STN).EQ."N") TYPN(I,STN) =  35  !  "N-Ph Sulfonamide -NHPh" 7 14.007 3 
        IF( TYP(I,STN).EQ."CA") TYPN(I,STN) =  38  !  "N-Ph Sulfonamide Cipso" 6 12.011 3 
        IF( TYP(I,STN).EQ."CA") TYPN(I,STN) =  38  !  "Benzoate C-COO-" 6 12.011 3  
        IF( TYP(I,STN).EQ."N") TYPN(I,STN) =  35  !  "N-Phenyl Urea N" 7 14.007 3 
        IF( TYP(I,STN).EQ."CA") TYPN(I,STN) =  38  !  "N-Phenyl Urea Cipso" 6 12.011 3 
        IF( TYP(I,STN).EQ."C") TYPN(I,STN) =  3  !  "Tertiary Amide -CO-NR2" 6 12.011 3 
        IF( TYP(I,STN).EQ."O") TYPN(I,STN) =  4  !  "Tertiary Amide -CO-NR2" 8 15.999 1 
        IF( TYP(I,STN).EQ."NM") TYPN(I,STN) =  97  !  "Tertiary Amide -CO-NR2" 7 14.007 3 
        IF( TYP(I,STN).EQ."CT") TYPN(I,STN) =  13  !  "Tertiary Amide -NRCH3" 6 12.011 4 
        IF( TYP(I,STN).EQ."CT") TYPN(I,STN) =  13  !  "Tertiary Amide -NRCH2R" 6 12.011 4 
        IF( TYP(I,STN).EQ."CT") TYPN(I,STN) =  13  !  "Tertiary Amide -NRCHR2" 6 12.011 4 
        IF( TYP(I,STN).EQ."CT") TYPN(I,STN) =  13  !  "Tertiary Amide -NRCR3" 6 12.011 4 
        IF( TYP(I,STN).EQ."HC") TYPN(I,STN) =  36  !  "Tertiary Amide H-C-N" 1 1.008 2 
        IF( TYP(I,STN).EQ."C") TYPN(I,STN) =  3  !  "Tertiary Formamide C=O" 6 12.011 3 
        IF( TYP(I,STN).EQ."O") TYPN(I,STN) =  4  !  "Tertiary Formamide C=O" 8 15.999 1 
        IF( TYP(I,STN).EQ."HC") TYPN(I,STN) =  36  !  "Tertiary Formamide H-C=O" 1 1.008 1 
        IF( TYP(I,STN).EQ."CT") TYPN(I,STN) =  13  !  "B2-Peptide CA" 6 12.011 4  
        IF( TYP(I,STN).EQ."CT") TYPN(I,STN) =  13  !  "B3-Peptide CA Main/N-Ter" 6 12.011 4 
        IF( TYP(I,STN).EQ."CT") TYPN(I,STN) =  13  !  "B3-Pep CB GLY Main/C-Ter" 6 12.011 4
        IF( TYP(I,STN).EQ."CT") TYPN(I,STN) =  13  !  "B3-Pep CB ALA Main/C-Ter" 6 12.011 4
        IF( TYP(I,STN).EQ."CT") TYPN(I,STN) =  13  !  "B3-Pep CB PRO Main/C-Ter" 6 12.011 4
        IF( TYP(I,STN).EQ."CT") TYPN(I,STN) =  13  !  "B3-Peptide CA C-Ter" 6 12.011 4 
        IF( TYP(I,STN).EQ."CT") TYPN(I,STN) =  13  !  "B3-Peptide CB ALA N-Ter" 6 12.011 4
        IF( TYP(I,STN).EQ."CT") TYPN(I,STN) =  13  !  "B3-Peptide CB GLY N-Ter" 6 12.011 4
        IF( TYP(I,STN).EQ."CT") TYPN(I,STN) =  13  !  "B3-Peptide CB PRO N-Ter" 6 12.011 4
        IF( TYP(I,STN).EQ."CT") TYPN(I,STN) =  13  !  "B3-Peptide CE PRO N-Ter" 6 12.011 4
        IF( TYP(I,STN).EQ."Si") TYPN(I,STN) =  98  !  "Alkyl Silane R4Si" 14 28.086 4 
        IF( TYP(I,STN).EQ."Si") TYPN(I,STN) =  98  !  "Alkyl Silane R3SiH" 14 28.086 4 
        IF( TYP(I,STN).EQ."Si") TYPN(I,STN) =  98  !  "Alkyl Silane R2SiH2" 14 28.086 4 
        IF( TYP(I,STN).EQ."Si") TYPN(I,STN) =  98  !  "Alkyl Silane RSiH3" 14 28.086 4 
        IF( TYP(I,STN).EQ."H") TYPN(I,STN) =  34  !  "Alkyl Silane H-C-Si" 1 1.008 1 
        IF( TYP(I,STN).EQ."CT") TYPN(I,STN) =  13  !  "Methyl Silane CH3-Si" 6 12.011 4 
        IF( TYP(I,STN).EQ."CT") TYPN(I,STN) =  13  !  "Alkyl Silane RCH2-Si" 6 12.011 4 
        IF( TYP(I,STN).EQ."CT") TYPN(I,STN) =  13  !  "Alkyl Silane R2CH-Si" 6 12.011 4 
        IF( TYP(I,STN).EQ."CT") TYPN(I,STN) =  13  !  "Alkyl Silane R3C-Si" 6 12.011 4 
        IF( TYP(I,STN).EQ."F") TYPN(I,STN) =  1  !  "Fluoride Ion (GBSA)" 9 18.998 0 
        IF( TYP(I,STN).EQ."Cl") TYPN(I,STN) =  30  !  "Chloride Ion (GBSA)" 17 35.453 0 
        IF( TYP(I,STN).EQ."Br") TYPN(I,STN) =  55  !  "Bromide Ion (GBSA)" 35 79.904 0 
        IF( TYP(I,STN).EQ."I") TYPN(I,STN) =  56  !  "Iodide Ion (GBSA)" 53 126.905 0 
        IF( TYP(I,STN).EQ."Li") TYPN(I,STN) =  58  !  "Lithium Ion (GBSA)" 3 6.941 0 
        IF( TYP(I,STN).EQ."Na") TYPN(I,STN) =  59  !  "Sodium Ion (GBSA)" 11 22.990 0 
        IF( TYP(I,STN).EQ."K") TYPN(I,STN) =  60  !  "Potassium Ion (GBSA)" 19 39.098 0 
        IF( TYP(I,STN).EQ."Rb") TYPN(I,STN) =  61  !  "Rubidium Ion (GBSA)" 37 85.468 0 
        IF( TYP(I,STN).EQ."Cs") TYPN(I,STN) =  62  !  "Cesium Ion (GBSA)" 55 132.905 0 
        IF( TYP(I,STN).EQ."Mg") TYPN(I,STN) =  63  !  "Magnesium Ion (GBSA)" 12 24.305 0 
        IF( TYP(I,STN).EQ."Ca") TYPN(I,STN) =  64  !  "Calcium Ion (GBSA)" 20 40.080 0 
        IF( TYP(I,STN).EQ."Sr") TYPN(I,STN) =  65  !  "Strontium Ion (GBSA)" 38 87.620 0 
        IF( TYP(I,STN).EQ."Ba") TYPN(I,STN) =  66  !  "Barium Ion (GBSA)" 56 137.330 0 
        IF( TYP(I,STN).EQ."CT") TYPN(I,STN) =  13  !  "Ammonium CH3-NR3+" 6 12.011 4  
        IF( TYP(I,STN).EQ."CT") TYPN(I,STN) =  13  !  "Ammonium RCH2-NR3+" 6 12.011 4  
        IF( TYP(I,STN).EQ."CT") TYPN(I,STN) =  13  !  "Ammonium R2CH-NR3+" 6 12.011 4  
        IF( TYP(I,STN).EQ."CT") TYPN(I,STN) =  13  !  "Ammonium R3C-NR3+" 6 12.011 4  
        IF( TYP(I,STN).EQ."HC") TYPN(I,STN) =  36  !  "Ammonium CH3-NR3+" 1 1.008 1  
        IF( TYP(I,STN).EQ."N3") TYPN(I,STN) =  43  !  "Anilinium Ar-NR3+" 7 14.007 4  
        IF( TYP(I,STN).EQ."CA") TYPN(I,STN) =  38  !  "Anilinium C-NR3+" 6 12.011 3  
        IF( TYP(I,STN).EQ."N3") TYPN(I,STN) =  43  !  "Anilinium Ar-NHR2+" 7 14.007 4  
        IF( TYP(I,STN).EQ."CA") TYPN(I,STN) =  38  !  "Anilinium C-NHR2+" 6 12.011 3  
        IF( TYP(I,STN).EQ."C#") TYPN(I,STN) =  99  !  "Triene R2-C= (mid C=C)" 6 12.011 3
        IF( TYP(I,STN).EQ."C#") TYPN(I,STN) =  99  !  "Triene RH-C= (mid C=C)" 6 12.011 3
        IF( TYP(I,STN).EQ."HC") TYPN(I,STN) =  36  !  "Allene/Ketene H-C=C=X" 1 1.008 1  
        IF( TYP(I,STN).EQ."CM") TYPN(I,STN) =  37  !  "Allene/Ketene H2C=C=X" 6 12.011 3  
        IF( TYP(I,STN).EQ."CM") TYPN(I,STN) =  37  !  "Allene/Ketene HRC=C=X" 6 12.011 3  
        IF( TYP(I,STN).EQ."CM") TYPN(I,STN) =  37  !  "Allene/Ketene R2C=C=X" 6 12.011 3  
        IF( TYP(I,STN).EQ."C:") TYPN(I,STN) =  100  !  "Allene =C=" 6 12.011 2  
        IF( TYP(I,STN).EQ."C:") TYPN(I,STN) =  100  !  "Ketene =C=" 6 12.011 2  
        IF( TYP(I,STN).EQ."O") TYPN(I,STN) =  4  !  "Ketene C=O" 8 15.999 1  
        IF( TYP(I,STN).EQ."CT") TYPN(I,STN) =  13  !  "N-Me-HIS CB" 6 12.011 4  
      ENDDO

      RETURN
      END SUBROUTINE  ATYPNOPLSAA 

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!     Version 5.0 06/29/2012 T. W. Kemper                      !
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!     Calculate dipole momment distribution for an index
!
      SUBROUTINE ATYPNAMOEBA(STN)
      USE structure
      USE const
      USE specify
      USE potential
!
      IMPLICIT none
!
      INTEGER :: STN,I
!
      DO I=1,NA(STN)
 
        IF( TYP(I,STN).EQ."He") TYPN(I,STN) =  1  !  "Helium Atom He" 2 4.003 0 
        IF( TYP(I,STN).EQ."Ne") TYPN(I,STN) =  2  !  "Neon Atom Ne" 10 20.179 0 
        IF( TYP(I,STN).EQ."Ar") TYPN(I,STN) =  3  !  "Argon Atom Ar" 18 39.948 0 
        IF( TYP(I,STN).EQ."Kr") TYPN(I,STN) =  4  !  "Krypton Atom Kr" 36 83.800 0 
        IF( TYP(I,STN).EQ."Xe") TYPN(I,STN) =  5  !  "Xenon Atom Xe" 54 131.290 0 
        IF( TYP(I,STN).EQ."Li+") TYPN(I,STN) =  6  !  "Lithium Ion Li+" 3 6.941 0 
        IF( TYP(I,STN).EQ."Na+") TYPN(I,STN) =  7  !  "Sodium Ion Na+" 11 22.990 0 
        IF( TYP(I,STN).EQ."K+") TYPN(I,STN) =  8  !  "Potassium Ion K+" 19 39.098 0 
        IF( TYP(I,STN).EQ."Rb+") TYPN(I,STN) =  9  !  "Rubidium Ion Rb+" 37 85.468 0 
        IF( TYP(I,STN).EQ."Cs+") TYPN(I,STN) =  10  !  "Cesium Ion Cs+" 55 132.905 0 
        IF( TYP(I,STN).EQ."Mg+") TYPN(I,STN) =  11  !  "Magnesium Ion Mg+2" 12 24.305 0 
        IF( TYP(I,STN).EQ."Ca+") TYPN(I,STN) =  12  !  "Calcium Ion Ca+2" 20 40.078 0 
        IF( TYP(I,STN).EQ."F-") TYPN(I,STN) =  13  !  "Fluoride Ion F-" 9 18.998 0 
        IF( TYP(I,STN).EQ."Cl-") TYPN(I,STN) =  14  !  "Chloride Ion Cl-" 17 35.453 0 
        IF( TYP(I,STN).EQ."Br-") TYPN(I,STN) =  15  !  "Bromide Ion Br-" 35 79.904 0 
        IF( TYP(I,STN).EQ."I-") TYPN(I,STN) =  16  !  "Iodide Ion I-" 53 126.904 0 
        IF( TYP(I,STN).EQ."C") TYPN(I,STN) =  17  !  "Cyanide Ion C" 6 12.001 1 
        IF( TYP(I,STN).EQ."N") TYPN(I,STN) =  18  !  "Cyanide Ion N" 7 14.007 1 
        IF( TYP(I,STN).EQ."B") TYPN(I,STN) =  19  !  "Tetrafluoroborate B" 5 10.810 4  
        IF( TYP(I,STN).EQ."F") TYPN(I,STN) =  20  !  "Tetrafluoroborate F" 9 18.998 1  
        IF( TYP(I,STN).EQ."P") TYPN(I,STN) =  21  !  "Hexafluorophosphate P" 15 30.974 6  
        IF( TYP(I,STN).EQ."F") TYPN(I,STN) =  22  !  "Hexafluorophosphate F" 9 18.998 1  
        IF( TYP(I,STN).EQ."N") TYPN(I,STN) =  23  !  "Dinitrogen N2" 7 14.007 1  
        IF( TYP(I,STN).EQ."C") TYPN(I,STN) =  24  !  "Methane CH4" 6 12.011 4  
        IF( TYP(I,STN).EQ."H") TYPN(I,STN) =  25  !  "Methane H4C" 1 1.008 1  
        IF( TYP(I,STN).EQ."C") TYPN(I,STN) =  26  !  "Ethane CH3" 6 12.011 4  
        IF( TYP(I,STN).EQ."H") TYPN(I,STN) =  27  !  "Ethane H3C" 1 1.008 1  
        IF( TYP(I,STN).EQ."C") TYPN(I,STN) =  26  !  "Alkane CH3-" 6 12.011 4  
        IF( TYP(I,STN).EQ."H") TYPN(I,STN) =  27  !  "Alkane H3C-" 1 1.008 1  
        IF( TYP(I,STN).EQ."C") TYPN(I,STN) =  28  !  "Alkane -CH2-" 6 12.011 4  
        IF( TYP(I,STN).EQ."H") TYPN(I,STN) =  29  !  "Alkane -H2C-" 1 1.008 1  
        IF( TYP(I,STN).EQ."C") TYPN(I,STN) =  30  !  "Alkane >CH-" 6 12.011 4  
        IF( TYP(I,STN).EQ."H") TYPN(I,STN) =  31  !  "Alkane -HC<" 1 1.008 1  
        IF( TYP(I,STN).EQ."C") TYPN(I,STN) =  32  !  "Alkane >C<" 6 12.011 4  
        IF( TYP(I,STN).EQ."O") TYPN(I,STN) =  33  !  "Water O" 8 15.999 2  
        IF( TYP(I,STN).EQ."H") TYPN(I,STN) =  34  !  "Water H" 1 1.008 1  
        IF( TYP(I,STN).EQ."O") TYPN(I,STN) =  35  !  "Methanol O" 8 15.999 2  
        IF( TYP(I,STN).EQ."H") TYPN(I,STN) =  36  !  "Methanol HO" 1 1.008 1  
        IF( TYP(I,STN).EQ."C") TYPN(I,STN) =  37  !  "Methanol CH3" 6 12.011 4  
        IF( TYP(I,STN).EQ."H") TYPN(I,STN) =  38  !  "Methanol H3C" 1 1.008 1  
        IF( TYP(I,STN).EQ."O") TYPN(I,STN) =  35  !  "Ethanol O" 8 15.999 2  
        IF( TYP(I,STN).EQ."H") TYPN(I,STN) =  36  !  "Ethanol HO" 1 1.008 1  
        IF( TYP(I,STN).EQ."C") TYPN(I,STN) =  39  !  "Ethanol CH2" 6 12.011 4  
        IF( TYP(I,STN).EQ."H") TYPN(I,STN) =  38  !  "Ethanol H2C" 1 1.008 1  
        IF( TYP(I,STN).EQ."C") TYPN(I,STN) =  39  !  "Ethanol CH3" 6 12.011 4  
        IF( TYP(I,STN).EQ."H") TYPN(I,STN) =  40  !  "Ethanol H3C" 1 1.008 1  
        IF( TYP(I,STN).EQ."C") TYPN(I,STN) =  39  !  "Propanol Me-CH2" 6 12.011 4  
        IF( TYP(I,STN).EQ."H") TYPN(I,STN) =  41  !  "Propanol Me-CH2" 1 1.008 1  
        IF( TYP(I,STN).EQ."C") TYPN(I,STN) =  39  !  "Propanol CH3" 6 12.011 4  
        IF( TYP(I,STN).EQ."H") TYPN(I,STN) =  40  !  "Propanol H3C" 1 1.008 1  
        IF( TYP(I,STN).EQ."O") TYPN(I,STN) =  35  !  "isoPropanol O" 8 15.999 2  
        IF( TYP(I,STN).EQ."H") TYPN(I,STN) =  36  !  "isoPropanol HO" 1 1.008 1  
        IF( TYP(I,STN).EQ."C") TYPN(I,STN) =  42  !  "isoPropanol >CH-" 6 12.011 4  
        IF( TYP(I,STN).EQ."H") TYPN(I,STN) =  38  !  "isoPropanol >CH-" 1 1.008 1  
        IF( TYP(I,STN).EQ."C") TYPN(I,STN) =  39  !  "isoPropanol CH3" 6 12.011 4  
        IF( TYP(I,STN).EQ."H") TYPN(I,STN) =  40  !  "isoPropanol H3C" 1 1.008 1  
        IF( TYP(I,STN).EQ."O") TYPN(I,STN) =  35  !  "Methyl Ether O" 8 15.999 2 
        IF( TYP(I,STN).EQ."C") TYPN(I,STN) =  39  !  "Methyl Ether CH3" 6 12.011 4 
        IF( TYP(I,STN).EQ."H") TYPN(I,STN) =  43  !  "Methyl Ether H3C" 1 1.008 1 
        IF( TYP(I,STN).EQ."N") TYPN(I,STN) =  44  !  "Ammonia N" 7 14.007 3  
        IF( TYP(I,STN).EQ."H") TYPN(I,STN) =  45  !  "Ammonia H3N" 1 1.008 1  
        IF( TYP(I,STN).EQ."N") TYPN(I,STN) =  46  !  "Ammonium Ion N+" 7 14.007 4 
        IF( TYP(I,STN).EQ."H") TYPN(I,STN) =  47  !  "Ammonium Ion H4N+" 1 1.008 1 
        IF( TYP(I,STN).EQ."N") TYPN(I,STN) =  48  !  "Methyl Amine N" 7 14.007 3 
        IF( TYP(I,STN).EQ."H") TYPN(I,STN) =  49  !  "Methyl Amine H2N" 1 1.008 1 
        IF( TYP(I,STN).EQ."C") TYPN(I,STN) =  39  !  "Methyl Amine CH3" 6 12.011 4 
        IF( TYP(I,STN).EQ."H") TYPN(I,STN) =  50  !  "Methyl Amine H3C" 1 1.008 1 
        IF( TYP(I,STN).EQ."N") TYPN(I,STN) =  48  !  "Ethyl Amine N" 7 14.007 3 
        IF( TYP(I,STN).EQ."H") TYPN(I,STN) =  49  !  "Ethyl Amine H2N" 1 1.008 1 
        IF( TYP(I,STN).EQ."C") TYPN(I,STN) =  39  !  "Ethyl Amine CH2" 6 12.011 4 
        IF( TYP(I,STN).EQ."H") TYPN(I,STN) =  50  !  "Ethyl Amine H2C" 1 1.008 1 
        IF( TYP(I,STN).EQ."C") TYPN(I,STN) =  39  !  "Ethyl Amine CH3" 6 12.011 4 
        IF( TYP(I,STN).EQ."H") TYPN(I,STN) =  40  !  "Ethyl Amine H3C" 1 1.008 1 
        IF( TYP(I,STN).EQ."C") TYPN(I,STN) =  39  !  "Propyl Amine CH3" 6 12.011 4 
        IF( TYP(I,STN).EQ."H") TYPN(I,STN) =  40  !  "Propyl Amine H3C" 1 1.008 1 
        IF( TYP(I,STN).EQ."C") TYPN(I,STN) =  39  !  "Propyl Amine Me-CH2" 6 12.011 4 
        IF( TYP(I,STN).EQ."H") TYPN(I,STN) =  41  !  "Propyl Amine Me-CH2" 1 1.008 1 
        IF( TYP(I,STN).EQ."N") TYPN(I,STN) =  48  !  "Dimethyl Amine N" 7 14.007 3 
        IF( TYP(I,STN).EQ."H") TYPN(I,STN) =  49  !  "Dimethyl Amine HN" 1 1.008 1 
        IF( TYP(I,STN).EQ."C") TYPN(I,STN) =  39  !  "Dimethyl Amine CH3" 6 12.011 4 
        IF( TYP(I,STN).EQ."H") TYPN(I,STN) =  50  !  "Dimethyl Amine H3C" 1 1.008 1 
        IF( TYP(I,STN).EQ."N") TYPN(I,STN) =  48  !  "Trimethyl Amine N" 7 14.007 3 
        IF( TYP(I,STN).EQ."C") TYPN(I,STN) =  39  !  "Trimethyl Amine CH3" 6 12.011 4 
        IF( TYP(I,STN).EQ."H") TYPN(I,STN) =  50  !  "Trimethyl Amine H3C" 1 1.008 1 
        IF( TYP(I,STN).EQ."N") TYPN(I,STN) =  48  !  "Pyrrolidine N" 7 14.007 3  
        IF( TYP(I,STN).EQ."H") TYPN(I,STN) =  49  !  "Pyrrolidine HN" 1 1.008 1  
        IF( TYP(I,STN).EQ."C") TYPN(I,STN) =  39  !  "Pyrrolidine C-CH2-C" 6 12.011 4  
        IF( TYP(I,STN).EQ."H") TYPN(I,STN) =  51  !  "Pyrrolidine C-CH2-C" 1 1.008 1  
        IF( TYP(I,STN).EQ."C") TYPN(I,STN) =  39  !  "Pyrrolidine CH2-N" 6 12.011 4  
        IF( TYP(I,STN).EQ."H") TYPN(I,STN) =  51  !  "Pyrrolidine H2C-N" 1 1.008 1  
        IF( TYP(I,STN).EQ."N") TYPN(I,STN) =  48  !  "NMePyrrolidine N" 7 14.007 3  
        IF( TYP(I,STN).EQ."C") TYPN(I,STN) =  39  !  "NMePyrrolidine CH2-N" 6 12.011 4  
        IF( TYP(I,STN).EQ."H") TYPN(I,STN) =  51  !  "NMePyrrolidine H2C-N" 1 1.008 1  
        IF( TYP(I,STN).EQ."C") TYPN(I,STN) =  39  !  "NMePyrrolidine CH2<" 6 12.011 4  
        IF( TYP(I,STN).EQ."H") TYPN(I,STN) =  51  !  "NMePyrrolidine H2C<" 1 1.008 1  
        IF( TYP(I,STN).EQ."C") TYPN(I,STN) =  39  !  "NMePyrrolidine CH3-N" 6 12.011 4  
        IF( TYP(I,STN).EQ."H") TYPN(I,STN) =  41  !  "NMePyrrolidine H3C-N" 1 1.008 1  
        IF( TYP(I,STN).EQ."C") TYPN(I,STN) =  52  !  "Formamide C=O" 6 12.011 3  
        IF( TYP(I,STN).EQ."H") TYPN(I,STN) =  53  !  "Formamide HCO" 1 1.008 1  
        IF( TYP(I,STN).EQ."O") TYPN(I,STN) =  54  !  "Formamide O" 8 15.999 1  
        IF( TYP(I,STN).EQ."N") TYPN(I,STN) =  55  !  "Formamide N" 7 14.007 3  
        IF( TYP(I,STN).EQ."H") TYPN(I,STN) =  56  !  "Formamide H2N" 1 1.008 1  
        IF( TYP(I,STN).EQ."C") TYPN(I,STN) =  52  !  "Acetamide C=O" 6 12.011 3  
        IF( TYP(I,STN).EQ."O") TYPN(I,STN) =  54  !  "Acetamide O" 8 15.999 1  
        IF( TYP(I,STN).EQ."N") TYPN(I,STN) =  55  !  "Acetamide N" 7 14.007 3  
        IF( TYP(I,STN).EQ."H") TYPN(I,STN) =  56  !  "Acetamide H2N" 1 1.008 1  
        IF( TYP(I,STN).EQ."C") TYPN(I,STN) =  39  !  "Acetamide CH3" 6 12.011 4  
        IF( TYP(I,STN).EQ."H") TYPN(I,STN) =  57  !  "Acetamide H3C" 1 1.008 1  
        IF( TYP(I,STN).EQ."C") TYPN(I,STN) =  39  !  "Propamide Me-CH2" 6 12.011 4  
        IF( TYP(I,STN).EQ."H") TYPN(I,STN) =  41  !  "Propamide Me-CH2" 1 1.008 1  
        IF( TYP(I,STN).EQ."C") TYPN(I,STN) =  39  !  "Propamide CH3" 6 12.011 4  
        IF( TYP(I,STN).EQ."H") TYPN(I,STN) =  40  !  "Propamide H3C" 1 1.008 1  
        IF( TYP(I,STN).EQ."C") TYPN(I,STN) =  52  !  "NMeFormamide C=O" 6 12.011 3  
        IF( TYP(I,STN).EQ."H") TYPN(I,STN) =  53  !  "NMeFormamide HCO" 1 1.008 1  
        IF( TYP(I,STN).EQ."O") TYPN(I,STN) =  54  !  "NMeFormamide O" 8 15.999 1  
        IF( TYP(I,STN).EQ."N") TYPN(I,STN) =  55  !  "NMeFormamide N" 7 14.007 3  
        IF( TYP(I,STN).EQ."H") TYPN(I,STN) =  56  !  "NMeFormamide HN" 1 1.008 1  
        IF( TYP(I,STN).EQ."C") TYPN(I,STN) =  39  !  "NMeFormamide CH3" 6 12.011 4  
        IF( TYP(I,STN).EQ."H") TYPN(I,STN) =  58  !  "NMeFormamide H3C" 1 1.008 1  
        IF( TYP(I,STN).EQ."C") TYPN(I,STN) =  39  !  "NEtFormamide CH2-N" 6 12.011 4  
        IF( TYP(I,STN).EQ."C") TYPN(I,STN) =  39  !  "NEtFormamide CH3-" 6 12.011 4  
        IF( TYP(I,STN).EQ."H") TYPN(I,STN) =  40  !  "NEtFormamide H3C" 1 1.008 1  
        IF( TYP(I,STN).EQ."C") TYPN(I,STN) =  52  !  "NMeAcetamide C=O" 6 12.011 3  
        IF( TYP(I,STN).EQ."O") TYPN(I,STN) =  54  !  "NMeAcetamide O" 8 15.999 1  
        IF( TYP(I,STN).EQ."N") TYPN(I,STN) =  55  !  "NMeAcetamide N" 7 14.007 3  
        IF( TYP(I,STN).EQ."H") TYPN(I,STN) =  56  !  "NMeAcetamide HN" 1 1.008 1  
        IF( TYP(I,STN).EQ."C") TYPN(I,STN) =  39  !  "NMeAcetamide CH3-N" 6 12.011 4  
        IF( TYP(I,STN).EQ."H") TYPN(I,STN) =  58  !  "NMeAcetamide H3C-N" 1 1.008 1  
        IF( TYP(I,STN).EQ."C") TYPN(I,STN) =  39  !  "NMeAcetamide CH3-C" 6 12.011 4  
        IF( TYP(I,STN).EQ."H") TYPN(I,STN) =  57  !  "NMeAcetamide H3C-C" 1 1.008 1  
        IF( TYP(I,STN).EQ."C") TYPN(I,STN) =  52  !  "DiMeFormamide C=O" 6 12.011 3  
        IF( TYP(I,STN).EQ."H") TYPN(I,STN) =  53  !  "DiMeFormamide HCO" 1 1.008 1  
        IF( TYP(I,STN).EQ."O") TYPN(I,STN) =  54  !  "DiMeFormamide O" 8 15.999 1  
        IF( TYP(I,STN).EQ."N") TYPN(I,STN) =  55  !  "DiMeFormamide N" 7 14.007 3  
        IF( TYP(I,STN).EQ."C") TYPN(I,STN) =  39  !  "DiMeFormamide CH3" 6 12.011 4  
        IF( TYP(I,STN).EQ."H") TYPN(I,STN) =  58  !  "DiMeFormamide H3C" 1 1.008 1  
        IF( TYP(I,STN).EQ."C") TYPN(I,STN) =  52  !  "DiMeAcetamide C=O" 6 12.011 3  
        IF( TYP(I,STN).EQ."O") TYPN(I,STN) =  54  !  "DiMeAcetamide O" 8 15.999 1  
        IF( TYP(I,STN).EQ."N") TYPN(I,STN) =  55  !  "DiMeAcetamide N" 7 14.007 3  
        IF( TYP(I,STN).EQ."C") TYPN(I,STN) =  39  !  "DiMeAcetamide CH3-N" 6 12.011 4  
        IF( TYP(I,STN).EQ."H") TYPN(I,STN) =  58  !  "DiMeAcetamide H3C-N" 1 1.008 1  
        IF( TYP(I,STN).EQ."C") TYPN(I,STN) =  39  !  "DiMeAcetamide CH3-C" 6 12.011 4  
        IF( TYP(I,STN).EQ."H") TYPN(I,STN) =  57  !  "DiMeAcetamide H3C-C" 1 1.008 1  
        IF( TYP(I,STN).EQ."O") TYPN(I,STN) =  59  !  "Formic Acid OH" 8 15.999 2 
        IF( TYP(I,STN).EQ."H") TYPN(I,STN) =  60  !  "Formic Acid HO" 1 1.008 1 
        IF( TYP(I,STN).EQ."C") TYPN(I,STN) =  61  !  "Formic Acid C=O" 6 12.011 3 
        IF( TYP(I,STN).EQ."O") TYPN(I,STN) =  54  !  "Formic Acid O=C" 8 15.999 1 
        IF( TYP(I,STN).EQ."H") TYPN(I,STN) =  62  !  "Formic Acid HC=O" 1 1.008 1 
        IF( TYP(I,STN).EQ."O") TYPN(I,STN) =  59  !  "Acetic Acid OH" 8 15.999 2 
        IF( TYP(I,STN).EQ."H") TYPN(I,STN) =  60  !  "Acetic Acid HO" 1 1.008 1 
        IF( TYP(I,STN).EQ."C") TYPN(I,STN) =  61  !  "Acetic Acid C=O" 6 12.011 3 
        IF( TYP(I,STN).EQ."O") TYPN(I,STN) =  54  !  "Acetic Acid O=C" 8 15.999 1 
        IF( TYP(I,STN).EQ."C") TYPN(I,STN) =  39  !  "Acetic Acid CH3" 6 12.011 4 
        IF( TYP(I,STN).EQ."H") TYPN(I,STN) =  63  !  "Acetic Acid H3C" 1 1.008 1 
        IF( TYP(I,STN).EQ."C") TYPN(I,STN) =  61  !  "Formaldehyde C=O" 6 12.011 3  
        IF( TYP(I,STN).EQ."O") TYPN(I,STN) =  54  !  "Formaldehyde O=C" 8 15.999 1  
        IF( TYP(I,STN).EQ."H") TYPN(I,STN) =  62  !  "Formaldehyde HC=O" 1 1.008 1  
        IF( TYP(I,STN).EQ."C") TYPN(I,STN) =  61  !  "Acetaldehyde C=O" 6 12.011 3  
        IF( TYP(I,STN).EQ."O") TYPN(I,STN) =  54  !  "Acetaldehyde O=C" 8 15.999 1  
        IF( TYP(I,STN).EQ."C") TYPN(I,STN) =  39  !  "Acetaldehyde CH3" 6 12.011 4  
        IF( TYP(I,STN).EQ."H") TYPN(I,STN) =  63  !  "Acetaldehyde H3C" 1 1.008 1  
        IF( TYP(I,STN).EQ."H") TYPN(I,STN) =  62  !  "Acetaldehyde HC=O" 1 1.008 1  
        IF( TYP(I,STN).EQ."S") TYPN(I,STN) =  64  !  "Hydrogen Sulfide S" 16 32.066 2 
        IF( TYP(I,STN).EQ."H") TYPN(I,STN) =  65  !  "Hydrogen Sulfide H2S" 1 1.008 1 
        IF( TYP(I,STN).EQ."S") TYPN(I,STN) =  66  !  "Methyl Sulfide S" 16 32.066 2 
        IF( TYP(I,STN).EQ."H") TYPN(I,STN) =  67  !  "Methyl Sulfide HS" 1 1.008 1 
        IF( TYP(I,STN).EQ."C") TYPN(I,STN) =  61  !  "Methyl Sulfide CH3" 6 12.011 4 
        IF( TYP(I,STN).EQ."H") TYPN(I,STN) =  68  !  "Methyl Sulfide H3C" 1 1.008 1 
        IF( TYP(I,STN).EQ."S") TYPN(I,STN) =  66  !  "Dimethyl Sulfide S" 16 32.066 2 
        IF( TYP(I,STN).EQ."C") TYPN(I,STN) =  61  !  "Dimethyl Sulfide CH3" 6 12.011 4 
        IF( TYP(I,STN).EQ."H") TYPN(I,STN) =  68  !  "Dimethyl Sulfide H3C" 1 1.008 1 
        IF( TYP(I,STN).EQ."S") TYPN(I,STN) =  66  !  "Dimethyl Disulfide S" 16 32.066 2 
        IF( TYP(I,STN).EQ."C") TYPN(I,STN) =  61  !  "Dimethyl Disulfide CH3" 6 12.011 4 
        IF( TYP(I,STN).EQ."H") TYPN(I,STN) =  68  !  "Dimethyl Disulfide H3C" 1 1.008 1 
        IF( TYP(I,STN).EQ."S") TYPN(I,STN) =  66  !  "Ethyl Sulfide S" 16 32.066 2 
        IF( TYP(I,STN).EQ."H") TYPN(I,STN) =  67  !  "Ethyl Sulfide HS" 1 1.008 1 
        IF( TYP(I,STN).EQ."C") TYPN(I,STN) =  61  !  "Ethyl Sulfide CH2" 6 12.011 4 
        IF( TYP(I,STN).EQ."H") TYPN(I,STN) =  68  !  "Ethyl Sulfide H2C" 1 1.008 1 
        IF( TYP(I,STN).EQ."C") TYPN(I,STN) =  69  !  "Ethyl Sulfide CH3" 6 12.011 4 
        IF( TYP(I,STN).EQ."H") TYPN(I,STN) =  70  !  "Ethyl Sulfide H3C" 1 1.008 1 
        IF( TYP(I,STN).EQ."S") TYPN(I,STN) =  66  !  "MeEt Sulfide S" 16 32.066 2 
        IF( TYP(I,STN).EQ."C") TYPN(I,STN) =  61  !  "MeEt Sulfide CH3-S" 6 12.011 4 
        IF( TYP(I,STN).EQ."H") TYPN(I,STN) =  68  !  "MeEt Sulfide H3C-S" 1 1.008 1 
        IF( TYP(I,STN).EQ."C") TYPN(I,STN) =  61  !  "MeEt Sulfide CH2-S" 6 12.011 4 
        IF( TYP(I,STN).EQ."H") TYPN(I,STN) =  68  !  "MeEt Sulfide CH2-S" 1 1.008 1 
        IF( TYP(I,STN).EQ."C") TYPN(I,STN) =  69  !  "MeEt Sulfide CH3-C" 6 12.011 4 
        IF( TYP(I,STN).EQ."H") TYPN(I,STN) =  70  !  "MeEt Sulfide H3C-C" 1 1.008 1 
        IF( TYP(I,STN).EQ."S") TYPN(I,STN) =  71  !  "Dimethyl Sulfoxide S=O" 16 32.066 3 
        IF( TYP(I,STN).EQ."O") TYPN(I,STN) =  72  !  "Dimethyl Sulfoxide S=O" 8 15.999 1 
        IF( TYP(I,STN).EQ."C") TYPN(I,STN) =  73  !  "Dimethyl Sulfoxide CH3" 6 12.011 4 
        IF( TYP(I,STN).EQ."H") TYPN(I,STN) =  74  !  "Dimethyl Sulfoxide H3C" 1 1.008 1 
        IF( TYP(I,STN).EQ."S") TYPN(I,STN) =  75  !  "Methyl Sulfonate SO3-" 16 32.066 4 
        IF( TYP(I,STN).EQ."O") TYPN(I,STN) =  76  !  "Methyl Sulfonate SO3-" 8 15.999 1 
        IF( TYP(I,STN).EQ."C") TYPN(I,STN) =  39  !  "Methyl Sulfonate CH3" 6 12.011 4 
        IF( TYP(I,STN).EQ."H") TYPN(I,STN) =  77  !  "Methyl Sulfonate H3C" 1 1.008 1 
        IF( TYP(I,STN).EQ."S") TYPN(I,STN) =  75  !  "Ethyl Sulfonate SO3-" 16 32.066 4 
        IF( TYP(I,STN).EQ."O") TYPN(I,STN) =  76  !  "Ethyl Sulfonate SO3-" 8 15.999 1 
        IF( TYP(I,STN).EQ."C") TYPN(I,STN) =  39  !  "Ethyl Sulfonate CH2" 6 12.011 4 
        IF( TYP(I,STN).EQ."H") TYPN(I,STN) =  77  !  "Ethyl Sulfonate H2C" 1 1.008 1 
        IF( TYP(I,STN).EQ."C") TYPN(I,STN) =  39  !  "Ethyl Sulfonate CH3" 6 12.011 4 
        IF( TYP(I,STN).EQ."H") TYPN(I,STN) =  40  !  "Ethyl Sulfonate H3C" 1 1.008 1 
        IF( TYP(I,STN).EQ."C") TYPN(I,STN) =  39  !  "Propyl Sulfonate Me-CH2" 6 12.011 4 
        IF( TYP(I,STN).EQ."H") TYPN(I,STN) =  41  !  "Propyl Sulfonate Me-CH2" 1 1.008 1 
        IF( TYP(I,STN).EQ."C") TYPN(I,STN) =  39  !  "Propyl Sulfonate CH3" 6 12.011 4 
        IF( TYP(I,STN).EQ."H") TYPN(I,STN) =  40  !  "Propyl Sulfonate H3C" 1 1.008 1 
        IF( TYP(I,STN).EQ."C") TYPN(I,STN) =  78  !  "Acetonitrile CN" 6 12.011 2  
        IF( TYP(I,STN).EQ."N") TYPN(I,STN) =  79  !  "Acetonitrile CN" 7 14.007 1  
        IF( TYP(I,STN).EQ."C") TYPN(I,STN) =  80  !  "Acetonitrile CH3" 6 12.011 4  
        IF( TYP(I,STN).EQ."H") TYPN(I,STN) =  81  !  "Acetonitrile H3C" 1 1.008 1  
        IF( TYP(I,STN).EQ."C") TYPN(I,STN) =  82  !  "Tricyanomethide CN" 6 12.011 2  
        IF( TYP(I,STN).EQ."N") TYPN(I,STN) =  83  !  "Tricyanomethide CN" 7 14.007 1  
        IF( TYP(I,STN).EQ."C") TYPN(I,STN) =  84  !  "Tricyanomethide >C-" 6 12.011 3  
        IF( TYP(I,STN).EQ."C") TYPN(I,STN) =  85  !  "Benzene C" 6 12.011 3  
        IF( TYP(I,STN).EQ."H") TYPN(I,STN) =  86  !  "Benzene HC" 1 1.008 1  
        IF( TYP(I,STN).EQ."C") TYPN(I,STN) =  87  !  "Ethylbenzene C1-CH2" 6 12.011 3  
        IF( TYP(I,STN).EQ."C") TYPN(I,STN) =  87  !  "Ethylbenzene C2" 6 12.011 3  
        IF( TYP(I,STN).EQ."C") TYPN(I,STN) =  87  !  "Ethylbenzene C3" 6 12.011 3  
        IF( TYP(I,STN).EQ."C") TYPN(I,STN) =  87  !  "Ethylbenzene C4" 6 12.011 3  
        IF( TYP(I,STN).EQ."H") TYPN(I,STN) =  88  !  "Ethylbenzene H2" 1 1.008 1  
        IF( TYP(I,STN).EQ."H") TYPN(I,STN) =  88  !  "Ethylbenzene H3" 1 1.008 1  
        IF( TYP(I,STN).EQ."H") TYPN(I,STN) =  88  !  "Ethylbenzene H4" 1 1.008 1  
        IF( TYP(I,STN).EQ."C") TYPN(I,STN) =  39  !  "Ethylbenzene CH2" 6 12.011 4  
        IF( TYP(I,STN).EQ."H") TYPN(I,STN) =  40  !  "Ethylbenzene H2C" 1 1.008 1  
        IF( TYP(I,STN).EQ."C") TYPN(I,STN) =  39  !  "Ethylbenzene CH3" 6 12.011 4  
        IF( TYP(I,STN).EQ."H") TYPN(I,STN) =  40  !  "Ethylbenzene H3C" 1 1.008 1  
        IF( TYP(I,STN).EQ."O") TYPN(I,STN) =  35  !  "Phenol HO" 8 15.999 2  
        IF( TYP(I,STN).EQ."H") TYPN(I,STN) =  36  !  "Phenol HO" 1 1.008 1  
        IF( TYP(I,STN).EQ."C") TYPN(I,STN) =  87  !  "Phenol C1-OH" 6 12.011 3  
        IF( TYP(I,STN).EQ."C") TYPN(I,STN) =  87  !  "Phenol C2" 6 12.011 3  
        IF( TYP(I,STN).EQ."C") TYPN(I,STN) =  87  !  "Phenol C3" 6 12.011 3  
        IF( TYP(I,STN).EQ."C") TYPN(I,STN) =  87  !  "Phenol C4" 6 12.011 3  
        IF( TYP(I,STN).EQ."H") TYPN(I,STN) =  88  !  "Phenol H2" 1 1.008 1  
        IF( TYP(I,STN).EQ."H") TYPN(I,STN) =  88  !  "Phenol H3" 1 1.008 1  
        IF( TYP(I,STN).EQ."H") TYPN(I,STN) =  88  !  "Phenol H4" 1 1.008 1  
        IF( TYP(I,STN).EQ."C") TYPN(I,STN) =  39  !  "Toluene CH3" 6 12.011 4  
        IF( TYP(I,STN).EQ."H") TYPN(I,STN) =  40  !  "Toluene H3C" 1 1.008 1  
        IF( TYP(I,STN).EQ."C") TYPN(I,STN) =  87  !  "Toluene C1-CH3" 6 12.011 3  
        IF( TYP(I,STN).EQ."C") TYPN(I,STN) =  87  !  "Toluene C2" 6 12.011 3  
        IF( TYP(I,STN).EQ."C") TYPN(I,STN) =  87  !  "Toluene C3" 6 12.011 3  
        IF( TYP(I,STN).EQ."C") TYPN(I,STN) =  87  !  "Toluene C4" 6 12.011 3  
        IF( TYP(I,STN).EQ."H") TYPN(I,STN) =  88  !  "Toluene H2" 1 1.008 1  
        IF( TYP(I,STN).EQ."H") TYPN(I,STN) =  88  !  "Toluene H3" 1 1.008 1  
        IF( TYP(I,STN).EQ."H") TYPN(I,STN) =  88  !  "Toluene H4" 1 1.008 1  
        IF( TYP(I,STN).EQ."O") TYPN(I,STN) =  35  !  "p-Cresol OH" 8 15.999 2  
        IF( TYP(I,STN).EQ."H") TYPN(I,STN) =  36  !  "p-Cresol HO" 1 1.008 1  
        IF( TYP(I,STN).EQ."C") TYPN(I,STN) =  39  !  "p-Cresol CH3" 6 12.011 4  
        IF( TYP(I,STN).EQ."H") TYPN(I,STN) =  50  !  "p-Cresol H3C" 1 1.008 1  
        IF( TYP(I,STN).EQ."C") TYPN(I,STN) =  87  !  "p-Cresol C1-CH3" 6 12.011 3  
        IF( TYP(I,STN).EQ."C") TYPN(I,STN) =  87  !  "p-Cresol C2" 6 12.011 3  
        IF( TYP(I,STN).EQ."C") TYPN(I,STN) =  87  !  "p-Cresol C3" 6 12.011 3  
        IF( TYP(I,STN).EQ."C") TYPN(I,STN) =  87  !  "p-Cresol C4-OH" 6 12.011 3  
        IF( TYP(I,STN).EQ."H") TYPN(I,STN) =  88  !  "p-Cresol H2" 1 1.008 1  
        IF( TYP(I,STN).EQ."H") TYPN(I,STN) =  88  !  "p-Cresol H3" 1 1.008 1  
        IF( TYP(I,STN).EQ."N") TYPN(I,STN) =  89  !  "Imidazole NH" 7 14.007 3  
        IF( TYP(I,STN).EQ."H") TYPN(I,STN) =  45  !  "Imidazole HN" 1 1.008 1  
        IF( TYP(I,STN).EQ."C") TYPN(I,STN) =  90  !  "Imidazole N-C-N" 6 12.011 3  
        IF( TYP(I,STN).EQ."H") TYPN(I,STN) =  91  !  "Imidazole HC" 1 1.008 1  
        IF( TYP(I,STN).EQ."N") TYPN(I,STN) =  89  !  "Imidazole N=C-" 7 14.007 2  
        IF( TYP(I,STN).EQ."C") TYPN(I,STN) =  90  !  "Imidazole C-N=C" 6 12.011 3  
        IF( TYP(I,STN).EQ."H") TYPN(I,STN) =  91  !  "Imidazole HC" 1 1.008 1  
        IF( TYP(I,STN).EQ."C") TYPN(I,STN) =  90  !  "Imidazole C-NH-" 6 12.011 3  
        IF( TYP(I,STN).EQ."H") TYPN(I,STN) =  91  !  "Imidazole HC" 1 1.008 1  
        IF( TYP(I,STN).EQ."C") TYPN(I,STN) =  39  !  "4-Ethylimidazole CH2" 6 12.011 4  
        IF( TYP(I,STN).EQ."H") TYPN(I,STN) =  40  !  "4-Ethylimidazole H2C" 1 1.008 1  
        IF( TYP(I,STN).EQ."C") TYPN(I,STN) =  39  !  "4-Ethylimidazole CH3" 6 12.011 4  
        IF( TYP(I,STN).EQ."H") TYPN(I,STN) =  40  !  "4-Ethylimidazole H3C" 1 1.008 1  
        IF( TYP(I,STN).EQ."N") TYPN(I,STN) =  89  !  "4-Ethylimidazole ND" 7 14.007 3  
        IF( TYP(I,STN).EQ."H") TYPN(I,STN) =  45  !  "4-Ethylimidazole HND" 1 1.008 1  
        IF( TYP(I,STN).EQ."C") TYPN(I,STN) =  90  !  "4-Ethylimidazole CE" 6 12.011 3  
        IF( TYP(I,STN).EQ."H") TYPN(I,STN) =  91  !  "4-Ethylimidazole HCE" 1 1.008 1  
        IF( TYP(I,STN).EQ."N") TYPN(I,STN) =  89  !  "4-Ethylimidazole NE" 7 14.007 3  
        IF( TYP(I,STN).EQ."H") TYPN(I,STN) =  45  !  "4-Ethylimidazole HNE" 1 1.008 1  
        IF( TYP(I,STN).EQ."C") TYPN(I,STN) =  90  !  "4-Ethylimidazole CD" 6 12.011 3  
        IF( TYP(I,STN).EQ."H") TYPN(I,STN) =  91  !  "4-Ethylimidazole HCD" 1 1.008 1  
        IF( TYP(I,STN).EQ."C") TYPN(I,STN) =  90  !  "4-Ethylimidazole CG" 6 12.011 3  
        IF( TYP(I,STN).EQ."N") TYPN(I,STN) =  89  !  "Indole N" 7 14.007 3  
        IF( TYP(I,STN).EQ."H") TYPN(I,STN) =  45  !  "Indole HN" 6 1.008 1  
        IF( TYP(I,STN).EQ."C") TYPN(I,STN) =  92  !  "Indole C2" 6 12.011 3  
        IF( TYP(I,STN).EQ."C") TYPN(I,STN) =  92  !  "Indole C3" 6 12.011 3  
        IF( TYP(I,STN).EQ."C") TYPN(I,STN) =  92  !  "Indole C4" 6 12.011 3  
        IF( TYP(I,STN).EQ."C") TYPN(I,STN) =  92  !  "Indole C5" 6 12.011 3  
        IF( TYP(I,STN).EQ."C") TYPN(I,STN) =  92  !  "Indole C6" 6 12.011 3  
        IF( TYP(I,STN).EQ."C") TYPN(I,STN) =  92  !  "Indole C7" 6 12.011 3  
        IF( TYP(I,STN).EQ."C") TYPN(I,STN) =  92  !  "Indole C3a" 6 12.011 3  
        IF( TYP(I,STN).EQ."C") TYPN(I,STN) =  92  !  "Indole C7a" 6 12.011 3  
        IF( TYP(I,STN).EQ."H") TYPN(I,STN) =  86  !  "Indole HC2" 1 1.008 1  
        IF( TYP(I,STN).EQ."H") TYPN(I,STN) =  86  !  "Indole HC3" 1 1.008 1  
        IF( TYP(I,STN).EQ."H") TYPN(I,STN) =  86  !  "Indole HC4" 1 1.008 1  
        IF( TYP(I,STN).EQ."H") TYPN(I,STN) =  86  !  "Indole HC5" 1 1.008 1  
        IF( TYP(I,STN).EQ."H") TYPN(I,STN) =  86  !  "Indole HC6" 1 1.008 1  
        IF( TYP(I,STN).EQ."H") TYPN(I,STN) =  86  !  "Indole HC7" 1 1.008 1  
        IF( TYP(I,STN).EQ."N") TYPN(I,STN) =  89  !  "3-Ethylindole N" 7 14.007 3  
        IF( TYP(I,STN).EQ."H") TYPN(I,STN) =  45  !  "3-Ethylindole HN" 1 1.008 1  
        IF( TYP(I,STN).EQ."C") TYPN(I,STN) =  92  !  "3-Ethylindole C2" 6 12.011 3  
        IF( TYP(I,STN).EQ."C") TYPN(I,STN) =  92  !  "3-Ethylindole C3" 6 12.011 3  
        IF( TYP(I,STN).EQ."C") TYPN(I,STN) =  92  !  "3-Ethylindole C4" 6 12.011 3  
        IF( TYP(I,STN).EQ."C") TYPN(I,STN) =  92  !  "3-Ethylindole C5" 6 12.011 3  
        IF( TYP(I,STN).EQ."C") TYPN(I,STN) =  92  !  "3-Ethylindole C6" 6 12.011 3  
        IF( TYP(I,STN).EQ."C") TYPN(I,STN) =  92  !  "3-Ethylindole C7" 6 12.011 3  
        IF( TYP(I,STN).EQ."C") TYPN(I,STN) =  92  !  "3-Ethylindole C3a" 6 12.011 3  
        IF( TYP(I,STN).EQ."C") TYPN(I,STN) =  92  !  "3-Ethylindole C7a" 6 12.011 3  
        IF( TYP(I,STN).EQ."H") TYPN(I,STN) =  86  !  "3-Ethylindole HC2" 1 1.008 1  
        IF( TYP(I,STN).EQ."H") TYPN(I,STN) =  86  !  "3-Ethylindole HC4" 1 1.008 1  
        IF( TYP(I,STN).EQ."H") TYPN(I,STN) =  86  !  "3-Ethylindole HC5" 1 1.008 1  
        IF( TYP(I,STN).EQ."H") TYPN(I,STN) =  86  !  "3-Ethylindole HC6" 1 1.008 1  
        IF( TYP(I,STN).EQ."H") TYPN(I,STN) =  86  !  "3-Ethylindole HC7" 1 1.008 1  
        IF( TYP(I,STN).EQ."C") TYPN(I,STN) =  93  !  "3-Ethylindole CH2" 6 12.011 4  
        IF( TYP(I,STN).EQ."H") TYPN(I,STN) =  40  !  "3-Ethylindole H2C" 1 1.008 1  
        IF( TYP(I,STN).EQ."C") TYPN(I,STN) =  93  !  "3-Ethylindole CH3" 6 12.011 4  
        IF( TYP(I,STN).EQ."H") TYPN(I,STN) =  40  !  "3-Ethylindole H3C" 1 1.008 1  
        IF( TYP(I,STN).EQ."N") TYPN(I,STN) =  89  !  "3-Formylindole N" 7 14.007 3  
        IF( TYP(I,STN).EQ."H") TYPN(I,STN) =  45  !  "3-Formylindole HN" 1 1.008 1  
        IF( TYP(I,STN).EQ."C") TYPN(I,STN) =  92  !  "3-Formylindole C2" 6 12.011 3  
        IF( TYP(I,STN).EQ."C") TYPN(I,STN) =  92  !  "3-Formylindole C3" 6 12.011 3  
        IF( TYP(I,STN).EQ."C") TYPN(I,STN) =  92  !  "3-Formylindole C4" 6 12.011 3  
        IF( TYP(I,STN).EQ."C") TYPN(I,STN) =  92  !  "3-Formylindole C5" 6 12.011 3  
        IF( TYP(I,STN).EQ."C") TYPN(I,STN) =  92  !  "3-Formylindole C6" 6 12.011 3  
        IF( TYP(I,STN).EQ."C") TYPN(I,STN) =  92  !  "3-Formylindole C7" 6 12.011 3  
        IF( TYP(I,STN).EQ."C") TYPN(I,STN) =  92  !  "3-Formylindole C3a" 6 12.011 3  
        IF( TYP(I,STN).EQ."C") TYPN(I,STN) =  92  !  "3-Formylindole C7a" 6 12.011 3  
        IF( TYP(I,STN).EQ."H") TYPN(I,STN) =  86  !  "3-Formylindole HC2" 1 1.008 1  
        IF( TYP(I,STN).EQ."H") TYPN(I,STN) =  86  !  "3-Formylindole HC4" 1 1.008 1  
        IF( TYP(I,STN).EQ."H") TYPN(I,STN) =  86  !  "3-Formylindole HC5" 1 1.008 1  
        IF( TYP(I,STN).EQ."H") TYPN(I,STN) =  86  !  "3-Formylindole HC6" 1 1.008 1  
        IF( TYP(I,STN).EQ."H") TYPN(I,STN) =  86  !  "3-Formylindole HC7" 1 1.008 1  
        IF( TYP(I,STN).EQ."C") TYPN(I,STN) =  52  !  "3-Formylindole C=O" 6 12.011 3  
        IF( TYP(I,STN).EQ."O") TYPN(I,STN) =  94  !  "3-Formylindole O=C" 8 15.999 1  
        IF( TYP(I,STN).EQ."H") TYPN(I,STN) =  62  !  "3-Formylindole HC=O" 1 1.008 1  
        IF( TYP(I,STN).EQ."N") TYPN(I,STN) =  95  !  "Benzamidine N" 7 14.007 3  
        IF( TYP(I,STN).EQ."H") TYPN(I,STN) =  96  !  "Benzamidine HN" 1 1.008 1  
        IF( TYP(I,STN).EQ."C") TYPN(I,STN) =  97  !  "Benzamidine N-C-N" 6 12.011 3  
        IF( TYP(I,STN).EQ."C") TYPN(I,STN) =  87  !  "Benzamidine C1-CN2" 6 12.011 3  
        IF( TYP(I,STN).EQ."C") TYPN(I,STN) =  87  !  "Benzamidine C2" 6 12.011 3  
        IF( TYP(I,STN).EQ."C") TYPN(I,STN) =  87  !  "Benzamidine C3" 6 12.011 3  
        IF( TYP(I,STN).EQ."C") TYPN(I,STN) =  87  !  "Benzamidine C4" 6 12.011 3  
        IF( TYP(I,STN).EQ."H") TYPN(I,STN) =  88  !  "Benzamidine H2" 1 1.008 1  
        IF( TYP(I,STN).EQ."H") TYPN(I,STN) =  88  !  "Benzamidine H3" 1 1.008 1  
        IF( TYP(I,STN).EQ."H") TYPN(I,STN) =  88  !  "Benzamidine H4" 1 1.008 1  
        IF( TYP(I,STN).EQ."N") TYPN(I,STN) =  98  !  "Pyridinium N" 7 14.007 3  
        IF( TYP(I,STN).EQ."C") TYPN(I,STN) =  85  !  "Pyridinium C2" 6 12.011 3  
        IF( TYP(I,STN).EQ."C") TYPN(I,STN) =  85  !  "Pyridinium C3" 6 12.011 3  
        IF( TYP(I,STN).EQ."C") TYPN(I,STN) =  85  !  "Pyridinium C4" 6 12.011 3  
        IF( TYP(I,STN).EQ."H") TYPN(I,STN) =  99  !  "Pyridinium HN" 1 1.008 1  
        IF( TYP(I,STN).EQ."H") TYPN(I,STN) =  86  !  "Pyridinium H2" 1 1.008 1  
        IF( TYP(I,STN).EQ."H") TYPN(I,STN) =  86  !  "Pyridinium H3" 1 1.008 1  
        IF( TYP(I,STN).EQ."H") TYPN(I,STN) =  86  !  "Pyridinium H4" 1 1.008 1  
      ENDDO

      RETURN
      END SUBROUTINE  ATYPNAMOEBA
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
