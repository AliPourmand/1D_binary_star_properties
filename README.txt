# 1D_binary_star_properties

The main purpose of this repository is to provide numerical tables describing binary properties. The description of the method with which the tables were obtained can be found in the accomponiing paper ``Properties of binary systems in a one-dimensional approximation". IF you use the tables in your research, we ask you to cite the paper.

The provided sample code demonstrates how to obtaine the 1D effective acceleration a star in a binary system feels (which is a function of the distance to the center of the star).


The effective acceleration at each local zone of the star is a function of the following parameters, which have been denoted in the way they appear in the gravity_table_use_example.f90 code:

1. M1: donor star's mass
2. M2: companion star's mass
3. a: the orbital separation
4. M_loc: the enclosed mass, up to the zone of interest of the donor star
5. r_loc: the distance from the center to the zone of interest of the donor star


These 5 parameters can be given as inputs, and the code will give the effective acceleration of that certain zone as an output.

There is an parameter called "flag_accel" which enables one to toggle between CGS units "flag_accel=1", and in the units of the star if it wasn't in a binary system "flag_accel=2".  
The folder "binary_tables" contains a list of 11 parameters that have been numerically calculated for 109 different mass ratios spanning from q =10^{-6}  all the way to q =10^+5. These parameters are a function of the equipotential shell's corresponding potential, which start from the shell passing x=0.05x_{L1}, all the way to the shell passing the outer lagrange point.

We have provided extra information as well, that could be useful for mass transfer calculations.


The full format of the table is listed below:

column 1: mass ratio of system (donor to companion)  
column 2: volume equivalent radius of shell divided by Eggleton Radius  
column 3: relative potential as defined in Mochnacki (1984)  eg, for potentials higher than psiL1: F=psiL1/psi  
column 4: potential of shell  
column 5: volume equivalent radius of shell divided by orbital separation  
column 6: avg acceleration on shell   
column 7: total area of shell  
column 8: avg acceleration on L plane cross section with shell  
column 9: area of L plane cross section with shell  
column 10: semi major axis of L plane cross section with shell  
column 11: semi minor axis of L plane cross section with shell  


The folder "binary_tables_compact_ss_units" (which should be unzipped for the code to work); only includes properties that are useful in the described subroutine, and for mass transfer calculations. The 3 properties in column 4, 6, and 8 are now provided in the units of a single point-mass star's gravitational acceleration at that zone.
The full format of these tables are listed below:

column 1: volume equivalent radius of shell divided by orbital separation  
column 2: relative potential of shell  
column 3: relative avg acceleration on shell   
column 4: relative avg acceleration on L plane cross section with shell  
column 5: area of L plane cross section with shell  

The file l1l2l3_props.dat includes properties of the three lagrange points for these 109 mass ratios. The properties are listed below:

column 1: mass ratio of system (donor to companion)  
column 2: volume equivalent radius of shell passing L1  
column 3: volume equivalent radius of shell passing L2  
column 4: volume equivalent radius of shell passing L3  
column 5: area of L plane cross section with shell passing L2  
column 6: area of L plane cross section with shell passing L3  
column 7: avg acceleration on L plane cross section with shell passing L2  
column 8: avg acceleration on L plane cross section with shell passing L3  
