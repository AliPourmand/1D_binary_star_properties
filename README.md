# 1D_binary_star_properties

The main purpose of this repository is to provide numerical tables describing binary properties. The description of the method with which the tables were obtained can be found in the accompanying paper ``Properties of binary systems in a one-dimensional approximation". If you use the tables in your research, we ask you to cite the paper. 

The provided sample code demonstrates how to obtaine the 1D effective acceleration a star in a binary system feels (which is a function of the distance to the center of the star).

The effective acceleration at each local zone of the star is a function of the following parameters, which have been denoted in the way they appear in the gravity_table_use_example.f90 code:

1. $M_1$: donor star's mass
2. $M_2$: companion star's mass
3. $a$: the orbital separation
4. $M_{loc}$: the enclosed mass, up to the zone of interest of the donor star
5. $r_{loc}$: the distance from the center to the zone of interest of the donor star


These 5 parameters can be given as inputs, and the code will give the effective acceleration of that certain zone as an output.

There is an parameter called "flag_accel" which enables one to toggle between CGS units "flag_accel=1", and in the units of the star if it wasn't in a binary system "flag_accel=2".  
The folder "binary_tables" contains a list of 11 parameters that have been numerically calculated for 109 different mass ratios spanning from $q =10^{-6}$  all the way to $q =10^{+5}$. These parameters are a function of the equipotential shell's corresponding potential, which start from the shell passing $x=0.05x_{L_1}$, all the way to the shell passing the outer lagrange point.

We have provided extra information as well, that could be useful for mass transfer calculations.


The full format of the table is listed below:

column 1: equipotential shell number

column 2: mass ratio of system (donor to companion)  

column 3: relative potential as defined in Mochnacki (1984)  eg, for potentials higher than psiL1: F=psiL1/psi  

column 4: potential of shell  

column 5: volume equivalent radius of shell divided by orbital separation  

column 6: avg acceleration on shell (3D surface integrations including the $L_1$ plane) 

column 7: inverse effective acceleration averaged on the equipotential shell (3D surface integrations including the $L_1$ plane) 

column 8: total area of shell (3D surface integrations including the $L_1$ plane) 

column 9: area of $L_1$ plane cross section with shell  

column 10: y intersection of $L_1$ plane cross section with shell  

column 11: z intersection of $L_1$ plane cross section with shell  

column 12: effective acceleration averaged over the intersection with $L_1$-plane (2D integrations)

column 13: effective acceleration averaged over $L_1$ cross-section area (2D integrations)

column 14: inverse effective acceleration averaged over $L_1$ cross-section area (2D integrations)

The folder "binary_tables_compact_ss_units" (which should be unzipped for the code to work); only includes properties that are useful in the described subroutine, and for mass transfer calculations. The 3 properties in column 4, 6, and 13 are now provided in the units of a single point-mass star's gravitational acceleration at that zone.
The full format of these tables are listed below:

column 1: volume equivalent radius of shell divided by orbital separation  
column 2: relative potential of shell  
column 3: relative avg acceleration on shell (including the $L_1$ plane)  
column 4: relative avg acceleration on L plane cross section with shell  
column 5: area of L plane cross section with shell  
column 6: total area of shell (including the $L_1$ plane)

Note that the provided code that uses the table above reconstructs the truncated (limited by $L_1$ plane) values using 2D values in L1 plane; see the paper. Should you need equipotentials for a star that has $\log_{10} q> 2.5$, you should reduce into the same format tables_without_l1

The file l1l2l3_props.dat includes properties of the three lagrange points for these 109 mass ratios. The properties are listed below:

column 1: mass ratio of system (donor to companion)  
column 2: volume equivalent radius of shell passing $L_1$  
column 3: volume equivalent radius of shell passing $L_2$  
column 4: volume equivalent radius of shell passing $L_3$  
column 5: area of L plane cross section with shell passing $L_2$  
column 6: area of L plane cross section with shell passing $L_3$  
column 7: avg acceleration on L plane cross section with shell passing $L_2$  
column 8: avg acceleration on L plane cross section with shell passing $L_3$  



The supplementary data tables:

tables_without_l1 is essentially the same as the main table, but 3D integrations were performed excluding $L_1$ plane 

l1_etax	provides integrations for the x-component of the effective acceleration on the $L_1$ plane

column 1: equipotential shell number

column 2: potential of shell

column 3: Area of $L_1$ plane within the potential

column 4: y intersection of $L_1$ plane cross section with the potential

column 5: z intersection of $L_1$ plane cross section with the potential

column 6: x-component of the effective acceleration averaged over $L_1$ cross-section area

column 7: x-component of the effective acceleration averaged over the intersection with $L_1$-plane


