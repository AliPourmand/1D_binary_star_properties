# 1D_binary_star_properties

The provided codes have been made with the purpose of obtaining the 1D effective acceleration a star in a binary system feels (which is a function of the distance to the center of the star).

The effective acceleration at each local zone of the star is a function of the following parameters, which have been denoted in the way they appear in the gravity_table_use_example.f90 code:

1. M1: donor star's mass
2. M2: companion star's mass
3. a: the orbital separation
4. M_loc: the enclosed mass, up to the zone of interest of the donor star
5. r_loc: the distance from the center to the zone of interest of the donor star


These 5 parameters can be given as inputs, and the code will give the effective acceleration of that certain zone as an output.

There is an parameter called "flag_accel" which enables one to toggle between CGS units "flag_accel=1", and in the units of the star if it wasn't in a binary system "flag_accel=2".

The folder "binary tables" contains a list of 11 parameters that have been numerically calculated for 109 different mass ratios spanning from q =1e-6  all the way to q =1e+5. These parameters are a function of the equipotential shell's corresponding potential. To find the effective acceleration we only need to store in memory two columns, 5 and 6.
We have provided extra information as well, that could be useful for mass transfer calculations.

The full format of the list is listed below:

!! 1 column: mass ratio of system (donor to companion)
!! 2 column: volume equivalent radius of shell divided by Eggleton Radius
!! 3 column: relative potential as defined in Mochnacki (1984)  eg, for potentials higher than psiL1: F=psiL1/psi
!! 4 column: potential of shell
!! 5 column: volume equivalent radius of shell divided by orbital separation
!! 6 column: avg acceleration on shell. 
!! 7 column: total area of shell
!! 8 column: avg acceleration on L plane cross section with shell
!! 9 column: area of L plane cross section with shell
!! 10 column: semi major axis of L plane cross section with shell
!! 11 column: semi minor axis of L plane cross section with shell
