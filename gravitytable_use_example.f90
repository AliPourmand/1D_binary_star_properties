
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  example of initialization, of the use of effective gravity
!!  and of the other useful routines
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      program modified_grav
        implicit none
        
        real*8 effective_gravity
        real*8 grav_final_real, M1, M2, M_loc,a,r_loc, r_max, r_l1
        integer flag_accel
        real*8 r_linner, r_louter

        integer i

        ! first star in solar masses
        M1=10.d0
        ! second star in solar masses
        M2=10.d0
        ! orbital separation in solar radii
        a=10.d0
        ! flag_accel = 1 returns acceleration in CGS units
        ! flag_accel = 2 returns acceleration in units of single star acceleration
        flag_accel=2


        ! this pre-reads tables as well as initializes L1/L2/L3 points
        call gravity_table_ini(M1,M2,r_linner,r_louter)


        ! r_lout is volume equivalent radius for the outer Lagrangian point
        r_louter= r_louter*a! here add the read of the outer L point from the table
        ! write(*,*) "r_lout=",r_louter
        ! r_lin is volume equivalent radius for the inner Lagrangian point
        r_linner= r_linner*a! here add the read of the L1 point from the table

        
        do i=1, 100
! this will cycle through point mass gravity at differenct distance approximately till it will fills L3           
        ! r_loc voluem equivalent radis
          r_loc=r_louter*(0.01d0*i)
        ! mass of the first star inside r_loc
          ! M_loc=(M1/M2)
           ! for experiments with the local mass ratio, here the case of a constant density star is provided. 
           ! note that the code will not allow to use a negative density gradient.
          ! if (r_loc.ge.r_lin) then
          M_loc = M1*(r_loc/r_louter)**3.d0
          ! M_loc = 0.5d0*M1
          ! write(*,*) "M_loc=", M_loc
          ! end if 
          grav_final_real = effective_gravity(M1,M2,a,r_loc,M_loc,flag_accel)
          write(*,*) "call for", i,  r_loc/a, a, grav_final_real
        end do
        
        stop
      end program modified_grav




