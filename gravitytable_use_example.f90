
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  example of initialization, of the use of effective gravity
!!  and of the other useful routines
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      program modified_grav
        implicit none
        
        real*8 effective_gravity
        real*8 grav_final_real, M1, M2, M_loc,a,r_loc
        integer flag_accel

        integer i

        M1=1.d0
        M2=5.d0
        a=10.d0
        M_loc=0.3d0*M1
        r_loc= 0.43*a  !!!!in the units of the orbital separaton a
        flag_accel=2

        call gravity_table_ini()

        do i=1, 50
! this will cycle through point mass gravity at differenct distance approximately till it will fill L3           
           r_loc=a*(0.0385*i)
           grav_final_real = effective_gravity(M1,M2,a,r_loc,M_loc,flag_accel)
           write(*,*) "call for", i,  r_loc, a, grav_final_real
        end do
        
        stop
      end program modified_grav




