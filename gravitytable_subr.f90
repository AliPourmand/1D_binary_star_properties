!!! #######################################################################
!!! ###  reading pre-calacuated table with effective acceleration       ###
!!! ###  by default, the files are assumed to be in ./binary directory  ###
!!! #######################################################################      

      subroutine gravity_table_ini

        implicit none
        
        include "gravitytable.inc"        
        
!!! internal variables                 
        integer :: j,l
        integer :: n_q  !! this is the number of table that is read
        
        real*8 :: r_min_array(c_gtable), lnq_table(c_gtable)
        character*255 :: file_name
                        
!!! dummy properties are those getting read from the tables        
        real*8 :: q_dummy, rfrac_dummy,F_dummy, xi_arr_dummy
        real*8 :: A_shell_dummy, g_Lplane_arr_dummy, A_tot_lplane_dummy
        real*8 :: y_lplane_dummy, z_lplane_dummy        
        
        flag_gtable_ini=0

        
        !!size gives length of array
        n_q=c_gtable

         do j=1,16

            lnq_table(j)=-6.d0+(j-1)*0.25d0
            q_table(j)=10.d0**lnq_table(j)

         end do

         do j=17,96

            lnq_table(j)=-2.d0+(j-16-1)*0.05d0
            q_table(j)=10.d0**lnq_table(j)

         end do

         do j=97,n_q

            lnq_table(j)=2.d0+(j-96-1)*0.25d0
            q_table(j)=10.d0**lnq_table(j)

         end do        
    

        do j=1,n_q
            if (lnq_table(j).lt.0) then 
               write(file_name,'(A31,F6.3,A4)') "binary_tables/binary_props_lnq="&
               ,lnq_table(j),".dat"
            end if
            if (lnq_table(j).ge.-0.001d0) then 
               write(file_name,'(A31,F5.3,A4)') "binary_tables/binary_props_lnq="&
               ,lnq_table(j),".dat"
            end if

            ! write (*,*) file_name
           
            open (unit = 23, file = file_name)

           do l=1,n_gtable

!! The complete format is as below. To find effective acceleration we only need to store in memory two columns, 5 and 6.
!! Extra information is provided for MT calcualtions.
!! 1 column: mass ratio of system (donor to companion)
!! 2 column: volume equivalent radius of shell divided by eggleton radius
!! 3 column: relative potential as defined in mochnacki  eg, for potentials higher than psiL1: F=psiL1/psi
!! 4 column: potential of shell
!! 5 column: volume equivalent radius of shell divided by orbital separation
!! 6 column: avg acceleration on shell. 
!! 7 column: total area of shell
!! 7 column: avg acceleration on L plane cross section with shell
!! 8 column: area of L plane cross section with shell
!! 9 column: semi major axis of L plane cross section with shell
!! 10 column: semi minor axis of L plane cross section with shell
            
              read (23,*) q_dummy, rfrac_dummy,F_dummy,&
                   xi_arr_dummy, table_rel_r(j,l), table_grav(j,l),&
                   A_shell_dummy, g_Lplane_arr_dummy, A_tot_lplane_dummy,&
                   y_lplane_dummy, z_lplane_dummy

           end do
           close (23)

           r_min_array(j)=table_rel_r(j,1)

        end do
        
        r_total_min=minval(r_min_array)
        
        flag_gtable_ini=1

      end subroutine gravity_table_ini
    


!!! ##################################################################################
!!! ###  uses pre-calacuated and pre-read tables with effective acceleration       ###
!!! ###  input is in solar units, returns effective acceleration                   ###      
!!! ###  M1 is the total mass (in Msun) of the star for which we find acceleration ###
!!! ###  M2 is the total mass (in Msun) of the companion star                      ###
!!! ###  a is the orbital separation (in Rsun)                                     ###
!!! ###  r_loc is radial coordinate (in Rsun) inside the "donor" star in 1D code   ###
!!! ###  M_loc is the mass coordinate inside the "donor" star (in Msun)            ###
!!! ###  flag_accel: 1: the function returns effective acceleration  in CGS units  ###
!!! ###              2 : returns in units of of "single star" acceleration         ###
!!! ##################################################################################      


       real*8  function effective_gravity(M1,M2,a,r_loc,M_loc, flag_accel)
  
         implicit none
  
         include "gravitytable.inc"        
        
!!! input and output variable
         real*8 :: M1, M2, M_loc
         real*8 :: r_loc, a
         integer flag_accel
           
!!! locally used variables

         real*8 grav_final
         real*8 :: q_local, q_total, r_rel_loc
          
         real*8 :: grav_found(2), r_rel_thresh(2)
          
         integer :: j,l,k, flag_use_analyt
         integer :: n_q, j_found
 
         real*8 :: table_f, g_sing, g_exterior
      
         !!scaling the input
         r_rel_loc=r_loc/a
         q_local=M_loc/M2
         q_total=M1/M2  

         flag_use_analyt=0

         
         !! Check if the tables are initialized and called for acceptable parameters
         
         if(flag_gtable_ini.ne.1) then
            write(*,*) "error, tables for effective gravity were not read"
            stop
         end if
         
         if(flag_accel.ne.1.and.flag_accel.ne.2) then
            write(*,*) flag_accel, "error, accelerations can be requested as 1 (in CGS units) or 2 (in relative units)"
            stop
         end if
         
         if (M_loc.le.0.d0.or.M_loc.gt.M1) then
            write (*,*) "error, local mass coordinate should be 0<M_local<=M1 donor star mass"
            stop
         end if
         
         if (r_loc.le.0.d0.or.r_loc.ge.a) then
            write (*,*) "error, radius of zone should be  0<=r<=a (orbital separation)", r_loc, a 
            stop
         end if

         if(flag_accel.eq.1) then
            effective_gravity = GT_accel*((M1/r_loc**2.d0) + &
                                 (-2.d0*(M1+M2)*r_loc/(3*a**3.d0)))
         else
            effective_gravity = 1+(-2.d0*(1.d0+q_total)*(r_rel_loc**3.d0)/(3.d0*q_total))

         end if

         if (q_total.gt.100000.d0) then
            write (*,*) "error, the acceptable mass ratio q=M1/M2 is 1e-6<=q_total<=1e5"
            stop
         end if

         if (q_total.lt. q_table(1)) then
!! return pre-assigned value as in a single star             
            write (*,*) "error, the acceptable mass ratio q=M1/M2 is 1e-6<=q_total<=1e5"
            stop
         end if


          !! Searching for the two mass ratio tables to interpolate 
    
         j_found =-1
         n_q=size(q_table)

         do j =1,n_q-1
            if (q_total.ge.q_table(j).and.&
                  q_total.le.q_table(j+1)) then
                j_found=j
               exit
            end if
         end do

         if (r_rel_loc.le.table_rel_r(j_found,1).or.&
            r_rel_loc.le.table_rel_r(j_found+1,1)) then
!! return pre-assigned value as in a single star             
            return
         end if 
         
          
          if(j_found.eq.-1) then
             write(*,*) j_found, q_total, "error, can not locate the table for this mass ratio"
          end if

          if (r_rel_loc.ge.table_rel_r(j_found,n_gtable)) then
             write (*,*)" Stars exceeds L3 equivalent volume, check input", r_rel_loc, table_rel_r(j_found,n_gtable)
             stop
          end if
          


          grav_found(1)=-1.
          grav_found(2)=-1. 
          do k=1,2  ! k is 1 and 2             
             do l=1,(n_gtable-1)                
                if (r_rel_loc.ge.&
                     table_rel_r(k+j_found-1,l).and.&
                     r_rel_loc.le.&
                     table_rel_r(k+j_found-1,l+1)) then
                   
                   grav_found(k) = linear_interpolant(table_rel_r(k+j_found-1,l),&
                        table_rel_r(k+j_found-1,l+1),table_grav(k+j_found-1,l),&
                        table_grav(k+j_found-1,l+1),r_rel_loc)
                   
                   exit                   
                end if
             end do
          end do

          if(grav_found(1).lt.0.or.grav_found(2).lt.0) then
             write(*,*) "did not find appropriate radius in the table", r_rel_loc, q_total
             stop
          end if
          
          grav_final = linear_interpolant(q_table(j_found),&
               q_table(j_found+1),grav_found(1),grav_found(2), q_total)

                                     
         g_exterior=2.d0*(q_local-q_total)&
         /((1.d0+q_total)*r_rel_loc**2.d0)


         if(flag_accel.eq.1) then
            effective_gravity=GT_accel*(M1+M2)*(grav_final+g_exterior)/(2.d0*a**2.d0)
         else
            effective_gravity=(grav_final+g_exterior)*((q_total+1.d0)*(r_loc/a)**2.d0)/2*q_total
         end if

         return      
       end function effective_gravity

!!##
!!## returns interpolated value y for x that is between x1 and x2
!!##       

       
       real*8 function linear_interpolant(x1,x2,y1,y2,x)
         implicit none
         real*8 :: x1, x2, y1, y2, x
         
         linear_interpolant=((x-x1)*(y2-y1)/(x2-x1))+y1
         
       end function linear_interpolant


