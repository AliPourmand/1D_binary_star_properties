!!! #######################################################################
!!! ###  reading pre-calacuated table with effective acceleration       ###
!!! ###  by default, the files are assumed to be in ./binary directory  ###
!!! #######################################################################      

      subroutine gravity_table_ini(M1,M2,r_lin_val,r_lout_val)

         implicit none
        
         include "gravitytable.inc"        
        
!!! internal variables                 
         integer :: j,l
         integer :: n_q  !! this is the number of table that is read
        
         real*8 :: lgq_table(c_gtable)

         character*255 :: file_name
                        
!!! dummy properties are those getting read from the tables        
         real*8 :: xi_arr_dummy
         real*8 :: q_dummy(c_gtable),r_l2(c_gtable)
         real*8 :: r_l3(c_gtable),a_lplane2(c_gtable)
         real*8 :: a_lplane3(c_gtable), g_lplane2(c_gtable), g_lplane3(c_gtable)
!!! inputs
         real*8 :: M1, M2, q_input
!!! output
         real*8 :: r_lin_val, r_lout_val

        
         flag_gtable_ini=0

         q_input = M1/M2   
         write(*,*) "q_tot_input=",q_input
         !!size gives length of array
         n_q=c_gtable

         do j=1,16

            lgq_table(j)=-6.d0+(j-1)*0.25d0
            q_table(j)=10.d0**lgq_table(j)

         end do

         do j=17,96

            lgq_table(j)=-2.d0+(j-16-1)*0.05d0
            q_table(j)=10.d0**lgq_table(j)

         end do

         do j=97,n_q

            lgq_table(j)=2.d0+(j-96-1)*0.25d0
            q_table(j)=10.d0**lgq_table(j)

         end do        
    

         do j=1,n_q
            if (lgq_table(j).lt.0) then 
               write(file_name,'(A49,F6.3,A4)') "binary_tables_compact_ss_units/binary_props_logq="&
               ,lgq_table(j),".dat"
            end if
            if (lgq_table(j).ge.-0.001d0) then 
               write(file_name,'(A49,F5.3,A4)') "binary_tables_compact_ss_units/binary_props_logq="&
               ,lgq_table(j),".dat"
            end if

            ! write (*,*) file_name
           
            open (unit = 23, file = file_name)

            do l=1,n_gtable
!!!! the following are being extracted for each mass ratio
! column 1: volume equivalent radius of shell divided by orbital separation  
! column 2: relative potential of shell  
! column 3: relative avg acceleration on shell   
! column 4: relative avg acceleration on L plane cross section with shell  
! column 5: area of L plane cross section with shell
! column 6: area of shell

            
               read (23,*) table_rel_r(j,l), xi_arr_dummy, table_grav(j,l),&
                   g_Lplane(j,l), A_Lplane(j,l), A_shell(j,l)
            end do
            close (23)


         end do
 

!!!!!!!!!!! reading the lagrange point properties from the file l1l2l3_props.dat

 

         write(file_name,'(A16)') "l1l2l3_props.dat"
           
         open (unit = 24, file = file_name)

         do j=1,n_q

! The file l1l2l3_props.dat includes properties of the three lagrange points for these 109 mass ratios. The properties are listed below:

! column 1: mass ratio of system (donor to companion)  
! column 2: volume equivalent radius of shell passing L1  
! column 3: volume equivalent radius of shell passing L2  
! column 4: volume equivalent radius of shell passing L3  
! column 5: area of L plane cross section with shell passing L2  
! column 6: area of L plane cross section with shell passing L3  
! column 7: avg acceleration on L plane cross section with shell passing L2  
! column 8: avg acceleration on L plane cross section with shell passing L3     
            
            read (24,*) q_dummy(j),r_l1(j),r_l2(j),r_l3(j),a_lplane2(j),&
                        a_lplane3(j), g_lplane2(j), g_lplane3(j)


!!!!!!!!!!!! outer lagrange point is l2 for q <= 1
            if (q_dummy(j).le.0.999d0) then
               r_lout(j) = r_l2(j)
            else
               r_lout(j) = r_l3(j)
            end if

            ! write (*,*) r_lout(j)


         end do

! find input mass ratio's corresponding rlin and rlout   
         
         do j =1,n_q-1
            if (q_input.ge.q_table(j).and.&
                  q_input.lt.q_table(j+1)) then
               
               r_lin_val = linear_interpolant(q_table(j),q_table(j+1),&
                  r_l1(j),r_l1(j+1),q_input)

               r_lout_val = linear_interpolant(q_table(j),q_table(j+1),&
                  r_lout(j),r_lout(j+1),q_input)

               write(*,*) "called for mass ratio =",q_input
    
               exit
            end if
         end do        




         close (24)



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

         real*8 grav_final, g_shell_minus_l1_1, g_shell_minus_l1_2
         real*8 :: q_local, q_total, r_rel_loc
          
         real*8 :: grav_found(2), r_rel_thresh(2)
          
         integer :: j,l,k, flag_use_analyt
         integer :: n_q, j_found, found_index, l_found(2)
 
         real*8 :: table_f, g_sing, g_exterior
         real*8 :: r_lin_val, r_lout_val

 
         n_q=size(q_table)

      
         !!scaling the input
         r_rel_loc=r_loc/a
         q_local=M_loc/M2
         q_total=M1/M2  

! find input mass ratio's corresponding rlin and rlout   

         do j =1,n_q-1
            if (q_total.ge.q_table(j).and.&
                  q_total.lt.q_table(j+1)) then
               
               found_index = j
               r_lin_val = linear_interpolant(q_table(j),q_table(j+1),&
                  r_l1(j),r_l1(j+1),q_total)

               r_lout_val = linear_interpolant(q_table(j),q_table(j+1),&
                  r_lout(j),r_lout(j+1),q_total)
                  
          
               exit
            end if
         end do  

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


!!!!!!!!!!!!!! the following equations are from the article, equation 23

         if(flag_accel.eq.1) then
            effective_gravity = GT_accel*(1.d0+(-2.d0*(1.d0+&
               q_total)*(r_rel_loc**3.d0)/(3.d0*q_local)))*M_loc/r_loc**2.d0
         else
            effective_gravity = (1.d0+(-2.d0*(1.d0+&
               q_total)*(r_rel_loc**3.d0)/(3.d0*q_local)))*q_local/q_total

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

         do j =1,n_q-1
            if (q_total.ge.q_table(j).and.&
                  q_total.lt.q_table(j+1)) then
               j_found=j
            
               exit
            end if
         end do



!!!!!!!!!!!!! uses analytical solution for values smaller than what is in the tables

         if (r_rel_loc.le.table_rel_r(j_found,1).or.&
            r_rel_loc.le.table_rel_r(j_found+1,1)) then
!! return pre-assigned value as in a single star     
            flag_use_analyt =1 
            write (*,*)"return"
            return
         end if 
         
          
         if(j_found.eq.-1) then
            write(*,*) j_found, q_total, "error, can not locate the table for this mass ratio"
         end if

         if (r_rel_loc.gt.r_lout_val.or.abs(r_rel_loc-r_lout_val).le.1.e-9) then
            write (*,*)" Stars exceeds outer Lagrange equivalent volume, check input", r_rel_loc !, table_rel_r(j_found,n_gtable) r_lout_val
            stop
         end if
          


         grav_found(1)=-1.
         grav_found(2)=-1. 
         do k=1,2  ! k is 1 and 2             
            do l=1,(n_gtable-1)                
               if (r_rel_loc.gt.&
                     table_rel_r(k+j_found-1,l).and.&
                     r_rel_loc.lt.&
                     table_rel_r(k+j_found-1,l+1)) then

                  l_found(k) = l
                  grav_found(k) = linear_interpolant(table_rel_r(k+j_found-1,l),&
                        table_rel_r(k+j_found-1,l+1),table_grav(k+j_found-1,l),&
                        table_grav(k+j_found-1,l+1),r_rel_loc)

                  if (l_found(k).gt.500) then 

                     if (log10(q_total).gt.2.5d0) then
                        write (*,*) "error, mass ratio q should be less than log(q)<2.5 for beyond L1 radii"
                        stop
                     end if 
                     !! if radius bigger than L1, L1plane gets reduced from g_bin

                     g_shell_minus_l1_1 = (table_grav(k+j_found-1,l)*A_shell(k+j_found-1,l)- &
                     g_Lplane(k+j_found-1,l)*A_Lplane(k+j_found-1,l))/(A_shell(k+j_found-1,l)-&
                        A_Lplane(k+j_found-1,l))

                     g_shell_minus_l1_2 = (table_grav(k+j_found-1,l+1)*A_shell(k+j_found-1,l+1)- &
                     g_Lplane(k+j_found-1,l+1)*A_Lplane(k+j_found-1,l+1))/(A_shell(k+j_found-1,l+1)-&
                        A_Lplane(k+j_found-1,l+1))

                     grav_found(k) = linear_interpolant(table_rel_r(k+j_found-1,l),&
                           table_rel_r(k+j_found-1,l+1),g_shell_minus_l1_1,&
                           g_shell_minus_l1_2,r_rel_loc)
                  end if 

                  exit                   
               end if
            end do
         end do

         if(grav_found(1).lt.0.or.grav_found(2).lt.0) then
             write(*,*) "did not find appropriate radius in the table", r_rel_loc, &
             table_rel_r(j_found,600),&
             l_found(1),&
             table_rel_r(j_found+1,l_found(2)),&
             l_found(2),&
              q_total
             stop
         end if


         grav_final = linear_interpolant(q_table(j_found),&
               q_table(j_found+1),grav_found(1),grav_found(2), q_total)

         !! equation 32 of article, in terms of singular gravity
                                     
         g_exterior=(q_local-q_total)/q_total


         if(flag_use_analyt.eq.0) then

            if(flag_accel.eq.1) then
               effective_gravity=GT_accel*(M1)*(grav_final+g_exterior)/(r_loc**2.d0)
            else
               effective_gravity=(grav_final+g_exterior)


            end if

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


