!!!	size of (each) precalculated gravity arrays.
      integer n_gtable
      parameter (n_gtable=600)

!!!   maximum number of precalculated gravity arrays, used to allocate pre-pread arrays
      integer c_gtable
      parameter (c_gtable=109)

!!! functions
      real*8 linear_interpolant


!!! external variables         
      integer flag_gtable_ini
      real*8 :: table_rel_r(c_gtable,n_gtable), table_grav(c_gtable,n_gtable), q_table(c_gtable)
      real*8 :: g_Lplane(c_gtable,n_gtable), A_Lplane(c_gtable,n_gtable), A_shell(c_gtable,n_gtable)      
      real*8 :: r_l1(c_gtable),r_lout(c_gtable)
      !real*8 :: r_lout_val, r_lin_val


        
      common /gravity_tables/ table_rel_r, A_Lplane, A_shell, g_Lplane, table_grav, q_table,r_l1, r_lout, flag_gtable_ini 
      ! , found_index
      

!
!     physical constants
!
      real*8    GT_gconst,GT_msun, GT_rsun, GT_accel

      parameter (GT_gconst = 6.6743e-8           ) ! gravity const
      parameter (GT_msun   = 1.989e33            ) ! Mass of the Sun
      parameter (GT_rsun   = 6.9598e10           ) ! Radius of the Sun
      parameter (GT_accel  = GT_gconst*GT_msun/GT_rsun/GT_rsun ) ! conversing acceleration in CGS units

      
