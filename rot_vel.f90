
!********************************************************************************!
!     angular velocity and correlation                                           !
!********************************************************************************!

!**********************************************************************************************!
!MODULE1: PRECISION CONTROL OF THE PROGRAM IS INCORPORATED HERE                                !
!**********************************************************************************************!
      module Precision_Control
      implicit none
      integer, parameter :: single = selected_real_kind(p=6,r=37)    !single precision
      integer, parameter :: double = selected_real_kind(p=15,r=37)  !double precision
      integer, parameter :: quad   = selected_real_kind(p=33,r=4931) !quadruple precision
      save
      end module

!***********************************************************************************************!
!MODULE2: PARAMETERS OF THE PROGRAM                                                             !
!***********************************************************************************************!
      module Parameter_Control
      Use Precision_Control

      integer::Ndmso , Natom , Nstep
      real(kind=double)::boxl(3) , dt


      end module
!***********************************************************************************************!
!MODULE3: ALL THE READ    VARIABLES                                                             !
!***********************************************************************************************!
      module Read_Variables
      Use Precision_Control

      real(kind=double),allocatable,dimension(:,:,:):: rs , ro , rc2
      
      end module
!***********************************************************************************************!
!MODULE3: ALL THE SIMULATION    VARIABLES                                                       !
!***********************************************************************************************!
      module Sim_Variables
      Use Precision_Control
      real(kind=double),allocatable,dimension(:,:)::phi , theta , psi
     
      end module
!***********************************************************************************************!
!PROGRAM MAIN                                                                                   !
!***********************************************************************************************!
      program main
      Use Precision_Control
      Use Parameter_Control
      Use Read_Variables
      Use Sim_Variables
      implicit none
      call INITIALIZE
      call READ_DATA
      call axis_angle
      end program main


!************************************************************************************************!
!INITIALIZE                                                                                      !
!************************************************************************************************!


      subroutine INITIALIZE
      Use Precision_Control
      Use Parameter_Control
      Use Read_Variables
      Use Sim_Variables
      implicit none
      Nstep = 8000
      Ndmso    =960
      Natom = 4*Ndmso
      dt = 0.01
      allocate(rs(Ndmso , Nstep ,3) , ro(Ndmso , Nstep , 3), rc2(Ndmso , Nstep ,3) )
      allocate(phi(Ndmso , Nstep) , theta(Ndmso , Nstep) , psi(Ndmso , Nstep))
open(unit = 1000 , file = "trjectory_500ps.gro" , status = "old")
open(unit = 300 , file ="angular_velocity.dat"  , action = "write")
      end subroutine INITIALIZE
!************************************************************************************************!
!READ DATA                                                                                       !
!************************************************************************************************!
      subroutine READ_DATA
      Use Precision_Control
      Use Parameter_Control
      Use Read_Variables
      Use Sim_Variables
      implicit none
      integer :: istep , ip
      character :: junk1 , junk2
      integer   :: junk3
     character(len =100) :: gro_format = "(A8 , A7 , I5 , F8.3 ,F8.3 , F8.3)"
      do istep = 1 , Nstep
         read(1000 , *)
         read(1000 , *)
         do ip = 1 , Ndmso
         
            read(1000 , gro_format) junk1 , junk2 , junk3 , rs(ip , istep , 1) , rs(ip , istep , 2), rs(ip , istep , 3)
            read(1000,  *) 
            read(1000 , *) junk1 , junk2 , junk3 , ro(ip , istep , 1) , ro(ip , istep , 2), ro(ip , istep , 3)
            read(1000 , *) junk1 , junk2 , junk3 , rc2(ip , istep , 1) , rc2(ip , istep , 2), rc2(ip , istep , 3)
         end do
         read(1000 , *) boxl(1) , boxl(2) , boxl(3)
      end do
      end subroutine READ_DATA



!****************************************************************************************************!
!moving ahead to calculate bond vector and its rotation in axis angle representation                 !
!****************************************************************************************************!
  

      subroutine axis_angle
      Use Precision_Control
      Use Parameter_Control
      Use Read_Variables
      Use Sim_Variables
      implicit none 
      real(kind=double)::cosangle , angle , ax , ay , az , rso_old(3) , rc2_new(3) , rso_new(3) , rso_ref(3)
      real(kind=double)::mod_a , mod_rso_old , mod_rso_new , arb_vec(3) , mod_rc2_new , mod_arbvec
      real(kind=double):: qx ,qy , qz , qw
      integer ::ii , jj , ip ,istep
      real(kind = double):: omgx, omgy , omgz , dphi , dtheta , dpsi
      rso_ref(1) =0
      rso_ref(2) =0
      rso_ref(3) =1
      rso_old(:) = rso_ref(:)
      do ii = 1 , Nstep
        do ip = 1 , Ndmso
           !rso_old(:) = rso_ref(:)
           rso_new(:) = ro(ip, ii , :) -rs(ip, ii , :) 
           !rso_old(:) = rso_old(:) - boxl(:)*anint(rso_old(:)/boxl(:))
           rso_new(:) = rso_new(:) - boxl(:)*anint(rso_new(:)/boxl(:))
           mod_rso_new = sqrt(rso_new(1)**2+rso_new(2)**2+rso_new(3)**2)
           rso_new = rso_new/mod_rso_new
           ax = rso_old(2)*rso_new(3) -rso_old(3)*rso_new(2)
           ay = -(rso_old(1)*rso_new(3) -rso_old(3)*rso_new(1) )
           az = rso_old(1)*rso_new(2) -rso_old(2)*rso_new(1)
           mod_a = sqrt(ax*ax +ay*ay +az*az)
           ax = ax/mod_a
           ay = ay/mod_a
           az = az/mod_a
           
           cosangle = dot_product(rso_old , rso_new)
           write(600 , *) ax
           
           angle    = acos(cosangle)
           if(cosangle .eq. 1) then
             qx = 0
             qy = 0
             qz = 0
             qw = 1
           else if(cosangle .eq. -1.0) then
             rc2_new(:) = rc2(ip ,ii ,:)-rs(ip, ii , :)
             rc2_new(:) = rc2_new(:) - boxl(:)*anint(rc2_new(:)/boxl(:))
             mod_rc2_new = sqrt(dot_product(rc2_new , rc2_new))
             arb_vec(1) = rc2_new(2)*rso_new(3) -rc2_new(3)*rso_new(2)
             arb_vec(2) = -(rc2_new(1)*rso_new(3) -rc2_new(3)*rso_new(1))
             arb_vec(3) = (rc2_new(1)*rso_new(2) -rc2_new(2)*rso_new(1))
             mod_arbvec = sqrt(dot_product(arb_vec , arb_vec))
             qx = arb_vec(1)/mod_arbvec
             qy = arb_vec(2)/mod_arbvec
             qz = arb_vec(3)/mod_arbvec
             qw = 0


           else
           qx = ax*sin(angle/2)
           qy = ay*sin(angle/2)
           qz = az*sin(angle/2)
           qw = cos(angle/2)
           end if
           write(400 , *) qx   , qy , qz , qw 
           phi(ip , ii) = atan2(2*(qw*qx+qy*qz) ,(1-2*(qx**2+qy**2)))
           theta (ip , ii)= asin(2*(qw*qy-qz*qx))
           psi (ip ,ii)  = atan2(2*(qw*qz+qx*qy),(1-2*(qy**2+qz**2)))
           write(500 ,*)phi(ip  , ii)    
              !forward difference to euler angles 
              
                                
        end do
        print*, "done for=" , ii
      end do   
      do istep = 1 , Nstep -1
       do ip = 1 ,Ndmso
         dphi = phi(ip , istep+1 )-phi(ip ,istep)
         dtheta = theta(ip , istep+1)-theta(ip , istep)
         dpsi   = psi(ip , istep+1)-psi(ip ,istep)
         dphi = dphi -3.14*anint(dphi/3.14)
         dtheta = dtheta -1.57*anint(dtheta/1.57)
         dpsi = dpsi -3.14*anint(dpsi/3.14)
         write(5000 , *)dphi , dtheta , dpsi
         dphi = dphi/dt
         dtheta = dtheta/dt
         dpsi = dpsi/dt
         omgx = dpsi*cos(theta(ip , istep))*cos(phi(ip , istep)) - dtheta*sin(phi(ip , istep))
         omgy = dphi -dpsi*sin(theta(ip , istep))
         omgz = dpsi*cos(theta(ip , istep))*sin(phi(ip , istep)) +dtheta*cos(phi(ip , istep)) 
         write(300 , *)  omgx , omgy , omgz
       end do 
         
         print*, "2nd part done for =" , istep
      end do
      end subroutine axis_angle
      


