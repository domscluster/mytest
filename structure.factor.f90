!*************************************************************************!
!                       precision control                                 !
!*************************************************************************!
module Precision_Control
implicit none
integer, parameter :: single = selected_real_kind(p=6,r=37)    !single precision
integer, parameter :: double = selected_real_kind(p=15,r=307)  !double precision
integer, parameter :: quad   = selected_real_kind(p=33,r=4931) !quadruple precision
save
end module
!*************************************************************************!
!                       simulation control                                !
!*************************************************************************!
module Sim_Control
Use Precision_Control
implicit none
integer          :: Nstep , Natom , max_k
real(kind=double)::boxlen
real(kind=double) , allocatable, dimension(: , :)::Rx , Ry , Rz
end module
!*************************************************************************!
!                       program main                                      !
!*************************************************************************!
program main
Use Precision_Control
Use Sim_Control
implicit none
integer :: istep
call INITIALIZE
call READ_FILE
call STRUCTURE_FACTOR
end program main
!*************************************************************************!
! INITIALIZATION OF ALL VARIABLES                                         !
!*************************************************************************!
subroutine INITIALIZE
Use Precision_Control
Use Sim_Control
implicit none
Nstep     = 50000
Natom     = 960

max_k     = 1

allocate(Rx(Natom, Nstep) , Ry(Natom , Nstep) , Rz(Natom, Nstep))
open(unit = 100 , file= "/home/milan/DMSO.ROTATIONAL.VELOCITY/program.old/trjectory_500ps.gro")

end subroutine INITIALIZE
!*************************************************************************!
!   READ COORDINATE FILE                                                  !
!*************************************************************************!
subroutine READ_FILE
Use Precision_Control
Use Sim_Control
implicit none
integer :: ii , istep
character(len =100) :: gro_format = "(A8 , A7 , I5 , F8.3 ,F8.3 , F8.3)"
character :: junk1 , junk2
integer   :: junk3


do istep = 1 , Nstep
   read(100 , *) 
   read(100 , *)
 do ii = 1 , Natom
    read(100 , gro_format)  junk1 , junk2 , junk3 , Rx(ii , istep) , Ry(ii , istep) , Rz(ii , istep)
    read(100 , *)
    read(100 , *) 
    read(100 , *)
 end do
 read(100 , *) boxlen , boxlen , boxlen
end do
print*, Rx(Natom , Nstep)
print*, "reading done"
end subroutine READ_FILE

!*************************************************************************!
!   STRUCTURE FACTOR CALCULATION                                          !
!*************************************************************************!
subroutine STRUCTURE_FACTOR
Use Precision_Control
Use Sim_Control
implicit none
real(kind=double):: fact
real(kind=double)::mod_k , min_k_x , min_k_y , min_k_z , mod_mink , mod_maxk
real(kind=double)::bin_size 
integer          :: ii , jj , kk , ip ,  istep ,jstep ,  iatom , jatom , tot_bin , bin_no , ibin
real(kind=double)::kx , ky , kz
real(kind = double) , allocatable , dimension(:) ::sum1 , sum_sf
real(kind=double)::  sum2 , sum3 , sum4 , sum5 ,sum45 ,  temp , static_sf , norm , term
real(kind=double)::rr(3)
real(kind=double), allocatable,dimension(:) :: sum_cos , sum_sin

integer ::counta , corr_step , count2
allocate(sum_cos(Nstep) ,sum_sin(Nstep)  )
fact = 2*3.14
tot_bin = 100
min_k_x = fact/boxlen
min_k_y = fact/boxlen
min_k_z = fact/boxlen
mod_mink = sqrt(min_k_x**2+min_k_y**2+min_k_z**2)
mod_maxk = sqrt(3*max_k*max_k*min_k_x*min_k_x)
bin_size = (mod_maxk - mod_mink)/tot_bin

corr_step = int(Nstep/2)
allocate (sum1(corr_step), sum_sf(corr_step))
sum1  = 0
counta = 0
! do ii = 1 , max_k
!   do jj = 1 , max_k
!    do kk = 1 , max_k
        ii =2
        jj =2
        kk =2
       counta = counta+1
           
       kx = ii*min_k_x
       ky = jj*min_k_y
       kz = kk*min_k_z
       mod_k = sqrt(kx*kx +ky*ky +kz*kz)
       print*, mod_k
!       stop
       sum1  = 0 
       sum_cos = 0
       sum_sin = 0

       do istep = 1 , Nstep
!          sum_cos = 0
!          sum_sin = 0
        do ip = 1 , Natom
           rx(ip , istep) =rx(ip , istep) - boxlen*anint(rx(ip , istep)/ boxlen)
           ry(ip , istep) =ry(ip , istep) - boxlen*anint(ry(ip , istep)/ boxlen)
           rz(ip , istep) =rz(ip , istep) - boxlen*anint(rz(ip , istep)/ boxlen)
           temp = kx*rx(ip , istep) +ky*ry(ip , istep)+kz*rz(ip , istep)
           sum_cos(istep) = sum_cos(istep) +cos(temp)
           sum_sin(istep) = sum_sin(istep) +sin(temp)
        end do
        !write(200 , *) sum_cos(istep)
       end do
       do istep = 1 , Nstep
       write(200 , *) sum_cos(istep)
       end do
       norm = 0
       sum_sf = 0
       do jstep = 1 , Nstep
!          print*, sum_cos(jstep)
          norm = norm+sum_cos(jstep)**2
       end do
       do istep = 1 , corr_step
          do jstep = 1 , Nstep
             if((istep+jstep) .gt. Nstep) exit
             term = sum_cos(jstep) * sum_cos(jstep+istep)
             sum_sf(istep) = sum_sf(istep)+term
          end do
          write(300, *) istep*0.01 , sum_sf(istep)/norm
       end do
!    end do
!   end do
!  end do
end subroutine STRUCTURE_FACTOR
