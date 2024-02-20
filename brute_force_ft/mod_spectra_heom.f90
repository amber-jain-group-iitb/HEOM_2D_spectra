Module mod_spectra_heom
!use matmul_lapack
!! Reproducing figure 5 of Strumpfer Schulten, JCTC 8, 2808 (2012)
!! Chen, Zheng, Shi, Tan, JCP 131, 094502 (2009)
implicit none
real*8, parameter :: clight=137.d0,av=6.0221367D23,hbar=1.d0!1.05457266D-34
real*8, parameter :: mass_h=1.007825d0,kb=1.d0!1.3806503d-23
real*8, parameter :: pascal_to_atm=9.86923267d-6,kcal_to_J=4184.d0
real*8, parameter :: wave_to_J=1.98644568d-23
real*8, parameter :: amu2kg=1.66053892d-27
real*8, parameter :: au2kg=9.10938291d-31,au2J=4.35974417d-18,au2s=2.418884326505d-17
complex*16,parameter :: iota = (0.d0,1.d0)
real*8 pi

!! HEOM variables
integer LL,KK
integer nn_tot
integer,allocatable :: nn(:,:,:)
integer,allocatable :: map_nplus(:,:,:),map_nneg(:,:,:),map_sum(:)
integer,allocatable :: zero(:,:)
complex*16,allocatable :: rho(:,:,:,:)
complex*16,allocatable :: c_matsubara(:),omg_matsubara(:)
real*8 tolerance
real*8 sum_c_over_nu

!! Output
complex*16 dip_mom_corr

!! 2d spectra
real*8 t1,t2,t3,dt_t1
real*8 t1_max,t3_max
integer nst_t1,nst_t3,nst_avg
complex*16,allocatable :: mu_t3(:,:,:)

!! System specific
integer nquant
real*8 gamma,eta,temperature
complex*16,allocatable :: Hamil_e(:,:),dip_moment_pos(:,:),dip_moment_neg(:,:)

!! Evolution
integer nsteps
real*8 dt,tim_tot,curr_time

!! Parallelization
integer iparallel,iflow,iproc,iwait,ifolder

!! Misc
real*8 tim_ind

contains
!---------------------------------------------------------- 
!---------------------------------------------------------- 

subroutine setup
  implicit none
  integer i,j
  character st_ch

  pi=dacos(-1.d0)
  open(10,file="spectra_heom.inp")
  read(10,*) iflow
  read(10,*) iproc
  read(10,*) iparallel
  read(10,*) iwait
  read(10,*) nquant
  read(10,*) LL
  read(10,*) KK
  read(10,*) t2
  read(10,*) t1_max
  read(10,*) t3_max
  read(10,*) nst_avg
  read(10,*) tolerance
  read(10,*) dt
  read(10,*) dt_t1
  read(10,*) tim_tot
  read(10,*) st_ch
  close(10)
  !----------------------------------------------------------

  if(st_ch.ne.'x') then
    write(6,*) "problem in reading input file"
    stop
  endif

  !---------------------------------------------------------- 

  dt=dt/au2s
  dt_t1=dt_t1/au2s
  t2=t2/au2s
  t1_max=t1_max/au2s
  t3_max=t3_max/au2s
  tim_tot=tim_tot/au2s

  !nn_tot=(factorial(LL+nquant*KK)/factorial(LL))/factorial(nquant*KK)
  call compute_nn_tot
  nsteps=nint(tim_tot/dt)
  nst_t1=nint(t1_max/dt_t1)+1
  nst_t3=nint(t3_max/dt)+1

  i=nst_t1/nst_avg+1
  j=nst_t3/nst_avg+1

  allocate(dip_moment_pos(nquant,nquant),dip_moment_neg(nquant,nquant))
  allocate(mu_t3(nquant,nquant,6))
  allocate(nn(nquant,0:KK,nn_tot),map_nplus(nquant,0:KK,nn_tot),map_nneg(nquant,0:KK,nn_tot))
  allocate(map_sum(0:LL),zero(nn_tot,6))
  allocate(rho(nquant,nquant,nn_tot,6))
  allocate(c_matsubara(0:KK),omg_matsubara(0:KK))
  allocate(Hamil_e(nquant,nquant))

  call compute_nn
  call compute_map

end subroutine setup
!---------------------------------------------------------- 

subroutine main
  implicit none
  integer clck_counts_beg, clck_counts_end, clck_rate
  real*8 t1_cpu,t2_cpu
  integer i,j,nst_t2
  complex*16 rho_sav(nquant,nquant,nn_tot,6)

  call system_clock ( clck_counts_beg, clck_rate )
  call cpu_time(t1_cpu)

  open(100,file="components_correlation.out")
  open(101,file="dip_mom_corr.out")
  open(102,file="spectra.out")
  call setup_parameters
  call init_cond
  call strike_lightning(0)  !! t=0
  rho_sav=rho
  t1=0.d0

  !! parallelizing over t1
!  if(iflow==2) then
!    open(10,file="ifolder.inp")
!    read(10,*) ifolder
!    close(10)
!    t1_max=t1_max/real(iparallel)
!    nst_t1=nint(t1_max/dt_t1)+1
!    if(ifolder>1) then
!      do i=1,ifolder-1
!        call evolve(nst_t1-1)
!        t1=t1+t1_max
!        call evolve(nst_avg)
!        t1=t1+dt_t1*nst_avg
!      enddo
!    endif
!  endif
!  rho_sav=rho
!
!  do i=1,nst_t1
!    call strike_lightning(1)  !! second laser hit
!    nst_t2=nint(t2/dt)
!    call evolve(nst_t2)    !! evolve for t2
!    call strike_lightning(2)  !! third laser hit
!    do j=1,nst_t3
!      t3=(j-1)*t3_max/real(nst_t3-1)
!      mu_t3=rho(:,:,1,:)
!      if(mod(i,nst_avg)==1.and.mod(j,nst_avg)==1) call compute_response(i,j) !! measurement
!      !write(100,'(5f15.7)') t1*au2s*1.d15,t3*au2s*1.d15,rho(2,1,1)!dip_mom_corr
!      call evolve(1)     !! evolve for t3
!    enddo
!    if(mod(i,nst_avg)==1)write(100,*)
!    if(mod(i,nst_avg)==1)write(101,*)
!
!    rho=rho_sav
!    call evolve(1)       !! evolve for t1
!    t1=t1+dt_t1!(i-1)*t1_max/real(nst_t1-1)
!    rho_sav=rho
!
!  enddo

  call fft_dip_mom_corr
  close(100);close(101)

  call system_clock ( clck_counts_end, clck_rate )
  call cpu_time(t2_cpu)
  write(6,*) "total time=",(clck_counts_end - clck_counts_beg) / real(clck_rate),t2_cpu-t1_cpu
  write(6,*) "tim_index=",tim_ind


  close(102)

end subroutine main
!---------------------------------------------------------- 

subroutine setup_parameters
  implicit none
  integer i,j,k
  real*8 tmp,omg,cc

  !JJ=1.0!*wave_to_J
  !gamma=50.d0/1.d-12*au2s
  !eta=2*198.d0/hbar*wave_to_J/au2J
  !gamma=4.47d0/1.d-12*au2s
  gamma=0.d0/100.d-15*au2s
  eta=120.d0/hbar*wave_to_J/au2J
  temperature=77.d0*1.38064852d-23/au2J

  omg_matsubara(0)=gamma 
  !c_matsubara(0)=eta*gamma/2.d0* (cotan(hbar*gamma/(2.d0*kb*temperature))-iota)
  c_matsubara(0)=eta*kb*temperature-iota*eta*gamma/2.d0
  do k=1,KK
    omg_matsubara(k)=2*k*pi*kb*temperature/hbar
    c_matsubara(k)=2*eta*gamma*kb*temperature/hbar * omg_matsubara(k)/(omg_matsubara(k)**2-gamma**2)
    !c_matsubara(k)=4*(k-1)*pi*eta*gamma/((2*(k-1)*pi)**2-(hbar*gamma/(kb*temperature))**2)
  enddo

  sum_c_over_nu=0.d0
  do k=KK+1,200
    omg=2*k*pi*kb*temperature/hbar
    cc=2*eta*gamma*kb*temperature/hbar * omg/(omg**2-gamma**2)
    sum_c_over_nu=sum_c_over_nu+cc/omg
  enddo

!write(6,*) eta*kb*temperature/gamma-sum(c_matsubara/omg_matsubara),sum_c_over_nu
!stop

  dip_moment_pos=0.d0
  dip_moment_pos(2,1)=1.d0
  dip_moment_pos(3,1)=-0.2d0
  dip_moment_neg=0.d0
  dip_moment_neg(1,2)=1.d0
  dip_moment_neg(1,3)=-0.2d0

  Hamil_e=0.d0
  Hamil_e(2,2)=-50.d0;Hamil_e(2,3)=100.d0
  Hamil_e(3,2)=100.d0;Hamil_e(3,3)=50.d0
  Hamil_e=Hamil_e*wave_to_J/au2J

write(6,*) "Parameters Set ..."

end subroutine setup_parameters
!-----------------------------------------------------------------  

subroutine init_cond
  implicit none
  integer i

  rho=0.d0 !! XXX !!
  rho(1,1,1,:)=1.d0
  !rho(:,:,1)=matmul_lap(dip_moment,rho(:,:,1))

  curr_time=0.d0
  dip_mom_corr=0.d0

  tim_ind=0.d0

write(6,*) "Intitial Conditions set ... "

end subroutine init_cond
!-----------------------------------------------------------------  

subroutine strike_lightning(iflag)
  implicit none
  integer,intent(in) :: iflag
  integer i

 !! dipole moment on the right --> 0
 !! dipole moment on the left  --> 1

  if(iflag==0) then   !! t=0
    call mult(rho(:,:,:,1),dip_moment_neg,0)
    call mult(rho(:,:,:,2),dip_moment_neg,0)
    call mult(rho(:,:,:,3),dip_moment_neg,0)
    call mult(rho(:,:,:,4),dip_moment_pos,1)
    call mult(rho(:,:,:,5),dip_moment_pos,1)
    call mult(rho(:,:,:,6),dip_moment_pos,1)
  endif

  if(iflag==1) then   !! t=t1
    call mult(rho(:,:,:,1),dip_moment_pos,1)
    call mult(rho(:,:,:,2),dip_moment_pos,0)
    call mult(rho(:,:,:,3),dip_moment_pos,1)
    call mult(rho(:,:,:,4),dip_moment_neg,0)
    call mult(rho(:,:,:,5),dip_moment_neg,1)
    call mult(rho(:,:,:,6),dip_moment_neg,0)
  endif

  if(iflag==2) then   !! t=t2
    call mult(rho(:,:,:,1),dip_moment_pos,0)
    call mult(rho(:,:,:,2),dip_moment_pos,1)
    call mult(rho(:,:,:,3),dip_moment_pos,1)
    call mult(rho(:,:,:,4),dip_moment_pos,0)
    call mult(rho(:,:,:,5),dip_moment_pos,1)
    call mult(rho(:,:,:,6),dip_moment_pos,1)
  endif

end subroutine strike_lightning
!-----------------------------------------------------------------  

subroutine mult(rho,mu,iflag)
  implicit none
  complex*16,intent(inout)::rho(nquant,nquant,nn_tot)
  complex*16,intent(in)::mu(nquant,nquant)
  integer,intent(in)::iflag
  integer i

  if(iflag==0) then !!  dipole moment on the right
    do i=1,nn_tot
      rho(:,:,i)=matmul(rho(:,:,i),mu)
    enddo
  endif

  if(iflag==1) then !!  dipole moment on the left
    do i=1,nn_tot
      rho(:,:,i)=matmul(mu,rho(:,:,i))
    enddo
  endif



end subroutine mult
!-----------------------------------------------------------------  

subroutine evolve(nsteps)
  implicit none
  integer,intent(in)::nsteps
  integer i

  do i=1,nsteps
    zero=0!if(mod(i,10)==1)call filter
    call rk4(rho,dt)
    curr_time=curr_time+dt
  enddo

end subroutine evolve
!-----------------------------------------------------------------  

subroutine compute_response(i,j)
  implicit none
  integer,intent(in) :: i,j
  integer k,l
  complex*16 mat(nquant,nquant),response(6)

  write(100,'(2f15.7$)') t1,t3
  do l=1,6
    mat=matmul(dip_moment_neg,mu_t3(:,:,l))

    response(l)=0.d0
    do k=1,nquant
      response(l)=response(l)+mat(k,k)
    enddo
    write(100,'(2f15.7$)') response(l)
  enddo
  write(100,*)

  write(101,'(6f15.7)') t1,t3,response(1)+response(2)-response(3),response(4)+response(5)-response(6)

end subroutine compute_response
!-----------------------------------------------------------------  

subroutine rk4(rho,dt)
  implicit none
  complex*16,intent(inout)::rho(nquant,nquant,nn_tot,6)
  real*8,intent(in)::dt
  complex*16,dimension(nquant,nquant,nn_tot) :: vec,k1,k2,k3,k4
  integer i

!write(6,*) "In RK4 ..."

  !call compute_deriv(rho,k1)
  !rho=rho+dt*k1

  do i=1,6
    vec=rho(:,:,:,i)
    call compute_deriv(vec,k1,zero(:,i))
    call compute_deriv(vec+0.5*dt*k1,k2,zero(:,i))
    call compute_deriv(vec+0.5*dt*k2,k3,zero(:,i))
    call compute_deriv(vec+dt*k3,k4,zero(:,i))

    rho(:,:,:,i)=rho(:,:,:,i)+dt/6.d0*(k1+2*k2+2*k3+k4)
  enddo

end subroutine rk4
!-----------------------------------------------------------------  

subroutine compute_deriv(rho,drho_dt,zero)
  implicit none
  integer,intent(in) :: zero(nn_tot)
  complex*16,intent(in)::rho(nquant,nquant,nn_tot)
  complex*16,intent(out)::drho_dt(nquant,nquant,nn_tot)
  complex*16 tmp(nquant,nquant)
  integer n
  integer tid,OMP_GET_THREAD_NUM

!  call omp_set_num_threads(10)

  !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(TID,tmp)
!  TID = OMP_GET_THREAD_NUM()
!  PRINT *, 'Hello World from thread = ', TID
  !$OMP DO SCHEDULE(STATIC)
  do n=1,nn_tot
    call compute_deriv_n(n,rho,tmp,zero)
    drho_dt(:,:,n)=tmp
  enddo
  !$OMP END DO NOWAIT
  !$OMP END PARALLEL

end subroutine compute_deriv
!-----------------------------------------------------------------  

subroutine compute_deriv_n(n,rho,drho_dt,zero)
  !! implicit none
  !! Eq. 15, JCP 131, 094502 (2009)
  integer,intent(in) :: n
  integer,intent(in) :: zero(nn_tot)
  complex*16,intent(in)::rho(nquant,nquant,nn_tot)
  complex*16,intent(out)::drho_dt(nquant,nquant)
  complex*16 tmp(nquant,nquant),mat_tmp(nquant,nquant),mat_tmp2(nquant,nquant)
  integer m,k,nplus,nminus
  integer nvec(nquant,0:KK)!,nvec_plus(nquant,KK),nvec_neg(nquant,KK)

    m=n
    call index(m,nvec,0)
    tmp=-iota*commute(Hamil_e,rho(:,:,n))/hbar

    if(zero(n)==0) then !! matrix at n is not filtered out

      do m=2,nquant
        do k=0,KK
          tmp=tmp-nvec(m,k)*omg_matsubara(k)*rho(:,:,n)
        enddo
      enddo

      !do m=2,nquant
      !  mat_tmp=0.d0;mat_tmp(m,m)=1.d0 !! matrix |m><m|
      !  tmp=tmp-(eta*kb*temperature/gamma-real(sum(c_matsubara/omg_matsubara)))*commute(mat_tmp,commute(mat_tmp,rho(:,:,n)))
      !  !tmp=tmp-(sum_c_over_nu)*commute(mat_tmp,commute(mat_tmp,rho(:,:,n)))
      !enddo
    endif

    do m=2,nquant
      mat_tmp=0.d0;mat_tmp(m,m)=1.d0 !! matrix |m><m|
      mat_tmp2=0.d0
      do k=0,KK
!        nvec_plus=nvec;nvec_plus(m,k)=nvec_plus(m,k)+1
!        call index(nplus,nvec_plus,1)
        nplus=map_nplus(m,k,n)

        if(nplus>0.and.nplus<=nn_tot) then
          if(zero(nplus)==0) tmp=tmp - iota*commute(mat_tmp,rho(:,:,nplus))!* sqrt((nvec(m,k)+1.d0)*abs(c_matsubara(k)))
            !mat_tmp2=mat_tmp2+rho(:,:,nplus)! * sqrt((nvec(m,k)+1.d0)*abs(c_matsubara(k)))
        endif
      enddo
      !tmp=tmp - iota*commute(mat_tmp,mat_tmp2)
    enddo

    do m=2,nquant
      mat_tmp=0.d0;mat_tmp(m,m)=1.d0 !! matrix |m><m|
      do k=0,KK
!        nvec_neg=nvec;nvec_neg(m,k)=nvec_neg(m,k)-1
!        call index(nminus,nvec_neg,1)
        nminus=map_nneg(m,k,n)
        if(nminus>0.and.nminus<=nn_tot) then
          if(zero(nminus)==0)then
            tmp=tmp-iota*nvec(m,k)*(c_matsubara(k)*matmul(mat_tmp,rho(:,:,nminus)) - dconjg(c_matsubara(k))*matmul(rho(:,:,nminus),mat_tmp) )!*sqrt(nvec(m,k)/c_matsubara(k))
          endif
        endif
      enddo
    enddo

    drho_dt=tmp
    !stop

end subroutine compute_deriv_n
!-----------------------------------------------------------------  

subroutine filter
  !! Shi, Chen, Nan, Xu, Yan, JCP 130, 084105 (2009)
  integer n,l

  zero=0
  do l=1,6
    do n=1,nn_tot
      if(maxval(abs(rho(:,:,n,l)))<tolerance) then
        rho(:,:,n,l)=0.d0
        zero(n,l)=1
      endif
    enddo
  enddo

!write(6,*) t1*au2s*1.d15,t3*au2s*1.d15,nn_tot-sum(zero),nn_tot

end subroutine filter
!-----------------------------------------------------------------  

subroutine index(n,nvec,iflag)
  implicit none
  integer,intent(inout)::n
  integer,intent(inout)::nvec(nquant,0:KK)
  integer,intent(in)::iflag
  integer m,k,ind,temp
  integer n_beg,n_end,l_sum
  integer clck_counts_beg,clck_counts_end, clck_rate
  real*8 t1,t2

!  call cpu_time(t1)

  if(iflag==0) then  !n is input, nvec is output
    nvec=nn(:,:,n)
  endif

  if(iflag==1) then  !nvec is input, n is output
    n=0
    l_sum=sum(nvec)
    if(l_sum<=LL) then
      n_beg=map_sum(l_sum)
      if(l_sum==LL) n_end=nn_tot
      if(l_sum<LL) n_end=map_sum(l_sum+1)-1
      do m=1,nn_tot!n_beg,n_end
        if(all(nvec==nn(:,:,m))) then
          n=m
          exit
        endif
      enddo
    endif
  endif
        
!  call cpu_time(t2)
!  tim_ind=tim_ind+t2-t1

end subroutine index
!-----------------------------------------------------------------  

subroutine binary(n,nvec,iflag,base)
  implicit none
  integer,intent(inout)::n
  integer,intent(inout)::nvec(nquant,0:KK)
  integer,intent(in)::iflag,base
  integer m,k,ind,temp

  if(iflag==0) then  !n is input, nvec is output
    temp = n-1
    nvec=0
    m=1;k=0
    DO WHILE(temp > 0)
      nvec(m,k)=mod(temp,base)
      k=k+1
      if(k>KK) then
        k=1;m=m+1
      endif
      temp = temp / base
    END DO
  endif

  if(iflag==1) then  !nvec is input, n is output
    n=0
    ind=0
    do m=1,nquant
      do k=0,KK
        if(nvec(m,k).ne.0)n=n+nvec(m,k)*base**ind
        ind=ind+1
      enddo
    enddo
    n=n+1
  endif

end subroutine binary
!-----------------------------------------------------------------  

subroutine compute_nn
  implicit none
  integer n,m,n_tot,tot
  integer nvec(nquant,0:KK)
  integer n_beg,n_end

!! This method works, but inefficient for large LL
!  nvec=0;nvec(nquant,KK)=LL
!  call binary(n_tot,nvec,1,LL+1)
!
!  tot=0
!  do n=1,n_tot
!    m=n
!    call binary(m,nvec,0,LL+1)
!    if(sum(nvec)<=LL) then
!      tot=tot+1
!      nn(:,:,tot)=nvec
!    endif
!  enddo

!! Alternate recursive method

  n_tot=0
  n_beg=0;n_end=0
  do n=0,LL
    call compute_nn_sum_L(n,n_beg,n_end)
    map_sum(n)=n_beg
write(6,*) n,n_beg,n_end
  enddo
write(6,*) n_end,nn_tot

end subroutine compute_nn
!-----------------------------------------------------------------  

subroutine compute_nn_sum_L(L,n_beg,n_end)
  implicit none
  integer,intent(in)::L
  integer,intent(inout)::n_beg,n_end!! at input, state/end point of entries for sum L-1, at output, for entries for sum L
  integer i,j,m,k,tot_L
  integer flag_new,nvec(nquant,0:KK)

  if(L==0) then
    nn(:,:,1)=0
    n_beg=1;n_end=1
  else
    tot_L=n_end
    do i=n_beg,n_end
      do m=1,nquant
        do k=0,KK
          
          nvec=nn(:,:,i);nvec(m,k)=nvec(m,k)+1
          flag_new=0
          if(tot_L>n_end) then
            do j=n_end+1,tot_L
              if(all(nvec==nn(:,:,j))) then
                flag_new=1;exit
              endif
            enddo
          endif

          if(flag_new==0) then
            tot_L=tot_L+1
            nn(:,:,tot_L)=nn(:,:,i)
            nn(m,k,tot_L)=nn(m,k,i)+1
!            map_nneg(m,k,tot_L)=i
!            map_nplus(m,k,i)=tot_L
          endif
        enddo
      enddo
    enddo
    n_beg=n_end+1
    n_end=tot_L
  endif

end subroutine compute_nn_sum_L
!-----------------------------------------------------------------  

subroutine compute_map
  implicit none
  integer n,m,k
  integer nvec(nquant,0:KK),nvec_plus(nquant,0:KK),nvec_neg(nquant,0:KK)
  integer nplus,nneg

  map_nneg=0
  map_nplus=0
  do n=1,nn_tot
    m=n
    call index(m,nvec,0)
    do m=1,nquant
      do k=0,KK
        nvec_plus=nvec;nvec_plus(m,k)=nvec_plus(m,k)+1
        call index(nplus,nvec_plus,1)
        map_nplus(m,k,n)=nplus
       
        if(nvec(m,k)>0) then
          nvec_neg=nvec;nvec_neg(m,k)=nvec_neg(m,k)-1
          call index(nneg,nvec_neg,1)
          map_nneg(m,k,n)=nneg
        endif
      enddo
    enddo
  enddo

write(6,*) "map completed ..."

end subroutine compute_map
!-----------------------------------------------------------------  

subroutine compute_nn_tot
  !nn_tot=(factorial(LL+nquant*KK)/factorial(LL))/factorial(nquant*KK)
  implicit none
  integer i
  real*8 tmp

  tmp=1.d0
  do i=1,nquant*(KK+1)
    tmp=tmp*(LL+i)/real(i)
  enddo
  nn_tot=nint(tmp)

!nn_tot=(factorial(LL+nquant*KK)/factorial(LL))/factorial(nquant*KK)

end subroutine compute_nn_tot
!-----------------------------------------------------------------  

function commute(mat1,mat2) result(mat3)
  complex*16,intent(in) :: mat1(:,:),mat2(:,:)
  complex*16,allocatable :: mat3(:,:)
  integer i1

  i1=size(mat1,1)
  allocate(mat3(i1,i1))

  mat3=matmul(mat1,mat2)-matmul(mat2,mat1)

end function commute
!-----------------------------------------------------------------  

integer function factorial(n)
  implicit none
  integer,intent(in) :: n
  integer i

  if(n==0) then
    factorial=1
  else
    factorial=1
    do i=1,n
      factorial=factorial*i
    enddo
  endif

end function factorial
!-----------------------------------------------------------------  

subroutine fft_dip_mom_corr
  implicit none
  integer i,j
  integer nst_w1,nst_w3
  real*8,allocatable :: freq_w1(:),freq_w3(:)
  real*8,allocatable :: tim_t1(:),tim_t3(:)
  real*8 re1,cc1,re2,cc2,w1,w3
  complex*16,allocatable :: rephase(:,:),non_rephase(:,:)
  complex*16,allocatable :: fft_rephase(:,:),fft_nonrephase(:,:),spectra(:,:)

  nst_t1=nst_t1/nst_avg+1
  nst_t3=nst_t3/nst_avg+1

  allocate(tim_t1(nst_t1),tim_t3(nst_t3),rephase(nst_t1,nst_t3),non_rephase(nst_t1,nst_t3))

  do i=1,nst_t1
    do j=1,nst_t3
      read(101,*)tim_t1(i),tim_t3(j),re1,cc1,re2,cc2
      rephase(i,j)=cmplx(re1,cc1)
      non_rephase(i,j)=cmplx(re2,cc2)
    enddo
    read(101,*)
  enddo

  call fft_2d(rephase,non_rephase,tim_t1,tim_t3,freq_w1,nst_w1,freq_w3,nst_w3,fft_rephase)
  !call fft_2d(non_rephase,tim_t1,tim_t3,freq_w1,nst_w1,freq_w3,nst_w3,fft_nonrephase)

  allocate(spectra(nst_w1,nst_w3))

  spectra=fft_rephase

!  do i=2,nst_w1
!    do j=1,nst_w3
!      spectra(i,j)=fft_nonrephase(i,j)+fft_rephase(nst_w1-i+2,j)
!    enddo
!  enddo

  do i=1,nst_w1
    w1=(hbar*freq_w1(i))*au2J/(wave_to_J)
!    if(w1>11000.d0.and.w1<12500.d0) then
      do j=1,nst_w3
        w3=(hbar*freq_w3(j))*au2J/(wave_to_J)
!        if(w3>11000.d0.and.w3<12500.d0) then
          write(102,'(8es15.5)')w1,w3,spectra(i,j)!,fft_rephase(i,j),fft_nonrephase(i,j)
!        endif
      enddo
      write(102,*)
      write(200,*)
!    endif
  enddo

end subroutine fft_dip_mom_corr
!-----------------------------------------------------------------  

subroutine fft_2d(rephase,non_rephase,tim_t1,tim_t3,freq_w1,nst_w1,freq_w3,nst_w3,dat_f_2d)
  complex*16,intent(in) :: rephase(nst_t1,nst_t3),non_rephase(nst_t1,nst_t3)
  real*8,intent(in) :: tim_t1(nst_t1),tim_t3(nst_t3)
  integer,intent(out) :: nst_w1,nst_w3
  real*8,allocatable,intent(out) :: freq_w1(:),freq_w3(:)
  complex*16,allocatable,intent(out) :: dat_f_2d(:,:)
  real*8 t1,t3,w1,w3
  real*8 re,cc
  complex*16 tmp
  integer i,j,k,l

  nst_w1=50
  nst_w3=50
  allocate(freq_w1(nst_w1),freq_w3(nst_w3),dat_f_2d(nst_w1,nst_w3))

  do i=1,nst_w1
    freq_w1(i)=-750.d0+1500.d0*(i-1)/real(nst_w1-1)
  enddo
  freq_w1=freq_w1*wave_to_J/(au2J*hbar)

  do i=1,nst_w3
    freq_w3(i)=-750.d0+1500.d0*(i-1)/real(nst_w3-1)
  enddo
  freq_w3=freq_w3*wave_to_J/(au2J*hbar)

  do k=1,nst_w1
    w1=freq_w1(k)
    do l=1,nst_w3
      w3=freq_w3(l)
      tmp=0.d0
      do i=1,nst_t1
        t1=tim_t1(i)
        do j=1,nst_t3
          t3=tim_t3(j)
          tmp=tmp+exp(iota*(w1*t1+w3*t3))*non_rephase(i,j)+exp(iota*(-w1*t1+w3*t3))*rephase(i,j)
        enddo
      enddo
      dat_f_2d(k,l)=tmp
    enddo
  enddo

end subroutine fft_2d
!-----------------------------------------------------------------  

End Module mod_spectra_heom
