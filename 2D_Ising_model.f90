! Program to generate a lattice and study its ferromagnetic properties

program ising
    IMPLICIT NONE 

    integer :: i,j,L,p,a,b,c,d,niter,time,mm,mn,N
    real :: r,E,M,mag,dE,Ei,Ef,u,h  

    real :: T=2.0,J_ising=1.0 ! assigning value to relevant parameter: k_BT=1

    integer,dimension(:,:),allocatable :: spin
    integer :: seed
    character(len=30) :: charac_b ! charac_b stores the name dump_pos

    seed=44859
    charac_b = 'store_config'

    print*,'Enter the number of lattice points in one dimension'
    read*,L
    print*,'Enter the number of iterations'
    read*,niter 

    allocate(spin(L,L))
    E=0.0  ! Instantaneous Energy of the lattice
    M=0.0  ! Instantaneous magnetization of the lattice
    N=L*L 

    call random_seed 

 ! Intialize your lattice
    open(71,file='initial_ising.dat')
    p=0
    do i=1,L
        do j=1,L
            call RANDOM_NUMBER(r)
            spin(j,i) = +1

            if(r<0.5)then
                spin(j,i) = -1
            else 
                spin(j,i) = +1
            end if 
            !WRITING DOWN INTIIAL CONFIGURATION 
            write(71,fmt='(4g10.8)') float(i),float(j),float(p),float(spin(j,i))
        end do 
    end do 
close(71)

! Calculate initial magnetization and energy
do i=1,L
    do j=1,L
        a=i+1;b=i-1;c=j+1;d=j-1  !Identifying four numbers of spin(i,j)

        if(i==L)a=1  ! PBC: periodic boundary conditions.
        if(i==1)b=L
        if(j==1)d=L
        if(j==L)c=1

        M=M+spin(i,j)
        E=E-J_ising*float((spin(j,i))*(spin(a,j)+spin(b,j)+spin(i,c)+spin(i,d)))
    end do 
end do

mag=M/(float(N))   ! N=L*L: magnetization (instantaneous) per spin 
E=E*0.5
print*,'intial energy E, E per spin =',E, E/float(N)
print*,'intial magnetization M, M per spin=',M,mag 

! INITIALIZATION COMPLETE
!____________________________________________________________________________
! Evolve it to reach equilibrium

open(10,file='ising_T2_L40_init_random.dat')
do time=1,niter ! loop over no of MCS 

    do mm=1,L
        do mn=1,L

            call RANDOM_NUMBER(r); i=int(r*float(L))+1  ! Choosing a lattice site
            call RANDOM_NUMBER(r); j=int(r*float(L))+1

            a=i+1;b=i-1;c=j+1;d=j-1  ! Identifying neighbours of spin(i,j)

            if(i==L)a=1; if(i==1)b=L; if(j==1)d=L; if(j==L)c=1  !PBC 
            Ei=-J_ising*float(spin(i,j)*(spin(a,j)+spin(b,j)+spin(i,c)+spin(i,d)))  ! BEFORE TRIAL FLIP

            spin(i,j) = -spin(i,j)  ! TRIAL FLIP

            Ef=-J_ising*float((spin(i,j))*(spin(a,j)+spin(b,j)+spin(i,c)+spin(i,d))) ! AFTER TRIAL FLIP

            dE = Ef-Ei  ! DIFFERENCE IN ENERGIES


            if(dE<=0.0)then
                E=E+dE  ! UPDATING ENERGY AND MAGNETIZATION OF LATTICE
                M=M+(2.0*float(spin(i,j)))
            else 
                u=exp(-dE/(T))
                call RANDOM_NUMBER(h)
                if(h<u)then
                    E=E+dE 
                    M=M+(2.0*float(spin(i,j)))  ! INSTANANEOUS MAG. of ENTIRE LATTICE
                else 
                    spin(i,j)=-spin(i,j) ! TRIAL FLIP NOT ACCEPTED; E and M NOT UPDATED 
                end if
            end if 
        end do
    end do 

    write(10,*)time,M/float(N),E/float(N)  ! writing down E and M with no. of iterations
end do 

close(10)

print*,'final energy E, E per spin =',E, E/float(N)
print*,'final magnetization M, M per spin=',M,mag 

End program ising 



            
