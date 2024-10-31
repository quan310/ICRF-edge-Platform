program monte_carlo
    
	use Field
    Use MPI
    implicit none

    Integer :: NArion,NArion1,Ndepos2D
    Integer :: Nsputter = 0
    Integer, parameter :: Nx=301, Ny=101, Nparmax = 5000000, NsputterMax = 8000000, N_E = 1000,Nspu=10, N_Eion=100, Nyy = 101,N_Aion=90  !123555        
    Integer, parameter :: Ny2D = 201
	real(8), parameter :: Xlength=3.0d-2,Ylength=1.0d-2,dt=2.7d-7	
	real(8), parameter :: Gap=1.0d-2
    real(8), parameter :: dy = Ylength/dble(Ny-1),VConvertFactor=dy/dt,regional1(1:Nx,1:Ny)=0.d0
    real(8) :: sputter(1:5,1:NsputterMax) = 0.d0			
	real(8) :: Arion(1:5,1:Nparmax) 
	real(8) :: Arionall(1:5,1:Nparmax)
	real(8) :: Arion1(1:5,1:Nparmax)						
	real(8) :: depos(1:5,1:NsputterMax)= 0.d0				
	real(8) :: depos2D(1:5,1:NsputterMax)= 0.d0
	real(8) :: Rho2D(1:Nx,1:Ny2D) = 0.d0
	real(8) :: Rho2Dall(1:Nx,1:Ny2D) = 0.d0	
	real(8) :: Rhodepos(Nx)=0.d0								
	real(8) :: Rhosputter(Ny)=0.d0								
	real(8) :: Rhoincident(Ny)=0.d0	
	real(8) :: Rhosputterall(Ny)=0.d0
	real(8) :: Rhoincidentall(Ny)=0.d0								
	real(8) :: YEDF(Nyy,N_Eion)=0.d0
	real(8) :: YEDFall(Nyy,N_Eion)=0.d0
	real(8) :: YADF(Nyy,N_Aion)=0.d0
	real(8) :: YADFall(Nyy,N_Aion)=0.d0    
	real(8):: IonEnergy
	real(8), parameter :: Mion = 1.66d-27
	Real(8),parameter :: JtoeV=1.6022d-19   
	integer :: i,j
	Character(len=30) :: name
	Character(len=30) :: dump
	character(len=100) :: filename
	integer :: fileunit

	Call MPI_INIT(ierr)
    Call MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr)
    Call MPI_COMM_SIZE(MPI_COMM_WORLD,NProcess,ierr) 
    dims(1)=NProcess
    Call MPI_CART_CREATE(MPI_COMM_WORLD,Ndim,dims,ParaPeriod,Reorder,SerialCOMM,ierr)
    Call MPI_COMM_RANK(SerialCOMM, myid2, ierr)
    
	Name='ArX'					!读取离子信息
    Call ParallelLoad(Name,NArion1,NParMax,Arion1(1:5,1:NParMax))
	Write(*,*) 'NArion1=',NArion1
	NArion = 1
	
	do i=1,NArion1									        
		if (Arion1(2,i)>0.0025) then						
				Arion(1,NArion) = Arion1(1,i)/dy
				Arion(2,NArion) = Arion1(2,i)/dy							
				Arion(3,NArion) = Arion1(3,i)
				Arion(4,NArion) = -Arion1(4,i)	
				Arion(5,NArion) = Arion1(5,i)			
				NArion = NArion + 1			
		end if		
	end do
	
	Write(*,*) 'NArion=',NArion
	
	call YEDFion(NArion,N_Eion,Nyy,YEDF(1:Nyy,1:N_Eion),Arion(1:5,1:NArion))		
	Call MPI_AllReduce(YEDF,YEDFall,Nyy*N_Eion,MPI_DOUBLE_PRECISION,MPI_SUM,SerialCOMM,ierr)	
	name='YEDF_ion'
	call Dump2D(Name,Nyy,N_Eion,YEDFall)	
	
	call YADFion(NArion,N_Aion,Nyy,YADF(1:Nyy,1:N_Aion),Arion(1:5,1:NArion))
	Call MPI_AllReduce(YADF,YADFall,Nyy*N_Aion,MPI_DOUBLE_PRECISION,MPI_SUM,SerialCOMM,ierr)
	name='YADF_ion'
	call Dump2D(Name,Nyy,N_Aion,YADFall)		
     
	Call DWeightingy(NArion,Arion(1:5,1:NParMax),Ny,Rhoincident)		
	Call MPI_AllReduce(Rhoincident,Rhoincidentall,Ny,MPI_DOUBLE_PRECISION,MPI_SUM,SerialCOMM,ierr)		
	if(myid==0) then
		name='Arrho'
		Call Dump1D(name,Ny,Rhoincidentall(1:Ny))			
	end if
		
	Ndepos2D = 0	
    Call sputteratom(NArion, Arion(1:5,1:NArion),Nsputter,NsputterMax,sputter(1:5,1:NsputterMax))	
	
	Name='sputter'													
    Call ParallelDump(Name,Nsputter,sputter(1:5,1:Nsputter))		
	Call SWeighting(Nsputter,sputter(1:5,1:Nsputter),Ny,Rhosputter)	 
	Call MPI_AllReduce(Rhosputter,Rhosputterall,Ny,MPI_DOUBLE_PRECISION,MPI_SUM,SerialCOMM,ierr)				
	if(myid==0) then	
		name='Rhosputter'
		Call Dump1D(name,Ny,Rhosputterall(1:Ny))		
	end if	
		
	Call Move(Nsputter,sputter(1:5,1:Nsputter),Nx,Ny2D,dy,dt,Nspu,Ndepos2D,depos2D(1:5,1:Nsputter*Nspu))
	call Weighting2D(Ndepos2D,depos2D(1:5,1:Ndepos2D),Nx,Ny2D,Rho2D)
	Call MPI_AllReduce(Rho2D,Rho2Dall,Nx*Ny2D,MPI_DOUBLE_PRECISION,MPI_SUM,SerialCOMM,ierr)
	name='Rho2D'
	call Dump2D(name,Nx,Ny2D,Rho2Dall)	
	
    Write(*,*) 'stop'
    Call MPI_FINALIZE(ierr)
    pause 0
	stop
	
end program monte_carlo


 