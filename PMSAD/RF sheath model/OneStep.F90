Module OneStepModule      
	Use FieldModule 
	Use Input
    Implicit none	
    Integer(4),parameter :: NSpecyMax=1_4
    Integer(4) :: NParticle=0  
    Type(ParticleBundle) ::  ParticleGlobal(0:NSpecyMax)      
    Integer(4),parameter :: NGasMax=2_4
    Integer(4) ::  Ngas 
    Type(Gas) ::  GasGlobal(1:NGasMax)     
    Integer(4) :: Period,Nx,Ny,Nxy,gap,Nelechight
    Real(8) :: Xmin,Xmax,Ymin,Ymax
    Real(8) :: dx,dt
    Type(PICParticlePara) ::  PICParaGlobal(0:NSpecyMax)
    Integer(4) :: tianjiaE(1:(InputNx),1:(InputNy))
	Integer(4) :: tianjiaI(1:(InputNx),1:(InputNy))
    Integer(4) :: rhoagain(1:InputNy)
    Integer(4) :: Nepergrid
	Integer(4) :: Nipergrid
    contains
	
    Subroutine Initilalization()
        Use RFSource
        Use Input
        Use GasModule
        Use InitilalizationModule
        Use ElectronModule
        Use HGasModule
        Use ParticleBoundary
        Use ParticleModule
        Use Parallel
        Use Move2D
        Implicit none
			 
        Integer(4) ::  i,NSpecy,Shift=0
        Logical :: alive
        character(len=30) name,Filename       
		Ngas=InputNGas
        Nx=InputNx
        Ny=InputNy
        Nxy=Nx*Ny
        dx=Inputdx
        dt=Inputdt       
		Xmin=0.d0
        XMax=dble(Nx-1)
        Ymin=0.d0
        YMax=dble(Ny-1)      
		gap=inputgap
        Nelechight=Nelech-1      
		Call DRandomInt(Myid)
        Call AborptionXYInit(Xmin,XMax,Ymin,YMax,Gap)
        Call AborptionXYInit1(Xmin,XMax,Ymin,YMax,Gap)
        Call AborptionXYInit2(Xmin,XMax,Ymin,YMax,Gap)
        Call RfVoltageXInitilalization(RFKind,RFFrequency,RfVol,dt,Period)
        Call RfVoltageYInitilalization(YFrequency,YVol,dt)
        Call XYFieldInit(Nx,Ny,Gap,dx,Nelechight)
        Call ParticleInit(ElectronGas,1_4,Electron,PICParaGlobal(0),ParticleGlobal(0),dx,dt,InitDensity,Nx,Ny,ParticlePerGrid,NProcess)   
		do i=1,Ngas
            select case (GasDefinition(i))
                Case(1_4) 
                    Call GasInit(GasPressure(i),GasTemperature(i),ArGas,GasGlobal(i),Shift) 
                    NSpecy=ArGas%NSpecy
                    Call ParticleInit(GasGlobal(i),NSpecy,ArSpecy,PICParaGlobal(Shift+1:Shift+NSpecy),ParticleGlobal(Shift+1:Shift+NSpecy),dx,dt,InitDensity,Nx,Ny,ParticlePerGrid,NProcess)                       
					Shift=Shift+NSpecy               
				Case(2_4) 
           End Select      
		End do		
		Nepergrid = ParticleGlobal(0)%Npar/(Nx*Ny)
		Nipergrid = ParticleGlobal(1)%Npar/(Nx*Ny)
        NParticle=Shift		
		BxField=2.0*cos(30*3.1415926/180.0) 
		ByField=0.0
		BzField=2.0*sin(30*3.1415926/180.0)
		Call MPI_Bcast(BxField,Nxy,mpi_double_precision,master, MPI_COMM_WORLD,Ierr)
        Call MPI_Bcast(ByField,Nxy,mpi_double_precision,master, MPI_COMM_WORLD,Ierr)
        Call MPI_Bcast(BzField,Nxy,mpi_double_precision,master, MPI_COMM_WORLD,Ierr)              				
        return      
	End Subroutine Initilalization
          
    Subroutine OneStep()
        Use WeightingModule
        Use Move2D
        Use ParticleBoundary 
        Use FileIO
        Use Diagnostics  
        Use ElectronModule
        Use HGasModule
        Use InitilalizationModule
        Implicit none
        
		Integer(4) :: i,j,k,N1,N2
        real(8) :: Timer1,Timer2,Line(inputNy)
        Type(Grid) :: TGrid
        Character(len=30) :: name        
		Call MPI_Barrier(MPI_COMM_WORLD,ierr)		
		Line(1:30)=1.0       
		Do i=31,Ny  
            Line(i)=1.0+DLog(Dble(i-1)/Dble(2-1))/DLog(Dble(Ny-1)/Dble(2-1))*(0.01d0-1.0)
        End do
		do i=1,inputNy
			tianjiaI(1:inputNx,i)= int(Line(i)*Nipergrid)
		end do		
		do i = ParticleGlobal(1)%NPar,1,-1			
			N1=Ceiling(ParticleGlobal(1)%PhaseSpace(i)%X)
            N2=Ceiling(ParticleGlobal(1)%PhaseSpace(i)%Y) 
            tianjiaI(N1,N2)= tianjiaI(N1,N2)-1				
		end do
        do i= 1,inputNx
            do j=1,inputNy
				if(tianjiaI(i,j) <0) then		
                    tianjiaI(i,j)=0
				end if                                               
			end do
		end do
		tianjiaI(1:50,1:inputNy)=0
		call add(GasGlobal(1),tianjiaI(1:inputNx,1:inputNy),ParticleGlobal(1),Xmin,Xmax,Ymin,Ymax,dx,k)		          
		Line(1:30)=1.0
        Do i=31,Ny  
            Line(i)=1.0+DLog(Dble(i-1)/Dble(2-1))/DLog(Dble(Ny-1)/Dble(2-1))*(0.01d0-1.0)
        End do
		do i=1,inputNy
			tianjiaE(1:inputNx,i)= int(Line(i)*Nipergrid)
		end do			
		do i = ParticleGlobal(0)%NPar,1,-1
			N1=Ceiling(ParticleGlobal(0)%PhaseSpace(i)%X)
            N2=Ceiling(ParticleGlobal(0)%PhaseSpace(i)%Y) 
            tianjiaE(N1,N2)= tianjiaE(N1,N2)-1			
		end do
        do i= 1,inputNx
            do j=1,inputNy
				if(tianjiaE(i,j) <0) then		
                    tianjiaE(i,j)=0
				end if                                               
			end do
		end do
		tianjiaE(1:50,1:inputNy)=0
		call add(ElectronGas,tianjiaE(1:inputNx,1:inputNy),ParticleGlobal(0),Xmin,Xmax,Ymin,Ymax,dx,k)		

		do i=0,NParticle           
			if(i==1) then
                Call AborptionXY1(ParticleGlobal(i),PICParaGlobal(i)%QFactor,PICParaGlobal(i)%QBottom,PICParaGlobal(i)%QTop,dx,rhoagain(1:Ny))			
            end if 	
            if(i==0) then
                Call AborptionXY2(ParticleGlobal(i),PICParaGlobal(i)%QFactor,PICParaGlobal(i)%QBottom,PICParaGlobal(i)%QTop,dx)
            end if 			
			call MoveB(ParticleGlobal(i),Nx,Ny,ExGlobal%Value,EyGlobal%Value,BxField,ByField,Bzfield,PICParaGlobal(i)%AccFactor,PICParaGlobal(i)%AccFactorB)
            if(i==1) then
                Call AborptionXY1(ParticleGlobal(i),PICParaGlobal(i)%QFactor,PICParaGlobal(i)%QBottom,PICParaGlobal(i)%QTop,dx,rhoagain(1:Ny))
				call againDione(ParticleGlobal(1),ParticleGlobal(0),Nx,Ny,rhoagain(1:Ny))
            end if  
            if(i==0) then
                Call AborptionXY2(ParticleGlobal(i),PICParaGlobal(i)%QFactor,PICParaGlobal(i)%QBottom,PICParaGlobal(i)%QTop,dx)
            end if 
            Call Weighting2D(ParticleGlobal(i),Nx,Ny, PICParaGlobal(i)%Value)			   
        end do
		Call MPI_Barrier(MPI_COMM_WORLD,ierr) 
		Call XYField(NParticle,PICParaGlobal(0:NParticle)) 		
        return           
    End  subroutine OneStep
       
    Subroutine OneStepDiag()
        Use WeightingModule
        Use Move2D
        Use ParticleBoundary 
        Use FileIO
        Use Diagnostics  
        Use input
        Use ElectronModule
        Use HGasModule
		Use InitilalizationModule
        Implicit none
			   
        Integer(4) :: i,j,N1,N2,k
        real(8) :: Timer1,Timer2,Nelechight,Xbottomgap ,line(inputNy) 
        Type(Grid) :: TGrid
        Nelechight=Nelech-1
        Xbottomgap=dble(Ny-1)-dble(Gap)

		Line(1:30)=1.0
        Do i=31,Ny  
            Line(i)=1.0+DLog(Dble(i-1)/Dble(2-1))/DLog(Dble(Ny-1)/Dble(2-1))*(0.01d0-1.0)
		End do

		do i=1,inputNy
			tianjiaI(1:inputNx,i)= int(Line(i)*Nipergrid)
		end do
		
		do i = ParticleGlobal(1)%NPar,1,-1
			N1=Ceiling(ParticleGlobal(1)%PhaseSpace(i)%X)
            N2=Ceiling(ParticleGlobal(1)%PhaseSpace(i)%Y) 
            tianjiaI(N1,N2)= tianjiaI(N1,N2)-1			
		end do

        do i= 1,inputNx
            do j=1,inputNy
				if(tianjiaI(i,j) <0) then		
                    tianjiaI(i,j)=0
				end if                                              
			end do
		end do

		tianjiaI(1:50,1:inputNy)=0

		call add(GasGlobal(1),tianjiaI(1:inputNx,1:inputNy),ParticleGlobal(1),Xmin,Xmax,Ymin,Ymax,dx,k)
		
		Line(1:30)=1.0
           
		Do i=31,Ny  
            Line(i)=1.0+DLog(Dble(i-1)/Dble(2-1))/DLog(Dble(Ny-1)/Dble(2-1))*(0.01d0-1.0)
        End do

		do i=1,inputNy
			tianjiaE(1:inputNx,i)= int(Line(i)*Nipergrid)
		end do
	
		do i = ParticleGlobal(0)%NPar,1,-1
			N1=Ceiling(ParticleGlobal(0)%PhaseSpace(i)%X)
            N2=Ceiling(ParticleGlobal(0)%PhaseSpace(i)%Y) 
            tianjiaE(N1,N2)= tianjiaE(N1,N2)-1			
		end do

        do i= 1,inputNx
            do j=1,inputNy
				if(tianjiaE(i,j) <0) then		
                    tianjiaE(i,j)=0
				end if                                               
			end do
		end do

		tianjiaE(1:50,1:inputNy)=0
		call add(ElectronGas,tianjiaE(1:inputNx,1:inputNy),ParticleGlobal(0),Xmin,Xmax,Ymin,Ymax,dx,k)
		
        do i=0,NParticle
		
            if(i==1) then
                Call AborptionXY1(ParticleGlobal(i),PICParaGlobal(i)%QFactor,PICParaGlobal(i)%QBottom,PICParaGlobal(i)%QTop,dx,rhoagain(1:Ny))
			end if 
			
			if(i==0) then
                Call AborptionXY2(ParticleGlobal(i),PICParaGlobal(i)%QFactor,PICParaGlobal(i)%QBottom,PICParaGlobal(i)%QTop,dx)
            end if 
			
			call MoveB(ParticleGlobal(i),Nx,Ny,ExGlobal%Value,EyGlobal%Value,BxField,ByField,Bzfield,PICParaGlobal(i)%AccFactor,PICParaGlobal(i)%AccFactorB)          
			Call DiagDumpParticle(ParticleGlobal(1)) 
			
            if(i==1) then
				Call AborptionXY1(ParticleGlobal(i),PICParaGlobal(i)%QFactor,PICParaGlobal(i)%QBottom,PICParaGlobal(i)%QTop,dx,rhoagain(1:Ny))                 
            end if 
 
            if(i==0) then
                Call AborptionXY2(ParticleGlobal(i),PICParaGlobal(i)%QFactor,PICParaGlobal(i)%QBottom,PICParaGlobal(i)%QTop,dx)
            end if 

            Call Weighting2D(ParticleGlobal(i),Nx,Ny, PICParaGlobal(i)%Value)
   
        end do

        Call MPI_Barrier(MPI_COMM_WORLD,ierr) 
        Call XYField(NParticle,PICParaGlobal(0:NParticle))        
 
		return           
    End  subroutine OneStepDiag  	   
	
End Module OneStepModule
