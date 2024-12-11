Program PIC_MCC_for_CCP

	Use OneStepModule
	Use FileIO 
	Use FieldModule 
	Use WeightingModule
	Use Diagnostics
	Implicit none
	
	Type(Grid) :: Temp
	Integer(4) :: i,j 
	Real(8) Ctime1,Ctime2 

	Call MPI_INIT(ierr)
	Call MPI_COMM_RANK(MPI_COMM_WORLD , myid , ierr)
	Call MPI_COMM_SIZE(MPI_COMM_WORLD , NProcess , ierr) 
	
	Call Initilalization()
	
	Write(*,*)  myid,'myid',ParticleGlobal%NPar,'init over'
	Temp%Nx=Nx
	Temp%Ny=Ny
	Temp%Name='erho'
	Call Weighting2D(ParticleGlobal(0),Nx,Ny, PICParaGlobal(0)%Value)
	Call MPI_ALLReduce(PICParaGlobal(0)%Value,Temp%Value,Nx*Ny,mpi_double_precision,MPI_SUM,MPI_COMM_WORLD,Ierr)
	If(Myid==Master) Call DumpGrid(Temp)		
	If(Myid==Master) then 
        
		Call DumpGrid(PhiGlobal)
        Call DumpGrid(ExGlobal)
        Call DumpGrid(EyGlobal)
        Call DumpGrid(RhoGlobal)
	
	End if
		
	do j=1, 1
		Ctime1=MPI_Wtime() 		
		do i=1,1
			Call OneStep()
		end do
		If(Myid==Master) then             
			Call DumpGrid(PhiGlobal)
			Call DumpGrid(ExGlobal)
			Call DumpGrid(EyGlobal)
			Call DumpGrid(RhoGlobal)			
		End if	  
		Ctime2=MPI_Wtime()  
	  
		If(myid==master) Write(*,*)  ParticleGlobal%NPar,'NParB' ,Ctime2-Ctime1,Period,i,j  
			if (myid==0)then
				Open(13,file="particle.dat",access='append') !指在最后面续写
				Write(13,*)  j,ParticleGlobal%NPar
				close(13)
			end if
		if(mod(j,3)==1) then
            
			Temp%Name='erho'
			Call Weighting2D(ParticleGlobal(0),Nx,Ny, PICParaGlobal(0)%Value)
			Call MPI_ALLReduce(PICParaGlobal(0)%Value,Temp%Value,Nx*Ny,mpi_double_precision,MPI_SUM,MPI_COMM_WORLD,Ierr)
			If(Myid==Master) Call DumpGrid(Temp)
    
			Temp%Name='irho' 
			Call Weighting2D(ParticleGlobal(1),Nx,Ny, PICParaGlobal(1)%Value)
			Call MPI_ALLReduce(PICParaGlobal(1)%Value,Temp%Value,Nx*Ny,mpi_double_precision,MPI_SUM,MPI_COMM_WORLD,Ierr) 
			If(Myid==Master) Call DumpGrid(Temp)
	   
			If(Myid==Master) then 
				Call DumpGrid(PhiGlobal)
				Call DumpGrid(ExGlobal)
				Call DumpGrid(EyGlobal)
				Call DumpGrid(RhoGlobal)
			End if
	  
			do i=0,NParticle   
				call ParallelDumpParticle(ParticleGlobal(i))
			end do
		end if
    End do

	Call DensityTemperatureInit(2,ParticleGlobal(0:1),PhiGlobal%Nx,PhiGlobal%Ny,PICParaGlobal(0:1))
	Call AvaPhiEInit(PhiGlobal,ExGlobal,EyGlobal)
	Call DiagDumpParticleInit(ParticleGlobal(1),Xmin,XMax,Ymin,YMax)
	Call CalPowerDensityInit(NParticle,ParticleGlobal(0:NParticle),Nx,Ny,dx,PICParaGlobal(0:1))
	Call WeightingEEDFInit() 

	do j=1, 0
		do i=1, 1 
			Call OneStepDiag ()
			Call DensityTemperature(2,ParticleGlobal(0:1))
			Call AvaPhiE(PhiGlobal,ExGlobal,EyGlobal)
			Call CalPowerDensity(NParticle,ParticleGlobal(0),ExGlobal,EyGlobal,PICParaGlobal(0:1))           
			Call WeightingEEDF(ParticleGlobal(0))
		end do
	end do

	Call DensityTemperatureFinal(2,PICParaGlobal(0:1)) 
	Call AvaPhiEFinal()
	Call DiagDumpParticleFinal()
	Call CalPowerDensityFinal()
	Call WeightingEEDFFinal()
   
	Call DiagDumpParticleFinal()
	Call MPI_Finalize(ierr)
	stop
end  Program  


               
                
