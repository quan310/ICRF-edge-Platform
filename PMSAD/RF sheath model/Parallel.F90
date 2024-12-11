Module Parallel
	Use MPI
	Use TypeModule
	Use ParticleModule
	implicit none  
	integer,parameter :: MAX_LOCAL_ROW=2048,MAX_LOCAL_COL=2048,MAX_PROC=128
	integer,parameter :: RHO_FIRST_PHASE=1001,RHO_SECOND_PHASE=1002,RHO_THIRD_PHASE=1003
	integer,private :: snds(MAX_PROC),disp(MAX_PROC)
	integer :: lrows(MAX_PROC),local_ii(MAX_PROC)
	Integer :: NProcess,Myid=0,ierr
	integer :: master=0
	integer,Private :: status(MPI_STATUS_SIZE)
	save snds,disp,lrows,local_ii,myid,NProcess
	contains
	
    subroutine ParallelDumpParticle(InputParticleBundle)
        implicit none
        Type(ParticleBundle),intent(in) :: InputParticleBundle
        Type(ParticleBundle) :: TempParticleBundle
        Integer(4) :: i,j,Dim,Nmax,NTemp
        Character(len=20) :: Filename
        TempParticleBundle%Name=InputParticleBundle%Name
        TempParticleBundle%XFactor=InputParticleBundle%XFactor
        TempParticleBundle%VFactor=InputParticleBundle%VFactor
        Call MPI_ALLGATHER(InputParticleBundle%NPar,1,mpi_integer,snds,1,mpi_integer, MPI_COMM_WORLD,Ierr)
        Call MPI_AllReduce(InputParticleBundle%NPar,Nmax,1,mpi_integer,Mpi_Sum,MPI_COMM_WORLD,Ierr)
        Dim=Sizeof(TempParticleBundle%PhaseSpace(1))/8_4
        Call MPI_Barrier(MPI_COMM_WORLD,ierr)
        
        If (myid/=Master) then
            NTemp=Dim*snds(myid+1)
            Call MPI_SEND(InputParticleBundle%PhaseSpace,NTemp,mpi_double_precision,master,2000+myid,MPI_COMM_WORLD,Ierr)
        end if 
        
        If (myid==Master) then
                Write(*,*) 'Parallel Saving ', trim(TempParticleBundle%Name), Nmax,' Please Wait...'
                Write(filename,*) trim(TempParticleBundle%Name),".dat"
                Open (10,file=filename)
                Write(10,*) Nmax
                TempParticleBundle=InputParticleBundle
                Call UpdatePositionBundle(TempParticleBundle,1_4)
                Call UpdateVelocityBundle(TempParticleBundle,1_4)
                do i=1,TempParticleBundle%NPar
                    Write(10,FMt="(<Dim>D23.15)")TempParticleBundle%PhaseSpace(i)
                end do
                do i=1,Nprocess-1 
					NTemp=Dim*snds(i+1)
					Call MPI_RECV(TempParticleBundle%PhaseSpace,NTemp,mpi_double_precision,i,2000+i,MPI_COMM_WORLD,status,Ierr)
                    TempParticleBundle%NPar=snds(i+1)
                    Call UpdatePositionBundle(TempParticleBundle,1_4)
                    Call UpdateVelocityBundle(TempParticleBundle,1_4)
                    Write(*,*) 'Myid=',i,'Segment Saving',snds(i+1)
                    do j=1,snds(i+1)
                        Write(10,FMt="(<Dim>D23.15)") TempParticleBundle%PhaseSpace(j)
                    end do
				end do
				Close(10) 
		End If
		return
    end subroutine ParallelDumpParticle
     
    subroutine ParallelLoadParticle(InputSpecy,InputParticleBundle,Logistat)
        implicit none
        Type(OneSpecy),intent(in) ::  InputSpecy
        Type(ParticleBundle),intent(out) :: InputParticleBundle
        Logical,intent(out) ::  LogiStat
        Type(ParticleBundle) :: TempParticleBundle
        Logical :: alive
        Character(len=20) :: Filename
        Integer(4) :: i,j,Dim,NMax,NDivide,NMod,NTemp       
        InputParticleBundle%Name=InputSpecy%Name
        Write(filename,*) trim(InputSpecy%Name),".dat"
        Write(*,*) InputSpecy%Name
        If (myid==Master) then
            Inquire(file=filename,exist=alive)
            If(alive) then
                Open (10,file=filename)
                Read(10,*) NMax 
                LogiStat=.False. 
                Write(*,*) 'Parallel Loading ', trim(InputParticleBundle%Name), NMax ,' Please Wait...' 
            else
                Write(*,*) 'Can not find the file for ', trim(InputSpecy%Name),' the particle will be randomly initilalized.'    
                LogiStat=.True.
            end if
		End If
		Call  MPI_BCAST(LogiStat,1,MPI_LOGICAL,master,MPI_COMM_WORLD,Ierr)   
		Call MPI_Barrier(MPI_COMM_WORLD,ierr) 
		If (LogiStat==.False.) Then            
            Call MPI_Bcast(Nmax,1,MPI_Integer,master, MPI_COMM_WORLD,Ierr)
            NDivide=NMax/Nprocess
            NMod=Mod(NMax,Nprocess)
            do i=1,Nprocess
                if (i<=NMod) then
                    snds(i)=NDivide+1
                else
                    snds(i)=NDivide
                end if
			end do 
			InputParticleBundle%NPar=snds(Myid+1)
			Dim=Sizeof(TempParticleBundle%PhaseSpace(1))/8_4
			If (Myid==0) then
                Do i=1,InputParticleBundle%NPar
                    Read(10,*) InputParticleBundle%PhaseSpace(i)
                End do                
                do i=1,Nprocess-1
                    do j=1,snds(i+1)
                        Read(10,*) TempParticleBundle%PhaseSpace(j)
                    End do
                    Write(*,*) 'Myid=',i,'Segment Loading',snds(i+1)
                    Call MPI_SEND(TempParticleBundle%PhaseSpace,Dim*snds(i+1),mpi_double_precision,i,1000+i,MPI_COMM_WORLD,Ierr)
                End do
                Close(10) 
			End If
			If (Myid/=0) then 
				Call MPI_RECV(InputParticleBundle%PhaseSpace,Dim*snds(myid+1),mpi_double_precision,master,1000+myid,MPI_COMM_WORLD,status,Ierr)
			End if
		End If    
		return
	end subroutine  ParallelLoadParticle
	
end module Parallel