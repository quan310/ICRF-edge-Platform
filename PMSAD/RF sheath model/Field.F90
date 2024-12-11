Module FieldModule
	Use TypeModule
	Use Constants
	Use FileIO
	Use Parallel 
	Implicit none	
	Type(GridSmall),private :: BottomGlobal,TopGlobal,LeftGlobal,RightGlobal,AreaGlobal
	Type(Grid) :: ExGlobal,EyGlobal
	Type(Grid) :: PhiGlobal,RhoGlobal,PhiLGlobal,PhiPGlobal
	Type(Grid),private :: LocalRho,LocalPhi
	Real(8),save :: QBottomGlobal,PhiZeroGlobal,QBottome,QBottomion
	Integer(4) :: specie
	contains
   
    subroutine  XYFieldInit(InputNx, InputNy,InputGap,Inputdx,Nelechight)
                Implicit none
                Integer(4),intent(in) :: NSpecy,InputNx, InputNy,InputGap,Nelechight
                Type(PICParticlePara),intent(inout) :: InputParticlePara(0:NSpecy)
                Real(8),intent(in) ::Inputdx
                Integer(4) :: i
                Integer(4),save ::Nx,Ny,Ntotal,gap2,Nelech
                Real(8),save :: dx,GeoFactor,checkerror               
				dx=Inputdx
                Nx=InputNx
                Ny=InputNy
                Ntotal=Nx*Ny
                gap2=InputGap
                GeoFactor=dx*dx/Epsilon
                Nelech=Nelechight				
                AreaGlobal%Name='Area' 
                AreaGlobal%Nx=Ny
                AreaGlobal%dx=dx               
				Call XYLineBoundaryInit(InputGap,AreaGlobal%Value,dx,Nx,Ny)				
                PhiLGlobal%Name='PhiL' 
                PhiLGlobal%Nx=Nx
                PhiLGlobal%Ny=Ny
                PhiLGlobal%dx=dx
                PhiLGlobal%Value=0.d0                
                PhiPGlobal%Name='PhiP' 
                PhiPGlobal%Nx=Nx
                PhiPGlobal%Ny=Ny
                PhiPGlobal%dx=dx
                PhiPGlobal%Value=0.d0               
                BottomGlobal%Name='Bottom' 
                BottomGlobal%Nx=Ny
                BottomGlobal%dx=dx                
                TopGlobal%Name='Top'
                TopGlobal%Nx=Ny
                TopGlobal%dx=dx                
                LeftGlobal%Name='Left'
                LeftGlobal%Nx=Nx
                LeftGlobal%dx=dx               
                RightGlobal%Name='Right'
                RightGlobal%Nx=Nx
                RightGlobal%dx=dx
                ExGlobal%Name='Ex' 
                ExGlobal%Nx=Nx
                ExGlobal%Ny=Ny
                ExGlobal%dx=dx                
                EyGlobal%Name='Ey' 
                EyGlobal%Nx=Nx
                EyGlobal%Ny=Ny
                EyGlobal%dx=dx                
                PhiGlobal%Name='Phi' 
                PhiGlobal%Nx=Nx
                PhiGlobal%Ny=Ny
                PhiGlobal%dx=dx                
                RhoGlobal%Name='Rho' 
                RhoGlobal%Nx=Nx
                RhoGlobal%Ny=Ny
                RhoGlobal%dx=dx                                                 
                LocalRho%Name='LocalRho'
                LocalRho%Nx=Nx-2
                LocalRho%Ny=Ny-1
                LocalRho%dx=dx                
                LocalPhi%Name='LocalPhi'
                LocalPhi%Nx=Nx-2
                LocalPhi%Ny=Ny-1
                LocalPhi%dx=dx
                return
                   
        Entry XYField(NSpecy,InputParticlePara)			
                RhoGlobal%Value(1:Ntotal)=0.d0               			
                BottomGlobal%Value=0.d0
                TopGlobal%Value=0.d0
                LeftGlobal%Value=0.d0
                RightGlobal%Value=0.d0
                QBottomGlobal=0.d0
                do i=0,NSpecy                    
					specie = i				
					Call  AddRhoChiB(Nx,Ny,RhoGlobal%Value,InputParticlePara(i)%Value,InputParticlePara(i)%RhoFactor,QBottomGlobal,QBottome,QBottomion,specie,InputParticlePara(i)%QBottom,dx)
				End do                
                If (NProcess>1) then
					Call ParallelRhoChiB(Nx,Ny,RhoGlobal%Value,QBottomGlobal,QBottome,QBottomion)
                    Call MPI_Barrier(MPI_COMM_WORLD,ierr) 
                End if
				RhoGlobal%Value = RhoGlobal%Value*GeoFactor              				
				call XYLineBoundary(Nx,Ny,BottomGlobal%Value,QBottome,QBottomion,TopGlobal%Value,LeftGlobal%Value,RightGlobal%Value)
                call SORB(Ny,Nx,RhoGlobal%Value,PhiLGlobal%Value,BottomGlobal%Value,TopGlobal%Value,LeftGlobal%Value,RightGlobal%Value,PhiGlobal%Value)				
				Call MPI_Barrier(MPI_COMM_WORLD,ierr)			           
                Call XYPhitoE2(PhiGlobal%Value,ExGlobal%Value,EyGlobal%Value,Nx,Ny,dx,BottomGlobal%Value,TopGlobal%Value,LeftGlobal%Value,RightGlobal%Value)
		return
    End subroutine XYFieldInit            
    subroutine AddRhoChiB(Nx,Ny,Rho,RhoOne,RhoFactor,QBottom,QBottome,QBottomion,specie,QBottomOne,dx)
            Implicit none
            Integer(4),intent(in) :: Nx,Ny,specie
            Real(8),intent(inout) :: Rho(Nx,Ny),RhoOne(Nx,Ny)
            Real(8),intent(in) :: RhoFactor,dx,QBottomOne
            Real(8),intent(inout) ::  QBottom,QBottome,QBottomion
            Integer(4)::i,j
			
            RhoOne=RhoOne*RhoFactor
            Rho=Rho+RhoOne			            			
            QBottom = QBottom+QBottomOne*RhoFactor*dx*dx     
			if (specie == 0) then
				QBottome = QBottome+QBottomOne*RhoFactor*dx*dx
			else
				QBottomion = QBottomion+QBottomOne*RhoFactor*dx*dx
			end if
            return
    End subroutine AddRhoChiB
    
	subroutine ParallelRhoChiB(Nx,Ny,Rho,QBottom,QBottome,QBottomion)
		Implicit none
        Integer(4),intent(in) :: Nx,Ny
        Real(8),intent(inout) :: Rho(Nx,Ny),QBottom,QBottome,QBottomion
        Integer(4) :: Nxy
        Real(8) :: Temp(Nx,Ny),Tbom
		Nxy=Nx*Ny
        Call MPI_ALLReduce(Rho,Temp,Nxy,mpi_double_precision,MPI_SUM,MPI_COMM_WORLD,Ierr)   
        Rho=Temp
        Call MPI_ALLReduce(QBottom,TBom,1,mpi_double_precision,MPI_SUM,MPI_COMM_WORLD,Ierr)    
        QBottom=TBom
        return 
    End  subroutine ParallelRhoChiB 
         
    subroutine XYLineBoundaryInit(InputGap,Area,dx,Nx,Ny)
            Use RFSource 
            Use Input,only: RFFrequency,Inputdt
            Implicit none			
            Integer(4),intent(in) :: Nx,Ny,InputGap
            Real(8),intent(out) :: Bottom(Ny),Top(Ny),Left(Nx),Right(Nx),Area(Ny)
            Integer(4),save :: Gap,Timer2=0,Timer3=0,Period
            Integer(4) :: i,NLeft
            Real(8),intent(in) :: dx
            Real(8) ::  VolX1,VolX2,VolY1,VolY2,Vsh(Ny,Int(1.d0/RFFrequency(1)/inputdt)),TopV,QBottomion,QBottome,TBome,TBomion,Boeee,Line(Ny)
            Real(8),save :: BottomV=0.0			
            Gap=InputGap
			Period=Int(1.d0/RFFrequency(1)/inputdt)	 
            Vsh=0.0			 
            Area(1)=0.25d0*PI*dx*dx			
            Do i=2,Ny
                Area(i)=2.d0*PI*(i-1)*dx*dx
            End do       
            Return 
			
			Entry XYLineBoundary(Nx,Ny,Bottom,QBottome,QBottomion,Top,Left,Right)         
            VolX1=1.d0
            VolX2=0.d0
            VolY1=0.d0
            VolY2=0.d0
			Left=0.0
			Right=0.0
            Bottom=VolX1
            TopV=600.0*DSin(2*PI*RFFrequency(1)*inputdt*Dble(Timer2))
			Timer2=Timer2+1
			Line(1:30)=1.0
            Do i=30,Ny  
                Line(i)=1.0+DLog(Dble(i-1)/Dble(2-1))/DLog(Dble(Ny-1)/Dble(2-1))*(0.d0-1.0)
            End do
			Do i=1,Ny
                Top(i)=Line(i)*TopV    
            End do 
			if(timer2==period) then
				Call MPI_ALLReduce(QBottome,TBome,1,mpi_double_precision,MPI_SUM,MPI_COMM_WORLD,Ierr)    
				QBottome=TBome
				Call MPI_ALLReduce(QBottomion,TBomion,1,mpi_double_precision,MPI_SUM,MPI_COMM_WORLD,Ierr)    
				QBottomion=TBomion
				timer3=timer3+1				
				if (timer3<0) then
					BottomV=0.0
				else				
					if(QBottomion==0.0) then
						else
							Boeee=(QBottomion+QBottome)/QBottomion
							if(Boeee<-1.0) then
								BottomV=BottomV-100.0
								else if(Boeee<-0.7) then
								BottomV=BottomV-50.0
								else if(Boeee<-0.6) then
								BottomV=BottomV-30.0
								else if(Boeee<-0.5) then
								BottomV=BottomV-10.0
								else if(Boeee<-0.4) then
								BottomV=BottomV-5.0
								else if(Boeee<-0.3) then
								BottomV=BottomV-3.0
								else if(Boeee<-0.2) then
								BottomV=BottomV-2.0
								else if(Boeee<-0.1) then
								BottomV=BottomV-1.0
								else if(Boeee<-0.07) then
								BottomV=BottomV-0.5
								else if(Boeee<-0.05) then
								BottomV=BottomV-0.1
								else if(Boeee>-0.05 .and. Boeee<0.05) then
								BottomV=BottomV-0.0
								else if(Boeee<0.07) then
								BottomV=BottomV+0.1
								else if(Boeee<0.1) then
								BottomV=BottomV+0.5
								else if(Boeee<0.2) then
								BottomV=BottomV+1.0
								else if(Boeee<0.3) then
								BottomV=BottomV+3.0
								else if(Boeee<0.4) then
								BottomV=BottomV+5.0
								else if(Boeee<0.5) then
								BottomV=BottomV+10.0
								else
								BottomV=BottomV+20.0
							end if
						end if
					end if				
				Timer2=0
				write(*,*)  '1111111111',QBottomion,QBottome,BottomV
				QBottomion=0.0
				QBottome=0.0
			end if
			Top=Top-BottomV           
			Do i=1,Ny             
                Bottom(i)=0.0                 
            End do    
        return
    End subroutine  XYLineBoundaryInit                                                                                 
    subroutine SORB(Ny,Nx,localrho,localphi,Bottom,Top,Left,Right,phit1)
		Implicit none		
        Integer(4),intent(in) :: Nx,Ny
        Real(8),intent(in) :: localrho(Nx,Ny)
		Real(8),intent(in) :: Bottom(Ny),Top(Ny),Left(Nx),Right(Nx)				
        Real(8),intent(out) ::localphi(Nx,Ny)		
        Integer(4) i,j,k,l,num
        real(8),parameter :: limit = 1e-4		
        real(8),parameter :: omg = 1.5				
        real(8) :: erro,erromax
        real(8) :: phi(Nx,Ny),phi_old(Nx,Ny),phit1(Nx,Ny) 
        character(40) :: name="phi"

        phi=0.d0   
		num=0
		do i=30,Ny 
            phi(1,i)=Bottom(i) 
        end do		
		do j=1,Ny
            phi(Nx,j)=Top(j)
        end do
        do k=1,5000
			erromax=0.d0
			phi_old=phi 
			num=num+1       
			do i=2,Nx-1
				do j=2,Ny-1 				
 			        phi(i,j)=(phi(i+1,j)+phi(i-1,j)+phi(i,j+1)+phi(i,j-1)+localrho(i,j))/4.0
                    phi(i,j)=(1.0-omg)*phi_old(i,j)+omg*phi(i,j)  					
				end do
			end do			
			do j=1,30					
				phi(1,j)=phi(2,j) 		 
			end do
			do j=31,Ny					
				phi(1,j)=Bottom(j)
			end do
			do j=1,Ny
				phi(Nx,j)=Top(j)
			end do
			do i=1,Nx
				phi(i,1)=phi(i,2)		
			end do
			do i=1,Nx
				phi(i,Ny)=0.0
			end do			
			phi_old=phi 
			do i=Nx-1,2,-1
				do j=Ny-1,2,-1                    
 			        phi(i,j)=(phi(i+1,j)+phi(i-1,j)+phi(i,j+1)+phi(i,j-1)+localrho(i,j))/4.0
                    phi(i,j)=(1.0-omg)*phi_old(i,j)+omg*phi(i,j)                            		
				end do
			end do 
			do j=1,30
				phi(1,j)=phi(2,j)  
			end do
			do j=31,Ny
				phi(1,j)=Bottom(j)
			end do
			do j=1,Ny
				phi(Nx,j)=Top(j)
			end do
			do i=1,Nx
				phi(i,1)=phi(i,2)			
			end do
			do i=1,Nx
				phi(i,Ny)=0.0
			end do
            erro=0.d0
			erromax=0.d0       
			do i = 2 , Nx-1
				do j = 2 , Ny-1               
					if (abs(phi(i,j))/=0.d0)   erro=abs(phi(i,j)-phi_old(i,j))/abs(phi(i,j))
					if (erromax<erro)   erromax=erro+erromax
				end do
			end do           
		end do		
		localphi=phi
		phit1=phi    
	end subroutine SORB 	
              
    Subroutine XYPhitoE2(PhiSum,Ex,Ey,Nx,Ny,dx,Bottom,Top,Left,Right)
        Implicit none
        Integer(4), intent(in) ::  Nx,Ny
        Real(8), intent(in) ::  dx,Bottom(Ny),Top(Ny),Left(Nx),Right(Nx)
        Real(8), intent(out) ::  Ex(Nx,Ny),Ey(Nx,Ny)
        Real(8), intent(inout) ::  PhiSum(Nx,Ny) 
        Real(8) :: PhiSum_Temp(Ny,Nx) 
        Integer(4) :: i,j
               
        PhiSum_Temp=TRANSPOSE(PhiSum) 
        do i=2,Ny-1
            do j=2,Nx-1
                Ex(j,i)=(PHISum(j-1,i)-PHISum(j+1,i))/2.d0/dx
            end do
        end do               
        do i=2,Nx-1
            do j=2,Ny-1
                Ey(i,j)=(PHISum_Temp(j-1,i)-PHISum_Temp(j+1,i))/2.d0/dx
            end do
        end do        
        do i=2,Nx-1
                   Ex(i,1)=Ex(i,2)
                   Ey(i,1)=0
                   Ex(i,Ny)=Ex(i,Ny-1)
                   Ey(i,Ny)=Ey(i,Ny-1) 
        enddo
		do i=2,Ny-1
                   Ex(1,i)=Ex(2,i)
                   Ey(1,i)=Ey(2,i)
                   Ex(Nx,i)=Ex(Nx-1,i)
                   Ey(Nx,i)=Ey(Nx-1,i)
        enddo    
        Ey(1,1)=(Ey(1,2)+Ey(2,1))/2.d0!0
        Ex(1,1)=(Ex(1,2)+Ex(2,1))/2.d0
        Ey(1,Ny)=(Ey(1,Ny-1)+Ey(2,Ny))/2.d0
        Ex(1,Ny)=(Ex(1,Ny-1)+Ex(2,Ny))/2.d0
        Ey(Nx,1)=(Ey(Nx,2)+Ey(Nx-1,1))/2.d0 !0!
        Ex(Nx,1)=(Ex(Nx,2)+Ex(Nx-1,1))/2.d0
        Ey(Nx,Ny)=(Ey(Nx-1,Ny)+Ey(Nx,Ny-1))/2.d0
        Ex(Nx,Ny)=(Ex(Nx-1,Ny)+Ex(Nx,Ny-1))/2.d0
		return
    end subroutine XYPHItoE2
	
end  Module FieldModule


              
              