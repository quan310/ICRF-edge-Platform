Module Field
	Use MPI
    Implicit none

	integer,parameter :: MAX_LOCAL_ROW=2048,MAX_LOCAL_COL=2048,MAX_PROC=128
	integer,parameter :: RHO_FIRST_PHASE=1001,RHO_SECOND_PHASE=1002,RHO_THIRD_PHASE=1003
	integer :: snds(MAX_PROC),disp(MAX_PROC)
	integer :: lrows(MAX_PROC),local_ii(MAX_PROC)
	Integer :: myid,ierr,myid2,NProcess
	integer :: master=0
	integer :: status(MPI_STATUS_SIZE)
	Integer :: SerialCOMM,Ndim=1
	Integer :: dims(1)  
	logical :: ParaPeriod(1)=.False.,Reorder=.True.  
	save snds,disp,lrows,local_ii,myid2,NProcess
    Contains
	
    subroutine ParallelDump(Name,NPar,Particle)
        implicit none
        Integer,intent(in) :: NPar
        integer :: i,j,NSum
		integer ::NMax
		real(8) :: Particle(5,NPar) 
        real(8),allocatable :: SumPar(:,:)
        character(len=30) name,Filename
        logical alive
        write(filename,*) trim(Name),".dat"
		
		Call MPI_ALLGATHER(NPar,1,mpi_integer,snds,1,mpi_integer, SerialCOMM,Ierr)        
		Nmax=0		
        do i = 0, Nprocess-1
            if  (Nmax<snds(i+1)) then
                Nmax=snds(i+1)
            end if
        end do
		
        Call MPI_Reduce(NPar,NSum,1,mpi_integer,mpi_sum,master, MPI_COMM_WORLD,Ierr)
		
        if (myid2/=0) then
            Call MPI_SEND(Particle(1:5,1:snds(myid2+1)),5*snds(myid2+1),mpi_double_precision,master,2000+myid2,SerialCOMM,Ierr)
        end if		
        if (myid2==0) then
            open (10,file=filename)			
            write(*,*) 'Nsum ',Name, NSum
            write(*,*) 'NPar ',Name,NPar,'myid2=',myid2
            write(10,*) NSum
			
            do j=1,NPar
                Write(10,FMt="(5D23.15)") Particle(1,j),Particle(2,j),Particle(3,j),Particle(4,j),Particle(5,j)
            end do 			
            allocate (SumPar(5,Nmax))			
            do i=1,Nprocess-1 
                Write(*,*) 'NPar ', Name,snds(i+1),'myid2=',i
                Call MPI_RECV(SumPar(1:5,1:snds(i+1)),5*snds(i+1),mpi_double_precision,i,2000+i,SerialCOMM,status,Ierr)
                
				do j=1,snds(i+1)
                    Write(10,FMt="(5D23.15)") SumPar(1,j),SumPar(2,j),SumPar(3,j),SumPar(4,j),SumPar(5,j)
                end do				
            end do			
            deallocate (SumPar) 
            close(10)
        end if
		
        return
    end subroutine ParallelDump

	subroutine ParallelLoad(Name,NPar,NParMax,Particle)
		implicit none
		Integer :: NPar,NParMax
		real(8) :: Particle(5,NParMax)
		real(8),allocatable :: SumPar(:,:)
		integer, allocatable :: snds(:)    
		integer :: i,j,NMax
		integer :: NDivide,NMod
		character(len=30) name,Filename
		logical alive
		write(filename,*) trim(Name),".dat"
		
		if (myid2 == 0) then
			inquire(file=filename, exist=alive)       
			if (alive) then
				open (10, file=filename)
				read(10,*) NMax
				Write(*,*) NMax, name, 'total number loaded.'
			else
				Write(*,*) "Warning: No ", trim(Name), " are loaded."
			end if
		end if   
		Call MPI_BCAST(NMax, 1, mpi_integer, 0, SerialCOMM, ierr)
		
		NDivide = NMax / (Nprocess-1)
		NMod = Mod(NMax, Nprocess-1)
		allocate(snds(Nprocess-1))
		snds(1) = 0

		do i = 1, Nprocess-1
			if (i <= NMod) then
				snds(i) = NDivide + 1
			else
				snds(i) = NDivide
			end if
		end do
		NPar = 0
		allocate(SumPar(5, NMax))
		if (myid2 == 0) then
			do i = 1, Nprocess-1
				do j = 1, snds(i)
					read(10, fmt="(5D23.15)") SumPar(1, j), SumPar(2, j), SumPar(3, j), SumPar(4, j), SumPar(5, j)
				end do
			end do
		end if
		Call MPI_BCast(SumPar, 5*NMax, mpi_double_precision, 0, SerialCOMM, ierr)
		do i = 1, snds(myid2)
			NPar = NPar + 1
			Particle(:, NPar) = SumPar(:, i)
		end do
		deallocate(SumPar)
		if (myid2 == 0) then
			close(10)
		end if
		return
	end subroutine ParallelLoad

    Subroutine Dump1D(Name,Nre,Res1)
        Implicit none
        Integer i,Nre
        Real(8) :: Res1(Nre)       
        Character(len=30) name,Filename
        Write(filename,*) trim(Name),".dat"
        Open (10,file=filename, status='REPLACE')
        Write(10,*) Nre
        
		do i=1,Nre
            Write(10,FMt="(D21.13)") Res1(i)                           
        End do
		
        Close(10)  
        Return
    End Subroutine Dump1D
    
    subroutine Dump2D(Name,Nx,Ny,Res)
        implicit none
        integer i,j,Nx,Ny
        real(8) :: Res(Nx,Ny)
        character(len=30) name,Filename
        logical alive
        write(filename,*) trim(Name),".dat"
        Write(*,*) "Saving ",trim(Name)," Please wait..."
        open (10,file=filename)
  
        do i=1,Ny
            Write(10,FMt="(D23.15\)")   (Res(j,i),j=1,Nx)
            Write(10,*)  
        end do
		
        close(10)
        Write(*,*) "Save ",trim(Name),"Complete!"  
        return
    end subroutine Dump2D
	
    subroutine sputteratom(NAr, Ar, Nsp, NsputterMax,sp)
    
		implicit none    
		real(8), parameter :: PI = 3.1415926535897932384626433832795		
		real(8), parameter :: MH = 1.66d-27, Mtarget = 3.06D-25
		real(8), parameter :: Zi = 1.0, Zt = 74.0				
		real(8), parameter :: Mi = 1.0079D0, Mt = 183.84D0		
		real(8), parameter :: Us = 8.79D0		
		real(8), parameter :: gamma = 4.0*Mi*Mt/((Mi+Mt)**2), Q_sput = 2.99D0		
		Real(8),parameter :: JtoeV=1.6022d-19			
		Real(8) :: sp(1:5,1:NsputterMax),Ar(1:5,1:NAr)									
		Real(8) :: Energy0,VxTemp,VyTemp,VzTemp,v,omige,ArX,ArY
		Real(8) :: V11z,V11y
		real(8) :: Ecut,epsilon_sput, alpha_star, K_sput, s_e, s_n, sp_rate
		real(8) :: r5, rr, r3, r4, mm, nn, r1, r2, bb, integral_prev, theta, phi
		real(8) :: E_step, dEnergy, Emax, ymax, v0, vx0, vy0, vz0, e_u, E_sput,vx1, vy1, vz1		
		real(8), dimension(100)  :: ETerpoint, E_diff, rect_area, rect_prob						
		Integer, intent(out) :: Nsp												
		Integer :: n, NAr , ii, kk, m, i, j, pp, n_samples, count_numrect, N_sput,NsputterMax
		Real(8), dimension(3) :: normal	
		Real(8) :: magnitude_V,magnitude_N,dotproduct,cos_theta,alpha,YE0,up,down,angle_para
		Real(8), parameter :: Q = 7.1574D-3,lambda = 1.3849d0,miu = 8.5638d-1, epsilon_L = 1.5136d-4,E_th = 4.8731d2
		Real(8) ::epsilon_1		
		Real(8), parameter :: f = 1.3323d0,b = 1.5132d-1,c = 1.0560d0,E_sp = 1.0D0
		Real(8) :: alpha_zero		
		
		call random_seed()
		
		normal = (/ 1.0, 0.0, 0.0 /)
		do n = 1, int(NAr)																		
			Energy0=0.0        
			VxTemp = Ar(3,n)														
			VyTemp = Ar(4,n)
			VzTemp = Ar(5,n)
			v=VxTemp*VxTemp+VyTemp*VyTemp+VzTemp*VzTemp											
			Energy0=MH*v/2.0/JtoeV												
			
			magnitude_V = sqrt(v)
			magnitude_N = sqrt(normal(1)**2+normal(2)**2+normal(3)**2)
			dotproduct = VxTemp*normal(1) + VyTemp*normal(2) + VzTemp*normal(3)
			cos_theta = dotproduct / (magnitude_V * magnitude_N)
			alpha = PI - acos(cos_theta)			
			
			if (Energy0 <= E_th) then																
			else																					
				Ecut = int(gamma*Energy0 - Us)				
				Write(*,*) 'Ecut=',Ecut
				if (Ecut<= 0.d0) then
					Ecut = 1
				end if
				
				epsilon_1 = epsilon_L*Energy0
				up = ((Energy0/E_th-1.0d0)**miu)*log(1.0d0+epsilon_1*1.2288d0)
				down = lambda + ((Energy0/E_th-1.0d0)**miu)*(epsilon_1+0.1728d0*sqrt(epsilon_1)+0.008d0*epsilon_1**0.1504d0)				
				YE0 = 0.5*Q*up/down			
				alpha_zero = PI - acos(sqrt(1.0d0 / (1.0d0 + Energy0 / E_sp)))
				angle_para = (cos((alpha*PI/alpha_zero/2.0d0)**c))**(-f)*exp(b*(1.0d0-1.0d0/cos((alpha*PI/alpha_zero/2.0d0)**c)))
				sp_rate = angle_para*YE0									

				call random_number(r5)									
				if (r5 < fraction(sp_rate)) then	
					N_sput = floor(sp_rate) + 1		
				else
					N_sput = floor(sp_rate)				
				end if

				ymax = 0.0D0
				E_step = 0.1D0
				m = (Ecut - 0.0D0) / E_step							
				
				do kk = 1, m
					dEnergy = E_step * kk		
					bb = y_func(dEnergy)							
					if (bb > ymax) then
						ymax = bb				
						Emax = dEnergy			
					end if			
				end do				
				eterpoint(1) = Emax		

				do								
					integral_prev = 0.0							
					do j = 1, m
						if ((e_step*j) > eterpoint(1)) exit		
						integral_prev = integral_prev + y_func(e_step*j) * e_step						
					end do					
					if (eterpoint(1) * ymax > 2.0 * integral_prev ) exit
						eterpoint(1) = eterpoint(1) + e_step						
				end do
								
				count_numrect = 1   
				do i = 2, 100
					eterpoint(i) = eterpoint(i-1) + 0.1							
					do
						integral_prev = 0.0							
						do j = 1, m      
							if (eterpoint(i-1) + e_step*j > eterpoint(i)) exit										
								integral_prev = integral_prev + y_func(eterpoint(i-1) + e_step*j) * e_step			
						end do						
						if ((eterpoint(i) - eterpoint(i-1)) * y_func(eterpoint(i-1)) > 2.0 * integral_prev ) exit	
							eterpoint(i) = eterpoint(i) + e_step													
					end do					
					if (eterpoint(i-1) + e_step*j > Ecut) exit					
						count_numrect = count_numrect + 1												
				end do
				      
				e_diff(1) = eterpoint(1)
				rect_area(1) = e_diff(1) * ymax				
				do i = 2, count_numrect
					e_diff(i) = eterpoint(i) - eterpoint(i-1)
					rect_area(i) = e_diff(i) * y_func(eterpoint(i-1))
				end do				   
				do i = 1, count_numrect
					rect_prob(i) = rect_area(i) / sum(rect_area)
				end do				
				do i = 2, count_numrect
					rect_prob(i) = rect_prob(i) + rect_prob(i-1)
				end do
				
				do ii = 1, N_sput        												
					outerloop: do       																		
						call random_number(rr)
			
						if ( rr < rect_prob(1)) then									     
							call random_number(r1)
							call random_number(r2)							
							e_u = ETerpoint(1) * r1
								
							if ((r2 * ymax) <= y_func(e_u)) then																
								E_sput = e_u
								exit outerloop
							end if								
						else															     
							innerloop: do i = 2, count_numrect							
								if ( rr < rect_prob(i)) then  
									call random_number(mm)
									call random_number(nn)
									e_u = ETerpoint(i-1) + E_diff(i) * nn
									if (mm * y_func(ETerpoint(i-1)) <= y_func(e_u)) then									
										E_sput = e_u																
										exit outerloop				
									end if
									exit
								end if
								
							end do innerloop							
						end if		
					end do outerloop
					          
					call random_number(r3)
					call random_number(r4)					
					theta = acos((r3)**0.1)										
					phi = 2.0*pi*r4										
					v0 = sqrt(2.0*E_sput*1.602d-19/Mtarget)			 					
					vx1 = 1.0*v0*cos(theta)  					
					vy1 = v0*sin(theta)*sin(phi)
					vz1 = v0*sin(theta)*cos(phi) 															
					Nsp = Nsp + 1
					sp(1,Nsp) = Ar(1,n)							
					sp(2,Nsp) = Ar(2,n)
					sp(3,Nsp) = vx1
					sp(4,Nsp) = vy1
					sp(5,Nsp) = vz1					   
				end do					
			end if
		end do	   
		contains
		
		function y_func(energy1)				!
			Implicit none
			real(8):: energy1
			real(8) :: y_func
			y_func= (1.0d0 - sqrt((us + energy1) / (gamma * energy0))) / (energy1**2 * (1.0d0 + us / energy1)**3)		
			RETURN
		end function y_func
		
	end subroutine sputteratom
		
	Subroutine Move(Nsp,sp,Nx,Ny,dy,dt,Nspu,Ndep,dep)		
		Implicit none		
		Integer,intent(in) :: Nsp,Nx,Ny,Nspu 
		Integer i,j,Yheight,Xheight
		Integer :: Ndep       
		real(8),intent(in)  :: dy,dt
		real(8),intent(in)  :: sp(1:5,1:Nsp)		
		Real(8) t, Xsp,Ysp		
		real(8) :: dep(1:5,1:Nsp*Nspu)		
		Yheight = dble(Ny-1)
		Xheight = dble(Nx-1)
		
		do j = 1,Nspu		
			do i = 1, Nsp			
				Xsp = sp(3,i)*dt*j/dy
				Ysp = sp(2,i)+sp(4,i)*dt*j/dy			
				if (Xsp > 0.d0 .and. Xsp <= Xheight) then       					
					if(Ysp > 0.d0 .and. Ysp <= Yheight) then 							
						Ndep = Ndep + 1
						dep(1,Ndep) = Xsp
						dep(2,Ndep) = Ysp
						dep(3,Ndep) = sp(3,i)
						dep(4,Ndep) = sp(4,i)
						dep(5,Ndep) = sp(5,i)						
					End If					
				End If			
			END DO
		End do
	  
		Return
    End Subroutine Move

	subroutine Weighting2D(Ndepos2D,depos2D,Nx,Ny,Rho)
        implicit none        
		Integer(4),intent(in) :: Nx,Ny,Ndepos2D
		real(8),intent(in) :: depos2D(1:5,1:Ndepos2D)
        real(8),intent(out) :: Rho(Nx,Ny)		
        real(8) :: X1,X2,Y1,Y2
        Integer(4) :: i,N1,N2
        Rho=0.d0		
        do i=1,Ndepos2D		
            N1=Ceiling(depos2D(1,i))
            X1=Dble(N1)-depos2D(1,i)
            X2=1.d0-X1			
            N2=Ceiling(depos2D(2,i)) 
            Y1=Dble(N2)-depos2D(2,i)
            Y2=1.d0-Y1			
            Rho(N1,N2)=Rho(N1,N2)+X1*Y1
            Rho(N1,N2+1)=Rho(N1,N2+1)+X1*Y2
            Rho(N1+1,N2)=Rho(N1+1,N2)+X2*Y1
            Rho(N1+1,N2+1)=Rho(N1+1,N2+1)+X2*Y2			
        end do
        return
	end subroutine Weighting2D
	
	Subroutine YADFion(Nsp,NFA,Ny,YADF,sp)	 
        Implicit none        
        real(8), parameter :: Mion = 1.66d-27	
        Integer :: j,N,NFA,Ny,NN     	
		Integer, intent(in) :: Nsp					
		real(8), intent(in) :: sp(5,Nsp)			
		Real(8) :: AnInterval = 1.d0
		Real(8) :: yInterval = 1.d0
        Real(8) :: YADF(Ny,NFA),sputterenergy,x1,x2,EE,s1,s2,yy
        real(8) :: Vx1,Vy1,Vz1,V
		Real(8) :: magnitude_V,magnitude_N,cos_theta,alpha
		REAL(8),dimension(3)::normal 
        real(8),parameter :: JtoeV=1.6022d-19
		real(8), parameter :: PI = 3.1415926535897932384626433832795
        
		normal = (/ 1.0, 0.0, 0.0 /)
		YADF=0.d0
		
        do j = 1,Nsp		
			Vx1 = sp(3,j)
			Vy1 = sp(4,j)
			Vz1 = sp(5,j)
			V = Vx1*Vx1+Vy1*Vy1+Vz1*Vz1			
			magnitude_V = sqrt(V)
			magnitude_N = sqrt(normal(1)**2+normal(2)**2+normal(3)**2)
			cos_theta = (Vx1*normal(1) + Vy1*normal(2) + Vz1*normal(3)) / (magnitude_V * magnitude_N)
			alpha = (PI - acos(cos_theta))*180/PI						
			EE = alpha/AnInterval
            N = Ceiling(EE)			
			yy = sp(2,j)/yInterval
			NN = Ceiling(yy)
			
			If( N> 0 .and. N <= NFA ) then					
				If( NN> 0 .and. NN <= Ny ) then
				
					X1 = dble(N)-EE
					X2 = 1.d0-X1									
					S1 = dble(NN)-yy								
					S2 = 1.d0-S1					
					YADF(NN,N) = YADF(NN,N)+S1*X1
					YADF(NN+1,N+1) = YADF(NN+1,N+1)+S2*X2
					YADF(NN,N+1) = YADF(NN,N+1)+S1*X2	
					YADF(NN+1,N) = YADF(NN+1,N)+S2*X1	
				
				end if			         
            End If
           
        End do
		
        Return
    End Subroutine YADFion

	Subroutine YEDFion(Nsp,NFE,Ny,YEDF,sp) 
        Implicit none        
        real(8), parameter :: Mion = 1.66d-27	!
        Integer :: j,N,NFE,Ny,NN     	
		Integer, intent(in) :: Nsp					
		real(8), intent(in) :: sp(5,Nsp)			
		Real(8) :: EnInterval = 10.d0
		Real(8) :: yInterval = 1.d0
        Real(8) :: YEDF(Ny,NFE),sputterenergy,x1,x2,EE,s1,s2,yy
        real(8) :: Vx1,Vy1,Vz1,V
        real(8),parameter :: JtoeV=1.6022d-19        
		YEDF=0.d0
		
        do j = 1,Nsp		
			Vx1 = sp(3,j)
			Vy1 = sp(4,j)
			Vz1 = sp(5,j)
			V = Vx1*Vx1+Vy1*Vy1+Vz1*Vz1
			sputterenergy = Mion*V/2.0/JtoeV						
			EE = sputterenergy/EnInterval
            N = Ceiling(EE)			
			yy = sp(2,j)/yInterval
			NN = Ceiling(yy)
			
			If( N> 0 .and. N <= NFE ) then					
				If( NN> 0 .and. NN <= Ny ) then				
					X1 = dble(N)-EE
					X2 = 1.d0-X1									
					S1 = dble(NN)-yy							
					S2 = 1.d0-S1					
					YEDF(NN,N) = YEDF(NN,N)+S1*X1
					YEDF(NN+1,N+1) = YEDF(NN+1,N+1)+S2*X2
					YEDF(NN,N+1) = YEDF(NN,N+1)+S1*X2	
					YEDF(NN+1,N) = YEDF(NN+1,N)+S2*X1	
				
				end if		         
            End If           
        End do		
        Return
    End Subroutine YEDFion
		
	Subroutine DWeightingy(Ndep,dep,Nyy,Rhodep)
		Implicit none
		Integer i,N
		Integer,intent(in) :: Nyy,Ndep
		Real(8),intent(in) :: dep(5,Ndep)		
		Real(8),intent(out) :: Rhodep(Nyy)
		Real(8) :: S1,S2
		Rhodep=0.d0
		
		do i=1,Ndep											
			N = Ceiling(dep(2,i))
			S1=dble(N)-dep(2,i)								
			S2=1.d0-S1			
			if(N<=0) then
			elseif(N>=dble(Nyy-1)) then
			else
				Rhodep(N)=Rhodep(N)+S1						
				Rhodep(N+1)=Rhodep(N+1)+S2
			end if			
		End do
		
		Return
    End Subroutine DWeightingy
	
	Subroutine SWeighting(Nsp,sp,Nyy,Rhosp)
        Implicit none
        Integer i,N
        Integer,intent(in) :: Nyy,Nsp
        Real(8),intent(in) :: sp(5,Nsp)
        Real(8),intent(out) :: Rhosp(Nyy)
        Real(8) :: S1,S2	
        Rhosp=0.d0
		
        do i=1,Nsp           
            N=Ceiling(sp(2,i))
			S1=dble(N)-sp(2,i)
            S2=1.d0-S1
			if(N<=0) then
				N=1
				S1=1.d0
				S2=0.d0
			elseif(N>=dble(Nyy-1)) then
				N=dble(Nyy-1)-1
				S1=0.d0
				S2=1.d0
			end if
            rhosp(N)=rhosp(N)+S1
            rhosp(N+1)=rhosp(N+1)+S2			
        End do
	
        Return
    End Subroutine SWeighting    
 
End Module Field