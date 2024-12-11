Module Constants
    Implicit none
    Real(8),parameter :: PI=3.141592653589793238D0
    Real(8),parameter :: a0=5.2918d-11
    Real(8),parameter :: lightspeed=3.0d8 
    Real(8),parameter :: kB=1.3807d-23
    Real(8),parameter :: exp=2.71828d0
    Real(8),parameter :: ElectronCharge=1.6022d-19
    Real(8),parameter :: ElectronMass=9.1095d-31
    Real(8),Parameter :: Epsilon=8.8542d-12
    Real(8),Parameter :: JtoeV=1.6022d-19
    Real(8),Parameter :: eVtoK=11605.d0
    Real(8),Parameter :: mTorrtoPa=0.13332d0
    Real(8)  :: R
    contains
    
	Subroutine DRandomInt(AA)
        USE IFPORT
        Implicit none
        Real(8),intent(out) :: RR
        Integer(4),intent(in) :: AA
        Integer(4) :: AA2
        Real(8) :: Temp 
        AA2=AA+1000
        Temp =DRAND(AA2)
    Return 
    Entry  DRandom(RR)
        RR=DRAND(0)
    Return
    End  subroutine DRandomInt 
	
End  Module Constants

Module TypeModule
	Implicit none
                Type  OneSpecy      
                       Character(len=16) Name
                       Integer(4) :: NAtom,Charge
                       Real(8) :: Mass
                End  Type  OneSpecy                
                Type GasPhysical    
                     Integer(4) ::  Model,NAtom,NSpecy
                     Real(8) ::  MGas,Radius,BetaMax 
                End Type  GasPhysical                                      
                Integer(4),parameter :: NParMax=4000000_4
                Integer(4),parameter :: NParMaxSmall=Int(0.05*dble(NParMax))              
                Type ParticleOne
                    Real(8) ::  X,Y,Vx,Vy,Vz,Ax,Ay
                EndType ParticleOne                  
                Type  ParticleBundle
                        Character(len=16) Name
                        Real(8) ::  XFactor,VFactor,Charge,Mass
                        Integer(4) :: NPar
                        Type(ParticleOne) ::  PhaseSpace(NParMax)
                End  Type  ParticleBundle
                Integer(4),parameter,private :: GridMax=129*1025
                Integer(4),parameter,private :: GridMaxSmall=3000_4 
                Type Grid  
					Character(len=16) Name
					Integer(4) :: Nx,Ny
					Real(8) :: Dx,Dy 
					Real(8) :: Value(GridMax)
                End  Type Grid
				Type GridSmall
                    Character(len=16) Name
                    Integer(4) :: Nx
                    Real(8) :: Dx
                    Real(8) :: Value(GridMaxSmall)
                End  Type GridSmall
                Type PICParticlePara
                    Integer(4) ::  SubCycle,Timer
                    Real(8) ::  RhoFactor,AccFactor,AccFactorB,Wq,QFactor
                    Real(8) ::  Value(GridMax)   ! For local rho and Electric field storage. 
                    Real(8) ::  Temp(GridMax)
                    Real(8) ::  Qbottom,QbottomTemp
                    Real(8) ::  Qtop,QtopTemp
                    Real(8) ::  BxfieldTemp(GridMax),ByfieldTemp(GridMax),BzfieldTemp(GridMax)
                    Real(8) ::  Omegax(GridMax),Omegay(GridMax),Omegaz(GridMax)
                    Real(8) ::  TransB11(GridMax),TransB12(GridMax),TransB13(GridMax),TransB21(GridMax),TransB22(GridMax),TransB23(GridMax),TransB31(GridMax),TransB32(GridMax),TransB33(GridMax)
                End Type PICParticlePara 
				Type Gas
                    Integer(4) ::  NSpecy,Shift
                    Real(8) :: MGas,TGas,PGas,BetaMax
                End Type  Gas                          
                                                                                                                                                                                                                                    
End Module TypeModule