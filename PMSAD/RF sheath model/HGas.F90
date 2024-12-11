Module HGasModule
    Use TypeModule
    Use Constants	
	Integer(4),Private,parameter :: NSpecyAr=1_4
	Type(GasPhysical),parameter :: ArGas=(GasPhysical(2,1,1,1.674d-27,0.d0,0.d0))    
	Type(OneSpecy),parameter :: ArSpecy(1:NSpecyAr)=(OneSpecy('Ar+',1,1,1.674d-27))                    
	contains 	
end module HGasModule

Module ElectronModule
	Use TypeModule
    Implicit none
    Type(OneSpecy),parameter :: Electron=(OneSpecy('Electron',1,-1,9.1095d-31))
    Type(Gas),parameter :: ElectronGas=(Gas(1,0,9.1095d-31,200.d0*11605.d0,0.d0,0.d0)) 
End Module ElectronModule
