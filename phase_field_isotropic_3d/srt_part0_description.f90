
!====================================================================
!          Program for phase field damage coupled with 
!          mechanical loading and hydrogen diffusion 
!          Mechanical model: standard Hooke's law elasticity 
!                            isotropic Von Mises plasticity
!          Hydrogen diffusion model: Fick's second law
!          Damage model: phase field damage model
!          by Nguyen Xuan Binh
!          binh.nguyen@aalto.fi
!          July 2024, Abaqus 2023
!          Model extended from the work of Emilio Martinez-Paneda
!          Paper: A phase field formulation for hydrogen assisted cracking
!          Computer Methods in Applied Mechanics and Engineering 342: 742-761 
!          (2018) doi: 10.1016/j.cma.2018.07.021
!          DO NOT DISTRIBUTE WITHOUT AUTHOR'S PERMISSION
!====================================================================

!     State variables  
!     statev(1:6): stress tensor 11, 22, 33, 12, 13, 23        
!     statev(7:12) : strain tensor 11, 22, 33, 12, 13, 23  
!     statev(13) : equivalent plastic strain PEEQ (eqplas)
!     statev(14) : hydrostatic stress (σm)
!     statev(15) : crack phase field (ϕ)
!     statev(16) : history variable field (H)
!     statev(17) : dislocation density rho_d (ρd)
!     statev(18) : total hydrogen concentration (C)
!     statev(19) : hydrogen concentration in lattice sites (CL)
!     statev(20) : hydrogen concentration in trap sites (CT)

!     Debug version
!     statev(1:6): stress tensor 11, 22, 33, 12, 13, 23        
!     statev(7:12) : strain tensor 11, 22, 33, 12, 13, 23  
!     statev(13) : crack phase field (ϕ)
!     statev(14) : hydrostatic stress (σm)
!     statev(15) : hydrogen concentration in lattice sites (CL)