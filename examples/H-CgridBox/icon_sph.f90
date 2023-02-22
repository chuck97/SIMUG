
!©2021 MPI-M, Carolin Mehlmann
                                                                                                                                                                                                                                                

!Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:                                                                                                                              
!1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.                                                                                                                                             
!2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
!3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
!THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSE
! AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.[7]                                               
  


!----------------------------
#include "omp_definitions.inc"
!----------------------------

MODULE mo_ice_test_2D
  !-------------------------------------------------------------------------
  !
  USE mo_kind,                ONLY: wp
  USE mo_parallel_config,     ONLY: nproma
  USE mo_sync,                ONLY: SYNC_C, SYNC_E, SYNC_V, sync_patch_array
  USE mo_model_domain,        ONLY: t_patch, t_patch_3D
  USE mo_grid_subset,         ONLY: t_subset_range, get_index_range
  USE mo_math_types,          ONLY: t_cartesian_coordinates

  
  IMPLICIT NONE

  PUBLIC  :: ice_test_2D

CONTAINS

  !===========Main Program Example 4.2 Sphere ======================  
  SUBROUTINE ice_test_2D( p_patch_3D, p_ice, p_os, p_as, atmos_fluxes, p_op_coeff, p_oce_sfc)

  
    TYPE(t_patch_3D), TARGET, INTENT(IN)     :: p_patch_3D
    TYPE(t_sea_ice),          INTENT(INOUT)  :: p_ice
    TYPE(t_hydro_ocean_state),INTENT(IN)     :: p_os
    TYPE (t_atmos_fluxes),    INTENT(IN)     :: atmos_fluxes
    TYPE(t_operator_coeff),   INTENT(IN)     :: p_op_coeff
    TYPE(t_atmos_for_ocean),  INTENT(INOUT)  :: p_as
    TYPE(t_subset_range),     POINTER        :: all_cells
    TYPE(t_subset_range), POINTER            :: all_edges
    TYPE(t_subset_range), POINTER            :: owned_cells
    TYPE(t_ocean_surface),    INTENT(INOUT)  :: p_oce_sfc
    TYPE(t_subset_range), POINTER :: edges_in_domain


    ! Local variables
    TYPE(t_patch), POINTER :: p_patch
    Real(wp):: boundary_cell_marker(nproma,n_zlev,p_patch_3D%p_patch_2D(1)%alloc_cell_blocks)
    REAL(wp)::  s11(nproma,n_zlev,p_patch_3D%p_patch_2d(1)%alloc_cell_blocks),&
         &s22(nproma,n_zlev,p_patch_3D%p_patch_2d(1)%alloc_cell_blocks),&
         &s21(nproma,n_zlev,p_patch_3D%p_patch_2d(1)%alloc_cell_blocks),&
         &s12(nproma,n_zlev,p_patch_3D%p_patch_2d(1)%alloc_cell_blocks),&
         &s31(nproma,n_zlev,p_patch_3D%p_patch_2d(1)%alloc_cell_blocks),&
         &s32(nproma,n_zlev,p_patch_3D%p_patch_2d(1)%alloc_cell_blocks),&
         &sigma_I(nproma,n_zlev,p_patch_3D%p_patch_2d(1)%alloc_cell_blocks),&
         &sigma_II(nproma,n_zlev,p_patch_3D%p_patch_2d(1)%alloc_cell_blocks),&
         &div_flux_A(nproma,n_zlev,p_patch_3D%p_patch_2d(1)%alloc_cell_blocks),&
         &div_flux_H(nproma,n_zlev,p_patch_3D%p_patch_2d(1)%alloc_cell_blocks),&
         &zeta_c(nproma,n_zlev,p_patch_3D%p_patch_2d(1)%alloc_cell_blocks),&
         &zeta_stabi(nproma,n_zlev,p_patch_3D%p_patch_2d(1)%alloc_cell_blocks),&
         &Delta_c(nproma,n_zlev,p_patch_3D%p_patch_2d(1)%alloc_cell_blocks),&
         &shear_c(nproma,n_zlev,p_patch_3D%p_patch_2d(1)%alloc_cell_blocks),&
         &cell_area_c(nproma,n_zlev,p_patch_3D%p_patch_2d(1)%alloc_cell_blocks),&
         &x1_c(nproma,n_zlev,p_patch_3D%p_patch_2d(1)%alloc_cell_blocks),&
         &x2_c(nproma,n_zlev,p_patch_3D%p_patch_2d(1)%alloc_cell_blocks),&
         &x3_c(nproma,n_zlev,p_patch_3D%p_patch_2d(1)%alloc_cell_blocks)
         
    

   REAL(wp),POINTER             ::  boundary_edge_marker(:,:) 
   REAL(wp),POINTER             ::  mass(:,:), atm_n(:,:), atm_t(:,:), u_old_n(:,:), u_old_t(:,:), u_ocean_n(:,:),&
        & u_ocean_t(:,:), u_ocean_x(:,:), u_ocean_y(:,:), Au_n(:,:), Au_t(:,:), h_e(:,:), grad_h(:,:),&
        &S_x(:,:),S_y(:,:), S_n(:,:),S_t(:,:), ice_x(:,:), ice_y(:,:), n1_s(:,:),n2_s(:,:),n3_s(:,:), ice_x_old(:,:), ice_y_old(:,:),&
        &Rn(:,:),Rt(:,:),Dn(:,:),Dt(:,:),Zn(:,:),Zt(:,:), zeta_e(:,:),wind_relativ(:,:),flux_A(:,:,:),flux_H(:,:,:),un(:,:,:)

    TYPE(t_cartesian_coordinates) :: p_tau_n_c(nproma,p_patch_3D%p_patch_2D(1)%alloc_cell_blocks)
    TYPE(t_cartesian_coordinates) :: ocean_c(nproma,p_patch_3D%p_patch_2D(1)%alloc_cell_blocks)

   INTEGER ::time_iter,outer_iter, time, edge_block_i,start_index,end_index,edge_index_i, cell_block, cell_index, neigbor,&
        &edge_block_1,edge_block_2,edge_block_3,edge_index_1,edge_index_2,edge_index_3, inner_iter&
       &,cell_index1,cell_index2,cell_block1,cell_block2,doy
   INTEGER vtkeverystep
   
       

   REAL(wp) :: a,x,y,z,x_0,y_0, nx,ny,tx,ty,Cda,Cdw, Cor, delta_min, Zeitschritt, alpha_evp, beta_evp,nix,niy,niz,tix,tiy,tiz&
        &,T,e,Delta, nx_cc,ny_cc, tx_cc,ty_cc, e1,e2,e3,s,R,nz,tz, nix1,nix2,nix3,niy2,niy3,niy1,O1,O2,O3,x123,grad_x,grad_y
   REAL(wp) ::  Oi, Oj, ei, e11, e22, e21, e12, e31,e32,P, Pstar, zeta,eta, x1,x2,x3, y1,y2,y3, tix1,tix2,tix3,tiy2,tiy3,tiy1,ux,uy,uz,&
        &sc_pn, l_sn, nix_1,niy_1,niz_1, tix_1,tiy_1,tiz_1,nix_l, niy_l, tix_l,tiy_l,U,V,W, U_ana, V_ana,W_ana, ice_z, u_rel,v_rel,w_rel,ocean_z,&
        &atop_n,atop_t,abot_n,abot_t,btop_n,btop_t,bbot_n,bbot_t,aaa,bbb,x_shift,y_shift,L,diff_n,diff_t,vdw,tp,U_xx,U_yy,atm_u,atm_v,mass_ice
   REAL(wp) :: alpha,vmax, ws,wx,wy,mx,my,rw,L_ref,dx,nn,nw,no,uo   
       
   !--------------------------------------------------------------------------------------------------
   
    p_patch => p_patch_3D%p_patch_2D(1)
    all_cells     => p_patch_3d%p_patch_2d(1)%cells%all
    all_edges =>     p_patch_3d%p_patch_2d(1)%edges%all
    owned_cells => p_patch%cells%owned
    edges_in_domain   => p_patch%edges%in_domain
     
    !===========initialize edge vector====================== 
    ALLOCATE(boundary_edge_marker(nproma,p_patch_3d%p_patch_2d(1)%nblks_e))
    ALLOCATE(mass(nproma,p_patch_3d%p_patch_2d(1)%nblks_e))
    ALLOCATE(atm_n(nproma,p_patch_3d%p_patch_2d(1)%nblks_e))
    ALLOCATE(atm_t(nproma,p_patch_3d%p_patch_2d(1)%nblks_e))
    ALLOCATE(u_old_n(nproma,p_patch_3d%p_patch_2d(1)%nblks_e))
    ALLOCATE(u_old_t(nproma,p_patch_3d%p_patch_2d(1)%nblks_e))
    ALLOCATE(wind_relativ(nproma,p_patch_3d%p_patch_2d(1)%nblks_e))
    ALLOCATE(u_ocean_n(nproma,p_patch_3d%p_patch_2d(1)%nblks_e))
    ALLOCATE(u_ocean_t(nproma,p_patch_3d%p_patch_2d(1)%nblks_e))
    ALLOCATE(u_ocean_x(nproma,p_patch_3d%p_patch_2d(1)%nblks_e))
    ALLOCATE(u_ocean_y(nproma,p_patch_3d%p_patch_2d(1)%nblks_e))
    ALLOCATE(ice_x(nproma,p_patch_3d%p_patch_2d(1)%nblks_e))
    ALLOCATE(ice_y(nproma,p_patch_3d%p_patch_2d(1)%nblks_e))
    ALLOCATE(ice_x_old(nproma,p_patch_3d%p_patch_2d(1)%nblks_e))
    ALLOCATE(ice_y_old(nproma,p_patch_3d%p_patch_2d(1)%nblks_e))
    ALLOCATE(Au_n(nproma,p_patch_3d%p_patch_2d(1)%nblks_e))
    ALLOCATE(Au_t(nproma,p_patch_3d%p_patch_2d(1)%nblks_e))
    ALLOCATE(h_e(nproma,p_patch_3d%p_patch_2d(1)%nblks_e))
    ALLOCATE(grad_h(nproma,p_patch_3d%p_patch_2d(1)%nblks_e))
    ALLOCATE(S_x(nproma,p_patch_3d%p_patch_2d(1)%nblks_e))
    ALLOCATE(S_y(nproma,p_patch_3d%p_patch_2d(1)%nblks_e))
    ALLOCATE(S_t(nproma,p_patch_3d%p_patch_2d(1)%nblks_e))
    ALLOCATE(S_n(nproma,p_patch_3d%p_patch_2d(1)%nblks_e))
    ALLOCATE(Rn(nproma,p_patch_3d%p_patch_2d(1)%nblks_e))
    ALLOCATE(Rt(nproma,p_patch_3d%p_patch_2d(1)%nblks_e))
    ALLOCATE(Dn(nproma,p_patch_3d%p_patch_2d(1)%nblks_e))
    ALLOCATE(Dt(nproma,p_patch_3d%p_patch_2d(1)%nblks_e))
    ALLOCATE(Zn(nproma,p_patch_3d%p_patch_2d(1)%nblks_e))
    ALLOCATE(Zt(nproma,p_patch_3d%p_patch_2d(1)%nblks_e))
    ALLOCATE(zeta_e(nproma,p_patch_3d%p_patch_2d(1)%nblks_e))
    ALLOCATE(flux_A(nproma,n_zlev,p_patch_3d%p_patch_2d(1)%nblks_e))
    ALLOCATE(flux_H(nproma,n_zlev,p_patch_3d%p_patch_2d(1)%nblks_e))
    ALLOCATE(un(nproma,n_zlev,p_patch_3d%p_patch_2d(1)%nblks_e))
    
    !===========set parameter======================
    p_ice%ice_iter=p_ice%ice_iter+1.0! Itertion count
    write(0,*) 'iter', p_ice%ice_iter
    mass_ice=900.0_wp ! ice density
    L_ref=0.2! resize the domain grid : AquaPlanet_Icos_0079km_mirror.nc
    dx=0.4
    L=dx*2
    Cdw=1026_wp*.0055_wp*L_ref!ocean_drag
    Cda=1.3_wp*.0012_wp*L_ref!athomsphere_drag
    Cor=0.0! Coriolis parameter
    Pstar=27500/(L_ref*L_ref)! ice strength parameter
    doy=0! count no. of cells
    Zeitschritt=1800!time step in seconds
    vtkeverystep = int(60*60*8/Zeitschritt)!output frequency 
    delta_min=0.000000002_wp! Hibler Delta
    
    
    !===========initialize cell vectors======================   
    DO cell_block = owned_cells%start_block, owned_cells%end_block
       CALL get_index_range(owned_cells, cell_block, start_index, end_index)
       DO cell_index = start_index, end_index
          s11(cell_index,1,cell_block)=0.0_wp
          s12(cell_index,1,cell_block)=0.0_wp
          s22(cell_index,1,cell_block)=0.0_wp
          s21(cell_index,1,cell_block)=0.0_wp
          s32(cell_index,1,cell_block)=0.0_wp
          s31(cell_index,1,cell_block)=0.0_wp
          sigma_I(cell_index,1,cell_block)=0.0_wp
          sigma_II(cell_index,1,cell_block)=0.0_wp
          zeta_c(cell_index,1,cell_block)=0.0_wp
          zeta_stabi(cell_index,1,cell_block)=0.0_wp
          div_flux_A(cell_index,1,cell_block)=0.0_wp
          p_ice%Delta(cell_index,cell_block)=0.0
          cell_area_c(cell_index,1,cell_block)=0.0
          x1_c(cell_index,1,cell_block)=0.0
          x2_c(cell_index,1,cell_block)=0.0
          x3_c(cell_index,1,cell_block)=0.0
          shear_c(cell_index,1,cell_block)=0.0
       ENDDO ! cell_index = start_index, end_index
    ENDDO !cell_block = owned_cells%start_block, owned_cells%end_block  

    CAll  interface_boundary_cell_marker(boundary_cell_marker, p_patch_3D, p_ice)
    Call interface_boundary_edge_marker(boundary_edge_marker,boundary_cell_marker, p_patch_3D, p_ice)
    Call cell_area(cell_area_c, p_patch_3D)
    Call  init_mass_matrix(mass,cell_area_c,p_patch_3D)
    p_ice%hi=0.0_wp
    p_ice%conc=0.0_wp

 !=========== initialize concentration and thickness====================== 
    DO cell_block = owned_cells%start_block, owned_cells%end_block
       CALL get_index_range(owned_cells, cell_block, start_index, end_index)
       DO cell_index = start_index, end_index
            if (boundary_cell_marker(cell_index,1,cell_block)>0.0) then
             doy=doy+1
             x= p_patch%cells%cartesian_center(cell_index,cell_block)%x(1)
             y= p_patch%cells%cartesian_center(cell_index,cell_block)%x(2)
             x_shift=(x+dx)/(2*L)
             y_shift=(y+dx)/(2*L)
             p_ice%hi(cell_index,1,cell_block)=2.0
             p_ice%conc(cell_index,1,cell_block)= (x+dx)/L
          endif
         ENDDO ! cell_index = start_index, end_index
    ENDDO !cell_block = owned_cells%start_block, owned_cells%end_block   
    write(0,*) doy/2, "number of cells"

 !=========== initialize concentration and thickness======================  
       DO cell_block = owned_cells%start_block, owned_cells%end_block
          CALL get_index_range(owned_cells, cell_block, start_index, end_index)
          DO cell_index = start_index, end_index
             DO neigbor=1,3 !no_primal_edges                                                                                                                                                                           
                edge_index_i =p_patch%cells%edge_idx(cell_index, cell_block, neigbor)
                edge_block_i = p_patch%cells%edge_blk(cell_index, cell_block, neigbor)
                !Caretesian coordinates at edge midpoiny
                x = p_patch%edges%cartesian_center(edge_index_i,edge_block_i)%x(1)
                y = p_patch%edges%cartesian_center(edge_index_i,edge_block_i)%x(2)
                z = p_patch%edges%cartesian_center(edge_index_i,edge_block_i)%x(3)
                !cartesian coordinates cell center
                x1_c(cell_index,1,cell_block)=x1_c(cell_index,1,cell_block) + 1.0_wp/3.0_wp * x
                x2_c(cell_index,1,cell_block)=x2_c(cell_index,1,cell_block) + 1.0_wp/3.0_wp * y
                x3_c(cell_index,1,cell_block)=x3_c(cell_index,1,cell_block) + 1.0_wp/3.0_wp * z
                !orientation of the normal vector
                Oi=p_patch%cells%edge_orientation(cell_index,cell_block,neigbor)
                !cartesian normal vector
                nix=p_patch%edges%primal_cart_normal(edge_index_i,edge_block_i)%x(1)
                niy=p_patch%edges%primal_cart_normal(edge_index_i,edge_block_i)%x(2)
                niz=p_patch%edges%primal_cart_normal(edge_index_i,edge_block_i)%x(3)
                !cartesian tangential vector
                tix=p_patch%edges%dual_cart_normal(edge_index_i,edge_block_i)%x(1)
                tiy=p_patch%edges%dual_cart_normal(edge_index_i,edge_block_i)%x(2)
                tiz=p_patch%edges%dual_cart_normal(edge_index_i,edge_block_i)%x(3)
                !wind forcing
                x_shift=x+dx
                y_shift=y+dx
                L=dx*2
                U_yy=(5.0_wp+(sin(2.0_wp*3.141*tp/4.0)-3)*sin(2.0_wp*3.141*y_shift/L)*sin(3.141*x_shift/L))/L_ref             
                U_xx=(5.0_wp+(sin(2.0_wp*3.141*tp/4.0)-3)*sin(2.0_wp*3.141*x_shift/L)*sin(3.141*y_shift/L))/L_ref
                atm_u=U_xx*dsqrt(U_xx*U_xx+U_yy*U_yy)*Cda
                atm_v=U_yy*dsqrt(U_xx*U_xx+U_yy*U_yy)*Cda
                atm_n(edge_index_i,edge_block_i)=atm_u*nix+atm_v*niy
                atm_t(edge_index_i,edge_block_i)=atm_u*tix+atm_v*tiy
                nn=sqrt(atm_n(edge_index_i,edge_block_i)*atm_n(edge_index_i,edge_block_i)+atm_t(edge_index_i,edge_block_i)*atm_t(edge_index_i,edge_block_i))
                nw=sqrt(atm_u*atm_u+atm_v*atm_v)
                atm_n(edge_index_i,edge_block_i)=atm_n(edge_index_i,edge_block_i)!/nn*nw
                atm_t(edge_index_i,edge_block_i)=atm_t(edge_index_i,edge_block_i)!/nn*nw
                !ocean forcing
                u_ocean_x(edge_index_i,edge_block_i)=0.1/L_ref*(2.0*y_shift-L)/(L)*boundary_edge_marker(edge_index_i,edge_block_i)                          
                u_ocean_y(edge_index_i,edge_block_i)=0.1/L_ref*(L-2.0*x_shift)/L*boundary_edge_marker(edge_index_i,edge_block_i)
                uo=sqrt(u_ocean_x(edge_index_i,edge_block_i)* u_ocean_x(edge_index_i,edge_block_i)+u_ocean_y(edge_index_i,edge_block_i)*u_ocean_y(edge_index_i,edge_block_i))
                u_ocean_n(edge_index_i,edge_block_i)=(u_ocean_x(edge_index_i,edge_block_i)*nix+&
                     &u_ocean_y(edge_index_i,edge_block_i)*niy)*boundary_edge_marker(edge_index_i,edge_block_i)
                u_ocean_t(edge_index_i,edge_block_i)=(u_ocean_x(edge_index_i,edge_block_i)*tix+&
                     & u_ocean_y(edge_index_i,edge_block_i)*tiy)*boundary_edge_marker(edge_index_i,edge_block_i)
                 no=sqrt(u_ocean_n(edge_index_i,edge_block_i)* u_ocean_n(edge_index_i,edge_block_i)+u_ocean_t(edge_index_i,edge_block_i)*u_ocean_t(edge_index_i,edge_block_i))                 
                !Initial sea ice velocity
                p_ice%vn_e(edge_index_i,edge_block_i)=0.0!
                p_ice%vt_e(edge_index_i,edge_block_i)=0.0!
                ! Initial value  CG
                ice_x(edge_index_i,edge_block_i)=0.0_wp
                ice_y(edge_index_i,edge_block_i)=0.0_wp
             ENDDO

          !cell center in cartesian coordinates
          x123 = sqrt(x1_c(cell_index,1,cell_block)*x1_c(cell_index,1,cell_block) + &
               &x2_c(cell_index,1,cell_block)*x2_c(cell_index,1,cell_block) + &
               &x3_c(cell_index,1,cell_block)*x3_c(cell_index,1,cell_block))
          !normalized to cell center
          x1_c(cell_index,1,cell_block)=x1_c(cell_index,1,cell_block)/x123
          x2_c(cell_index,1,cell_block)=x2_c(cell_index,1,cell_block)/x123
          x3_c(cell_index,1,cell_block)=x3_c(cell_index,1,cell_block)/x123
       ENDDO
    ENDDO
    
    !Initialize stress tensor at edges  
    Au_n=0.0_wp
    Au_t=0.0_wp
    !!Initialize stabilization at edges
    S_n=0.0_wp
    S_t=0.0_wp

    !=========== time loop ======================         
    DO time_iter=1,48 ! time steps
       tp=(Zeitschritt*time_iter)/(3600.0_wp*24.0_wp)! time in days
       write(0,*) "Zeit/day", tp, time_iter
       !sea ice vloxity from the old time step
        DO edge_block_i = all_edges%start_block, all_edges%end_block
          CALL get_index_range(all_edges, edge_block_i, start_index, end_index)
          DO edge_index_i =  start_index, end_index
             ice_x_old(edge_index_i,edge_block_i)=ice_x(edge_index_i,edge_block_i)*boundary_edge_marker(edge_index_i,edge_block_i)
             ice_y_old(edge_index_i,edge_block_i)=ice_y(edge_index_i,edge_block_i)*boundary_edge_marker(edge_index_i,edge_block_i)
          ENDDO
       ENDDO
       !subiteration per time step Picard iteration
       DO outer_iter=1,25       
          DO edge_block_i = all_edges%start_block, all_edges%end_block
             CALL get_index_range(all_edges, edge_block_i, start_index, end_index)
             DO edge_index_i =  start_index, end_index      
                !Update sea ice velocity
                p_ice%vn_e(edge_index_i,edge_block_i)=ice_x(edge_index_i,edge_block_i)*boundary_edge_marker(edge_index_i,edge_block_i)
                p_ice%vt_e(edge_index_i,edge_block_i)=ice_y(edge_index_i,edge_block_i)*boundary_edge_marker(edge_index_i,edge_block_i)
             ENDDO
          ENDDO
       
          DO cell_block = owned_cells%start_block, owned_cells%end_block
             CALL get_index_range(owned_cells, cell_block, start_index, end_index)
             DO cell_index = start_index, end_index
                x = p_patch%cells%cartesian_center(cell_index,cell_block)%x(1)
                y = p_patch%cells%cartesian_center(cell_index,cell_block)%x(2)
                z = p_patch%cells%cartesian_center(cell_index,cell_block)%x(3)
                !Update strain rate tensor and VP rheology
                if (boundary_cell_marker(cell_index,1,cell_block)>0.0) then
                   !Calculates strain rate tensor
                   Call compute_2Dvector_grad(cell_index,cell_block,e11,e12,e21,e22,cell_area_c(cell_index,1,cell_block),&
                        &x1_c(cell_index,1,cell_block),x2_c(cell_index,1,cell_block),x3_c(cell_index,1,cell_block),p_patch_3D, p_ice)
                   Delta=dsqrt(delta_min*delta_min+ e12*e12 + 1.25_wp*(e11*e11+e22*e22)+1.5_wp*e11*e22)!delta Hibler
                   Delta_c(cell_index,1,cell_block)=Delta! save Delta
                   shear_c(cell_index,1,cell_block)=dsqrt(Delta_min*Delta_min+(e11-e22)*(e11-e22)+(e12+e21)*(e12+e21))!shear stress
                   P=p_ice%hi(cell_index,1,cell_block)*Pstar*exp(-20.0_wp*(1.0_wp-p_ice%conc(cell_index,1,cell_block)))!ice strength
                   zeta_c(cell_index,1,cell_block)=(p_ice%hi(cell_index,1,cell_block)*P/(2*Delta))*boundary_cell_marker(cell_index,1,cell_block)!viscosity
                   !Update stresses VP rheology Hibler
                   Call VP( e11,e12,e21,e22,s11(cell_index,1,cell_block),s12(cell_index,1,cell_block),&
                        &s21(cell_index,1,cell_block),s22(cell_index,1,cell_block),zeta_c(cell_index,1,cell_block),P)
                endif
             ENDDO ! cell_index = start_index, end_index                                                                                                                                                             
          ENDDO !cell_block = owned_cells%start_block, owned_cells%end_block  
       
          DO edge_block_i = all_edges%start_block, all_edges%end_block
             CALL get_index_range(all_edges, edge_block_i, start_index, end_index)
             DO edge_index_i =  start_index, end_index
                cell_index1 = p_patch%edges%cell_idx(edge_index_i,edge_block_i,1)
                cell_block1 = p_patch%edges%cell_blk(edge_index_i,edge_block_i,1)
                cell_index2 = p_patch%edges%cell_idx(edge_index_i,edge_block_i,2)
                cell_block2 = p_patch%edges%cell_blk(edge_index_i,edge_block_i,2)
                !average zeta for stabilization and thickness for momentum equation
                zeta_e(edge_index_i,edge_block_i)=0.5_wp*boundary_edge_marker(edge_index_i,edge_block_i)*(zeta_c(cell_index1,1,cell_block1)+zeta_c(cell_index2,1,cell_block2))
                h_e(edge_index_i,edge_block_i)=0.5_wp*boundary_edge_marker(edge_index_i,edge_block_i)*(p_ice%hi(cell_index1,1,cell_block1)+p_ice%hi(cell_index2,1,cell_block2))
             ENDDO
          ENDDO
          Au_n=0.0_wp
          Au_t=0.0_wp
          !compute div(sigma) at edges 
          CALL compute_2D_vector_Laplace(Au_n,Au_t,boundary_cell_marker,s11,s12,s21,s22,s31,s32,cell_area_c,x1_c,x2_c,x3_c,p_patch_3D, p_ice)
          S_x=0.0_wp
          S_y=0.0_wp
          Call  Stabilization_sum(S_x,S_y,boundary_cell_marker,x1_c,x2_c,x3_c,p_ice,p_patch_3D)
          S_n=0.0_wp
          S_t=0.0_wp
          Call  Stabilization(S_x,S_y,S_n,S_t,boundary_cell_marker,x1_c,x2_c,x3_c,zeta_e,p_patch_3D)
      
          DO edge_block_i = all_edges%start_block, all_edges%end_block
             CALL get_index_range(all_edges, edge_block_i, start_index, end_index)
             DO edge_index_i =  start_index, end_index
                Rn(edge_index_i,edge_block_i)=0.0_wp
                Rt(edge_index_i,edge_block_i)=0.0_wp
                if(boundary_edge_marker(edge_index_i,edge_block_i)>0.0_wp) then
                   x = p_patch%edges%cartesian_center(edge_index_i,edge_block_i)%x(1)
                   y = p_patch%edges%cartesian_center(edge_index_i,edge_block_i)%x(2)
                   z = p_patch%edges%cartesian_center(edge_index_i,edge_block_i)%x(3)
                   nix=p_patch%edges%primal_cart_normal(edge_index_i,edge_block_i)%x(1)
                   niy=p_patch%edges%primal_cart_normal(edge_index_i,edge_block_i)%x(2)
                   niz=p_patch%edges%primal_cart_normal(edge_index_i,edge_block_i)%x(3)
                   tix=p_patch%edges%dual_cart_normal(edge_index_i,edge_block_i)%x(1)
                   tiy=p_patch%edges%dual_cart_normal(edge_index_i,edge_block_i)%x(2)
                   tiz=p_patch%edges%dual_cart_normal(edge_index_i,edge_block_i)%x(3)
                   x_shift=x+dx
                   y_shift=y+dx
                   L=dx*2
                   !vwd ocean stress
                   diff_n= u_ocean_n(edge_index_i,edge_block_i)-ice_x(edge_index_i,edge_block_i)
                   diff_t= u_ocean_t(edge_index_i,edge_block_i)-ice_y(edge_index_i,edge_block_i)
                   wind_relativ(edge_index_i,edge_block_i)= sqrt(diff_n*diff_n+diff_t*diff_t)
                   vdw=wind_relativ(edge_index_i,edge_block_i)*boundary_edge_marker(edge_index_i,edge_block_i)
                   U_yy=(5.0_wp+(sin(2.0_wp*3.141*tp/4.0)-3)*sin(2.0_wp*3.141*y_shift/L)*sin(3.141*x_shift/L))/L_ref
                   U_xx=(5.0_wp+(sin(2.0_wp*3.141*tp/4.0)-3)*sin(2.0_wp*3.141*x_shift/L)*sin(3.141*y_shift/L))/L_ref
                   atm_u=U_xx*dsqrt(U_xx*U_xx+U_yy*U_yy)*Cda
                   atm_v=U_yy*dsqrt(U_xx*U_xx+U_yy*U_yy)*Cda
                   atm_n(edge_index_i,edge_block_i)=(atm_u*nix+atm_v*niy)*boundary_edge_marker(edge_index_i,edge_block_i)
                   atm_t(edge_index_i,edge_block_i)=(atm_u*tix+atm_v*tiy)*boundary_edge_marker(edge_index_i,edge_block_i)
                   nn=sqrt(atm_n(edge_index_i,edge_block_i)*atm_n(edge_index_i,edge_block_i)+atm_t(edge_index_i,edge_block_i)*atm_t(edge_index_i,edge_block_i))
                   nw=sqrt(atm_u*atm_u+atm_v*atm_v)

                   mass_ice=900*h_e(edge_index_i,edge_block_i)
                   u_old_n(edge_index_i,edge_block_i)=&
                        &ice_x(edge_index_i,edge_block_i)&
                        &+Cdw*vdw*Zeitschritt/mass_ice*ice_x(edge_index_i,edge_block_i)&
                        &+Zeitschritt/mass_ice*(Au_n(edge_index_i,edge_block_i)/mass(edge_index_i,edge_block_i)&!*zeta_e(edge_index_i,edge_block_i)&
                        &-S_n(edge_index_i,edge_block_i)/mass(edge_index_i,edge_block_i))&
                        &-Cor/Zeitschritt*ice_y(edge_index_i,edge_block_i)
                   u_old_t(edge_index_i,edge_block_i)=&
                        &ice_y(edge_index_i,edge_block_i)&
                        &+Cdw*vdw*Zeitschritt/mass_ice*ice_y(edge_index_i,edge_block_i)&
                        &+Zeitschritt/mass_ice*(Au_t(edge_index_i,edge_block_i)/mass(edge_index_i,edge_block_i)&!*zeta_e(edge_index_i,edge_block_i)&
                        &-S_t(edge_index_i,edge_block_i)/mass(edge_index_i,edge_block_i))&
                        &+Cor/Zeitschritt*ice_x(edge_index_i,edge_block_i)
             
                ! Righthand side and initial value CG                                                                                               
                   Rn(edge_index_i,edge_block_i)=(Zeitschritt/mass_ice*Cdw*vdw*u_ocean_n(edge_index_i,edge_block_i)+Zeitschritt/mass_ice*atm_n(edge_index_i,edge_block_i)&
                        &+ice_x_old(edge_index_i,edge_block_i)-&
                        &u_old_n(edge_index_i,edge_block_i))&
                        &-Cor/Zeitschritt*u_ocean_t(edge_index_i,edge_block_i)
                   Rt(edge_index_i,edge_block_i)=(Zeitschritt/mass_ice*Cdw*vdw*u_ocean_t(edge_index_i,edge_block_i)+Zeitschritt/mass_ice*atm_t(edge_index_i,edge_block_i)&
                        &+ice_y_old(edge_index_i,edge_block_i)-&
                        &u_old_t(edge_index_i,edge_block_i))&
                        &+Cor/Zeitschritt*u_ocean_n(edge_index_i,edge_block_i)
                   if ( h_e(edge_index_i,edge_block_i)==0.0_wp) then
                      write(*,*) 'assert message h_e',h_e(edge_index_i,edge_block_i)
                      call exit(1)
                   endif
                endif
                Dn(edge_index_i,edge_block_i)=Rn(edge_index_i,edge_block_i)
                Dt(edge_index_i,edge_block_i)=Rt(edge_index_i,edge_block_i)
             ENDDO
          ENDDO
          write(0,*) "outer_iter" ,maxval(ice_x)*L_ref  !output x-coordinate velocity
          
          !CG method
          DO inner_iter=1,400
             ! CG: Prepare for z = A d     
             DO edge_block_i = all_edges%start_block, all_edges%end_block
                CALL get_index_range(all_edges, edge_block_i, start_index, end_index)
                DO edge_index_i =  start_index, end_index
                   p_ice%vn_e(edge_index_i,edge_block_i)=Dn(edge_index_i,edge_block_i)*boundary_edge_marker(edge_index_i,edge_block_i)
                   p_ice%vt_e(edge_index_i,edge_block_i)=Dt(edge_index_i,edge_block_i)*boundary_edge_marker(edge_index_i,edge_block_i)
                ENDDO
             ENDDO
             DO cell_block = owned_cells%start_block, owned_cells%end_block
                CALL get_index_range(owned_cells, cell_block, start_index, end_index)
                DO cell_index = start_index, end_index
                   x = p_patch%cells%cartesian_center(cell_index,cell_block)%x(1)
                   y = p_patch%cells%cartesian_center(cell_index,cell_block)%x(2)
                   z = p_patch%cells%cartesian_center(cell_index,cell_block)%x(3)
                   if (boundary_cell_marker(cell_index,1,cell_block)>0.0) then
                      Call compute_2Dvector_grad(cell_index,cell_block,e11,e12,e21,e22,cell_area_c(cell_index,1,cell_block),&
                           &x1_c(cell_index,1,cell_block),x2_c(cell_index,1,cell_block),x3_c(cell_index,1,cell_block),p_patch_3D, p_ice) 
                      P=p_ice%hi(cell_index,1,cell_block)*Pstar*exp(-20.0_wp*(1.0_wp-p_ice%conc(cell_index,1,cell_block)))
                      Call  VP( e11,e12,e21,e22,s11(cell_index,1,cell_block),s12(cell_index,1,cell_block),&
                           &s21(cell_index,1,cell_block),s22(cell_index,1,cell_block),zeta_c(cell_index,1,cell_block),P)
                   endif
                ENDDO ! cell_index = start_index, end_index
             ENDDO !cell_block = owned_cells%start_block, owned_cells%end_block
             Au_n=0.0_wp
             Au_t=0.0_wp
             CALL compute_2D_vector_Laplace(Au_n,Au_t,boundary_cell_marker,s11,s12,s21,s22,s31,s32,cell_area_c,x1_c,x2_c,x3_c,p_patch_3D, p_ice)
             S_x=0.0_wp
             S_y=0.0_wp
             Call  Stabilization_sum(S_x,S_y,boundary_cell_marker,x1_c,x2_c,x3_c,p_ice,p_patch_3D)
             S_n=0.0_wp
             S_t=0.0_wp
             Call  Stabilization(S_x,S_y,S_n,S_t,boundary_cell_marker,x1_c,x2_c,x3_c,zeta_e,p_patch_3D)
             ! CG set z of:  z = A d
             Zn = 0.0_wp
             Zt = 0.0_wp
             DO edge_block_i = all_edges%start_block, all_edges%end_block
                CALL get_index_range(all_edges, edge_block_i, start_index, end_index)
                DO edge_index_i =  start_index, end_index
                   vdw=wind_relativ(edge_index_i,edge_block_i)
                   Zn(edge_index_i,edge_block_i)=0.0_wp
                   Zt(edge_index_i,edge_block_i)=0.0_wp
                   if(boundary_edge_marker(edge_index_i,edge_block_i)>0.0_wp) then
                      mass_ice=900*h_e(edge_index_i,edge_block_i)
                      Zn(edge_index_i,edge_block_i)=&
                           & Dn(edge_index_i,edge_block_i)&
                           &+Zeitschritt/mass_ice*Cdw*vdw*Dn(edge_index_i,edge_block_i)&
                           &+Zeitschritt/mass_ice*(Au_n(edge_index_i,edge_block_i)/mass(edge_index_i,edge_block_i)&!*zeta_e(edge_index_i,edge_block_i)&
                           &-S_n(edge_index_i,edge_block_i)/mass(edge_index_i,edge_block_i))&
                           &-Cor/Zeitschritt*Dt(edge_index_i,edge_block_i)
                      Zt(edge_index_i,edge_block_i)=&
                           & Dt(edge_index_i,edge_block_i)&
                           &+Zeitschritt/mass_ice*Cdw*vdw*Dt(edge_index_i,edge_block_i)&
                           &+Zeitschritt/mass_ice*(Au_t(edge_index_i,edge_block_i)/mass(edge_index_i,edge_block_i)&!*zeta_e(edge_index_i,edge_block_i)&
                           &-S_t(edge_index_i,edge_block_i)/mass(edge_index_i,edge_block_i))&
                           &+Cor/Zeitschritt*Dn(edge_index_i,edge_block_i)
                   endif
                ENDDO
             ENDDO
             ! alpha = (r,r)/(d,z)
             Call       CG_scalarp(atop_n,atop_t,Rn,Rn,Rt,Rt,p_patch_3D)
             Call       CG_scalarp(abot_n,abot_t,Dn,Zn,Dt,Zt,p_patch_3D)
             aaa = (atop_n+atop_t)/(abot_n+abot_t)
             ! x = x + aaa * D
             Call       CG_addvector(aaa,ice_x,Dn,p_patch_3D)
             Call      CG_addvector(aaa,ice_y,Dt,p_patch_3D)
             ! R = R - aaa * Z
             Call      CG_addvector(-aaa,Rn,Zn,p_patch_3D)
             Call      CG_addvector(-aaa,Rt,Zt,p_patch_3D)
             ! beta = (r,r)/atop
             bbot_n = atop_n
             bbot_t = atop_t
             Call   CG_scalarp(btop_n,btop_t,Rn,Rn,Rt,Rt,p_patch_3D)
             bbb = (btop_n+btop_t)/(atop_n+atop_t)
             !output error
             write (0,*) "Fehler ", (btop_n+btop_t), inner_iter, outer_iter, time_iter
             ! d = r + bbb * d
             Call  CG_addvector2(bbb,Rn,Dn,p_patch_3D)
             Call  CG_addvector2(bbb,Rt,Dt,p_patch_3D)
          ENDDO !inner iter
       ENDDO !outer iter
      
       DO edge_block_i = all_edges%start_block, all_edges%end_block
          CALL get_index_range(all_edges, edge_block_i, start_index, end_index)
          DO edge_index_i =  start_index, end_index
             p_ice%vn_e(edge_index_i,edge_block_i)=ice_x(edge_index_i,edge_block_i)*boundary_edge_marker(edge_index_i,edge_block_i)
          ENDDO
       ENDDO
    ENDDO! time iter

    !Visualization
    DO edge_block_i = all_edges%start_block, all_edges%end_block
       CALL get_index_range(all_edges, edge_block_i, start_index, end_index)
       DO edge_index_i =  start_index, end_index
          p_ice%vn_e(edge_index_i,edge_block_i)=ice_x(edge_index_i,edge_block_i)
          p_ice%vt_e(edge_index_i,edge_block_i)=ice_y(edge_index_i,edge_block_i)
       ENDDO
    ENDDO
    
      ice_x=0.0_wp
      ice_y=0.0_wp
      
      DO cell_block = owned_cells%start_block, owned_cells%end_block
         CALL get_index_range(owned_cells, cell_block, start_index, end_index)
         DO cell_index = start_index, end_index
            DO neigbor=1,3 !no_primal_edges
               edge_index_i = p_patch%cells%edge_idx(cell_index, cell_block, neigbor)
               edge_block_i = p_patch%cells%edge_blk(cell_index, cell_block, neigbor)
               x = p_patch%edges%cartesian_center(edge_index_i,edge_block_i)%x(1)
               y = p_patch%edges%cartesian_center(edge_index_i,edge_block_i)%x(2)
               z = p_patch%edges%cartesian_center(edge_index_i,edge_block_i)%x(3)
               nix=p_patch%edges%primal_cart_normal(edge_index_i,edge_block_i)%x(1)
               niy=p_patch%edges%primal_cart_normal(edge_index_i,edge_block_i)%x(2)
               niz=p_patch%edges%primal_cart_normal(edge_index_i,edge_block_i)%x(3)
               tix=p_patch%edges%dual_cart_normal(edge_index_i,edge_block_i)%x(1)
               tiy=p_patch%edges%dual_cart_normal(edge_index_i,edge_block_i)%x(2)
               tiz=p_patch%edges%dual_cart_normal(edge_index_i,edge_block_i)%x(3)
                ice_x(edge_index_i,edge_block_i)=p_ice%vn_e(edge_index_i,edge_block_i)*nix+&
                     & p_ice%vt_e(edge_index_i,edge_block_i)*tix
                ice_y(edge_index_i,edge_block_i)=p_ice%vn_e(edge_index_i,edge_block_i)*niy+&
                     & p_ice%vt_e(edge_index_i,edge_block_i)*tiy
             ENDDO
          ENDDO
       ENDDO
       
       DO cell_block = owned_cells%start_block, owned_cells%end_block
          CALL get_index_range(owned_cells, cell_block, start_index, end_index)
          DO cell_index = start_index, end_index
             p_ice%Delta(cell_index,cell_block) =Delta_c(cell_index,1,cell_block)*boundary_cell_marker(cell_index,1,cell_block)
          ENDDO ! cell_index = start_index, end_index
       ENDDO
       CALL visu_vtk('ice_drift_sphere',p_ice%ice_iter,ice_x,ice_y,boundary_cell_marker,p_ice,p_patch_3D)      
  END SUBROUTINE ice_test_2D

     SUBROUTINE VP( e11,e12,e21,e22,s11,s12,s21,s22,zeta,P)
      REAL(wp), INTENT(in) :: e11
      REAL(wp), INTENT(in) :: e12
      REAL(wp), INTENT(in) :: e21
      REAL(wp), INTENT(in) :: e22

      REAL(wp), INTENT(out) :: s12
      REAL(wp), INTENT(out) :: s21
      REAL(wp), INTENT(out) :: s11
      REAL(wp), INTENT(out) :: s22
      REAL(wp), INTENT(in) :: zeta
      REAL(wp), INTENT(in) :: P
      REAL :: eta,Delta
      eta=0.25_wp*zeta
      s11=zeta*e11+(e11+e22)*(zeta-eta)-0.5_wp*P                                                                                                                               
      s22=zeta*e22+(e11+e22)*(zeta-eta)-0.5_wp*P                                                                                                                               
      s12=zeta*e12
      s21=zeta*e21
    END SUBROUTINE VP


         
        

   SUBROUTINE interface_boundary_cell_marker(boundary_cell_marker, p_patch_3D, p_ice)
     TYPE(t_patch_3D ),TARGET, INTENT(IN)   :: p_patch_3D
     TYPE(t_patch),  POINTER :: p_patch
     TYPE(t_subset_range), POINTER :: all_cells
     TYPE(t_sea_ice),          INTENT(INOUT)  :: p_ice
    REAL(wp),TARGET,INTENT(inout) :: boundary_cell_marker(nproma,n_zlev,p_patch_3D%p_patch_2D(1)%alloc_cell_blocks)  
    Real(wp):: x,y,z
    INTEGER :: cell_block,start_index,end_index,cell_index
    !-------------------------------------------------------------------
    p_patch         => p_patch_3D%p_patch_2D(1)
    all_cells        =>p_patch%cells%all
    DO cell_block = all_cells%start_block, all_cells%end_block
       CALL get_index_range(all_cells, cell_block, start_index, end_index)
       DO cell_index =  start_index, end_index
           x= p_patch%cells%cartesian_center(cell_index,cell_block)%x(1)
           y= p_patch%cells%cartesian_center(cell_index,cell_block)%x(2)
           z= p_patch%cells%cartesian_center(cell_index,cell_block)%x(3)
           boundary_cell_marker(cell_index,1,cell_block)=0.0_wp
         if(x > -0.4 .AND. x< 0.4 .AND. y>-0.4 .AND. y< 0.4 .AND. z>0.0)then
            boundary_cell_marker(cell_index,1,cell_block)=1.0_wp
         endif

      END DO
    END DO
  END SUBROUTINE interface_boundary_cell_marker


SUBROUTINE interface_boundary_edge_marker(boundary_edge_marker,boundary_cell_marker, p_patch_3D, p_ice)
  TYPE(t_patch_3D ),TARGET, INTENT(IN)   :: p_patch_3D
  TYPE(t_subset_range), POINTER            :: all_edges                                   
  TYPE(t_patch),  POINTER :: p_patch
  TYPE(t_sea_ice),          INTENT(INOUT)  :: p_ice
  REAL(wp),TARGET,INTENT(in) :: boundary_cell_marker(nproma,n_zlev,p_patch_3D%p_patch_2D(1)%alloc_cell_blocks)
  REAL(wp),TARGET, INTENT(inout) :: boundary_edge_marker(nproma,p_patch_3d%p_patch_2d(1)%nblks_e)
  INTEGER :: cell_block1,cell_block2,start_index,end_index,cell_index1,cell_index2, edge_index_i,&
       &edge_block_i,cell_index_1,cell_index_2, doy
  p_patch         => p_patch_3D%p_patch_2D(1)
  all_edges => p_patch%edges%all

  DO edge_block_i = all_edges%start_block, all_edges%end_block
     CALL get_index_range(all_edges, edge_block_i, start_index, end_index)
     DO edge_index_i =  start_index, end_index
        cell_index1 = p_patch%edges%cell_idx(edge_index_i,edge_block_i,1)
        cell_block1 = p_patch%edges%cell_blk(edge_index_i,edge_block_i,1)
        cell_index2 = p_patch%edges%cell_idx(edge_index_i,edge_block_i,2)
        cell_block2 = p_patch%edges%cell_blk(edge_index_i,edge_block_i,2)
        boundary_edge_marker(edge_index_i,edge_block_i)=0.0_wp
        doy=boundary_cell_marker(cell_index1,1,cell_block1)+boundary_cell_marker(cell_index2,1,cell_block2)
        if( doy>1.0 ) then
           boundary_edge_marker(edge_index_i,edge_block_i)=1.0_wp
        endif
     ENDDO
  ENDDO
END SUBROUTINE interface_boundary_edge_marker

SUBROUTINE cell_area(cell_area_c,p_patch_3D)
  TYPE(t_patch_3D ),TARGET, INTENT(IN)   :: p_patch_3D
  TYPE(t_patch),  POINTER :: p_patch
  TYPE(t_subset_range), POINTER :: all_cells
  TYPE(t_subset_range), POINTER :: owned_cells
  REAL(wp),TARGET,INTENT(inout) :: cell_area_c(nproma,n_zlev,p_patch_3D%p_patch_2D(1)%alloc_cell_blocks)  
  INTEGER :: cell_block,start_index,end_index,cell_index,&
       &edge_index_1,edge_block_1,edge_index_2,edge_block_2,edge_index_3,edge_block_3
  Real(wp) :: s,e1,e2,e3,n1x,n2x,n3x,n1y,n2y,n3y,n1z,n2z,n3z
  p_patch         => p_patch_3D%p_patch_2D(1)
  all_cells => p_patch%cells%all
  owned_cells => p_patch%cells%owned
  DO cell_block = owned_cells%start_block, owned_cells%end_block
     CALL get_index_range(owned_cells, cell_block, start_index, end_index)
     DO cell_index = start_index, end_index
        edge_index_1 = p_patch%cells%edge_idx(cell_index, cell_block, 1)
        edge_block_1 = p_patch%cells%edge_blk(cell_index, cell_block, 1)
        edge_index_2 = p_patch%cells%edge_idx(cell_index, cell_block, 2)
        edge_block_2 = p_patch%cells%edge_blk(cell_index, cell_block, 2)
        edge_index_3 = p_patch%cells%edge_idx(cell_index, cell_block, 3)
        edge_block_3 = p_patch%cells%edge_blk(cell_index, cell_block, 3)
        e1=p_patch%edges%primal_edge_length(edge_index_1,edge_block_1)
        e2=p_patch%edges%primal_edge_length(edge_index_2,edge_block_2)
        e3=p_patch%edges%primal_edge_length(edge_index_3,edge_block_3)
        s=(e1+e2+e3)*0.5_wp
        cell_area_c(cell_index,1,cell_block)=sqrt(s*(s-e1)*(s-e2)*(s-e3)) 
     ENDDO ! cell_index = start_index, end_index                                                                                                          
  ENDDO !cell_block = owned_cells%start_block, owned_cells%end_block
END SUBROUTINE cell_area

    
SUBROUTINE init_mass_matrix(mass,cell_area_c,p_patch_3D)
  TYPE(t_patch_3D ),TARGET, INTENT(IN)   :: p_patch_3D
  TYPE(t_patch),  POINTER :: p_patch
  TYPE(t_subset_range), POINTER :: all_cells
  TYPE(t_subset_range), POINTER :: owned_cells
  REAL(wp),TARGET,INTENT(in) :: cell_area_c(nproma,n_zlev,p_patch_3D%p_patch_2D(1)%alloc_cell_blocks)  
  INTEGER :: cell_block,start_index,end_index,cell_index, edge_index_i,edge_block_i, neigbor
  REAL(wp), TARGET, INTENT(inout) :: mass(nproma,p_patch_3d%p_patch_2d(1)%nblks_e)
  p_patch         => p_patch_3D%p_patch_2D(1)
  all_cells => p_patch%cells%all
  owned_cells => p_patch%cells%owned
  DO cell_block = owned_cells%start_block, owned_cells%end_block
     CALL get_index_range(owned_cells, cell_block, start_index, end_index)
     DO cell_index = start_index, end_index
        DO neigbor=1,3!p_patch%cells%num_edges(cell_index,cell_block)!no_primal_edges
           edge_index_i = p_patch%cells%edge_idx(cell_index, cell_block, neigbor)
           edge_block_i = p_patch%cells%edge_blk(cell_index, cell_block, neigbor)
           mass(edge_index_i,edge_block_i)=mass(edge_index_i,edge_block_i)&
                &+cell_area_c(cell_index,1,cell_block)/3.0_wp
        ENDDO !neigbor=1,patch_2D%num_edges i
     ENDDO ! cell_index = start_index, end_index
  ENDDO !cell_block = owned_cells%start_block, owned_cells%end_block    
END SUBROUTINE init_mass_matrix

   

SUBROUTINE compute_2Dvector_grad(cell_index,cell_block,uxx,uxy,uyx,uyy,cell_area,x1_c,x2_c,x3_c,p_patch_3D, p_ice)
  TYPE(t_patch_3D ),TARGET, INTENT(IN)     :: p_patch_3D
  TYPE(t_patch),  POINTER                  :: p_patch
  TYPE(t_sea_ice),          INTENT(INOUT)  :: p_ice
  TYPE(t_subset_range), POINTER            :: all_edges
  TYPE(t_subset_range),     POINTER        :: all_cells
  TYPE(t_subset_range), POINTER            :: owned_cells
  INTEGER, INTENT(in) :: cell_index
  INTEGER, INTENT(in) :: cell_block
  REAL(wp), INTENT(in) :: cell_area
  REAL(wp), INTENT(in) :: x1_c
  REAL(wp), INTENT(in) :: x2_c
  REAL(wp), INTENT(in) :: x3_c
  REAL(wp), INTENT(out) :: uxx,uyy, uxy,uyx
  INTEGER :: edge_block_1,edge_index_1,&
       &edge_index_i, edge_block_i, neigbor,edge_index_j, edge_block_j, neigbor_j
  REAL(wp) :: ei,h_ei,Pi,Oi,l_sn,sc_pn,l_st,tix_1,tiy_1,tiz_1,tjx_1,tjy_1,tjz_1,correct_n,correct_t,&
       &nix_l,niy_l,tix_l, tiy_l,nix,niz,niy,tix,tiy,tiz,O1,nix_1,niy_1,niz_1,njx_1,njy_1,njz_1,&
       & U,V,Z,uxx_l,uxy_l,uyy_l,uyx_l,uzx_l,uzy_l, tix_k,tiy_k,tiz_k,nix_k,niy_k,niz_k, Ux,Vx,Uy,Vy,Uz,Vz 
    p_patch   => p_patch_3D%p_patch_2D(1)
    all_cells => p_patch%cells%all
    owned_cells => p_patch%cells%owned
    all_edges => p_patch%edges%all
    uxx_l=0.0_wp
    uxy_l=0.0_wp
    uyx_l=0.0_wp
    uyy_l=0.0_wp
    uzx_l=0.0_wp
    uzy_l=0.0_wp

    DO neigbor=1,3 !no_primal_edges                                                                                                                                                                               
       edge_index_i =p_patch%cells%edge_idx(cell_index, cell_block, neigbor)
       edge_block_i =p_patch%cells%edge_blk(cell_index, cell_block, neigbor)
       ei=p_patch%edges%primal_edge_length(edge_index_i,edge_block_i)
       h_ei = 2.0_wp*cell_area/ei
       Pi=2.0_wp/h_ei
       Oi=p_patch%cells%edge_orientation(cell_index,cell_block,neigbor)
       nix=p_patch%edges%primal_cart_normal(edge_index_i,edge_block_i)%x(1)
       niy=p_patch%edges%primal_cart_normal(edge_index_i,edge_block_i)%x(2)
       niz=p_patch%edges%primal_cart_normal(edge_index_i,edge_block_i)%x(3)
       tix=p_patch%edges%dual_cart_normal(edge_index_i,edge_block_i)%x(1)
       tiy=p_patch%edges%dual_cart_normal(edge_index_i,edge_block_i)%x(2)
       tiz=p_patch%edges%dual_cart_normal(edge_index_i,edge_block_i)%x(3)
       nix=nix*Oi
       niy=niy*Oi
       niz=niz*Oi
       tix=tix*Oi
       tiy=tiy*Oi
       tiz=tiz*Oi
       if(neigbor==1) then
          O1=Oi
          sc_pn=nix*x1_c + niy*x2_c + niz*x3_c
          nix_1=nix-sc_pn*x1_c
          niy_1=niy-sc_pn*x2_c
          niz_1=niz-sc_pn*x3_c
          l_sn=sqrt(nix_1*nix_1 + niy_1*niy_1 + niz_1*niz_1)
          nix_1=nix_1/l_sn
          niy_1=niy_1/l_sn
          niz_1=niz_1/l_sn
          sc_pn=tix*x1_c + tiy*x2_c + tiz*x3_c
          tix_1=tix-sc_pn*x1_c
          tiy_1=tiy-sc_pn*x2_c
          tiz_1=tiz-sc_pn*x3_c
          l_sn=sqrt(tix_1*tix_1 + tiy_1*tiy_1 + tiz_1*tiz_1)
          tix_1=tix_1/l_sn
          tiy_1=tiy_1/l_sn
          tiz_1=tiz_1/l_sn
       endif
       !Projection in plane     
       sc_pn=nix*x1_c + niy*x2_c + niz*x3_c
       nix=nix-sc_pn*x1_c
       niy=niy-sc_pn*x2_c
       niz=niz-sc_pn*x3_c
       l_sn=sqrt(nix*nix + niy*niy + niz*niz)
       nix=nix/l_sn
       niy=niy/l_sn
       niz=niz/l_sn
       
       sc_pn=tix*x1_c + tiy*x2_c + tiz*x3_c
       tix=tix-sc_pn*x1_c
       tiy=tiy-sc_pn*x2_c
       tiz=tiz-sc_pn*x3_c          
       l_sn=sqrt(tix*tix + tiy*tiy + tiz*tiz)
       tix=tix/l_sn
       tiy=tiy/l_sn
       tiz=tiz/l_sn
       !Projection end
       nix_l =   tix_1*nix + tiy_1*niy + tiz_1*niz
       niy_l = - nix_1*nix - niy_1*niy - niz_1*niz
       tix_l =   tix_1*tix + tiy_1*tiy + tiz_1*tiz
       tiy_l = - nix_1*tix - niy_1*tiy - niz_1*tiz
       U=(p_ice%vn_e(edge_index_i,edge_block_i)*nix_l + p_ice%vt_e(edge_index_i,edge_block_i)*tix_l)*Oi          
       V=(p_ice%vn_e(edge_index_i,edge_block_i)*niy_l + p_ice%vt_e(edge_index_i,edge_block_i)*tiy_l)*Oi
       ! 0.5*(nabal u +nabla u^T)
       uxx_l=uxx_l+1.0_wp*nix_l*Pi*U
       uxy_l=uxy_l+0.5_wp*niy_l*Pi*U+0.5_wp*nix_l*Pi*V
       uyx_l=uxy_l
       uyy_l=uyy_l+1.0_wp*niy_l*Pi*V
    ENDDO !neigbor=1,patch_2D%num_edges i
    uxx=uxx_l
    uxy=uxy_l
    uyx=uyx_l
    uyy=uyy_l
  END SUBROUTINE compute_2Dvector_grad

  SUBROUTINE compute_2D_vector_Laplace(Au_n,Au_t,boundary_cell_marker,s11,s12,s21,s22,s31,s32,cell_area_c,x1_c,x2_c,x3_c,p_patch_3D, p_ice)
    TYPE(t_patch_3D ),TARGET, INTENT(IN)     :: p_patch_3D
    TYPE(t_patch),  POINTER                  :: p_patch
    TYPE(t_sea_ice),          INTENT(INOUT)  :: p_ice
    TYPE(t_subset_range), POINTER            :: all_edges
    TYPE(t_subset_range),     POINTER        :: all_cells
    TYPE(t_subset_range), POINTER            :: owned_cells
    REAL(wp),TARGET, INTENT(in) :: x1_c(nproma,n_zlev,p_patch_3D%p_patch_2D(1)%alloc_cell_blocks) 
    REAL(wp),TARGET, INTENT(in) :: x2_c(nproma,n_zlev,p_patch_3D%p_patch_2D(1)%alloc_cell_blocks) 
    REAL(wp),TARGET, INTENT(in) :: x3_c(nproma,n_zlev,p_patch_3D%p_patch_2D(1)%alloc_cell_blocks)
    REAL(wp),TARGET,INTENT(in) :: boundary_cell_marker(nproma,n_zlev,p_patch_3D%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp),TARGET,INTENT(in) :: s11(nproma,n_zlev,p_patch_3D%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp),TARGET,INTENT(in) :: s12(nproma,n_zlev,p_patch_3D%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp),TARGET,INTENT(in) :: s22(nproma,n_zlev,p_patch_3D%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp),TARGET,INTENT(in) :: s21(nproma,n_zlev,p_patch_3D%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp),TARGET,INTENT(in) ::s31(nproma,n_zlev,p_patch_3D%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp),TARGET,INTENT(in) :: s32(nproma,n_zlev,p_patch_3D%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp), TARGET, INTENT(out) ::  Au_n(nproma,p_patch_3d%p_patch_2d(1)%nblks_e)
    REAL(wp), TARGET, INTENT(out) ::  Au_t(nproma,p_patch_3d%p_patch_2d(1)%nblks_e)
    REAL(wp),TARGET,INTENT(in) :: cell_area_c(nproma,n_zlev,p_patch_3D%p_patch_2D(1)%alloc_cell_blocks)    
    INTEGER :: edge_block_1,edge_index_1,cell_index,cell_block,start_index,end_index,&
         &edge_index_i, edge_block_i, neigbor,edge_index_j, edge_block_j, neigbor_j
    
    REAL(wp) :: ei,h_ei,Pi,Oi,l_sn,sc_pn,tix_1,tiy_1,tiz_1,tjx_1,tjy_1,tjz_1,&
         &nix_l,niy_l,nix,niz,niy,tix,tiy,tiz,tix_l, tiy_l,O1,nix_1,niy_1,niz_1,njx_1,njy_1,njz_1,&
         &cell_area,x1,x2,x3
    
    
    p_patch   => p_patch_3D%p_patch_2D(1)
    all_cells => p_patch%cells%all
    owned_cells => p_patch%cells%owned
    all_edges => p_patch%edges%all
    
    DO cell_block = owned_cells%start_block, owned_cells%end_block
       CALL get_index_range(owned_cells, cell_block, start_index, end_index)
       DO cell_index = start_index, end_index

          if (boundary_cell_marker(cell_index,1,cell_block)>0.0) then
             cell_area = cell_area_c(cell_index,1,cell_block)
             x1=x1_c(cell_index,1,cell_block)
             x2=x2_c(cell_index,1,cell_block)
             x3=x3_c(cell_index,1,cell_block)
             DO neigbor=1,3 !no_primal_edges
                edge_index_i =p_patch%cells%edge_idx(cell_index, cell_block, neigbor)
                edge_block_i =p_patch%cells%edge_blk(cell_index, cell_block, neigbor)
                ei=p_patch%edges%primal_edge_length(edge_index_i,edge_block_i)
                h_ei = 2.0_wp*cell_area/ei
                Pi=2.0_wp/h_ei
                Oi=p_patch%cells%edge_orientation(cell_index,cell_block,neigbor)
                nix=p_patch%edges%primal_cart_normal(edge_index_i,edge_block_i)%x(1)*Oi
                niy=p_patch%edges%primal_cart_normal(edge_index_i,edge_block_i)%x(2)*Oi
                niz=p_patch%edges%primal_cart_normal(edge_index_i,edge_block_i)%x(3)*Oi
                tix=p_patch%edges%dual_cart_normal(edge_index_i,edge_block_i)%x(1)*Oi
                tiy=p_patch%edges%dual_cart_normal(edge_index_i,edge_block_i)%x(2)*Oi
                tiz=p_patch%edges%dual_cart_normal(edge_index_i,edge_block_i)%x(3)*Oi
                if(neigbor==1) then
                   O1=Oi
                   sc_pn=nix*x1 + niy*x2 + niz*x3
                   nix_1=nix-sc_pn*x1
                   niy_1=niy-sc_pn*x2
                   niz_1=niz-sc_pn*x3
                   l_sn=sqrt(nix_1*nix_1 + niy_1*niy_1 + niz_1*niz_1)
                   nix_1=nix_1/l_sn
                   niy_1=niy_1/l_sn
                   niz_1=niz_1/l_sn
                   sc_pn=tix*x1 + tiy*x2 + tiz*x3
                   tix_1=tix-sc_pn*x1
                   tiy_1=tiy-sc_pn*x2
                   tiz_1=tiz-sc_pn*x3
                   l_sn=sqrt(tix_1*tix_1 + tiy_1*tiy_1 + tiz_1*tiz_1)
                   tix_1=tix_1/l_sn
                   tiy_1=tiy_1/l_sn
                   tiz_1=tiz_1/l_sn
                endif
                sc_pn=nix*x1 + niy*x2 + niz*x3
                nix=nix-sc_pn*x1
                niy=niy-sc_pn*x2
                niz=niz-sc_pn*x3
                l_sn=sqrt(nix*nix + niy*niy + niz*niz)
                nix=nix/l_sn
                niy=niy/l_sn
                niz=niz/l_sn
                sc_pn=tix*x1 + tiy*x2 + tiz*x3
                tix=tix-sc_pn*x1
                tiy=tiy-sc_pn*x2
                tiz=tiz-sc_pn*x3
                l_sn=sqrt(tix*tix + tiy*tiy + tiz*tiz)
                tix=tix/l_sn
                tiy=tiy/l_sn
                tiz=tiz/l_sn
             
             nix_l =   tix_1*nix + tiy_1*niy + tiz_1*niz
             niy_l = - nix_1*nix - niy_1*niy - niz_1*niz
             tix_l =   tix_1*tix + tiy_1*tiy + tiz_1*tiz
             tiy_l = - nix_1*tix - niy_1*tiy - niz_1*tiz
             Au_n(edge_index_i,edge_block_i)=Au_n(edge_index_i,edge_block_i)+&
                  &Pi*cell_area*Oi*(&
                  & s11(cell_index,1,cell_block)*nix_l*nix_l+s12(cell_index,1,cell_block)*nix_l*niy_l&
                  &+s21(cell_index,1,cell_block)*niy_l*nix_l+s22(cell_index,1,cell_block)*niy_l*niy_l)
             Au_t(edge_index_i,edge_block_i)=Au_t(edge_index_i,edge_block_i)+&
                  &Pi*cell_area*Oi*(&
                  & s11(cell_index,1,cell_block)*tix_l*nix_l+s12(cell_index,1,cell_block)*tix_l*niy_l&
                  &+s21(cell_index,1,cell_block)*tiy_l*nix_l+s22(cell_index,1,cell_block)*tiy_l*niy_l)
          ENDDO !neigbor=1,patch_2D%num_edges i
       endif
       ENDDO
    ENDDO
  END SUBROUTINE compute_2D_vector_Laplace


 
  SUBROUTINE CG_addvector(aaa,ice_x,Dn,p_patch_3D)
    TYPE(t_patch_3D ),TARGET, INTENT(IN)     :: p_patch_3D
    TYPE(t_patch),  POINTER                  :: p_patch
    TYPE(t_subset_range), POINTER            :: all_edges
    TYPE(t_subset_range),     POINTER        :: all_cells
    TYPE(t_subset_range), POINTER            :: owned_cells
    REAL(wp), INTENT(in) :: aaa
    REAL(wp), TARGET, INTENT(inout) ::  ice_x(nproma,p_patch_3d%p_patch_2d(1)%nblks_e)
    REAL(wp), TARGET, INTENT(in) ::  Dn(nproma,p_patch_3d%p_patch_2d(1)%nblks_e)
    INTEGER :: cell_index,edge_block_i,start_index,end_index,&
         &edge_index_i
    REAL(wp) :: ei,h_ei,Pi,Oi,l_sn,sc_pn,tix_1,tiy_1,tiz_1,tjx_1,tjy_1,tjz_1,&
         &nix_l,niy_l,nix,niz,niy,tix,tiy,tiz,tix_l, tiy_l,O1,nix_1,niy_1,niz_1,njx_1,njy_1,njz_1,&
         &cell_area,x1,x2,x3,R
    p_patch   => p_patch_3D%p_patch_2D(1)
    all_cells => p_patch%cells%all
    owned_cells => p_patch%cells%owned
    all_edges => p_patch%edges%all
    DO edge_block_i = all_edges%start_block, all_edges%end_block
       CALL get_index_range(all_edges, edge_block_i, start_index, end_index)
       DO edge_index_i =  start_index, end_index
          ice_x(edge_index_i,edge_block_i)=ice_x(edge_index_i,edge_block_i)+&
               &aaa*Dn(edge_index_i,edge_block_i)            
       ENDDO
    ENDDO
  END SUBROUTINE CG_addvector
! d = r + bbb * d                                                                                                                                                                                                                                                      
       
 SUBROUTINE CG_addvector2(bbb,Rn,Dn,p_patch_3D)
    TYPE(t_patch_3D ),TARGET, INTENT(IN)     :: p_patch_3D
    TYPE(t_patch),  POINTER                  :: p_patch
    TYPE(t_subset_range), POINTER            :: all_edges
    TYPE(t_subset_range),     POINTER        :: all_cells
    TYPE(t_subset_range), POINTER            :: owned_cells
    REAL(wp), INTENT(in) :: bbb
    REAL(wp), TARGET, INTENT(in) ::  Rn(nproma,p_patch_3d%p_patch_2d(1)%nblks_e)
    REAL(wp), TARGET, INTENT(inout) ::  Dn(nproma,p_patch_3d%p_patch_2d(1)%nblks_e)
    INTEGER :: cell_index,edge_block_i,start_index,end_index,&
         &edge_index_i
    REAL(wp) :: ei,h_ei,Pi,Oi,l_sn,sc_pn,tix_1,tiy_1,tiz_1,tjx_1,tjy_1,tjz_1,&
         &nix_l,niy_l,nix,niz,niy,tix,tiy,tiz,tix_l, tiy_l,O1,nix_1,niy_1,niz_1,njx_1,njy_1,njz_1,&
         &cell_area,x1,x2,x3,R
    p_patch   => p_patch_3D%p_patch_2D(1)
    all_cells => p_patch%cells%all
    owned_cells => p_patch%cells%owned
    all_edges => p_patch%edges%all
    DO edge_block_i = all_edges%start_block, all_edges%end_block
       CALL get_index_range(all_edges, edge_block_i, start_index, end_index)
       DO edge_index_i =  start_index, end_index
          Dn(edge_index_i,edge_block_i)=Rn(edge_index_i,edge_block_i)+&
                  &bbb*Dn(edge_index_i,edge_block_i)
       ENDDO
    ENDDO
  END SUBROUTINE CG_addvector2
  

SUBROUTINE CG_scalarp(alpha_n,alpha_t,s1n,s2n,s1t,s2t,p_patch_3D)
    TYPE(t_patch_3D ),TARGET, INTENT(IN)     :: p_patch_3D
    TYPE(t_patch),  POINTER                  :: p_patch
    TYPE(t_subset_range), POINTER            :: all_edges
    TYPE(t_subset_range),     POINTER        :: all_cells
    TYPE(t_subset_range), POINTER            :: owned_cells

    REAL(wp), INTENT(out) :: alpha_n, alpha_t
    REAL(wp), TARGET, INTENT(in) ::  s1n(nproma,p_patch_3d%p_patch_2d(1)%nblks_e)
    REAL(wp), TARGET, INTENT(in) ::  s2n(nproma,p_patch_3d%p_patch_2d(1)%nblks_e)
    REAL(wp), TARGET, INTENT(in) ::  s1t(nproma,p_patch_3d%p_patch_2d(1)%nblks_e)
    REAL(wp), TARGET, INTENT(in) ::  s2t(nproma,p_patch_3d%p_patch_2d(1)%nblks_e)
    INTEGER :: cell_index,edge_block_i,start_index,end_index,&
         &edge_index_i
    REAL(wp) :: ei,h_ei,Pi,Oi,l_sn,sc_pn,tix_1,tiy_1,tiz_1,tjx_1,tjy_1,tjz_1,&
         &nix_l,niy_l,nix,niz,niy,tix,tiy,tiz,tix_l, tiy_l,O1,nix_1,niy_1,niz_1,njx_1,njy_1,njz_1,&
         &cell_area,x1,x2,x3,R
    R=6.371229e6_wp
    alpha_n = 0.0_wp
    alpha_t = 0.0_wp
    p_patch   => p_patch_3D%p_patch_2D(1)
    all_cells => p_patch%cells%all
    owned_cells => p_patch%cells%owned
    all_edges => p_patch%edges%all
    DO edge_block_i = all_edges%start_block, all_edges%end_block
       CALL get_index_range(all_edges, edge_block_i, start_index, end_index)
       DO edge_index_i =  start_index, end_index
          alpha_n=alpha_n+&
               &s1n(edge_index_i,edge_block_i)*s2n(edge_index_i,edge_block_i)
          alpha_t=alpha_t+&
               &s1t(edge_index_i,edge_block_i)*s2t(edge_index_i,edge_block_i)      
       ENDDO
    ENDDO
  END SUBROUTINE CG_scalarp
  
  SUBROUTINE Stabilization_sum(S_x,S_y,boundary_cell_marker,x1_c,x2_c,x3_c,p_ice,p_patch_3D)
    TYPE(t_patch_3D), TARGET, INTENT(IN)     :: p_patch_3D
    TYPE(t_sea_ice),          INTENT(INOUT)  :: p_ice
    TYPE(t_patch),  POINTER                  :: p_patch
   
    REAL(wp), TARGET, INTENT(inout) ::  S_x(nproma,p_patch_3d%p_patch_2d(1)%nblks_e)
    REAL(wp), TARGET, INTENT(inout) ::  S_y(nproma,p_patch_3d%p_patch_2d(1)%nblks_e)
    REAL(wp),TARGET,INTENT(in) :: boundary_cell_marker(nproma,n_zlev,p_patch_3D%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp),TARGET, INTENT(in) :: x1_c(nproma,n_zlev,p_patch_3D%p_patch_2D(1)%alloc_cell_blocks) 
    REAL(wp),TARGET, INTENT(in) :: x2_c(nproma,n_zlev,p_patch_3D%p_patch_2D(1)%alloc_cell_blocks) 
    REAL(wp),TARGET, INTENT(in) :: x3_c(nproma,n_zlev,p_patch_3D%p_patch_2D(1)%alloc_cell_blocks)
    REAL(wp) ::nx1,ny1,nx2,ny2,nx3,ny3,tx1,ty1,tx2,ty2,tx3,ty3
    INTEGER :: cell_block, start_index,  end_index, cell_index,&
        &edge_index_1,edge_index_2,edge_index_3,edge_block_1,edge_block_2,edge_block_3
    INTEGER :: edge_index_i, edge_block_i, neigbor
    REAL(wp) :: Oi,l_sn,sc_pn,tix_1,tiy_1,tiz_1,&
        &nix,niz,niy,tix,tiy,tiz,O1,O2,O3,nix_1,niy_1,niz_1,&
        &cell_area,x1,x2,x3,  Sx1,Sx2,Sx3,Sy1,Sy2,Sy3, a1,a2,a3,z1,z2,z3
    TYPE(t_subset_range), POINTER :: owned_cells
    TYPE(t_subset_range), POINTER :: all_edges
    p_patch   => p_patch_3D%p_patch_2D(1)
    owned_cells =>p_patch%cells%owned
    all_edges => p_patch%edges%all
    DO cell_block = owned_cells%start_block, owned_cells%end_block
       CALL get_index_range(owned_cells, cell_block, start_index, end_index)
       DO cell_index = start_index, end_index
          if(boundary_cell_marker(cell_index,1,cell_block)>0.0)then
             x1=x1_c(cell_index,1,cell_block)
             x2=x2_c(cell_index,1,cell_block)
             x3=x3_c(cell_index,1,cell_block)
             DO neigbor=1,3 !no_primal_edges
                edge_index_i =p_patch%cells%edge_idx(cell_index, cell_block, neigbor)
                edge_block_i =p_patch%cells%edge_blk(cell_index, cell_block, neigbor)
                Oi=p_patch%cells%edge_orientation(cell_index,cell_block,neigbor)
                nix=p_patch%edges%primal_cart_normal(edge_index_i,edge_block_i)%x(1)*Oi
                niy=p_patch%edges%primal_cart_normal(edge_index_i,edge_block_i)%x(2)*Oi
                niz=p_patch%edges%primal_cart_normal(edge_index_i,edge_block_i)%x(3)*Oi
                tix=p_patch%edges%dual_cart_normal(edge_index_i,edge_block_i)%x(1)*Oi
                tiy=p_patch%edges%dual_cart_normal(edge_index_i,edge_block_i)%x(2)*Oi
                tiz=p_patch%edges%dual_cart_normal(edge_index_i,edge_block_i)%x(3)*Oi
                if(neigbor==1) then
                   O1=Oi
                   sc_pn=nix*x1 + niy*x2 + niz*x3
                   nix_1=nix-sc_pn*x1
                   niy_1=niy-sc_pn*x2
                   niz_1=niz-sc_pn*x3
                   l_sn=sqrt(nix_1*nix_1 + niy_1*niy_1 + niz_1*niz_1)
                   nix_1=nix_1/l_sn
                   niy_1=niy_1/l_sn
                   niz_1=niz_1/l_sn
                   sc_pn=tix*x1 + tiy*x2 + tiz*x3
                   tix_1=tix-sc_pn*x1
                   tiy_1=tiy-sc_pn*x2
                   tiz_1=tiz-sc_pn*x3
                   l_sn=sqrt(tix_1*tix_1 + tiy_1*tiy_1 + tiz_1*tiz_1)
                   tix_1=tix_1/l_sn
                   tiy_1=tiy_1/l_sn
                   tiz_1=tiz_1/l_sn
                endif
                sc_pn=nix*x1 + niy*x2 + niz*x3
                nix=nix-sc_pn*x1
                niy=niy-sc_pn*x2
                niz=niz-sc_pn*x3
                l_sn=sqrt(nix*nix + niy*niy + niz*niz)
                nix=nix/l_sn
                niy=niy/l_sn
                niz=niz/l_sn
                sc_pn=tix*x1 + tiy*x2 + tiz*x3
                tix=tix-sc_pn*x1
                tiy=tiy-sc_pn*x2
                tiz=tiz-sc_pn*x3
                l_sn=sqrt(tix*tix + tiy*tiy + tiz*tiz)
                tix=tix/l_sn
                tiy=tiy/l_sn
                tiz=tiz/l_sn
                if(neigbor==1) then
                   nx1 =(   tix_1*nix + tiy_1*niy + tiz_1*niz)*Oi
                   ny1 =( - nix_1*nix - niy_1*niy - niz_1*niz)*Oi
                   tx1 = (  tix_1*tix + tiy_1*tiy + tiz_1*tiz)*Oi
                   ty1 = (- nix_1*tix - niy_1*tiy - niz_1*tiz)*Oi
                endif
                if(neigbor==2) then
                   O2=Oi
                   nx2 =  ( tix_1*nix + tiy_1*niy + tiz_1*niz)*Oi
                   ny2 = (- nix_1*nix - niy_1*niy - niz_1*niz)*Oi
                   tx2 =  ( tix_1*tix + tiy_1*tiy + tiz_1*tiz)*Oi
                   ty2 = (- nix_1*tix - niy_1*tiy - niz_1*tiz)*Oi
                endif
                if(neigbor==3) then
                   O3=Oi
                   nx3 =  ( tix_1*nix + tiy_1*niy + tiz_1*niz)*Oi
                   ny3 = (- nix_1*nix - niy_1*niy - niz_1*niz)*Oi
                   tx3 =  ( tix_1*tix + tiy_1*tiy + tiz_1*tiz)*Oi
                   ty3 = (- nix_1*tix - niy_1*tiy - niz_1*tiz)*Oi
                endif
             ENDDO !neigbor=1,patch_2D%num_edges i     
             edge_index_1 = p_patch%cells%edge_idx(cell_index, cell_block, 1)
             edge_block_1 = p_patch%cells%edge_blk(cell_index, cell_block, 1)
             edge_index_2 = p_patch%cells%edge_idx(cell_index, cell_block, 2)
             edge_block_2 = p_patch%cells%edge_blk(cell_index, cell_block, 2)
             edge_index_3 = p_patch%cells%edge_idx(cell_index, cell_block, 3)
             edge_block_3 = p_patch%cells%edge_blk(cell_index, cell_block, 3)
             z1= p_patch%edges%cartesian_center(edge_index_1,edge_block_1)%x(3)
             z2= p_patch%edges%cartesian_center(edge_index_2,edge_block_2)%x(3)
             z3= p_patch%edges%cartesian_center(edge_index_3,edge_block_3)%x(3)
             a1=1.0_wp
             a2=1.0_wp
             a3=1.0_wp
             if((-1.e8<z1 < 1.e-8) .and. (x3>0.0)) then
                a1=-1.0_wp
             endif
             if((-1.e8 <z2 < 1.e-8) .and. (x3>0.0)) then
                a2=-1.0_wp
             endif
             if((-1.e8 <z3 < 1.e-8) .and. (x3>0.0)) then
                a3=-1.0_wp
             endif
             Sx1 = p_ice%vn_e(edge_index_2,edge_block_2)*(nx2)-p_ice%vn_e(edge_index_3,edge_block_3)*(nx3)+&
                  &p_ice%vt_e(edge_index_2,edge_block_2)*(tx2)-p_ice%vt_e(edge_index_3,edge_block_3)*(tx3)
             Sx2 = p_ice%vn_e(edge_index_3,edge_block_3)*(nx3)-p_ice%vn_e(edge_index_1,edge_block_1)*(nx1)+&
                  &p_ice%vt_e(edge_index_3,edge_block_3)*(tx3)-p_ice%vt_e(edge_index_1,edge_block_1)*(tx1)
             Sx3 = p_ice%vn_e(edge_index_1,edge_block_1)*(nx1)-p_ice%vn_e(edge_index_2,edge_block_2)*(nx2)+&
                  &p_ice%vt_e(edge_index_1,edge_block_1)*(tx1)-p_ice%vt_e(edge_index_2,edge_block_2)*(tx2) 
             Sy1 = p_ice%vn_e(edge_index_2,edge_block_2)*(ny2)-p_ice%vn_e(edge_index_3,edge_block_3)*(ny3)+&
                  &p_ice%vt_e(edge_index_2,edge_block_2)*(ty2)-p_ice%vt_e(edge_index_3,edge_block_3)*(ty3)
             Sy2 = p_ice%vn_e(edge_index_3,edge_block_3)*(ny3)-p_ice%vn_e(edge_index_1,edge_block_1)*(ny1)+&
                  &p_ice%vt_e(edge_index_3,edge_block_3)*(ty3)-p_ice%vt_e(edge_index_1,edge_block_1)*(ty1)
             Sy3 = p_ice%vn_e(edge_index_1,edge_block_1)*(ny1)-p_ice%vn_e(edge_index_2,edge_block_2)*(ny2)+&
                  &p_ice%vt_e(edge_index_1,edge_block_1)*(ty1)-p_ice%vt_e(edge_index_2,edge_block_2)*(ty2)                            
             S_x(edge_index_1,edge_block_1)=S_x(edge_index_1,edge_block_1)+ (Sx1 * nx1 + Sy1 * ny1)*a1
             S_x(edge_index_2,edge_block_2)=S_x(edge_index_2,edge_block_2)+ (Sx2 * nx2 + Sy2 * ny2)*a2
             S_x(edge_index_3,edge_block_3)=S_x(edge_index_3,edge_block_3)+ (Sx3 * nx3 + Sy3 * ny3)*a3
             S_y(edge_index_1,edge_block_1)=S_y(edge_index_1,edge_block_1)+ (Sx1 * tx1 + Sy1 * ty1)*a1
             S_y(edge_index_2,edge_block_2)=S_y(edge_index_2,edge_block_2)+ (Sx2 * tx2 + Sy2 * ty2)*a2
             S_y(edge_index_3,edge_block_3)=S_y(edge_index_3,edge_block_3)+ (Sx3 * tx3 + Sy3 * ty3)*a3
          endif
       ENDDO ! cell_index = start_index, end_index
    ENDDO !cell_block = owned_cells%start_block, owned_cells%end_block
  END SUBROUTINE Stabilization_sum

  
  SUBROUTINE Stabilization(S_x,S_y,S_n,S_t,boundary_cell_marker,x1_c,x2_c,x3_c,zeta_e,p_patch_3D)
    TYPE(t_patch_3D), TARGET, INTENT(IN)     :: p_patch_3D
    TYPE(t_patch),  POINTER                  :: p_patch
    REAL(wp), TARGET, INTENT(inout) ::  S_n(nproma,p_patch_3d%p_patch_2d(1)%nblks_e)
    REAL(wp), TARGET, INTENT(inout) ::  S_t(nproma,p_patch_3d%p_patch_2d(1)%nblks_e)
    REAL(wp), TARGET, INTENT(in) ::  S_x(nproma,p_patch_3d%p_patch_2d(1)%nblks_e)
    REAL(wp), TARGET, INTENT(in) ::  S_y(nproma,p_patch_3d%p_patch_2d(1)%nblks_e)
    REAL(wp),TARGET, INTENT(in) :: boundary_cell_marker(nproma,n_zlev,p_patch_3D%p_patch_2D(1)%alloc_cell_blocks) 
    REAL(wp),TARGET, INTENT(in) :: x1_c(nproma,n_zlev,p_patch_3D%p_patch_2D(1)%alloc_cell_blocks) 
    REAL(wp),TARGET, INTENT(in) :: x2_c(nproma,n_zlev,p_patch_3D%p_patch_2D(1)%alloc_cell_blocks) 
    REAL(wp),TARGET, INTENT(in) :: x3_c(nproma,n_zlev,p_patch_3D%p_patch_2D(1)%alloc_cell_blocks)
    REAL(wp), TARGET, INTENT(in) ::  zeta_e(nproma,p_patch_3d%p_patch_2d(1)%nblks_e) 
    REAL(wp) ::nx1,ny1,nx2,ny2,nx3,ny3,tx1,ty1,tx2,ty2,tx3,ty3,ei,alpha,O2,O3
    INTEGER :: cell_block, start_index,  end_index, cell_index,&
         &edge_index_1,edge_index_2,edge_index_3,edge_block_1,edge_block_2,edge_block_3
    INTEGER :: edge_index_i, edge_block_i, neigbor
    REAL(wp) :: Oi,l_sn,sc_pn,tix_1,tiy_1,tiz_1,&
         &nix,niz,niy,tix,tiy,tiz,O1,nix_1,niy_1,niz_1,&
         &cell_area,x1,x2,x3, Sx1,Sx2,Sx3,Sy1,Sy2,Sy3, a1,a2,a3,z1,z2,z3, zeta_1,zeta_2,zeta_3 
    TYPE(t_subset_range), POINTER :: owned_cells
    TYPE(t_subset_range), POINTER :: all_edges
    p_patch   => p_patch_3D%p_patch_2D(1)
    owned_cells => p_patch%cells%owned
    all_edges => p_patch%edges%all
    DO cell_block = owned_cells%start_block, owned_cells%end_block
       CALL get_index_range(owned_cells, cell_block, start_index, end_index)
       DO cell_index = start_index, end_index
          if(boundary_cell_marker(cell_index,1,cell_block)>0.0) then 
          x1=x1_c(cell_index,1,cell_block)
          x2=x2_c(cell_index,1,cell_block)
          x3=x3_c(cell_index,1,cell_block)
          DO neigbor=1,3 !no_primal_edges
             edge_index_i =p_patch%cells%edge_idx(cell_index, cell_block, neigbor)
             edge_block_i =p_patch%cells%edge_blk(cell_index, cell_block, neigbor)
             Oi=p_patch%cells%edge_orientation(cell_index,cell_block,neigbor)
             nix=p_patch%edges%primal_cart_normal(edge_index_i,edge_block_i)%x(1)*Oi
             niy=p_patch%edges%primal_cart_normal(edge_index_i,edge_block_i)%x(2)*Oi
             niz=p_patch%edges%primal_cart_normal(edge_index_i,edge_block_i)%x(3)*Oi
             tix=p_patch%edges%dual_cart_normal(edge_index_i,edge_block_i)%x(1)*Oi
             tiy=p_patch%edges%dual_cart_normal(edge_index_i,edge_block_i)%x(2)*Oi
             tiz=p_patch%edges%dual_cart_normal(edge_index_i,edge_block_i)%x(3)*Oi
             if(neigbor==1) then
                O1=Oi
                sc_pn=nix*x1 + niy*x2 + niz*x3
                nix_1=nix-sc_pn*x1
                niy_1=niy-sc_pn*x2
                niz_1=niz-sc_pn*x3
                l_sn=sqrt(nix_1*nix_1 + niy_1*niy_1 + niz_1*niz_1)
                nix_1=nix_1/l_sn
                niy_1=niy_1/l_sn
                niz_1=niz_1/l_sn
                sc_pn=tix*x1 + tiy*x2 + tiz*x3
                tix_1=tix-sc_pn*x1
                tiy_1=tiy-sc_pn*x2
                tiz_1=tiz-sc_pn*x3
                l_sn=sqrt(tix_1*tix_1 + tiy_1*tiy_1 + tiz_1*tiz_1)
                tix_1=tix_1/l_sn
                tiy_1=tiy_1/l_sn
                tiz_1=tiz_1/l_sn
             endif
             sc_pn=nix*x1 + niy*x2 + niz*x3
             nix=nix-sc_pn*x1
             niy=niy-sc_pn*x2
             niz=niz-sc_pn*x3
             l_sn=sqrt(nix*nix + niy*niy + niz*niz)
             nix=nix/l_sn
             niy=niy/l_sn
             niz=niz/l_sn
             sc_pn=tix*x1 + tiy*x2 + tiz*x3
             tix=tix-sc_pn*x1
             tiy=tiy-sc_pn*x2
             tiz=tiz-sc_pn*x3
             l_sn=sqrt(tix*tix + tiy*tiy + tiz*tiz)
             tix=tix/l_sn
             tiy=tiy/l_sn
             tiz=tiz/l_sn
             if(neigbor==1) then
                nx1 = (   tix_1*nix + tiy_1*niy + tiz_1*niz)*Oi
                ny1 = ( - nix_1*nix - niy_1*niy - niz_1*niz)*Oi
                
                tx1 =   ( tix_1*tix + tiy_1*tiy + tiz_1*tiz)*Oi
                ty1 = ( - nix_1*tix - niy_1*tiy - niz_1*tiz)*Oi
             endif
             if(neigbor==2) then
                O2=Oi
                nx2 =   ( tix_1*nix + tiy_1*niy + tiz_1*niz)*Oi
                ny2 = ( - nix_1*nix - niy_1*niy - niz_1*niz)*Oi
                
                tx2 =   ( tix_1*tix + tiy_1*tiy + tiz_1*tiz)*Oi
                ty2 = ( - nix_1*tix - niy_1*tiy - niz_1*tiz)*Oi
             endif
             if(neigbor==3) then
                O3=Oi
                nx3 =   ( tix_1*nix + tiy_1*niy + tiz_1*niz)*Oi
                ny3 = ( - nix_1*nix - niy_1*niy - niz_1*niz)*Oi
                
                tx3 =   ( tix_1*tix + tiy_1*tiy + tiz_1*tiz)*Oi
                ty3 = ( - nix_1*tix - niy_1*tiy - niz_1*tiz)*Oi
             endif
          ENDDO !neigbor=1,patch_2D%num_edges i          
          edge_index_1 = p_patch%cells%edge_idx(cell_index, cell_block, 1)
          edge_block_1 = p_patch%cells%edge_blk(cell_index, cell_block, 1)
          edge_index_2 = p_patch%cells%edge_idx(cell_index, cell_block, 2)
          edge_block_2 = p_patch%cells%edge_blk(cell_index, cell_block, 2)
          edge_index_3 = p_patch%cells%edge_idx(cell_index, cell_block, 3)
          edge_block_3 = p_patch%cells%edge_blk(cell_index, cell_block, 3)
          ei=p_patch%edges%primal_edge_length(edge_index_1,edge_block_1)
          alpha=-1.0_wp/ei 
          Sx1 = S_x(edge_index_1,edge_block_1)*nx1 + S_y(edge_index_1,edge_block_1)*tx1
          Sx2 = S_x(edge_index_2,edge_block_2)*nx2 + S_y(edge_index_2,edge_block_2)*tx2
          Sx3 = S_x(edge_index_3,edge_block_3)*nx3 + S_y(edge_index_3,edge_block_3)*tx3
          Sy1 = S_x(edge_index_1,edge_block_1)*ny1 + S_y(edge_index_1,edge_block_1)*ty1
          Sy2 = S_x(edge_index_2,edge_block_2)*ny2 + S_y(edge_index_2,edge_block_2)*ty2
          Sy3 = S_x(edge_index_3,edge_block_3)*ny3 + S_y(edge_index_3,edge_block_3)*ty3
          z1= p_patch%edges%cartesian_center(edge_index_1,edge_block_1)%x(3)
          z2= p_patch%edges%cartesian_center(edge_index_2,edge_block_2)%x(3)
          z3= p_patch%edges%cartesian_center(edge_index_3,edge_block_3)%x(3)
          zeta_1=zeta_e(edge_index_1,edge_block_1)
          zeta_2=zeta_e(edge_index_2,edge_block_2)
          zeta_3=zeta_e(edge_index_3,edge_block_3)
          a1=1.0_wp
          a2=1.0_wp
          a3=1.0_wp
          if((-1.e8<z1 < 1.e-8) .and. (x3>0.0)) then
             a1=-1.0_wp
          endif
          
          if((-1.e8 <z2 < 1.e-8) .and. (x3>0.0)) then
             a2=-1.0_wp
          endif
          if((-1.e8 <z3 < 1.e-8) .and. (x3>0.0)) then
             a3=-1.0_wp
          endif
          S_n(edge_index_1,edge_block_1)=S_n(edge_index_1,edge_block_1)+ei/3.0_wp*alpha*&
               &a1*((Sx3*nx1+Sy3*ny1)*zeta_3 -(Sx2*nx1+Sy2*ny1)*zeta_2)
          S_t(edge_index_1,edge_block_1)=S_t(edge_index_1,edge_block_1)+ei/3.0_wp*alpha*&
               &a1*((Sx3*tx1+Sy3*ty1)*zeta_3 -(Sx2*tx1+Sy2*ty1)*zeta_2)
          S_n(edge_index_2,edge_block_2)=S_n(edge_index_2,edge_block_2)+ei/3.0_wp*alpha*&
               &a2*((Sx1*nx2+Sy1*ny2)*zeta_1 -(Sx3*nx2+Sy3*ny2)*zeta_3)
          S_t(edge_index_2,edge_block_2)=S_t(edge_index_2,edge_block_2)+ei/3.0_wp*alpha*&
               &a2*((Sx1*tx2+Sy1*ty2)*zeta_1 - (Sx3*tx2+Sy3*ty2)*zeta_3)
          S_n(edge_index_3,edge_block_3)=S_n(edge_index_3,edge_block_3)+ei/3.0_wp*alpha*&
               &a3*((Sx2*nx3+Sy2*ny3)*zeta_2 - (Sx1*nx3+Sy1*ny3)*zeta_1)
          S_t(edge_index_3,edge_block_3)=S_t(edge_index_3,edge_block_3)+ei/3.0_wp*alpha*&
               &a3*((Sx2*tx3+Sy2*ty3)*zeta_2 - (Sx1*tx3+Sy1*ty3)*zeta_1)
       endif
    ENDDO ! cell_index = start_index, end_index
 ENDDO !cell_block = owned_cells%start_block, owned_cells%end_block
END SUBROUTINE Stabilization

SUBROUTINE visu_vtk(dname,nummer,ice_x,ice_y,boundary_cell_marker,p_ice,p_patch_3D)    !initial value
  TYPE(t_patch_3D), TARGET, INTENT(IN)     :: p_patch_3D
  TYPE(t_patch),  POINTER                  :: p_patch
  TYPE(t_sea_ice),          INTENT(INOUT)  :: p_ice
  character *(*) dname
  INTEGER, INTENT(in) :: nummer
  REAL(wp), TARGET, INTENT(inout) :: ice_x(nproma,p_patch_3d%p_patch_2d(1)%nblks_e)
  REAL(wp), TARGET, INTENT(inout) :: ice_y(nproma,p_patch_3d%p_patch_2d(1)%nblks_e)
  REAL(wp), TARGET, INTENT(in) :: boundary_cell_marker(nproma,1,p_patch_3d%p_patch_2d(1)%nblks_e)
  TYPE(t_subset_range), POINTER :: all_cells
  TYPE(t_subset_range), POINTER :: all_edges
  TYPE(t_subset_range), POINTER :: verts_in_domain
  INTEGER :: cell_block,start_index,end_index,cell_index, neigbor, edge_index_i, edge_block_i,vert_block_1, vert_index_1,&
       &  edge_index_1, edge_block_1, edge_index_2, edge_block_2, edge_index_3, edge_block_3 
  REAL(wp) :: x,y,z
  character(len=8) :: fmt
  character(len=5) :: x1
  character(len=100) :: filename
  fmt = '(I5.5)'
  p_patch   => p_patch_3D%p_patch_2D(1)
  all_cells => p_patch%cells%all
  all_edges => p_patch%edges%all
  verts_in_domain   => p_patch%verts%in_domain

  write (x1,fmt) nummer 
  filename=trim(dname)//trim(x1)//'.dat'
  write(0,*) filename
  open (unit=20,file=trim(filename))
  write(20,*)"edges"

    DO cell_block = all_cells%start_block, all_cells%end_block
       CALL get_index_range(all_cells, cell_block, start_index, end_index)
       DO cell_index =  start_index, end_index
          DO neigbor=1, 3!no_primal_edges
             edge_index_i =p_patch%cells%edge_idx(cell_index, cell_block, neigbor)
             edge_block_i =p_patch%cells%edge_blk(cell_index, cell_block, neigbor)                                                                       
             x = p_patch%edges%cartesian_center(edge_index_i,edge_block_i)%x(1)
             y = p_patch%edges%cartesian_center(edge_index_i,edge_block_i)%x(2)
             z = p_patch%edges%cartesian_center(edge_index_i,edge_block_i)%x(3) 
             write(20,*)x,y,z
          ENDDO
       END DO
    END DO
    write(20,*)"edgevector"
    write(20,*)"velocity"
    DO cell_block = all_cells%start_block, all_cells%end_block
       CALL get_index_range(all_cells, cell_block, start_index, end_index)
       DO cell_index =  start_index, end_index
          DO neigbor=1, 3!no_primal_edges
             edge_index_i = p_patch%cells%edge_idx(cell_index, cell_block, neigbor)
             edge_block_i = p_patch%cells%edge_blk(cell_index, cell_block, neigbor)
             write(20,*) ice_x(edge_index_i,edge_block_i), ice_y(edge_index_i,edge_block_i)
          ENDDO
       END DO
    END DO

    write(20,*)"cellvalue"
    write(20,*)"Delta_c"
    DO cell_block = all_cells%start_block, all_cells%end_block
       CALL get_index_range(all_cells, cell_block, start_index, end_index)
       DO cell_index =  start_index, end_index
          write(20,*) p_ice%Delta(cell_index,cell_block)
       END DO
    END DO
 
    write(20,*)"cellvalue"
    write(20,*)"h"
    DO cell_block = all_cells%start_block, all_cells%end_block                                                                                                                                             
       CALL get_index_range(all_cells, cell_block, start_index, end_index)                                                                                                                                 
       DO cell_index =  start_index, end_index                                                                                                                                                             
             write(20,*) p_ice%hi(cell_index,1,cell_block)
       END DO
    END DO
    
    write(20,*)"cellvalue"
    write(20,*)"conc"
    DO cell_block = all_cells%start_block, all_cells%end_block
       CALL get_index_range(all_cells, cell_block, start_index, end_index)
       DO cell_index =  start_index, end_index
          write(20,*) p_ice%conc(cell_index,1,cell_block)
       END DO
    END DO
    close(unit=20)
  END SUBROUTINE visu_vtk
  
  

END MODULE mo_ice_test_2D