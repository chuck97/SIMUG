!©2021 MPI-M, Carolin Mehlmann                    
!Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
!1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
!2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
!3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
!THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYR
!IGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
! AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.[7]                                                 



  SUBROUTINE test_div_onTorus()
    ! Grid: OceanBasin_512x512_2km.nc

    Call Ice_drift_EVP !Example 4.2

  END SUBROUTINE test_div_onTorus

!============= main routine to compute Example 4.2=========================
 SUBROUTINE  Ice_Drift_EVP
    TYPE(t_subset_range), POINTER :: all_edges
    TYPE(t_subset_range), POINTER :: all_cells
    TYPE(t_subset_range), POINTER :: owned_cells


    INTEGER :: cell_block, start_index,  end_index, cell_index,&
         &edge_block_i,edge_block_j, neigbor, edge_index_i,edge_index_j, neigbor_j, cell_block1,&
         &cell_block2, cell_index1, cell_index2,edge_index_1,edge_index_2,edge_block_1,edge_block_2

    INTEGER outer_iter, inner_iter

    INTEGER vtkeverystep

    REAL(wp) ::  x, y, u, v, Zeitschritt,k,cell_area,x_0,y_0,P_0, Cda,Cdw

    REAL(wp) ::  L_ref,Oi, Oj, ei, e11, e22, e21, e12, P, Delta,eta,zeta,e, Delta_min,&
         &T,Pstar,Cor,L_x,L_y,tP,&
         &alpha_evp, beta_evp


    REAL(wp),POINTER             ::  Au_n(:,:,:), Au_t(:,:,:),&
         &atm_n(:,:,:), atm_t(:,:,:), u_n(:,:,:), u_t(:,:,:), u_x(:,:,:), u_y(:,:,:),&
         & mass(:,:,:),boundary_edge_marker(:,:,:),&
         &atm_u(:,:,:), atm_v(:,:,:),u_ocean_n(:,:,:), u_ocean_t(:,:,:),h_e(:,:,:),&
         &u_ocean_x(:,:,:), u_ocean_y(:,:,:),A_e(:,:,:),u_old_n(:,:,:), u_old_t(:,:,:),&
         &S_x(:,:,:),S_y(:,:,:), S_n(:,:,:),S_t(:,:,:), zeta_e(:,:,:)

    REAL(wp)::  s11(nproma,n_zlev,patch_3D%p_patch_2d(1)%alloc_cell_blocks),&
         &s22(nproma,n_zlev,patch_3D%p_patch_2d(1)%alloc_cell_blocks),&
         &s21(nproma,n_zlev,patch_3D%p_patch_2d(1)%alloc_cell_blocks),&
         &s12(nproma,n_zlev,patch_3D%p_patch_2d(1)%alloc_cell_blocks),&
         &sigma_I(nproma,n_zlev,patch_3D%p_patch_2d(1)%alloc_cell_blocks),&
         &sigma_II(nproma,n_zlev,patch_3D%p_patch_2d(1)%alloc_cell_blocks),&
         &boundary_cell_marker(nproma,n_zlev,patch_3D%p_patch_2d(1)%alloc_cell_blocks),&
         &A(nproma,n_zlev,patch_3D%p_patch_2d(1)%alloc_cell_blocks),&
         &H(nproma,n_zlev,patch_3D%p_patch_2d(1)%alloc_cell_blocks),&
         & P_e(nproma,n_zlev,patch_3D%p_patch_2d(1)%alloc_cell_blocks),&
         &zeta_c(nproma,n_zlev,patch_3D%p_patch_2d(1)%alloc_cell_blocks),&
         &Delta_c(nproma,n_zlev,patch_3D%p_patch_2d(1)%alloc_cell_blocks),&
         &shear_c(nproma,n_zlev,patch_3D%p_patch_2d(1)%alloc_cell_blocks)

    all_edges => patch_2d%edges%all
    all_cells => patch_2d%cells%all
    owned_cells => patch_2d%cells%owned
    verts_in_domain   => patch_2D%verts%in_domain

    !============= initalize edge vectors=========================
    ALLOCATE(mass(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e))
    ALLOCATE(atm_n(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e))
    ALLOCATE(atm_t(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e))
    ALLOCATE(u_n(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e))
    ALLOCATE(u_t(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e))
    ALLOCATE(u_old_n(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e))
    ALLOCATE(u_old_t(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e))
    ALLOCATE(u_ocean_n(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e))
    ALLOCATE(u_ocean_t(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e))
    ALLOCATE(u_ocean_x(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e))
    ALLOCATE(u_ocean_y(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e))
    ALLOCATE(u_x(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e))
    ALLOCATE(u_y(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e))
    ALLOCATE(atm_u(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e))
    ALLOCATE(atm_v(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e))
    ALLOCATE(Au_n(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e))
    ALLOCATE(Au_t(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e))
    ALLOCATE(S_x(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e))
    ALLOCATE(S_y(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e))
    ALLOCATE(S_n(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e))
    ALLOCATE(S_t(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e))
    ALLOCATE(h_e(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e))
    ALLOCATE(A_e(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e))
    ALLOCATE(zeta_e(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e))
    ALLOCATE(boundary_edge_marker(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e))

  !============= initalize cell vectors=========================
    Delta_c(nproma,1:n_zlev,patch_3D%p_patch_2d(1)%alloc_cell_blocks)=0.0_wp
    s11(nproma,1:n_zlev,patch_3D%p_patch_2d(1)%alloc_cell_blocks)=0.0_wp
    s22(nproma,1:n_zlev,patch_3D%p_patch_2d(1)%alloc_cell_blocks)=0.0_wp
    s21(nproma,1:n_zlev,patch_3D%p_patch_2d(1)%alloc_cell_blocks)=0.0_wp
    s12(nproma,1:n_zlev,patch_3D%p_patch_2d(1)%alloc_cell_blocks)=0.0_wp
    sigma_I(nproma,1:n_zlev,patch_3D%p_patch_2d(1)%alloc_cell_blocks)=0.0_wp
    sigma_II(nproma,1:n_zlev,patch_3D%p_patch_2d(1)%alloc_cell_blocks)=0.0_wp
    shear_c(nproma,1:n_zlev,patch_3D%p_patch_2d(1)%alloc_cell_blocks)=0.0_wp
  !============= initalize parameters=========================
    Cda=1.3_wp*.0012_wp
    Cdw=1026_wp*.0055_wp
    Cor=0.0
    k=0.0_wp! Initial time
    Delta_min=0.000000002_wp ! Hibler Delta_min
    e=2.0_wp! Hibler ellipse ratio
    L_ref=2.0
    L_x=511971.206305485_wp*L_ref
    L_y=512658.206708388_wp*L_ref
    Zeitschritt=300! timestep in seconds
    vtkeverystep = int(60*60*24/Zeitschritt)! set daily output
    T=0.0 ! count no. of cells
    alpha_evp=800!Parameter mEVP solver
    beta_evp=800!Parameter mEVP solver
    !============= initalize vectors=========================
    CALL init_boundary_cell_marker(boundary_cell_marker)! mark the cells with ice
    CALL init_mass_matrix(mass) ! Mass matrix for finite element approach  
    CALL initial_cell_values_EVP_drift(A,H,L_x) ! Inialize cell vectors
    CALL initial_edge_values_EVP_drift(u_n,u_t, u_x,u_y, u_ocean_x,u_ocean_y, boundary_edge_marker, L_x, L_y) ! Inizialize egde vectors

  DO cell_block = owned_cells%start_block, owned_cells%end_block
       CALL get_index_range(owned_cells, cell_block, start_index, end_index)
       DO cell_index = start_index, end_index
          s11(cell_index,1,cell_block)=0.0_wp
          s12(cell_index,1,cell_block)=0.0_wp
          s22(cell_index,1,cell_block)=0.0_wp
          s21(cell_index,1,cell_block)=0.0_wp
          T=T+1
          sigma_I(cell_index,1,cell_block)=0.0_wp
          sigma_II(cell_index,1,cell_block)=0.0_wp
       ENDDO ! cell_index = start_index, end_index
    ENDDO !cell_block = owned_cells%start_block, owned_cells%end_block
    write(0,*) T, "number of cell"

          Au_n=0.0_wp
          Au_t=0.0_wp
          !TIME_LOOP
          DO outer_iter=1,288 ! one day
             write(0,*) "======== Time step ===========", outer_iter, "time in days:", k/(60*60*24)
             k=k+Zeitschritt ! time in s
             Call wind_ocean_ice_drift(atm_n,atm_t,u_ocean_n,u_ocean_t,u_old_n,u_old_t,u_x,u_y,u_n,u_t,u_ocean_x,u_ocean_y,Cda,Cdw,L_x,L_y,k)! initialize egde components 
             Call H_A_interpolate_e(h_e,A_e,H,A,boundary_cell_marker, boundary_edge_marker)

       DO inner_iter=1,100! Looops of mEVP solver
          DO cell_block = owned_cells%start_block, owned_cells%end_block
             CALL get_index_range(owned_cells, cell_block, start_index, end_index)
             DO cell_index = start_index, end_index
                cell_area = patch_2D%cells%area(cell_index,cell_block)
                Call compute_strain_rate_tensor(cell_index,cell_block, cell_area, u_n,u_t, e11,e12,e21,e22,h_e,A_e,P_e)
                Pstar=27500._wp*exp(-20.0_wp*(1.0_wp-A(cell_index,1,cell_block)))/(L_ref*L_ref)
                P_0=Pstar*H(cell_index,1,cell_block)
                Delta=dsqrt(Delta_min*Delta_min + e12*e12 + 1.25_wp*(e11*e11+e22*e22)+1.5_wp*e11*e22)
                zeta=P_0/(2.0_wp*Delta)*1/(L_ref*L_ref)
                P=P_0/L_ref
                eta=0.25_wp*zeta
                if (boundary_cell_marker(cell_index,1,cell_block)==1) then
                   Delta_c(cell_index,1,cell_block)=Delta
                   shear_c(cell_index,1,cell_block)=dsqrt(Delta_min*Delta_min+(e11-e22)*(e11-e22)+(e12+e21)*(e12+e21))
                endif
                zeta_c(cell_index,1,cell_block)=0.0_wp
                Call  mEVP_drift(alpha_evp, e11,e12,e21,e22,zeta,P,P_0,sigma_I(cell_index,1,cell_block),sigma_II(cell_index,1,cell_block)&
                 &,s11(cell_index,1,cell_block),s12(cell_index,1,cell_block),&
                 &s21(cell_index,1,cell_block),s22(cell_index,1,cell_block),Delta,&
                zeta_c(cell_index,1,cell_block),boundary_cell_marker(cell_index,1,cell_block),A(cell_index,1,cell_block))
             ENDDO ! cell_index = start_index, end_index
          ENDDO !cell_block = owned_cells%start_block, owned_cells%end_block

          DO edge_block_i = all_edges%start_block, all_edges%end_block
              CALL get_index_range(all_edges, edge_block_i, start_index, end_index)
              DO edge_index_i =  start_index, end_index
                 cell_index1 = patch_2D%edges%cell_idx(edge_index_i,edge_block_i,1)
                 cell_block1 = patch_2D%edges%cell_blk(edge_index_i,edge_block_i,1)
                 cell_index2 = patch_2D%edges%cell_idx(edge_index_i,edge_block_i,2)
                 cell_block2 = patch_2D%edges%cell_blk(edge_index_i,edge_block_i,2)
                 zeta_e(edge_index_i,1,edge_block_i)=0.5_wp*(zeta_c(cell_index1,1,cell_block1)+zeta_c(cell_index2,1,cell_block2))
              END DO
           END DO

          Au_n=0.0_wp
          Au_t=0.0_wp


          call  compute_sigma(Au_n,Au_t,s11,s12,s22,s21)
          S_x=0.0_wp
          S_y=0.0_wp

          call Stabilization_sum(S_x,S_y,u_n,u_t)

          DO edge_block_i = all_edges%start_block, all_edges%end_block
             CALL get_index_range(all_edges, edge_block_i, start_index, end_index)
             DO edge_index_i =  start_index, end_index
                S_x(edge_index_i,1,edge_block_i)=S_x(edge_index_i,1,edge_block_i)*boundary_edge_marker(edge_index_i,1,edge_block_i)
                S_y(edge_index_i,1,edge_block_i)=S_y(edge_index_i,1,edge_block_i)*boundary_edge_marker(edge_index_i,1,edge_block_i)
             ENDDO
          ENDDO

          S_t=0.0_wp
          S_n=0.0_wp

          call Stabilization(S_x,S_y,S_n,S_t,zeta_c,zeta_e)
          call  solve_mEVP_drift(mass,h_e,A_e,boundary_edge_marker,atm_n,atm_t,&
             &u_old_n,u_old_t,u_ocean_x,u_ocean_y,S_t,S_n,Au_n,Au_t,u_x,u_y,u_n,u_t,beta_evp,Cdw,Zeitschritt,Cor)
       ENDDO ! inner iter mEVP

         if (modulo(outer_iter,vtkeverystep)==0) then
            CALL visu_vtk('Ice_drift_',int(outer_iter/vtkeverystep+0.00001_wp),A,H,Delta_c,shear_c,u_x,u_y, boundary_cell_marker)
         endif

          write(0,*) maxval(u_x)
       ENDDO !ENDTIME LOOP
     END SUBROUTINE Ice_drift_EVP



!============= visualize fields as vtk=========================  
   SUBROUTINE visu_vtk(dname,nummer,A,H,Delta_c,shear_c,u_x,u_y, boundary_cell_marker)    !initial value
    character *(*) dname

    INTEGER, INTENT(in) :: nummer

    REAL(wp),TARGET,INTENT(inout) :: Delta_c(nproma,n_zlev,patch_3D%p_patch_2d(1)%alloc_cell_blocks)
     REAL(wp),TARGET,INTENT(inout) :: shear_c(nproma,n_zlev,patch_3D%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp),TARGET,INTENT(inout) :: A(nproma,n_zlev,patch_3D%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp),TARGET,INTENT(inout) :: H(nproma,n_zlev,patch_3D%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp), TARGET, INTENT(inout) :: u_x(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e)
    REAL(wp), TARGET, INTENT(inout) :: u_y(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e)
    REAL(wp),TARGET,INTENT(in) :: boundary_cell_marker(nproma,n_zlev,patch_3D%p_patch_2d(1)%alloc_cell_blocks)

    TYPE(t_subset_range), POINTER :: all_cells
    TYPE(t_subset_range), POINTER :: all_edges
    INTEGER :: cell_block,start_index,end_index,cell_index, neigbor, edge_index_i, edge_block_i,vert_block_1, vert_index_1 
    REAL(wp) :: x,y,x_0,y_0

    character(len=8) :: fmt
    character(len=5) :: x1

    character(len=100) :: filename

    fmt = '(I5.5)'

    all_cells => patch_2d%cells%all
    all_edges => patch_2d%edges%all

    verts_in_domain   => patch_2D%verts%in_domain

    write (x1,fmt) nummer 
    filename=trim(dname)//trim(x1)//'.dat'
    write(0,*) filename

    open (unit=20,file=trim(filename))

    write(20,*)"edges"
    DO cell_block = all_cells%start_block, all_cells%end_block                                                                                                                                             
       CALL get_index_range(all_cells, cell_block, start_index, end_index)                                                                                                                                 
       DO cell_index =  start_index, end_index                                                                                                                                                             
          if (boundary_cell_marker(cell_index,1,cell_block)==1.0) then


             DO neigbor=1, 3!no_primal_edges                                                                                                                                                                 
                edge_index_i = patch_2D%cells%edge_idx(cell_index, cell_block, neigbor)                                                                                                                       
                edge_block_i = patch_2D%cells%edge_blk(cell_index, cell_block, neigbor)                                                                                                                       

                vert_index_1 = patch_2d %cells%vertex_idx(cell_index, cell_block,neigbor)                                                                                                                     
                vert_block_1 = patch_2d %cells%vertex_blk(cell_index, cell_block,neigbor)                                                                                                                     

                x_0 = patch_2d%edges%cartesian_center(edge_index_i,edge_block_i)%x(1)

                y_0 = patch_2d%edges%cartesian_center(edge_index_i,edge_block_i)%x(2)

                write(20,*)cell_block,x_0,y_0
             ENDDO

          endif

       END DO
    END DO

    write(20,*)"edgevector"
    write(20,*)"velocity"
    DO cell_block = all_cells%start_block, all_cells%end_block                                                                                                                                             
       CALL get_index_range(all_cells, cell_block, start_index, end_index)                                                                                                                                 
       DO cell_index =  start_index, end_index                                                                                                                                                             
          if (boundary_cell_marker(cell_index,1,cell_block)==1.0) then



             DO neigbor=1, 3!no_primal_edges                                                                                                                                                                 
                edge_index_i = patch_2D%cells%edge_idx(cell_index, cell_block, neigbor)                                                                                                                       
                edge_block_i = patch_2D%cells%edge_blk(cell_index, cell_block, neigbor)                                                                                                                       

                vert_index_1 = patch_2d %cells%vertex_idx(cell_index, cell_block,neigbor)                                                                                                                     
                vert_block_1 = patch_2d %cells%vertex_blk(cell_index, cell_block,neigbor)                                                                                                                     
                write(20,*) u_x(edge_index_i,1,edge_block_i), u_y(edge_index_i,1,edge_block_i)
             ENDDO
          endif

       END DO
    END DO

    write(20,*)"cellvalue"
    write(20,*)"h"
    DO cell_block = all_cells%start_block, all_cells%end_block                                                                                                                                             
       CALL get_index_range(all_cells, cell_block, start_index, end_index)                                                                                                                                 
       DO cell_index =  start_index, end_index                                                                                                                                                             
          if (boundary_cell_marker(cell_index,1,cell_block)==1.0) then
             write(20,*)H(cell_index,1,cell_block)
          endif
       END DO
    END DO

    write(20,*)"cellvalue"
    write(20,*)"a"
    DO cell_block = all_cells%start_block, all_cells%end_block                                                                                                                                             
       CALL get_index_range(all_cells, cell_block, start_index, end_index)                                                                                                                                 
       DO cell_index =  start_index, end_index                                                                                                                                                             
          if (boundary_cell_marker(cell_index,1,cell_block)==1.0) then
             write(20,*)A(cell_index,1,cell_block)
          endif
       END DO
    END DO

    write(20,*)"cellvalue"
    write(20,*)"Delta_c"
    DO cell_block = all_cells%start_block, all_cells%end_block
       CALL get_index_range(all_cells, cell_block, start_index, end_index)
       DO cell_index =  start_index, end_index
          if (boundary_cell_marker(cell_index,1,cell_block)==1.0) then
             write(20,*)Delta_c(cell_index,1,cell_block)
          endif
       END DO
    END DO


    write(20,*)"cellvalue"
    write(20,*)"shear_c"
    DO cell_block = all_cells%start_block, all_cells%end_block
       CALL get_index_range(all_cells, cell_block, start_index, end_index)
       DO cell_index =  start_index, end_index
          if (boundary_cell_marker(cell_index,1,cell_block)==1.0) then
             write(20,*)shear_c(cell_index,1,cell_block)
          endif
       END DO
    END DO


    
    close(unit=20)

  END SUBROUTINE visu_vtk
  
  SUBROUTINE init_boundary_cell_marker(boundary_cell_marker)
    REAL(wp),TARGET,INTENT(inout) :: boundary_cell_marker(nproma,n_zlev,patch_3D%p_patch_2d(1)%alloc_cell_blocks)

    TYPE(t_subset_range), POINTER :: all_cells
    INTEGER :: cell_block,start_index,end_index,cell_index

    all_cells => patch_2d%cells%all

    DO cell_block = all_cells%start_block, all_cells%end_block
       CALL get_index_range(all_cells, cell_block, start_index, end_index) 
       DO cell_index =  start_index, end_index 
          boundary_cell_marker(cell_index,1,cell_block)=patch_3d%p_patch_1d(1)%dolic_c(cell_index,cell_block)*0.5_wp! cell w ice=1.0_wp cell without ice=0.0_wp
       END DO
    END DO
  END SUBROUTINE init_boundary_cell_marker

   SUBROUTINE init_mass_matrix(mass)
    REAL(wp), TARGET, INTENT(inout) :: mass(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e)

    TYPE(t_subset_range), POINTER :: all_cells
    TYPE(t_subset_range), POINTER :: owned_cells

    INTEGER :: cell_block,start_index,end_index,cell_index, edge_index_i,edge_block_i, neigbor

    all_cells => patch_2d%cells%all
    owned_cells => patch_2d%cells%owned

    DO cell_block = owned_cells%start_block, owned_cells%end_block
       CALL get_index_range(owned_cells, cell_block, start_index, end_index)
       DO cell_index = start_index, end_index
          DO neigbor=1, patch_2D%cells%num_edges(cell_index,cell_block)!no_primal_edges

             edge_index_i = patch_2D%cells%edge_idx(cell_index, cell_block, neigbor)
             edge_block_i = patch_2D%cells%edge_blk(cell_index, cell_block, neigbor)

             mass(edge_index_i,1,edge_block_i)=mass(edge_index_i,1,edge_block_i)&
                  &+patch_2D%cells%area(cell_index,cell_block)/3.0_wp
          ENDDO !neigbor=1,patch_2D%num_edges i                                                                                                         \

       ENDDO ! cell_index = start_index, end_index       

    ENDDO !cell_block = owned_cells%start_block, owned_cells%end_block    
  END SUBROUTINE init_mass_matrix

   SUBROUTINE initial_cell_values_EVP_drift(A,H,L_x)                                                                                                                                                                                                                                                                 
    REAL(wp),TARGET,INTENT(inout) :: A(nproma,n_zlev,patch_3D%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp),TARGET,INTENT(inout) :: H(nproma,n_zlev,patch_3D%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp), INTENT(in) :: L_x

    TYPE(t_subset_range), POINTER :: all_cells
    INTEGER :: cell_block,start_index,end_index,cell_index
    REAL(wp) :: x,y,x_0,y_0

    all_cells => patch_2d%cells%all

    DO cell_block = all_cells%start_block, all_cells%end_block
       CALL get_index_range(all_cells, cell_block, start_index, end_index)
       DO cell_index =  start_index, end_index
          x = patch_2d%cells%cartesian_center(cell_index,cell_block)%x(1)
          y = patch_2d%cells%cartesian_center(cell_index,cell_block)%x(2)
          x_0=(x+253985.715628112_wp)*2
          A(cell_index,1,cell_block)=(x_0/L_x)
          H(cell_index,1,cell_block)=2.0                                                                                                                                                                                                                                                                                  
       END DO
    END DO
  END SUBROUTINE initial_cell_values_EVP_drift

  SUBROUTINE initial_edge_values_EVP_drift(u_n,u_t, u_x,u_y, u_ocean_x,u_ocean_y, boundary_edge_marker, L_x, L_y)
    REAL(wp), TARGET, INTENT(inout) :: u_n(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e)
    REAL(wp), TARGET, INTENT(inout) :: u_t(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e)
    REAL(wp), TARGET, INTENT(inout) :: u_x(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e)
    REAL(wp), TARGET, INTENT(inout) :: u_y(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e)
    REAL(wp), TARGET, INTENT(inout) :: u_ocean_x(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e)
    REAL(wp), TARGET, INTENT(inout) :: u_ocean_y(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e)
    REAL(wp), TARGET, INTENT(inout) :: boundary_edge_marker(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e)
    REAL(wp), INTENT(in) :: L_x
    REAL(wp), INTENT(in) :: L_y
    TYPE(t_subset_range), POINTER :: all_edges
    
    INTEGER :: edge_block_i,start_index,end_index,edge_index_i
    REAL(wp) :: x,y ,x_0,y_0, nx,ny,tx,ty

    all_edges => patch_2d%edges%all

    DO edge_block_i = all_edges%start_block, all_edges%end_block
       CALL get_index_range(all_edges, edge_block_i, start_index, end_index)
       DO edge_index_i =  start_index, end_index
          x = patch_2d%edges%cartesian_center(edge_index_i,edge_block_i)%x(1)
          y = patch_2d%edges%cartesian_center(edge_index_i,edge_block_i)%x(2)
          nx=patch_2d%edges%primal_normal(edge_index_i,edge_block_i)%v1
          ny=patch_2d%edges%primal_normal(edge_index_i,edge_block_i)%v2
          tx=patch_2d%edges%dual_normal(edge_index_i,edge_block_i)%v1
          ty=patch_2d%edges%dual_normal(edge_index_i,edge_block_i)%v2
          
           x_0=(x+253985.715628112_wp)*2
           y_0=(y+256329.103354194_wp)*2
       
          boundary_edge_marker(edge_index_i,1,edge_block_i)= patch_3d%p_patch_1d(1)%dolic_e(edge_index_i,edge_block_i)*0.5_wp ! inizialize edges with ice
          u_x(edge_index_i,1,edge_block_i)=0.0_wp!initial sea ice velocity x component                                                                                                                                                                                                                                  
          u_y(edge_index_i,1,edge_block_i)=0.0_wp!initial sea ice velocity y component                                                                                                                                                                                         
          u_n(edge_index_i,1,edge_block_i)= u_x(edge_index_i,1,edge_block_i)*nx+u_y(edge_index_i,1,edge_block_i)*ny
          u_t(edge_index_i,1,edge_block_i)= u_x(edge_index_i,1,edge_block_i)*tx+u_y(edge_index_i,1,edge_block_i)*ty
          u_ocean_x(edge_index_i,1,edge_block_i)= 0.1_wp*(2.0*y_0-L_y)/L_y!!initial ocean x component
          u_ocean_y(edge_index_i,1,edge_block_i)=-0.1_wp*(2.0*x_0-L_x)/L_x!!initial ocean y component

       END DO
    END DO

  END SUBROUTINE initial_edge_values_EVP_drift


   SUBROUTINE wind_ocean_ice_drift(atm_n,atm_t,u_ocean_n,u_ocean_t,u_old_n,u_old_t,u_x,u_y,u_n,u_t,u_ocean_x,u_ocean_y,Cda,Cdw,L_x,L_y,k)
    REAL(wp), TARGET, INTENT(inout) :: atm_n(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e)
    REAL(wp), TARGET, INTENT(inout) :: atm_t(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e)
    REAL(wp), TARGET, INTENT(inout) :: u_ocean_n(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e)
    REAL(wp), TARGET, INTENT(inout) :: u_ocean_t(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e)
    REAL(wp), TARGET, INTENT(inout) :: u_old_n(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e)
    REAL(wp), TARGET, INTENT(inout) :: u_old_t(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e)
    REAL(wp), TARGET, INTENT(in) :: u_ocean_x(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e)
    REAL(wp), TARGET, INTENT(in) :: u_ocean_y(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e)
    REAL(wp), TARGET, INTENT(in) :: u_x(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e)
    REAL(wp), TARGET, INTENT(in) :: u_y(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e)
    REAL(wp), TARGET, INTENT(in) :: u_n(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e)
    REAL(wp), TARGET, INTENT(in) :: u_t(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e)
    REAL(wp), INTENT(in) :: Cda
    REAL(wp), INTENT(in) :: Cdw
    REAL(wp), INTENT(in) :: L_x
    REAL(wp), INTENT(in) :: L_y
    REAL(wp), INTENT(in) :: k
    
    TYPE(t_subset_range), POINTER :: all_edges
    INTEGER :: edge_block_i,start_index,end_index,edge_index_i

    REAL(wp) :: x,y ,x_0,y_0, diff_x, diff_y, vwd, nx, ny,tx,ty,U_xx,U_yy, atm_u,atm_v, tP,vmax,ws,alpha,mx,my,wx,rw,s,wy

    all_edges => patch_2d%edges%all

    DO edge_block_i = all_edges%start_block, all_edges%end_block
       CALL get_index_range(all_edges, edge_block_i, start_index, end_index)
       DO edge_index_i =  start_index, end_index
          x = patch_2d%edges%cartesian_center(edge_index_i,edge_block_i)%x(1)
          y = patch_2d%edges%cartesian_center(edge_index_i,edge_block_i)%x(2)
          nx=patch_2d%edges%primal_normal(edge_index_i,edge_block_i)%v1
          ny=patch_2d%edges%primal_normal(edge_index_i,edge_block_i)%v2
          tx=patch_2d%edges%dual_normal(edge_index_i,edge_block_i)%v1
          ty=patch_2d%edges%dual_normal(edge_index_i,edge_block_i)%v2

          ! Update velcoity
          u_old_n(edge_index_i,1,edge_block_i)=u_n(edge_index_i,1,edge_block_i)
          u_old_t(edge_index_i,1,edge_block_i)=u_t(edge_index_i,1,edge_block_i)
          
          !athmosphere stress
          tP=k/(3600.0_wp*24.0_wp) ! time in days
          x_0=(x+253985.715628112_wp)*2
          y_0=(y+256329.103354194_wp)*2
         
          U_yy=(5.0_wp+(sin(2.0_wp*3.141*tp/4.0)-3)*sin(2.0_wp*3.141*y_0/L_y)*sin(3.141*x_0/L_x))
          U_xx=(5.0_wp+(sin(2.0_wp*3.141*tp/4.0)-3)*sin(2.0_wp*3.141*x_0/L_x)*sin(3.141*y_0/L_y))
          atm_u=U_xx*dsqrt(U_xx*U_xx+U_yy*U_yy)*Cda
          atm_v=U_yy*dsqrt(U_xx*U_xx+U_yy*U_yy)*Cda
          atm_n(edge_index_i,1,edge_block_i)=atm_u*nx+atm_v*ny
          atm_t(edge_index_i,1,edge_block_i)=atm_u*tx+atm_v*ty

          !Ocean stress
          diff_x=u_x(edge_index_i,1,edge_block_i)-u_ocean_x(edge_index_i,1,edge_block_i)
          diff_y=u_y(edge_index_i,1,edge_block_i)-u_ocean_y(edge_index_i,1,edge_block_i)
          vwd=dsqrt(diff_x*diff_x+diff_y*diff_y)
          diff_x=-u_ocean_x(edge_index_i,1,edge_block_i)
          diff_y=-u_ocean_y(edge_index_i,1,edge_block_i)
          u_ocean_n(edge_index_i,1,edge_block_i)=Cdw*vwd*(diff_x*nx&
               &+diff_y*ny)
          u_ocean_t(edge_index_i,1,edge_block_i)=Cdw*vwd*(diff_x*tx&
               &+diff_y*ty)
       ENDDO
    ENDDO
  END SUBROUTINE wind_ocean_ice_drift



  

  SUBROUTINE H_A_interpolate_e(h_e,A_e,H,A,boundary_cell_marker, boundary_edge_marker)
    REAL(wp), TARGET, INTENT(out) :: h_e(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e)
    REAL(wp), TARGET, INTENT(out) :: A_e(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e)
    REAL(wp),TARGET,INTENT(in) :: boundary_cell_marker(nproma,n_zlev,patch_3D%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp), TARGET, INTENT(in) :: boundary_edge_marker(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e)
    REAL(wp),TARGET,INTENT(in) :: H(nproma,n_zlev,patch_3D%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp),TARGET,INTENT(in) :: A(nproma,n_zlev,patch_3D%p_patch_2d(1)%alloc_cell_blocks)
    TYPE(t_subset_range), POINTER :: all_edges
    INTEGER :: edge_block_i,start_index,end_index,edge_index_i, cell_index1, cell_index2, cell_block1, cell_block2, &
         &vert_block, vert_index

    REAL(wp) :: x,y ,x_0,y_0, diff_x, diff_y, vwd, nx, ny,tx,ty,U_xx,U_yy, atm_u,atm_v, factor

    all_edges => patch_2d%edges%all

    DO edge_block_i = all_edges%start_block, all_edges%end_block
       CALL get_index_range(all_edges, edge_block_i, start_index, end_index)
       DO edge_index_i =  start_index, end_index
          cell_index1 = patch_2D%edges%cell_idx(edge_index_i,edge_block_i,1)
          cell_block1 = patch_2D%edges%cell_blk(edge_index_i,edge_block_i,1)
          cell_index2 = patch_2D%edges%cell_idx(edge_index_i,edge_block_i,2)
          cell_block2 = patch_2D%edges%cell_blk(edge_index_i,edge_block_i,2)
          A_e(edge_index_i,1,edge_block_i)=0.0_wp
          h_e(edge_index_i,1,edge_block_i)=0.0_wp 

          factor=boundary_cell_marker(cell_index1,1,cell_block1)+boundary_cell_marker(cell_index2,1,cell_block2)

          if(boundary_edge_marker(edge_index_i,1,edge_block_i)>0.0_wp)then
             A_e(edge_index_i,1,edge_block_i)=(A(cell_index1,1,cell_block1)*boundary_cell_marker(cell_index1,1,cell_block1)&
                  +A(cell_index2,1,cell_block2)*boundary_cell_marker(cell_index2,1,cell_block2))/factor
             h_e(edge_index_i,1,edge_block_i)=(H(cell_index1,1,cell_block1)*boundary_cell_marker(cell_index1,1,cell_block1)&
                  +H(cell_index2,1,cell_block2)*boundary_cell_marker(cell_index2,1,cell_block2))/factor                                
          end if
       ENDDO
    ENDDO

  END SUBROUTINE H_A_interpolate_e
  
 SUBROUTINE compute_strain_rate_tensor(cell_index,cell_block, cell_area, u_n,u_t, e11,e12,e21,e22,h_e,A_e,P_e)
    INTEGER, INTENT(in) :: cell_index
    INTEGER, INTENT(in) :: cell_block
    REAL(wp), INTENT(in) :: cell_area
    REAL(wp), TARGET, INTENT(in) :: u_n(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e) 
    REAL(wp), TARGET, INTENT(in) :: u_t(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e)
    REAL(wp), TARGET, INTENT(in) :: h_e(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e) 
    REAL(wp), TARGET, INTENT(in) :: A_e(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e)
    REAL(wp), INTENT(out) :: e11
    REAL(wp), INTENT(out) :: e12
    REAL(wp), INTENT(out) :: e21
    REAL(wp), INTENT(out) :: e22

    REAL(wp),TARGET,INTENT(inout) :: P_e(nproma,n_zlev,patch_3D%p_patch_2d(1)%alloc_cell_blocks)



    INTEGER :: neigbor, edge_block_i,edge_index_i, max_nbgh

    REAL(wp) :: ei,h_ei,Pi, nix,niy, tix,tiy, Oi, Pstar

    P_e(cell_index,1,cell_block)=0.0_wp
    e11=0.0_wp
    e22=0.0_wp
    e12=0.0_wp
    e21=0.0_wp

    max_nbgh=patch_2D%cells%num_edges(cell_index,cell_block)

    DO neigbor=1,max_nbgh !no_primal_edges                                                                      
       edge_index_i = patch_2D%cells%edge_idx(cell_index, cell_block, neigbor)
       edge_block_i = patch_2D%cells%edge_blk(cell_index, cell_block, neigbor)

       ei=patch_2D%edges%primal_edge_length(edge_index_i,edge_block_i)

       h_ei = 2.0_wp*cell_area/ei
       Pi=2.0_wp/h_ei

       nix=patch_2d%edges%primal_normal(edge_index_i,edge_block_i)%v1
       niy=patch_2d%edges%primal_normal(edge_index_i,edge_block_i)%v2
       tix=patch_2d%edges%dual_normal(edge_index_i,edge_block_i)%v1
       tiy=patch_2d%edges%dual_normal(edge_index_i,edge_block_i)%v2



       Oi=patch_2D%cells%edge_orientation(cell_index,cell_block,neigbor)

       e11=e11+1.0_wp*nix*nix*Pi*u_n(edge_index_i,1,edge_block_i)*Oi+&
            &1.0_wp*tix*nix*Pi*u_t(edge_index_i,1,edge_block_i)*Oi

       e22= e22+1.0_wp*niy*niy*Pi*u_n(edge_index_i,1,edge_block_i)*Oi+&
            &1.0_wp*tiy*niy*Pi*u_t(edge_index_i,1,edge_block_i)*Oi

       e12=e12+1.0_wp*nix*niy*Pi*u_n(edge_index_i,1,edge_block_i)*Oi+&
            &(0.5_wp*tix*niy+0.5_wp*nix*tiy)*Pi*u_t(edge_index_i,1,edge_block_i)*Oi!

       e21=e21+1.0_wp*nix*niy*Pi*u_n(edge_index_i,1,edge_block_i)*Oi+&
            &(0.5_wp*tix*niy+0.5_wp*nix*tiy)*Pi*u_t(edge_index_i,1,edge_block_i)*Oi

    ENDDO !neigbor=1,patch_2D%num_edges i


  END SUBROUTINE compute_strain_rate_tensor
  
  SUBROUTINE mEVP_drift( alpha_evp, e11,e12,e21,e22,zeta,P,P_0,sigma_I,sigma_II,s11,s12,s21,s22,Delta,zeta_c,marker_c,A)
    REAL(wp), INTENT(in) :: e11
    REAL(wp), INTENT(in) :: e12
    REAL(wp), INTENT(in) :: e21
    REAL(wp), INTENT(in) :: e22
    REAL(wp), INTENT(in) :: zeta
    REAL(wp), INTENT(in) :: P,P_0
    REAL(wp), INTENT(in) :: alpha_evp
     REAL(wp), INTENT(in) :: marker_c
     REAL(wp), INTENT(in) :: A
     REAL(wp), INTENT(inout) :: sigma_I
     REAL(wp), INTENT(inout) :: sigma_II
     REAL(wp), INTENT(inout) :: s12
     REAL(wp), INTENT(inout) :: s21
     REAL(wp), INTENT(inout) :: s11
     REAL(wp), INTENT(inout) :: s22
     REAL(wp), INTENT(inout) :: Delta
     REAL(wp), INTENT(inout) :: zeta_c
     REAL(wp) :: B_evp, C_evp
     TYPE(t_subset_range), POINTER :: all_edges
     INTEGER :: edge_block_i,start_index,end_index,edge_index_i

    
     B_evp=(alpha_evp-1.0_wp)/alpha_evp                                                                                                                   
     C_evp=1.0_wp/alpha_evp                                                                                                                                           

     sigma_I= B_evp*sigma_I+C_evp*&
          &(2.0_wp*zeta*(e11+e22)-P)

     sigma_II=B_evp*sigma_II+C_evp*&
          &0.5_wp*zeta*(e11-e22)

     s12=B_evp*s12+C_evp*0.5_wp*zeta*e12
     s21=B_evp*s21+C_evp*0.5_wp*zeta*e21
     
     s11=0.5_wp*(sigma_I+sigma_II)
     s22=0.5_wp*(sigma_I-sigma_II)

     if (marker_c==1.0) then
        zeta_c=P_0/Delta
     endif
   END SUBROUTINE mEVP_drift


    SUBROUTINE solve_mEVP_drift(mass,h_e,A_e,boundary_edge_marker,atm_n,atm_t,&
                   &u_old_n,u_old_t,u_ocean_x,u_ocean_y,S_t,S_n,Au_n,Au_t,u_x,u_y,u_n,u_t,beta_evp,Cdw,Zeitschritt,Cor)
    REAL(wp), TARGET, INTENT(in) :: mass(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e)
    REAL(wp), TARGET, INTENT(in) :: h_e(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e)            
    REAL(wp), TARGET, INTENT(in) :: A_e(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e)    
    REAL(wp), TARGET, INTENT(in) :: boundary_edge_marker(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e)
    REAL(wp), TARGET, INTENT(in) :: atm_n(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e)
    REAL(wp), TARGET, INTENT(in) :: atm_t(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e)
    REAL(wp), TARGET, INTENT(in) :: u_old_n(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e)
    REAL(wp), TARGET, INTENT(in) :: u_old_t(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e)
    REAL(wp), TARGET, INTENT(in) :: u_ocean_x(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e)
    REAL(wp), TARGET, INTENT(in) :: u_ocean_y(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e)
    REAL(wp), TARGET, INTENT(in) :: S_t(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e)
    REAL(wp), TARGET, INTENT(in) :: S_n(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e)
    REAL(wp), TARGET, INTENT(in) :: Au_t(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e)
    REAL(wp), TARGET, INTENT(in) :: Au_n(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e)
    REAL(wp), TARGET, INTENT(inout) :: u_x(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e)
    REAL(wp), TARGET, INTENT(inout) :: u_y(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e)
    REAL(wp), TARGET, INTENT(inout) :: u_n(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e)
    REAL(wp), TARGET, INTENT(inout) :: u_t(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e)
    REAL(wp) :: u_ocean_n(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e)
    REAL(wp) :: u_ocean_t(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e)
    REAL(wp), INTENT(in) :: beta_evp
    REAL(wp), INTENT(in) :: Cdw
    REAL(wp), INTENT(in) :: Zeitschritt
    REAL(wp), INTENT(in) :: Cor
    INTEGER :: edge_block_i,start_index,end_index,edge_index_i

    REAL(wp) :: x,y, diff_x, diff_y, vwd, nx, ny,tx,ty,mass_ice,C_imp

    TYPE(t_subset_range), POINTER :: all_edges
    
    all_edges => patch_2d%edges%all

    DO edge_block_i = all_edges%start_block, all_edges%end_block
       CALL get_index_range(all_edges, edge_block_i, start_index, end_index)
       DO edge_index_i =  start_index, end_index
          x = patch_2d%edges%cartesian_center(edge_index_i,edge_block_i)%x(1)
          y = patch_2d%edges%cartesian_center(edge_index_i,edge_block_i)%x(2)
          nx=patch_2d%edges%primal_normal(edge_index_i,edge_block_i)%v1
          ny=patch_2d%edges%primal_normal(edge_index_i,edge_block_i)%v2
          tx=patch_2d%edges%dual_normal(edge_index_i,edge_block_i)%v1
          ty=patch_2d%edges%dual_normal(edge_index_i,edge_block_i)%v2
                if(boundary_edge_marker(edge_index_i,1,edge_block_i)>0.0_wp) then
                   mass_ice=900.0_wp*h_e(edge_index_i,1,edge_block_i)
                   diff_x=u_x(edge_index_i,1,edge_block_i)-u_ocean_x(edge_index_i,1,edge_block_i)
                   diff_y=u_y(edge_index_i,1,edge_block_i)-u_ocean_y(edge_index_i,1,edge_block_i)
                   vwd=dsqrt(diff_x*diff_x+diff_y*diff_y)
                   
                   u_ocean_n(edge_index_i,1,edge_block_i)=&
                   &Cdw*vwd*(u_ocean_x(edge_index_i,1,edge_block_i)*nx&
                   &+u_ocean_y(edge_index_i,1,edge_block_i)*ny)
                   
                   u_ocean_t(edge_index_i,1,edge_block_i)=&
                    &Cdw*vwd*(u_ocean_x(edge_index_i,1,edge_block_i)*tx&
                    &+u_ocean_y(edge_index_i,1,edge_block_i)*ty)

                   C_imp=1.0_wp/(beta_evp+1.0_wp+vwd*Cdw*Zeitschritt/mass_ice)

                   u_n(edge_index_i,1,edge_block_i)=C_imp*u_old_n(edge_index_i,1,edge_block_i)+&
                        &C_imp*beta_evp*u_n(edge_index_i,1,edge_block_i)+&
                        &Zeitschritt/mass_ice*C_imp*&
                        &(atm_n(edge_index_i,1,edge_block_i)+u_ocean_n(edge_index_i,1,edge_block_i)&
                        &-Au_n(edge_index_i,1,edge_block_i)/mass(edge_index_i,1,edge_block_i))&
                        &+C_imp*Zeitschritt*Cor*(diff_y*nx-diff_x*ny)&
                        &+C_imp*Zeitschritt/mass_ice*S_n(edge_index_i,1,edge_block_i)/mass(edge_index_i,1,edge_block_i)

                   u_t(edge_index_i,1,edge_block_i)= C_imp*u_old_t(edge_index_i,1,edge_block_i)+&
                        &C_imp*beta_evp*u_t(edge_index_i,1,edge_block_i)+&
                        &C_imp*Zeitschritt/mass_ice&
                        &(atm_t(edge_index_i,1,edge_block_i)+u_ocean_t(edge_index_i,1,edge_block_i)&
                        -Au_t(edge_index_i,1,edge_block_i)/mass(edge_index_i,1,edge_block_i))&
                        &+C_imp*Zeitschritt*Cor*(diff_y*tx-diff_x*ty)&
                       &+C_imp*Zeitschritt/mass_ice*S_t(edge_index_i,1,edge_block_i)/mass(edge_index_i,1,edge_block_i)

                   if ( h_e(edge_index_i,1,edge_block_i)==0.0_wp) then
                      write(*,*) 'assert message h_e',h_e(edge_index_i,1,edge_block_i)
                      call exit(1)
                   endif
                end if

                u_n(edge_index_i,1,edge_block_i)=u_n(edge_index_i,1,edge_block_i)*boundary_edge_marker(edge_index_i,1,edge_block_i)
                u_t(edge_index_i,1,edge_block_i)=u_t(edge_index_i,1,edge_block_i)*boundary_edge_marker(edge_index_i,1,edge_block_i)
                u_x(edge_index_i,1,edge_block_i)=nx*u_n(edge_index_i,1,edge_block_i)+tx*u_t(edge_index_i,1,edge_block_i)
                u_y(edge_index_i,1,edge_block_i)=ny*u_n(edge_index_i,1,edge_block_i)+ty*u_t(edge_index_i,1,edge_block_i)
             ENDDO
          ENDDO

        END SUBROUTINE solve_mEVP_drift


       SUBROUTINE compute_sigma(Au_n,Au_t,s11,s12,s22,s21)
         REAL(wp), TARGET, INTENT(out) ::  Au_n(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e)
         REAL(wp), TARGET, INTENT(out) ::  Au_t(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e)
         REAL(wp),TARGET,INTENT(in) :: s11(nproma,n_zlev,patch_3D%p_patch_2d(1)%alloc_cell_blocks)
         REAL(wp),TARGET,INTENT(in) :: s12(nproma,n_zlev,patch_3D%p_patch_2d(1)%alloc_cell_blocks)
         REAL(wp),TARGET,INTENT(in) :: s22(nproma,n_zlev,patch_3D%p_patch_2d(1)%alloc_cell_blocks)
         REAL(wp),TARGET,INTENT(in) :: s21(nproma,n_zlev,patch_3D%p_patch_2d(1)%alloc_cell_blocks)
         REAL(wp) ::nix,niy,tix,tiy,Oi,cell_area, Pi,ei,h_ei
         INTEGER :: cell_block, start_index,  end_index, cell_index,&
              &edge_index_i,edge_block_i, neigbor
         TYPE(t_subset_range), POINTER :: owned_cells
         TYPE(t_subset_range), POINTER :: all_edges
         owned_cells => patch_2d%cells%owned
         all_edges => patch_2d%edges%all
       DO cell_block = owned_cells%start_block, owned_cells%end_block
             CALL get_index_range(owned_cells, cell_block, start_index, end_index)
             DO cell_index = start_index, end_index
                cell_area = patch_2D%cells%area(cell_index,cell_block)
                DO neigbor=1, patch_2D%cells%num_edges(cell_index,cell_block)!no_primal_edges                                                                                                               
                   edge_index_i = patch_2D%cells%edge_idx(cell_index, cell_block, neigbor)
                   edge_block_i = patch_2D%cells%edge_blk(cell_index, cell_block, neigbor)
                   ei=patch_2D%edges%primal_edge_length(edge_index_i,edge_block_i)
                   h_ei = 2.0_wp*cell_area/patch_2D%edges%primal_edge_length(edge_index_i,edge_block_i)
                   Pi=2.0_wp/h_ei
                   nix=patch_2d%edges%primal_normal(edge_index_i,edge_block_i)%v1
                   niy=patch_2d%edges%primal_normal(edge_index_i,edge_block_i)%v2
                   tix=patch_2d%edges%dual_normal(edge_index_i,edge_block_i)%v1
                   tiy=patch_2d%edges%dual_normal(edge_index_i,edge_block_i)%v2
                   Oi=patch_2D%cells%edge_orientation(cell_index,cell_block,neigbor)

                   Au_n(edge_index_i,1,edge_block_i)=Au_n(edge_index_i,1,edge_block_i)+&
                        &Pi*Oi*cell_area*(&
                        &s11(cell_index,1,cell_block)*nix*nix+s12(cell_index,1,cell_block)*nix*niy&
                        &+s21(cell_index,1,cell_block)*niy*nix+s22(cell_index,1,cell_block)*niy*niy)

                   Au_t(edge_index_i,1,edge_block_i)=Au_t(edge_index_i,1,edge_block_i)+&
                        &Pi*Oi*cell_area*(&
                        &0.5_wp*(s11(cell_index,1,cell_block)*tix*nix+s12(cell_index,1,cell_block)*tix*niy&
                        &+s21(cell_index,1,cell_block)*tiy*nix+s22(cell_index,1,cell_block)*tiy*niy)+&
                        &0.5_wp*(s11(cell_index,1,cell_block)*tix*nix+s12(cell_index,1,cell_block)*tiy*nix&
                        &+s21(cell_index,1,cell_block)*tix*niy+s22(cell_index,1,cell_block)*tiy*niy))

                ENDDO !neigbor=1,patch_2D%num_edges i                                                                                                                                                       
             ENDDO ! cell_index = start_index, end_index                                                                                                                                                    
          ENDDO !cell_block = owned_cells%start_block, owned_cells%end_block
          
        END SUBROUTINE compute_sigma
        
   SUBROUTINE Stabilization_sum(S_x,S_y,u_n,u_t)

    REAL(wp), TARGET, INTENT(inout) ::  S_x(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e)
    REAL(wp), TARGET, INTENT(inout) ::  S_y(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e)
    REAL(wp), TARGET, INTENT(in) ::  u_n(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e)
    REAL(wp), TARGET, INTENT(in) ::  u_t(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e)

    REAL(wp) ::nx1,ny1,nx2,ny2,nx3,ny3,tx1,ty1,tx2,ty2,tx3,ty3

    INTEGER :: cell_block, start_index,  end_index, cell_index,&
         &edge_index_1,edge_index_2,edge_index_3,edge_block_1,edge_block_2,edge_block_3

    TYPE(t_subset_range), POINTER :: owned_cells
    TYPE(t_subset_range), POINTER :: all_edges

    owned_cells => patch_2d%cells%owned
    all_edges => patch_2d%edges%all



    DO cell_block = owned_cells%start_block, owned_cells%end_block
       CALL get_index_range(owned_cells, cell_block, start_index, end_index)
       DO cell_index = start_index, end_index


          edge_index_1 = patch_2D%cells%edge_idx(cell_index, cell_block, 1)
          edge_block_1 = patch_2D%cells%edge_blk(cell_index, cell_block, 1)

          edge_index_2 = patch_2D%cells%edge_idx(cell_index, cell_block, 2)
          edge_block_2 = patch_2D%cells%edge_blk(cell_index, cell_block, 2)

          edge_index_3 = patch_2D%cells%edge_idx(cell_index, cell_block, 3)
          edge_block_3 = patch_2D%cells%edge_blk(cell_index, cell_block, 3)

          nx1=patch_2d%edges%primal_normal(edge_index_1,edge_block_1)%v1
          ny1=patch_2d%edges%primal_normal(edge_index_1,edge_block_1)%v2
          tx1=patch_2d%edges%dual_normal(edge_index_1,edge_block_1)%v1
          ty1=patch_2d%edges%dual_normal(edge_index_1,edge_block_1)%v2

          nx2=patch_2d%edges%primal_normal(edge_index_2,edge_block_2)%v1
          ny2=patch_2d%edges%primal_normal(edge_index_2,edge_block_2)%v2
          tx2=patch_2d%edges%dual_normal(edge_index_2,edge_block_2)%v1
          ty2=patch_2d%edges%dual_normal(edge_index_2,edge_block_2)%v2


          nx3=patch_2d%edges%primal_normal(edge_index_3,edge_block_3)%v1
          ny3=patch_2d%edges%primal_normal(edge_index_3,edge_block_3)%v2
          tx3=patch_2d%edges%dual_normal(edge_index_3,edge_block_3)%v1
          ty3=patch_2d%edges%dual_normal(edge_index_3,edge_block_3)%v2


          S_x(edge_index_1,1,edge_block_1)=S_x(edge_index_1,1,edge_block_1)+&
               &u_n(edge_index_2,1,edge_block_2)*(nx2)-u_n(edge_index_3,1,edge_block_3)*(nx3)+&
               &u_t(edge_index_2,1,edge_block_2)*(tx2)-u_t(edge_index_3,1,edge_block_3)*(tx3)

          S_x(edge_index_2,1,edge_block_2)=S_x(edge_index_2,1,edge_block_2)+&
               &u_n(edge_index_3,1,edge_block_3)*(nx3)-u_n(edge_index_1,1,edge_block_1)*(nx1)+&
               &u_t(edge_index_3,1,edge_block_3)*(tx3)-u_t(edge_index_1,1,edge_block_1)*(tx1)

          S_x(edge_index_3,1,edge_block_3)=S_x(edge_index_3,1,edge_block_3)+&
               &u_n(edge_index_1,1,edge_block_1)*(nx1)-u_n(edge_index_2,1,edge_block_2)*(nx2)+&
               &u_t(edge_index_1,1,edge_block_1)*(tx1)-u_t(edge_index_2,1,edge_block_2)*(tx2) 

          S_y(edge_index_1,1,edge_block_1)=S_y(edge_index_1,1,edge_block_1)+&
               &u_n(edge_index_2,1,edge_block_2)*(ny2)-u_n(edge_index_3,1,edge_block_3)*(ny3)+&
               &u_t(edge_index_2,1,edge_block_2)*(ty2)-u_t(edge_index_3,1,edge_block_3)*(ty3)

          S_y(edge_index_2,1,edge_block_2)=S_y(edge_index_2,1,edge_block_2)+&
               &u_n(edge_index_3,1,edge_block_3)*(ny3)-u_n(edge_index_1,1,edge_block_1)*(ny1)+&
               &u_t(edge_index_3,1,edge_block_3)*(ty3)-u_t(edge_index_1,1,edge_block_1)*(ty1)

          S_y(edge_index_3,1,edge_block_3)=S_y(edge_index_3,1,edge_block_3)+&
               &u_n(edge_index_1,1,edge_block_1)*(ny1)-u_n(edge_index_2,1,edge_block_2)*(ny2)+&
               &u_t(edge_index_1,1,edge_block_1)*(ty1)-u_t(edge_index_2,1,edge_block_2)*(ty2)                            





       ENDDO ! cell_index = start_index, end_index
    ENDDO !cell_block = owned_cells%start_block, owned_cells%end_block

  END SUBROUTINE Stabilization_sum



  SUBROUTINE Stabilization(S_x,S_y,S_n,S_t,zeta_c,zeta_e)

    REAL(wp), TARGET, INTENT(inout) ::  S_n(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e)
    REAL(wp), TARGET, INTENT(inout) ::  S_t(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e)
    REAL(wp), TARGET, INTENT(in) ::  S_x(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e)
    REAL(wp), TARGET, INTENT(in) ::  S_y(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e)
     REAL(wp), TARGET, INTENT(in) ::  zeta_e(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e)

    REAL(wp),TARGET,INTENT(in) :: zeta_c(nproma,n_zlev,patch_3D%p_patch_2d(1)%alloc_cell_blocks)

    REAL(wp) ::nx1,ny1,nx2,ny2,nx3,ny3,tx1,ty1,tx2,ty2,tx3,ty3,ei,alpha

    INTEGER :: cell_block, start_index,  end_index, cell_index,&
         &edge_index_1,edge_index_2,edge_index_3,edge_block_1,edge_block_2,edge_block_3

    TYPE(t_subset_range), POINTER :: owned_cells
    TYPE(t_subset_range), POINTER :: all_edges

    owned_cells => patch_2d%cells%owned
    all_edges => patch_2d%edges%all



    DO cell_block = owned_cells%start_block, owned_cells%end_block
       CALL get_index_range(owned_cells, cell_block, start_index, end_index)
       DO cell_index = start_index, end_index



          edge_index_1 = patch_2D%cells%edge_idx(cell_index, cell_block, 1)
          edge_block_1 = patch_2D%cells%edge_blk(cell_index, cell_block, 1)

          edge_index_2 = patch_2D%cells%edge_idx(cell_index, cell_block, 2)
          edge_block_2 = patch_2D%cells%edge_blk(cell_index, cell_block, 2)

          edge_index_3 = patch_2D%cells%edge_idx(cell_index, cell_block, 3)
          edge_block_3 = patch_2D%cells%edge_blk(cell_index, cell_block, 3)

          nx1=patch_2d%edges%primal_normal(edge_index_1,edge_block_1)%v1
          ny1=patch_2d%edges%primal_normal(edge_index_1,edge_block_1)%v2
          tx1=patch_2d%edges%dual_normal(edge_index_1,edge_block_1)%v1
          ty1=patch_2d%edges%dual_normal(edge_index_1,edge_block_1)%v2

          nx2=patch_2d%edges%primal_normal(edge_index_2,edge_block_2)%v1
          ny2=patch_2d%edges%primal_normal(edge_index_2,edge_block_2)%v2
          tx2=patch_2d%edges%dual_normal(edge_index_2,edge_block_2)%v1
          ty2=patch_2d%edges%dual_normal(edge_index_2,edge_block_2)%v2


          nx3=patch_2d%edges%primal_normal(edge_index_3,edge_block_3)%v1
          ny3=patch_2d%edges%primal_normal(edge_index_3,edge_block_3)%v2
          tx3=patch_2d%edges%dual_normal(edge_index_3,edge_block_3)%v1
          ty3=patch_2d%edges%dual_normal(edge_index_3,edge_block_3)%v2

          ei=patch_2D%edges%primal_edge_length(edge_index_1,edge_block_1)

          alpha=-1.0/ei

          S_n(edge_index_1,1,edge_block_1)=S_n(edge_index_1,1,edge_block_1)+ei/3.0_wp*alpha*&
               &(zeta_e(edge_index_3,1,edge_block_3)*S_x(edge_index_3,1,edge_block_3)*nx1+&
               &zeta_e(edge_index_3,1,edge_block_3)*S_y(edge_index_3,1,edge_block_3)*ny1&
               &-zeta_e(edge_index_2,1,edge_block_2)*S_x(edge_index_2,1,edge_block_2)*nx1&
               &-zeta_e(edge_index_2,1,edge_block_2)*S_y(edge_index_2,1,edge_block_2)*ny1)
          
          S_t(edge_index_1,1,edge_block_1)=S_t(edge_index_1,1,edge_block_1)+ei/3.0_wp*alpha*&
               &(zeta_e(edge_index_3,1,edge_block_3)*S_x(edge_index_3,1,edge_block_3)*tx1+&
               &zeta_e(edge_index_3,1,edge_block_3)*S_y(edge_index_3,1,edge_block_3)*ty1&
               &-zeta_e(edge_index_2,1,edge_block_2)*S_x(edge_index_2,1,edge_block_2)*tx1&
               &-zeta_e(edge_index_2,1,edge_block_2)*S_y(edge_index_2,1,edge_block_2)*ty1)

          S_n(edge_index_2,1,edge_block_2)=S_n(edge_index_2,1,edge_block_2)+ei/3.0_wp*alpha*&
               &(zeta_e(edge_index_1,1,edge_block_1)*S_x(edge_index_1,1,edge_block_1)*nx2+&
               &zeta_e(edge_index_1,1,edge_block_1)*S_y(edge_index_1,1,edge_block_1)*ny2&
               &-zeta_e(edge_index_3,1,edge_block_3)*S_x(edge_index_3,1,edge_block_3)*nx2&
               &-zeta_e(edge_index_3,1,edge_block_3)*S_y(edge_index_3,1,edge_block_3)*ny2)

          S_t(edge_index_2,1,edge_block_2)=S_t(edge_index_2,1,edge_block_2)+ei/3.0_wp*alpha*&
               &(zeta_e(edge_index_1,1,edge_block_1)*S_x(edge_index_1,1,edge_block_1)*tx2+&
               &zeta_e(edge_index_1,1,edge_block_1)*S_y(edge_index_1,1,edge_block_1)*ty2&
               &-zeta_e(edge_index_3,1,edge_block_3)*S_x(edge_index_3,1,edge_block_3)*tx2&
               &-zeta_e(edge_index_3,1,edge_block_3)*S_y(edge_index_3,1,edge_block_3)*ty2)                          


          S_n(edge_index_3,1,edge_block_3)=S_n(edge_index_3,1,edge_block_3)+ei/3.0_wp*alpha*&
               &(zeta_e(edge_index_2,1,edge_block_2)*S_x(edge_index_2,1,edge_block_2)*nx3+&
               &zeta_e(edge_index_2,1,edge_block_2)*S_y(edge_index_2,1,edge_block_2)*ny3&
               &-zeta_e(edge_index_1,1,edge_block_1)*S_x(edge_index_1,1,edge_block_1)*nx3&
               &-zeta_e(edge_index_1,1,edge_block_1)*S_y(edge_index_1,1,edge_block_1)*ny3)

          S_t(edge_index_3,1,edge_block_3)=S_t(edge_index_3,1,edge_block_3)+ei/3.0_wp*alpha*&
               &(zeta_e(edge_index_2,1,edge_block_2)*S_x(edge_index_2,1,edge_block_2)*tx3+&
               &zeta_e(edge_index_2,1,edge_block_2)*S_y(edge_index_2,1,edge_block_2)*ty3&
               &-zeta_e(edge_index_1,1,edge_block_1)*S_x(edge_index_1,1,edge_block_1)*tx3&
               &-zeta_e(edge_index_1,1,edge_block_1)*S_y(edge_index_1,1,edge_block_1)*ty3)                         



       ENDDO ! cell_index = start_index, end_index
    ENDDO !cell_block = owned_cells%start_block, owned_cells%end_block


  END SUBROUTINE Stabilization




!============= main routine to compute Example 4.3=========================

 SUBROUTINE  Ice_Drift_EVP

   
    TYPE(t_subset_range), POINTER :: all_edges
    TYPE(t_subset_range), POINTER :: all_cells
    TYPE(t_subset_range), POINTER :: owned_cells


    INTEGER :: cell_block, start_index,  end_index, cell_index,&
         &edge_block_i,edge_block_j, neigbor, edge_index_i,edge_index_j, neigbor_j, cell_block1,&
         &cell_block2, cell_index1, cell_index2,edge_index_1,edge_index_2,edge_block_1,edge_block_2

    INTEGER outer_iter, inner_iter

    INTEGER vtkeverystep

    REAL(wp) ::  x, y, u, v, Zeitschritt,k,cell_area,x_0,y_0,P_0, Cda,Cdw

    REAL(wp) ::  L_ref,Oi, Oj, ei, e11, e22, e21, e12, P, Delta,eta,zeta,e, Delta_min,&
         &T,Pstar,Cor,L_x,L_y,tP,&
         &alpha_evp, beta_evp


    REAL(wp),POINTER             ::  Au_n(:,:,:), Au_t(:,:,:),&
         &atm_n(:,:,:), atm_t(:,:,:), u_n(:,:,:), u_t(:,:,:), u_x(:,:,:), u_y(:,:,:),&
         & mass(:,:,:),boundary_edge_marker(:,:,:),&
         &atm_u(:,:,:), atm_v(:,:,:),u_ocean_n(:,:,:), u_ocean_t(:,:,:),h_e(:,:,:),&
         &u_ocean_x(:,:,:), u_ocean_y(:,:,:),A_e(:,:,:),u_old_n(:,:,:), u_old_t(:,:,:),&
         &S_x(:,:,:),S_y(:,:,:), S_n(:,:,:),S_t(:,:,:), zeta_e(:,:,:)

    REAL(wp)::  s11(nproma,n_zlev,patch_3D%p_patch_2d(1)%alloc_cell_blocks),&
         &s22(nproma,n_zlev,patch_3D%p_patch_2d(1)%alloc_cell_blocks),&
         &s21(nproma,n_zlev,patch_3D%p_patch_2d(1)%alloc_cell_blocks),&
         &s12(nproma,n_zlev,patch_3D%p_patch_2d(1)%alloc_cell_blocks),&
         &sigma_I(nproma,n_zlev,patch_3D%p_patch_2d(1)%alloc_cell_blocks),&
         &sigma_II(nproma,n_zlev,patch_3D%p_patch_2d(1)%alloc_cell_blocks),&
         &boundary_cell_marker(nproma,n_zlev,patch_3D%p_patch_2d(1)%alloc_cell_blocks),&
         &A(nproma,n_zlev,patch_3D%p_patch_2d(1)%alloc_cell_blocks),&
         &H(nproma,n_zlev,patch_3D%p_patch_2d(1)%alloc_cell_blocks),&
         & P_e(nproma,n_zlev,patch_3D%p_patch_2d(1)%alloc_cell_blocks),&
         &zeta_c(nproma,n_zlev,patch_3D%p_patch_2d(1)%alloc_cell_blocks),&
         &Delta_c(nproma,n_zlev,patch_3D%p_patch_2d(1)%alloc_cell_blocks),&
         &shear_c(nproma,n_zlev,patch_3D%p_patch_2d(1)%alloc_cell_blocks)


    !-----------------------------------------------------------------------
    all_edges => patch_2d%edges%all
    all_cells => patch_2d%cells%all
    owned_cells => patch_2d%cells%owned
    verts_in_domain   => patch_2D%verts%in_domain

    !============= initalize edge vectors========================= 

    ALLOCATE(mass(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e))


    ALLOCATE(atm_n(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e))
    ALLOCATE(atm_t(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e))

    ALLOCATE(u_n(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e))
    ALLOCATE(u_t(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e))
    ALLOCATE(u_old_n(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e))
    ALLOCATE(u_old_t(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e))


    ALLOCATE(u_ocean_n(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e))
    ALLOCATE(u_ocean_t(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e))

    ALLOCATE(u_ocean_x(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e))
    ALLOCATE(u_ocean_y(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e))

    ALLOCATE(u_x(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e))
    ALLOCATE(u_y(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e))

    ALLOCATE(atm_u(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e))
    ALLOCATE(atm_v(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e))

    ALLOCATE(Au_n(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e))
    ALLOCATE(Au_t(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e))



    ALLOCATE(S_x(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e))
    ALLOCATE(S_y(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e))

    ALLOCATE(S_n(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e))
    ALLOCATE(S_t(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e))


    ALLOCATE(h_e(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e))
    ALLOCATE(A_e(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e))
    ALLOCATE(zeta_e(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e))
    ALLOCATE(boundary_edge_marker(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e))

  !============= initalize cell vectors========================= 

    Delta_c(nproma,1:n_zlev,patch_3D%p_patch_2d(1)%alloc_cell_blocks)=0.0_wp
    s11(nproma,1:n_zlev,patch_3D%p_patch_2d(1)%alloc_cell_blocks)=0.0_wp
    s22(nproma,1:n_zlev,patch_3D%p_patch_2d(1)%alloc_cell_blocks)=0.0_wp
    s21(nproma,1:n_zlev,patch_3D%p_patch_2d(1)%alloc_cell_blocks)=0.0_wp
    s12(nproma,1:n_zlev,patch_3D%p_patch_2d(1)%alloc_cell_blocks)=0.0_wp
    sigma_I(nproma,1:n_zlev,patch_3D%p_patch_2d(1)%alloc_cell_blocks)=0.0_wp
    sigma_II(nproma,1:n_zlev,patch_3D%p_patch_2d(1)%alloc_cell_blocks)=0.0_wp
    shear_c(nproma,1:n_zlev,patch_3D%p_patch_2d(1)%alloc_cell_blocks)=0.0_wp


  !============= initalize parameters========================= 
    Cda=1.3_wp*.0012_wp
    Cdw=1026_wp*.0055_wp
    Cor=0.0

    k=0.0_wp! Initial time
    Delta_min=0.000000002_wp ! Hibler Delta_min
    e=2.0_wp! Hibler ellipse ratio
    
    L_ref=2.0
    L_x=511971.206305485_wp*L_ref
    L_y=512658.206708388_wp*L_ref
      
    Zeitschritt=300! timestep in seconds
    
    vtkeverystep = int(60*60*24/Zeitschritt)! set daily output                                                                                                                                                                                                                                                                        
    T=0.0 ! count no. of cells
    alpha_evp=800!Parameter mEVP solver
    beta_evp=800!Parameter mEVP solver

    !============= initalize vectors=========================

    CALL init_boundary_cell_marker(boundary_cell_marker)! mark the cells with ice
    CALL init_mass_matrix(mass) ! Mass matrix for finite element approach
    CALL initial_cell_values_EVP_drift(A,H,L_x) ! Inialize cell vectors
    CALL initial_edge_values_EVP_drift(u_n,u_t, u_x,u_y, u_ocean_x,u_ocean_y, boundary_edge_marker, L_x, L_y) ! Inizialize egde vectors


  DO cell_block = owned_cells%start_block, owned_cells%end_block
       CALL get_index_range(owned_cells, cell_block, start_index, end_index)
       DO cell_index = start_index, end_index
          s11(cell_index,1,cell_block)=0.0_wp
          s12(cell_index,1,cell_block)=0.0_wp
          s22(cell_index,1,cell_block)=0.0_wp
          s21(cell_index,1,cell_block)=0.0_wp
          T=T+1
          sigma_I(cell_index,1,cell_block)=0.0_wp
          sigma_II(cell_index,1,cell_block)=0.0_wp
       ENDDO ! cell_index = start_index, end_index
    ENDDO !cell_block = owned_cells%start_block, owned_cells%end_block
    write(0,*) T, "number of cell"
    
          Au_n=0.0_wp
          Au_t=0.0_wp


          !TIME_LOOP
          DO outer_iter=1,288 ! one day 
             write(0,*) "======== Time step ===========", outer_iter, "time in days:", k/(60*60*24)
   
             k=k+Zeitschritt ! time in s
             Call wind_ocean_ice_drift(atm_n,atm_t,u_ocean_n,u_ocean_t,u_old_n,u_old_t,u_x,u_y,u_n,u_t,u_ocean_x,u_ocean_y,Cda,Cdw,L_x,L_y,k)! initialize egde components
             Call H_A_interpolate_e(h_e,A_e,H,A,boundary_cell_marker, boundary_edge_marker)
                                                                                                                                                                                                                                                                                                          
       DO inner_iter=1,100! Looops of mEVP solver                                                                                                                                                                                                                                         
          DO cell_block = owned_cells%start_block, owned_cells%end_block
             CALL get_index_range(owned_cells, cell_block, start_index, end_index)
             DO cell_index = start_index, end_index
                cell_area = patch_2D%cells%area(cell_index,cell_block)

                Call compute_strain_rate_tensor(cell_index,cell_block, cell_area, u_n,u_t, e11,e12,e21,e22,h_e,A_e,P_e)
      
                Pstar=27500._wp*exp(-20.0_wp*(1.0_wp-A(cell_index,1,cell_block)))/(L_ref*L_ref)
                P_0=Pstar*H(cell_index,1,cell_block)
                Delta=dsqrt(Delta_min*Delta_min + e12*e12 + 1.25_wp*(e11*e11+e22*e22)+1.5_wp*e11*e22)
                zeta=P_0/(2.0_wp*Delta)*1/(L_ref*L_ref)
                P=P_0/L_ref
                eta=0.25_wp*zeta

                if (boundary_cell_marker(cell_index,1,cell_block)==1) then
                   Delta_c(cell_index,1,cell_block)=Delta
                   shear_c(cell_index,1,cell_block)=dsqrt(Delta_min*Delta_min+(e11-e22)*(e11-e22)+(e12+e21)*(e12+e21))
                endif

                zeta_c(cell_index,1,cell_block)=0.0_wp
             
                
                Call  mEVP_drift(alpha_evp, e11,e12,e21,e22,zeta,P,P_0,sigma_I(cell_index,1,cell_block),sigma_II(cell_index,1,cell_block)&
                 &,s11(cell_index,1,cell_block),s12(cell_index,1,cell_block),&
                 &s21(cell_index,1,cell_block),s22(cell_index,1,cell_block),Delta,&
                zeta_c(cell_index,1,cell_block),boundary_cell_marker(cell_index,1,cell_block),A(cell_index,1,cell_block))


             ENDDO ! cell_index = start_index, end_index
          ENDDO !cell_block = owned_cells%start_block, owned_cells%end_block

          DO edge_block_i = all_edges%start_block, all_edges%end_block
              CALL get_index_range(all_edges, edge_block_i, start_index, end_index)
              DO edge_index_i =  start_index, end_index
                 cell_index1 = patch_2D%edges%cell_idx(edge_index_i,edge_block_i,1)
                 cell_block1 = patch_2D%edges%cell_blk(edge_index_i,edge_block_i,1)
                 cell_index2 = patch_2D%edges%cell_idx(edge_index_i,edge_block_i,2)
                 cell_block2 = patch_2D%edges%cell_blk(edge_index_i,edge_block_i,2)
                 zeta_e(edge_index_i,1,edge_block_i)=0.5_wp*(zeta_c(cell_index1,1,cell_block1)+zeta_c(cell_index2,1,cell_block2))
              END DO
           END DO

          Au_n=0.0_wp
          Au_t=0.0_wp


          call  compute_sigma(Au_n,Au_t,s11,s12,s22,s21)
          S_x=0.0_wp
          S_y=0.0_wp

          call Stabilization_sum(S_x,S_y,u_n,u_t)

          DO edge_block_i = all_edges%start_block, all_edges%end_block
             CALL get_index_range(all_edges, edge_block_i, start_index, end_index)
             DO edge_index_i =  start_index, end_index
                S_x(edge_index_i,1,edge_block_i)=S_x(edge_index_i,1,edge_block_i)*boundary_edge_marker(edge_index_i,1,edge_block_i)
                S_y(edge_index_i,1,edge_block_i)=S_y(edge_index_i,1,edge_block_i)*boundary_edge_marker(edge_index_i,1,edge_block_i)
             ENDDO
          ENDDO

          S_t=0.0_wp
          S_n=0.0_wp

          call Stabilization(S_x,S_y,S_n,S_t,zeta_c,zeta_e)
          call  solve_mEVP_drift(mass,h_e,A_e,boundary_edge_marker,atm_n,atm_t,&
             &u_old_n,u_old_t,u_ocean_x,u_ocean_y,S_t,S_n,Au_n,Au_t,u_x,u_y,u_n,u_t,beta_evp,Cdw,Zeitschritt,Cor)
       ENDDO ! inner iter mEVP                                                                                                                                                                                                                                                                                                                                                                                                                                                                

         if (modulo(outer_iter,vtkeverystep)==0) then
            CALL visu_vtk('Ice_drift_',int(outer_iter/vtkeverystep+0.00001_wp),A,H,Delta_c,shear_c,u_x,u_y, boundary_cell_marker)
         endif

          write(0,*) maxval(u_x)

       ENDDO !ENDTIME LOOP

     END SUBROUTINE Ice_drift_EVP
