! CODES FOR CM:


module aux_module_CM

  use redistribute2D
  implicit none

contains


    subroutine cm_update_fluc(fm,fm2,gm,gm2,fp,fp2,gp,gp2,gradQ,gradQ2,ii,jj,N_cells,&
      qold,qold2,qnew,qnew2,aux,auxu,type_supper,type_sunder,dist_to_ecen_up,&
      dist_to_ecen_down,order,un_area_ij,up_area_ij,area_supper,area_sunder,&
      intersections,wall_height,dt,dx,dy,mx,my,ecen_un_x,ecen_un_y,&
      ecen_up_x,ecen_up_y,cen_grid_up,cen_grid_down,lengths_sunder,lengths_supper,&
      mthlim)

      ! This subroutine updates qnew and qnew2 such that at the small cells the cell merging
      ! method is applied, using qold and qold2.
      ! Order can be 1 (where no reconstruction at cuts are made) or 2 (with reconstruction)

      implicit none
      integer :: N_cells,type_supper(:),type_sunder(:),mx,my,ilat3,mthlim(3)
      integer :: order, i,j,i_0,j_0,ii(:),jj(:),ixy,itop,ibot,ilat,ilat2,is
      real(8) fp(3,-1:mx+2,-1:my+2)
      real(8) fm(3,-1:mx+2,-1:my+2)
      real(8) gp(3,-1:mx+2,-1:my+2)
      real(8) gm(3,-1:mx+2,-1:my+2)
      real(8) fp2(3,-1:mx+2,-1:my+2)
      real(8) fm2(3,-1:mx+2,-1:my+2)
      real(8) gp2(3,-1:mx+2,-1:my+2)
      real(8) gm2(3,-1:mx+2,-1:my+2)
      real(8) :: gradQ(3,2,-1:mx+2,-1:my+2),gradQ2(3,2,-1:mx+2,-1:my+2)
      real(8) :: qold(3,-1:mx+2,-1:my+2), qold2(3,-1:mx+2,-1:my+2)
      real(8) :: qnew(3,-1:mx+2,-1:my+2), qnew2(3,-1:mx+2,-1:my+2)
      real(8) :: aux(1,-1:mx+2,-1:my+2),auxu(1,-1:mx+2,-1:my+2)
      real(8) :: cen_grid_up(2,-1:mx+2,-1:my+2), cen_grid_down(2,-1:mx+2,-1:my+2)
      real(8) :: ecen_un_x(:,:,:),ecen_un_y(:,:,:)
      real(8) :: ecen_up_x(:,:,:),ecen_up_y(:,:,:)
      real(8) :: dist_to_ecen_up(:,:,:),dist_to_ecen_down(:,:,:)
      real(8) :: intersections(:,:),n_vec(2),t_vec(2),n_vec2(2),x,y
      real(8) :: dt,dx,dy,lengths_sunder(:,:),lengths_supper(:,:)
      real(8) :: cm_avg2(3),q1(3),q2(3),q3(3),q4(3),q5(3),q6(3)
      real(8) :: q1l(3),q2l(3),q3l(3),q4l(3),q5l(3),q6l(3),cm_avg(3)
      real(8) :: aux_un,aux_up,wall_height,temp_dist(2)
      real(8) :: un_area_ij(-1:mx+2,-1:my+2),up_area_ij(-1:mx+2,-1:my+2),apdq_l(3,3)
      real(8) :: amdq(3,3),apdq(3,3),amdq_wall(3),apdq_wall(3),amdq_l(3,3)
      real(8) :: amdq_r(3,3),apdq_r(3,3),s(3,3),fwave(3,3,3),area_sunder(:),area_supper(:)
      real(8) :: amdq_c(3),apdq_c(3)
        do i=1,N_cells
          i_0 = ii(i)
          j_0 = jj(i)
          x = intersections(1,i+1)- intersections(1,i)
          y = intersections(2,i+1) - intersections(2,i)


          ! correct via CM method the UPPER small cells and affected cells first. (there are some repeated calculations this way, but not much)
          select case(type_supper(i))
          case(4)
            ! get the old merged average  & update the old merged average:
            cm_avg2 = (1/(area_supper(i)+up_area_ij(i_0,j_0+1)))* &
        (area_supper(i)*qold2(:,i_0,j_0)+up_area_ij(i_0,j_0+1)*qold2(:,i_0,j_0+1))

             if (order .eq. 2) then
             ! get the fluctuations at each cell edge:
              ! get the cell edge value via gradient:
              ! barrier edge at indexed small cell
               q1 = qold2(:,i_0,j_0) + gradQ2(:,1,i_0,j_0)*dist_to_ecen_up(1,1,i) &
                 + gradQ2(:,2,i_0,j_0)*dist_to_ecen_up(2,1,i)   ! at the barrier center
              aux_up = auxu(1,i_0,j_0)

               q1l = qold(:,i_0,j_0) + gradQ(:,1,i_0,j_0)*dist_to_ecen_down(1,1,i) &
                 + gradQ(:,2,i_0,j_0)*dist_to_ecen_down(2,1,i)
               aux_un = aux(1,i_0,j_0)
             else
               q1 = qold2(:,i_0,j_0)
               aux_up = auxu(1,i_0,j_0)
               q1l = qold(:,i_0,j_0)
               aux_un = aux(1,i_0,j_0)
             end if
             call rotated_redist_fwave(q1,q1l,aux_up,aux_un, &
                  x,y,amdq_wall,apdq_wall,2,wall_height)
             cm_avg2 = cm_avg2 - dt/dx *lengths_supper(1,i)*(1/(area_supper(i)+&
               up_area_ij(i_0,j_0+1))) *amdq_wall

              ! at the edge 2:
              if (order .eq. 2) then
                call riemann_solver(qold2(:,i_0-2:i_0,j_0),qold2(:,i_0-1:i_0+1,j_0),&
                   auxu(1,i_0-2:i_0,j_0),auxu(1,i_0-1:i_0+1,j_0),s,fwave,&
                    amdq,apdq,1,2,dt/dx,mthlim(3))
                apdq_c = apdq(:,2)
              else
                apdq_c = fp2(:,i_0,j_0)
              end if

              cm_avg2 = cm_avg2 - dt/dx * lengths_supper(3,i)*(1/(area_supper(i)+&
                  up_area_ij(i_0,j_0+1)))* apdq_c

              ! at the merging neighbor cell edges:
              ! at the barrier edge and moving clockwise
              x = intersections(1,i+2)- intersections(1,i+1)
              y = intersections(2,i+2) - intersections(2,i+1)
            if (order .eq. 2) then
              q1 = qold2(:,i_0,j_0+1) + gradQ2(:,1,i_0,j_0+1)*dist_to_ecen_up(1,1,i+1)&
                 + gradQ2(:,2,i_0,j_0+1)*dist_to_ecen_up(2,1,i+1)
              aux_up = auxu(1,i_0,j_0+1)
              q1l = qold(:,i_0,j_0+1) + gradQ(:,1,i_0,j_0+1)*dist_to_ecen_down(1,1,i+1)&
                + gradQ(:,2,i_0,j_0+1) * dist_to_ecen_down(2,1,i+1)
              aux_un = aux(1,i_0,j_0+1)
            else
              q1 = qold2(:,i_0,j_0+1)
              q1l = qold(:,i_0,j_0+1)
              aux_up = auxu(1,i_0,j_0+1)
              aux_un = aux(1,i_0,j_0+1)
            end if
            call rotated_redist_fwave(q1,q1l,aux_up,aux_un, &
                 x,y,amdq_wall,apdq_wall,2,wall_height)
            cm_avg2 = cm_avg2 - dt/dx *lengths_supper(1,i+1)*(1/(area_supper(i)+&
              up_area_ij(i_0,j_0+1))) *amdq_wall

              ! left edge of merging neighbor cell:
            if (order .eq. 2) then
              call riemann_solver(qold2(:,i_0-2:i_0,j_0+1),qold2(:,i_0-1:i_0+1,j_0+1),&
               auxu(1,i_0-2:i_0,j_0+1),auxu(1,i_0-1:i_0+1,j_0+1),s,fwave,&
                  amdq,apdq,1,2,dt/dx,mthlim(3))
              apdq_c = apdq(:,2)
            else
              apdq_c = fp2(:,i_0,j_0+1)
            end if
            cm_avg2 = cm_avg2 - dt/dx * (1/(area_supper(i)+&
                up_area_ij(i_0,j_0+1)))* apdq_c

             ! upper edge of merging neighbor cell:
           if (order .eq. 2) then
             call riemann_solver(qold2(:,i_0,j_0:j_0+2),qold2(:,i_0,j_0+1:j_0+3),&
               auxu(1,i_0,j_0:j_0+2),auxu(1,i_0,j_0+1:j_0+3),s,fwave,&
               amdq,apdq,2,2,dt/dx,mthlim)
               amdq_c = amdq(:,2)
            else
              amdq_c = gm2(:,i_0,j_0+2)
            end if
           cm_avg2 = cm_avg2 - dt/dx * (1/(area_supper(i)+&
               up_area_ij(i_0,j_0+1)))* amdq_c

             ! right edge of merging neighbor cell:
           if (order .eq. 2) then
              call riemann_solver(qold2(:,i_0-1:i_0+1,j_0+1),qold2(:,i_0:i_0+2,j_0+2),&
              auxu(1,i_0-1:i_0+1,j_0+1),auxu(1,i_0:i_0+2,j_0+1),s,fwave,&
                amdq,apdq,1,2,dt/dx,mthlim)
              amdq_c = amdq(:,2)
            else
              amdq_c = fm2(:,i_0+1,j_0+1)
            end if
            cm_avg2 = cm_avg2 - dt/dx * (lengths_supper(2,i+1))* (1/(area_supper(i)+&
                up_area_ij(i_0,j_0+1)))* amdq_c

            ! update the corresponding affected cells:
            qnew2(:,i_0,j_0) =cm_avg2!(area_supper(i))/(area_supper(i)+area_supper(i+1))*cm_avg2
            qnew2(:,i_0,j_0+1)= cm_avg2!(area_supper(i+1)/(area_supper(i)+area_supper(i+1)))*cm_avg2

          case(3,6)
            if (type_supper(i).eq.3) then
              ixy = 2
            else if (type_supper(i).eq.6) then
              ixy = 1
            end if
            aux_un = aux(1,i_0,j_0)
            aux_up = auxu(1,i_0,j_0)
            if (order .eq. 2) then
              q1l = qold2(:,i_0,j_0) + gradQ2(:,1,i_0,j_0)*dist_to_ecen_up(1,1,i)&
                 + gradQ2(:,2,i_0,j_0)*dist_to_ecen_up(2,1,i)
              q1 = qold(:,i_0,j_0) + gradQ(:,1,i_0,j_0)*dist_to_ecen_down(1,1,i)&
                 + gradQ(:,2,i_0,j_0)*dist_to_ecen_down(2,1,i)
            else
              q1l = qold2(:,i_0,j_0)
              q1 = qold(:,i_0,j_0)
            end if
              call rotated_redist_fwave(q1l,q1,aux_up,aux_un, &
                    x,y,amdq_wall,apdq_wall,ixy,wall_height)

          ! the left edge
          if (order .eq. 2) then
           call riemann_solver(qold2(:,i_0-2:i_0,j_0),qold2(:,i_0-1:i_0+1,j_0),&
             auxu(1,i_0-2:i_0,j_0),auxu(1,i_0-1:i_0+1,j_0),s,fwave,&
            amdq_l,apdq_l,1,2,dt/dx,mthlim)
            apdq_c = apdq_l(:,2)
          ! the right edge
            call riemann_solver(qold2(:,i_0-1:i_0+1,j_0),qold2(:,i_0:i_0+2,j_0),&
              auxu(1,i_0-1:i_0+1,j_0),auxu(1,i_0:i_0+2,j_0),s,fwave, &
               amdq_r,apdq_r,1,2,dt/dx,mthlim)
               amdq_c = amdq_r(:,2)
            else
              amdq_c = fm2(:,i_0+1,j_0)
              apdq_c = fp2(:,i_0,j_0)
            end if
            if (area_supper(i)>0.5d0) then
              cm_avg2 = qold2(:,i_0,j_0)

               cm_avg2 = cm_avg2 - dt/dx * (lengths_supper(2,i)) * &
                 (1/area_supper(i)) * amdq_c

              cm_avg2 = cm_avg2 - dt/dx * (lengths_supper(4,i)) * &
                 (1/area_supper(i)) * apdq_c
            ! the upper edge
                if (order .eq. 2) then
                  call riemann_solver(qold2(:,i_0,j_0-1:j_0+1),qold2(:,i_0,j_0:j_0+2),&
                    auxu(:,i_0,j_0-1:j_0+1),auxu(:,i_0,j_0:j_0+2),s,fwave,&
                    amdq,apdq,2,2,dt/dx,mthlim)
                    amdq_c = amdq(:,2)
                else
                    amdq_c = gm2(:,i_0,j_0+1)
                end if
              cm_avg2 = cm_avg2 - dt/dx * (1/area_supper(i))*amdq_c
            ! the barrier edge
             cm_avg2 = cm_avg2 - dt/dx *lengths_supper(1,i)*(1/(area_supper(i)))&
                    *amdq_wall

            ! update the corresponding affected cells:
              qnew2(:,i_0,j_0) = cm_avg2
          else
              cm_avg2 = (area_supper(i)*qold2(:,i_0,j_0) + qold2(:,i_0,j_0+1)) &
               / ( 1+ area_supper(i))
               ! barrier edge:
               cm_avg2 = cm_avg2 - dt/dx *lengths_supper(1,i)*(1/(1+ area_supper(i)))&
                      *amdq_wall
              ! left edge:
              cm_avg2 = cm_avg2 - dt/dx * (lengths_supper(2,i)) * &
                (1/( 1+ area_supper(i)) ) * amdq_c
              ! right edge:
               cm_avg2 = cm_avg2 - dt/dx * (lengths_supper(4,i)) * &
                  (1/( 1+ area_supper(i)) ) * apdq_c
              ! upper left edge:
              if (order .eq. 2) then
                call riemann_solver(qold2(:,i_0-2:i_0,j_0+1),qold2(:,i_0-1:i_0+1,j_0+1),&
                  auxu(1,i_0-2:i_0,j_0+1),auxu(1,i_0-1:i_0+1,j_0+1),s,fwave,amdq,apdq,1,&
                   2, dt/dx,mthlim)
                apdq_c =apdq(:,2)
              else
                apdq_c = fp2(:,i_0,j_0+1)
              end if
              cm_avg2 = cm_avg2 - dt/dx* (1/( 1+ area_supper(i)) )*apdq_c
              ! upper upper edge:
              if (order .eq. 2) then
                call riemann_solver(qold2(:,i_0,j_0:j_0+2),qold2(:,i_0,j_0+1:j_0+3),&
                  auxu(1,i_0,j_0:j_0+2),auxu(1,i_0,j_0+1:j_0+3),s,fwave,amdq,apdq,2,&
                  2,dt/dx,mthlim)
                amdq_c =amdq(:,2)
              else
                  amdq_c = gm2(:,i_0,j_0+2)
              end if
              cm_avg2 = cm_avg2 - dt/dx*(1/(1+area_supper(i)))*amdq_c
              !upper right edge :
              if (order .eq. 2) then
                call riemann_solver(qold2(:,i_0-1:i_0+1,j_0+1),qold2(:,i_0:i_0+2,j_0+1),&
                  auxu(:,i_0-1:i_0+1,j_0+1),auxu(1,i_0:i_0+2,j_0+1),s,fwave,amdq,apdq,1,&
                  2, dt/dx,mthlim)
                amdq_c = amdq(:,2)
              else
                  amdq_c = fm2(:,i_0+1,j_0+1)
              end if
              cm_avg2 = cm_avg2 - dt/dx*(1/(1+area_supper(i)))*amdq_c

              ! update the corresponding affected cells:
              qnew2(:,i_0,j_0) =cm_avg2!(area_supper(i))/(area_supper(i)+1)*cm_avg2
              qnew2(:,i_0,j_0+1)=cm_avg2!1/(area_supper(i)+1)*cm_avg2

            end if
          case(8)
            ! get the old merged average  & update the old merged average:
              cm_avg2 = (1/(area_supper(i)+up_area_ij(i_0,j_0-1)))* &
        (area_supper(i)*qold2(:,i_0,j_0)+up_area_ij(i_0,j_0-1)*qold2(:,i_0,j_0-1))
              aux_up = auxu(1,i_0,j_0)
              aux_un = aux(1,i_0,j_0)

             ! get the fluctuations at each cell edge:
              ! get the cell edge value via gradient:
              ! barrier edge at indexed small cell
              if (order.eq.2) then
               q1 = qold2(:,i_0,j_0) + gradQ2(:,1,i_0,j_0)*dist_to_ecen_up(1,1,i) &
                 + gradQ2(:,2,i_0,j_0)*dist_to_ecen_up(2,1,i)   ! at the barrier center
               q1l = qold(:,i_0,j_0)+ gradQ(:,1,i_0,j_0)*dist_to_ecen_down(1,1,i) &
                 + gradQ(:,2,i_0,j_0)*dist_to_ecen_down(2,1,i)
             else
               q1 = qold2(:,i_0,j_0)
               q1l = qold(:,i_0,j_0)
             end if
             call rotated_redist_fwave(q1,q1l,aux_up,aux_un, &
                  x,y,amdq_wall,apdq_wall,1,wall_height)
             cm_avg2 = cm_avg2 - dt/dx *lengths_supper(1,i)*(1/(area_supper(i)+&
               up_area_ij(i_0,j_0-1))) *amdq_wall

              ! at the right:
              if (order .eq. 2) then
                call riemann_solver(qold2(:,i_0-1:i_0+1,j_0),qold2(:,i_0:i_0+2,j_0),&
                  auxu(1,i_0-1:i_0+1,j_0),auxu(1,i_0:i_0+2,j_0),s,fwave,&
                  amdq,apdq,1,2,dt/dx,mthlim)
                amdq_c = amdq(:,2)
              else
                  amdq_c = fm2(:,i_0+1,j_0)
              end if

              cm_avg2 = cm_avg2 - dt/dx *(1/(area_supper(i)+&
                  up_area_ij(i_0,j_0-1)))* amdq_c

              ! at the merging neighbor cell edges:
              ! at the barrier edge and moving clockwise
              aux_up = auxu(1,i_0,j_0+1)
              aux_un = aux(1,i_0,j_0+1)
              x = intersections(1,i+2)- intersections(1,i+1)
              y = intersections(2,i+2) - intersections(2,i+1)
              if (order .eq. 2) then
                q1 = qold2(:,i_0,j_0-1) + gradQ2(:,1,i_0,j_0-1)*dist_to_ecen_up(1,1,i+1)&
                   + gradQ2(:,2,i_0,j_0-1)*dist_to_ecen_up(2,1,i+1)
                q1l = qold(:,i_0,j_0-1) + gradQ(:,1,i_0,j_0-1)*dist_to_ecen_down(1,1,i+1)&
                  + gradQ(:,2,i_0,j_0-1) * dist_to_ecen_down(2,1,i+1)
              else
                  q1l = qold2(:,i_0,j_0-1)
                  q1 = qold(:,i_0,j_0-1)
              end if
              call rotated_redist_fwave(q1,q1l,aux_un,aux_up, &
                   x,y,apdq_wall,amdq_wall,1,wall_height)
              cm_avg2 = cm_avg2 - dt/dx *lengths_supper(1,i+1)*(1/(area_supper(i)+&
                up_area_ij(i_0,j_0-1))) *amdq_wall

              ! left edge of merging neighbor cell:
              if (order .eq. 2 ) then
                call riemann_solver(qold2(:,i_0-2:i_0,j_0),qold2(:,i_0-1:i_0+1,j_0),&
                  auxu(1,i_0-2:i_0,j_0),auxu(1,i_0-1:i_0+1,j_0),s,fwave,&
                    amdq,apdq,1,2,dt/dx,mthlim)
                apdq_c = apdq(:,2)
              else
                apdq_c = fp2(:,i_0,j_0)
              end if
              cm_avg2 = cm_avg2 - dt/dx * (lengths_supper(5,i))* (1/(area_supper(i)+&
                  up_area_ij(i_0,j_0-1)))* apdq_c

             ! upper edge of merging neighbor cell:
             if (order .eq. 2) then
               call riemann_solver(qold2(:,i_0,j_0-1:j_0+1),qold2(:,i_0,j_0:j_0+2),&
                 auxu(1,i_0,j_0-1:j_0+1),auxu(1,i_0,j_0:j_0+2),s,fwave,&
                 amdq,apdq,2,2,dt/dx,mthlim)
              amdq_c = amdq(:,2)
             else
               amdq_c = gm2(:,i_0,j_0+1)
             end if
             cm_avg2 = cm_avg2 - dt/dx * (1/(area_supper(i)+&
               up_area_ij(i_0,j_0-1)))* amdq_c

             ! right edge of merging neighbor cell:
             if (order .eq. 2) then
              call riemann_solver(qold2(:,i_0-1:i_0+1,j_0-1),qold2(:,i_0:i_0+2,j_0-1),&
                auxu(1,i_0-1:i_0+1,j_0-1),auxu(1,i_0:i_0+2,j_0-1),s,fwave,&
                amdq,apdq,1,2,dt/dx,mthlim)
              amdq_c = amdq(:,2)
            else
              amdq_c = fm2(:,i_0+1,j_0-1)
            end if
            cm_avg2 = cm_avg2 - dt/dx * lengths_supper(2,i+1)*(1/(area_supper(i)+&
                up_area_ij(i_0,j_0-1)))* amdq_c

            ! update the corresponding affected cells:
            qnew2(:,i_0,j_0) = cm_avg2!(area_supper(i))/(area_supper(i)+area_supper(i-1))*cm_avg2
            qnew2(:,i_0,j_0-1)= cm_avg2!(area_supper(i-1)/(area_supper(i)+area_supper(i-1)))*cm_avg2

          case(2,7)
            similar to case(3,6) but just merge with i_0+/-1,j_0
            if (type_supper(i).eq.2) then
              ixy = 2
              is = 1
              itop = 2
              ibot = 4
              ilat = i_0-1
              ilat2 = 1  ! for lateral edge center index for 0.5d0 case
              ilat3 = i_0-2  ! the lateral lateral edges for smaller than 0.5d0 case
            else if (type_supper(i).eq.7) then
              ixy = 1
              is = -1
              itop = 4
              ibot = 2
              ilat = i_0+1
              ilat2 = 2  ! for lateral edge center index for 0.5d0 case
              ilat3 = i_0+2  ! the lateral lateral edges for smaller than 0.5d0 case
            end if
            ! the barrier edge
            aux_un = aux(1,i_0,j_0)
            aux_up = auxu(1,i_0,j_0)
            if (order .eq. 2) then
              q1l = qold2(:,i_0,j_0) + gradQ2(:,1,i_0,j_0)*dist_to_ecen_up(1,1,i)&
                 + gradQ2(:,2,i_0,j_0)*dist_to_ecen_up(2,1,i)
              q1 = qold(:,i_0,j_0) + gradQ(:,1,i_0,j_0)*dist_to_ecen_down(1,1,i)&
                 + gradQ(:,2,i_0,j_0)*dist_to_ecen_down(2,1,i)
            else
              q1l = qold2(:,i_0,j_0)
              q1 = qold(:,i_0,j_0)
            end if
            call rotated_redist_fwave(q1l,q1,aux_up,aux_un, &
                x,y,amdq_wall,apdq_wall,ixy,wall_height)
                ! the bottom edge
            if (order .eq. 2) then
                 call riemann_solver(qold2(:,i_0,j_0-2:j_0),qold2(:,i_0,j_0-1:j_0+1),&
                   auxu(1,i_0,j_0-2:j_0),auxu(1,i_0,j_0-1:j_0+1),s,fwave,&
                  amdq_l,apdq_l,2,2,dt/dx,mthlim)
                  apdq_c = apdq_l(:,2)
            else
                  apdq_c = gp2(:,i_0,j_0)
            end if
                  ! the top edge
            if (order .eq. 2) then
                  call riemann_solver(qold2(:,i_0,j_0-1:j_0+1),qold2(:,i_0,j_0:j_0+2),&
                    auxu(1,i_0,j_0-1:j_0+1),auxu(1,i_0,j_0:j_0+2),s,fwave, &
                     amdq_r,apdq_r,2,2,dt/dx,mthlim)
                     amdq_c = amdq_r(:,2)
            else
                  amdq_c = gm2(:,i_0,j_0+1)
            end if
            if (area_supper(i)>0.5d0) then
                cm_avg2 = qold2(:,i_0,j_0)

                 cm_avg2 = cm_avg2 - dt/dx * (lengths_supper(ibot,i)) * &
                   (1/area_supper(i)) * apdq_c

                cm_avg2 = cm_avg2 - dt/dx * (lengths_supper(itop,i)) * &
                   (1/area_supper(i)) * amdq_c
            ! the lateral edge
                if (order .eq. 2) then
                  if(type_supper(i) .eq. 2) then
                    call riemann_solver(qold2(:,ilat-1:i_0,j_0),qold2(:,ilat:i_0+1,j_0),&
                       auxu(1,ilat-1:i_0,j_0),auxu(1,ilat:i_0+1,j_0),s,fwave,&
                      amdq,apdq,1,2,dt/dx,mthlim )
                    apdq_c = apdq(:,2)
                  else
                    call riemann_solver(qold2(:,i_0-1:ilat,j_0),qold2(:,i_0:ilat+1,j_0),&
                       auxu(1,i0-1:ilat,j_0),auxu(1,i_0:ilat+1,j_0),s,fwave,&
                      amdq,apdq,1,2,dt/dx,mthlim )
                    apdq_c = amdq(:,2)
                  end if
                else
                    if (type_supper(i) .eq. 2) then
                      apdq_c = fp2(:,i_0,j_0)
                    else
                      apdq_c = fm2(:,i_0+1,j_0)
                    end if
                end if
              cm_avg2 = cm_avg2 - dt/dx * (1/area_supper(i))*apdq_c
            ! the barrier edge
               cm_avg2 = cm_avg2 - dt/dx *lengths_supper(1,i)*(1/(area_supper(i)))&
                      *amdq_wall

            ! update the corresponding affected cells:
              qnew2(:,i_0,j_0) = cm_avg2
          else
              cm_avg2 = (area_supper(i)*qold2(:,i_0,j_0) + qold2(:,ilat,j_0)) &
               / ( 1+ area_supper(i))
               ! barrier edge:
               cm_avg2 = cm_avg2 - dt/dx *lengths_supper(1,i)*(1/(1+ area_supper(i)))&
                      *amdq_wall
              ! bottom edge:
              cm_avg2 = cm_avg2 - dt/dx * (lengths_supper(ibot,i)) * &
                (1/( 1+ area_supper(i)) ) * amdq_c
              ! top edge:
               cm_avg2 = cm_avg2 - dt/dx * (lengths_supper(itop,i)) * &
                  (1/( 1+ area_supper(i)) ) * apdq_c
              ! lateral top edge:
              if (order .eq. 2) then
                call riemann_solver(qold2(:,ilat,j_0-1:j_0+1),qold2(:,ilat,j_0:j_0+2),&
                  auxu(1,ilat,j_0-1:j_0+1),auxu(1,ilat,j_0:j_0+2),s,fwave,amdq,apdq,2,&
                  2,dt/dx,mthlim)
                  amdq_c =amdq(:,2)
              else
                amdq_c = gm2(:,ilat,j_0+1)
              end if
              cm_avg2 = cm_avg2 - dt/dx* (1/( 1+ area_supper(i)) )*amdq_c

              ! lateral lateral edge:
              if (order .eq. 2) then
                if (type_supper(i) .eq. 2) then
                  call riemann_solver(qold2(:,ilat3-1:ilat,j_0),qold2(:,ilat3:ilat+1,j_0),&
                    auxu(1,ilat3-1:ilat,j_0),auxu(1,ilat3:ilat+1,j_0),s,fwave,amdq,apdq,1,&
                    2,dt/dx,mthlim)
                  apdq_c = apdq(:,2)
                else
                  call riemann_solver(qold2(:,ilat-1:ilat3,j_0),qold2(:,ilat:ilat3+1,j_0),&
                    auxu(1,ilat-1:ilat3,j_0),auxu(1,ilat:ilat3+1,j_0),s,fwave,amdq,apdq,1,&
                    2,dt/dx,mthlim)
                  apdq_c = amdq(:,2)
                end if
              else
                if (type_supper(i) .eq. 2) then
                  apdq_c = fp2(:,ilat,j_0)
                else
                  apdq_c = fm2(:,ilat3,j_0)
                end if
              end if
              cm_avg2 = cm_avg2 - dt/dx*(1/(1+area_supper(i)))*apdq_c
              !lateral bottom edge :
              if (order .eq. 2) then
                call riemann_solver(qold2(:,ilat,j_0-2:j_0),qold2(:,ilat,j_0-1:j_0+1),&
                  auxu(1,ilat,j_0-2:j_0),auxu(1,ilat,j_0-1:j_0+1),s,fwave,amdq,apdq,2,2,&
                  dt/dx,mthlim)
                apdq_c = apdq(:,2)
              else
                apdq_c = gp2(:,ilat,j_0)
              end if
              cm_avg2 = cm_avg2 - dt/dx*(1/(1+area_supper(i)))*apdq_c

              ! update the corresponding affected cells:
              qnew2(:,i_0,j_0) =  cm_avg2!(area_supper(i))/(area_supper(i)+1)*cm_avg2
              qnew2(:,ilat,j_0)=  cm_avg2!(1)/(area_supper(i)+1)*cm_avg2

              ! do the same for lower cell?
            end if
          end select

!!!!!!!!!!!!!!!!!!!!!!!!   ! FOR LOWER CELLS  !!!!!!!!!!!!!!!!!!!
          select case(type_sunder(i))
          case(1)
            ! get the old merged average  & update the old merged average:
            cm_avg2 = (1/(area_sunder(i)+un_area_ij(i_0,j_0-1)))* &
        (area_sunder(i)*qold(:,i_0,j_0)+un_area_ij(i_0,j_0-1)*qold(:,i_0,j_0-1))
            aux_up = auxu(1,i_0,j_0)
            aux_un = aux(1,i_0,j_0)

             ! get the fluctuations at each cell edge:
              ! get the cell edge value via gradient:
              ! barrier edge at indexed small cell
            if (order.eq. 2) then
               q1 = qold2(:,i_0,j_0) + gradQ2(:,1,i_0,j_0)*dist_to_ecen_up(1,1,i) &
                 + gradQ2(:,2,i_0,j_0)*dist_to_ecen_up(2,1,i)   ! at the barrier center
               q1l = qold(:,i_0,j_0) + gradQ(:,1,i_0,j_0)*dist_to_ecen_down(1,1,i) &
                 + gradQ(:,2,i_0,j_0)*dist_to_ecen_down(2,1,i)
            else
                q1 =qold2(:,i_0,j_0)
                q1l =qold(:,i_0,j_0)
            end if
           call rotated_redist_fwave(q1,q1l,aux_up,aux_un, &
                x,y,amdq_wall,apdq_wall,2,wall_height)
           cm_avg2 = cm_avg2 - dt/dx *lengths_sunder(1,i)*(1/(area_sunder(i)+&
             un_area_ij(i_0,j_0-1))) *apdq_wall

              ! at the edge 2:
              if (order .eq. 2) then
                call riemann_solver(qold(:,i_0-1:i_0+1,j_0),qold(:,i_0:i_0+2,j_0),&
                  auxu(1,i_0-1:i_0+1,j_0),aux(1,i_0:i_0+2,j_0),s,fwave,&
                    amdq,apdq,1,2,dt/dx,mthlim)
                amdq_c = amdq(:,2)
              else
                amdq_c = fm(:,i_0+1,j_0)
              end if
              cm_avg2 = cm_avg2 - dt/dx * lengths_sunder(2,i)*(1/(area_sunder(i)+&
                  un_area_ij(i_0,j_0-1)))* amdq_c

              ! at the merging neighbor cell edges:
              ! at the barrier edge and moving clockwise
              aux_up = auxu(1,i_0,j_0-1)
              aux_un = aux(1,i_0,j_0-1)
              x = intersections(1,i)- intersections(1,i-1)
              y = intersections(2,i) - intersections(2,i-1)
              if (order .eq. 2) then
                q1 = qold2(:,i_0,j_0-1) + gradQ2(:,1,i_0,j_0-1)*dist_to_ecen_up(1,1,i-1)&
                   + gradQ2(:,2,i_0,j_0-1)*dist_to_ecen_up(2,1,i-1)
                q1l = qold(:,i_0,j_0-1) + gradQ(:,1,i_0,j_0-1)*dist_to_ecen_down(1,1,i-1)&
                  + gradQ(:,2,i_0,j_0-1) * dist_to_ecen_down(2,1,i-1)
              else
                q1 = qold2(:,i_0,j_0-1)
                q1l = qold(:,i_0,j_0-1)
              end if
              call rotated_redist_fwave(q1,q1l,aux_up,aux_un, &
                   x,y,amdq_wall,apdq_wall,2,wall_height)
              cm_avg2 = cm_avg2 - dt/dx *lengths_sunder(1,i-1)*(1/(area_sunder(i)+&
                un_area_ij(i_0,j_0-1))) *apdq_wall

              ! left edge of merging neighbor cell:
              if (order .eq. 2) then
                call riemann_solver(qold(:,i_0-2:i_0,j_0-1),qold(:,i_0-1:i_0+1,j_0-1),&
                  aux(1,i_0-2:i_0,j_0-1),aux(1,i_0-1:i_0+1,j_0-1),s,fwave,&
                    amdq,apdq,1,2,dt/dx,mthlim)
                apdq_c =apdq(:,2)
              else
                apdq_c = fp(:,i_0,j_0-1)
              end if
              cm_avg2 = cm_avg2 - dt/dx *lengths_sunder(5,i-1)*(1/(area_sunder(i)+&
                  un_area_ij(i_0,j_0-1)))* apdq_c

             ! lower edge of merging neighbor cell:
             if (order .eq.  2) then
               call riemann_solver(qold(:,i_0,j_0-3:j_0-1),qold(:,i_0,j_0-2:j_0),&
                 aux(1,i_0,j_0-3:j_0-1),aux(1,i_0,j_0-2:j_0),s,fwave,&
                 amdq,apdq,2,2,dt/dx,mthlim)
              apdq_c =apdq(:,2)
             else
               apdq_c = gp(:,i_0,j_0-1)
             end if
           cm_avg2 = cm_avg2 - dt/dx * (1/(area_sunder(i)+&
               un_area_ij(i_0,j_0-1)))* apdq_c

             ! right edge of merging neighbor cell:
             if (order .eq. 2) then
              call riemann_solver(qold(:,i_0-1:i_0+1,j_0-1),qold(:,i_0:i_0+2,j_0-1),&
              aux(1,i_0-1:i_0+1,j_0-1),aux(1,i_0:i_0+2,j_0-1),s,fwave,&
                amdq,apdq,1,2,dt/dx,mthlim )
              amdq_c = amdq(:,2)
            else
              amdq_c=fm(:,i_0+1,j_0-1)
            end if
              cm_avg2 = cm_avg2 - dt/dx * (1/(area_sunder(i)+&
                  un_area_ij(i_0,j_0-1)))* amdq_c

            ! update the corresponding affected cells:
            qnew(:,i_0,j_0) = cm_avg2!(area_sunder(i))/(area_sunder(i)+area_sunder(i-1))*cm_avg2
            qnew(:,i_0,j_0-1)= cm_avg2!(area_sunder(i-1)/(area_sunder(i)+area_sunder(i-1)))*cm_avg2

          case(3,6)
            if (type_sunder(i).eq.3) then
              ixy = 2
            else if (type_sunder(i).eq.6) then
              ixy = 1
            end if
            aux_un = aux(1,i_0,j_0)
            aux_up = auxu(1,i_0,j_0)
            if (order .eq. 2) then
              q1l = qold2(:,i_0,j_0) + gradQ2(:,1,i_0,j_0)*dist_to_ecen_up(1,1,i)&
                 + gradQ2(:,2,i_0,j_0)*dist_to_ecen_up(2,1,i)
              q1 = qold(:,i_0,j_0) + gradQ(:,1,i_0,j_0)*dist_to_ecen_down(1,1,i)&
                 + gradQ(:,2,i_0,j_0)*dist_to_ecen_down(2,1,i)
            else
                q1l = qold2(:,i_0,j_0)
                q1 = qold(:,i_0,j_0)
            end if
             call rotated_redist_fwave(q1l,q1,aux_up,aux_un, &
                  x,y,amdq_wall,apdq_wall,ixy,wall_height)
                ! the left edge
              if (order .eq. 2) then
                 call riemann_solver(qold(:,i_0-2:i_0,j_0),qold(:,i_0-1:i_0+1,j_0),&
                   aux(1,i_0-2:i_0,j_0),aux(1,i_0-1:i_0+1,j_0),s,fwave,&
                  amdq_l,apdq_l,1,2,dt/dx,mthlim)
                  apdq_c = apdq_l(:,2)
              else
                  apdq_c = fp(:,i_0,j_0)
              end if
                  ! the right edge
              if (order .eq. 2) then
                  call riemann_solver(qold(:,i_0-1:i_0+1,j_0),qold(:,i_0:i_0+2,j_0),&
                     aux(:,i_0-1:i_0+1,j_0),aux(:,i_0:i_0+2,j_0),s,fwave, &
                     amdq_r,apdq_r,1,2,dt/dx,mthlim)
                  amdq_c = amdq_r(:,2)
              else
                    amdq_c = fm(:,i_0+1,j_0)
              end if
            if (area_sunder(i)>0.5d0) then
              cm_avg2 = qold(:,i_0,j_0)

               cm_avg2 = cm_avg2 - dt/dx * (lengths_sunder(2,i)) * &
                 (1/area_sunder(i)) * amdq_c

              cm_avg2 = cm_avg2 - dt/dx * (lengths_sunder(4,i)) * &
                 (1/area_sunder(i)) * apdq_c
            ! the lower edge
                if (order .eq. 2) then
                  call riemann_solver(qold(:,i_0,j_0-2:j_0),qold(:,i_0,j_0-1:j_0+1),&
                    aux(:,i_0,j_0-2:j_0),aux(:,i_0,j_0-1:j_0+1),s,fwave,&
                    amdq,apdq,2,2,dt/dx,mthlim)
                  apdq_c =apdq(:,2)
                else
                  apdq_c = gp(:,i_0,j_0)
                end if
              cm_avg2 = cm_avg2 - dt/dx * (1/area_sunder(i))*apdq_c
            ! the barrier edge
               cm_avg2 = cm_avg2 - dt/dx *lengths_sunder(1,i)*(1/(area_sunder(i)))&
                      *apdq_wall

            ! update the corresponding affected cells:
              qnew(:,i_0,j_0) = cm_avg2
            else
              cm_avg2 = (area_sunder(i)*qold(:,i_0,j_0) + qold(:,i_0,j_0-1)) &
               / ( 1+ area_sunder(i))
               ! barrier edge:
               cm_avg2 = cm_avg2 - dt/dx *lengths_sunder(1,i)*(1/(1+ area_sunder(i)))&
                      *apdq_wall
              ! left edge:
              cm_avg2 = cm_avg2 - dt/dx * (lengths_sunder(2,i)) * &
                (1/( 1+ area_sunder(i)) ) * amdq_c
              ! right edge:
               cm_avg2 = cm_avg2 - dt/dx * (lengths_sunder(4,i)) * &
                  (1/( 1+ area_sunder(i)) ) * apdq_c
              ! lower left edge:
              if (order .eq. 2) then
                call riemann_solver(qold(:,i_0-2:i_0,j_0-1),qold(:,i_0-1:i_0+1,j_0-1),&
                  aux(1,i_0-2:i_0,j_0-1),aux(1,i_0-1:i_0+1,j_0-1),s,fwave,amdq,apdq,1,&
                  2,dt/dx,mthlim)
                apdq_c = apdq(:,2)
              else
                apdq_c = fp(:,i_0,j_0-1)
              end if
              cm_avg2 = cm_avg2 - dt/dx* (1/( 1+ area_sunder(i)) )*apdq_c
              ! lower lower edge:
              if (order .eq. 2) then
                  call riemann_solver(qold(:,i_0,j_0-3:j_0-1),qold(:,i_0,j_0-2:j_0),&
                    aux(1,i_0,j_0-3:j_0-1),aux(1,i_0,j_0-2:j_0),s,fwave,amdq,apdq,2,&
                    2,dt/dx,mthlim)
               apdq_c =apdq(:,2)
              else
                apdq_c = gp(:,i_0,j_0-1)
              end if
              cm_avg2 = cm_avg2 - dt/dx*(1/(1+area_sunder(i)))*apdq_c
              !lower right edge :
              if (order .eq. 2) then
                call riemann_solver(qold(:,i_0-1:i_0+1,j_0-1),qold(:,i_0:i_0+2,j_0-1),&
                  aux(1,i_0-1:i_0+1,j_0-1),aux(1,i_0:i_0+2,j_0-1),s,fwave,amdq,apdq,1,&
                  2,dt/dx,mthlim)
                amdq_c = amdq(:,2)
              else
                amdq_c = fm(:,i_0+1,j_0)
              end if
              cm_avg2 = cm_avg2 - dt/dx*(1/(1+area_sunder(i)))*amdq_c

              ! update the corresponding affected cells:
              qnew(:,i_0,j_0) = cm_avg2!(area_sunder(i))/(area_sunder(i)+1)*cm_avg2
              qnew(:,i_0,j_0-1)= cm_avg2!(1)/(area_sunder(i)+1)*cm_avg2

            end if
          case(5)
            ! get the old merged average  & update the old merged average:
            cm_avg2 = (1/(area_sunder(i)+un_area_ij(i_0,j_0+1)))* &
        (area_sunder(i)*qold(:,i_0,j_0)+un_area_ij(i_0,j_0+1)*qold(:,i_0,j_0+1))
            aux_up = auxu(1,i_0,j_0)
            aux_un = aux(1,i_0,j_0)

             ! get the fluctuations at each cell edge:
              ! get the cell edge value via gradient:
              ! barrier edge at indexed small cell
            if (order .eq. 2) then
             q1 = qold2(:,i_0,j_0) + gradQ2(:,1,i_0,j_0)*dist_to_ecen_up(1,1,i) &
               + gradQ2(:,2,i_0,j_0)*dist_to_ecen_up(2,1,i)   ! at the barrier center

             q1l = qold(:,i_0,j_0) + gradQ(:,1,i_0,j_0)*dist_to_ecen_down(1,1,i) &
               + gradQ(:,2,i_0,j_0)*dist_to_ecen_down(2,1,i)
            else
              q1 = qold2(:,i_0,j_0)
              q1l = qold(:,i_0,j_0)
            end if
           call rotated_redist_fwave(q1,q1l,aux_up,aux_un, &
                x,y,amdq_wall,apdq_wall,1,wall_height)
           cm_avg2 = cm_avg2 - dt/dx *lengths_sunder(1,i)*(1/(area_sunder(i)+&
             un_area_ij(i_0,j_0+1))) *apdq_wall

              ! at the left:
            if (order .eq. 2) then
              call riemann_solver(qold(:,i_0-2:i_0,j_0),qold(:,i_0-1:i_0+1,j_0),&
                aux(1,i_0-2:i_0,j_0),aux(1,i_0-1:i_0+1,j_0),s,fwave,&
                  amdq,apdq,1,2,dt/dx,mthlim)
              apdq_c =apdq(:,2)
            else
              apdq_c = fp(:,i_0,j_0)
            end if
            cm_avg2 = cm_avg2 - dt/dx * (1/(area_sunder(i)+&
                un_area_ij(i_0,j_0+1)))* apdq_c

              ! at the merging neighbor cell edges:
              ! at the barrier edge and moving clockwise
            aux_up = aux(1,i_0,j_0+1)
            aux_un = auxu(1,i_0,j_0+1)
            x = intersections(1,i)- intersections(1,i-1)
            y = intersections(2,i) - intersections(2,i-1)

            if (order .eq. 2) then
              q1 = qold(:,i_0,j_0+1) + gradQ(:,1,i_0,j_0+1)*dist_to_ecen_down(1,1,i-1)&
                 + gradQ(:,2,i_0,j_0+1)*dist_to_ecen_down(2,1,i-1)
              q1l = qold2(:,i_0,j_0+1) + gradQ2(:,1,i_0,j_0+1)*dist_to_ecen_up(1,1,i-1)&
                + gradQ2(:,2,i_0,j_0+1) * dist_to_ecen_up(2,1,i-1)
            else
              q1 = qold(:,i_0,j_0+1)
              q1l = qold2(:,i_0,j_0+1)
            end if

            call rotated_redist_fwave(q1l,q1,aux_un,aux_up, &
                 x,y,amdq_wall,apdq_wall,1,wall_height)
            cm_avg2 = cm_avg2 - dt/dx *lengths_sunder(1,i-1)*(1/(area_sunder(i)+&
              un_area_ij(i_0,j_0+1))) *apdq_wall

              ! left edge of merging neighbor cell:
              if (order .eq. 2) then
                call riemann_solver(qold(:,i_0-2:i_0,j_0+1),qold(:,i_0-1:i_0+1,j_0+1),&
                 aux(1,i_0-2:i_0,j_0+1),aux(1,i_0-1:i_0+1,j_0+1),s,fwave,&
                    amdq,apdq,1,2,dt/dx,mthlim)
                apdq_c =apdq(:,2)
              else
                apdq_c = fp(:,i_0,j_0+1)
              end if
              cm_avg2 = cm_avg2 - dt/dx * lengths_sunder(3,i-1)* (1/(area_sunder(i)+&
                  un_area_ij(i_0,j_0+1)))* apdq_c

             ! lower edge of merging neighbor cell:
             if (order .eq. 2) then
               call riemann_solver(qold(:,i_0,j_0-2:j_0),qold(:,i_0,j_0-1:j_0+1),&
                aux(1,i_0,j_0-2:j_0),aux(1,i_0,j_0-1:j_0+1),s,fwave,&
                 amdq,apdq,2,2,dt/dx,mthlim)
              apdq_c =apdq(:,2)
             else
               apdq_c = gp(:,i_0,j_0)
             end if
             cm_avg2 = cm_avg2 - dt/dx * (1/(area_sunder(i)+&
               un_area_ij(i_0,j_0+1)))* apdq_c

             ! right edge of merging neighbor cell:
             if (order .eq. 2) then
              call riemann_solver(qold(:,i_0-1:i_0+1,j_0),qold(:,i_0:i_0+2,j_0),&
                aux(1,i_0-1:i_0+1,j_0),aux(1,i_0:i_0+2,j_0),s,fwave,&
                amdq,apdq,1,2,dt/dx,mthlim)
              amdq_c = amdq(:,2)
            else
              amdq_c = fm(:,i_0+1,j_0)
            end if
              cm_avg2 = cm_avg2 - dt/dx * (lengths_sunder(2,i))*(1/(area_sunder(i)+&
                  un_area_ij(i_0,j_0+1)))* amdq_c

            ! update the corresponding affected cells:
            qnew(:,i_0,j_0) = cm_avg2!(area_sunder(i))/(area_sunder(i)+area_sunder(i+1))*cm_avg2
            qnew(:,i_0,j_0+1)=  cm_avg2!(area_sunder(i+1)/(area_sunder(i)+area_sunder(i+1)))*cm_avg2

          case(2,7)
            ! similar to case(3,6) but just merge with i_0+/-1,j_0
            if (type_sunder(i).eq.2) then
              ixy = 2
              is = -1
              itop = 2
              ibot = 4
              ilat = i_0+1
              ilat2 = 1  ! for lateral edge center index for 0.5d0 case
              ilat3 = i_0+2  ! the lateral lateral edges for smaller than 0.5d0 case
            else if (type_sunder(i).eq.7) then
              ixy = 1
              is = 1
              itop = 4
              ibot = 2
              ilat = i_0-1
              ilat2 = 2  ! for lateral edge center index for 0.5d0 case
              ilat3 = i_0-2  ! the lateral lateral edges for smaller than 0.5d0 case
            end if
            ! the barrier edge
            aux_un = aux(1,i_0,j_0)
            aux_up = auxu(1,i_0,j_0)
            if (order .eq. 2) then
              q1l = qold2(:,i_0,j_0) + gradQ2(:,1,i_0,j_0)*dist_to_ecen_up(1,1,i)&
                 + gradQ2(:,2,i_0,j_0)*dist_to_ecen_up(2,1,i)
              q1 = qold(:,i_0,j_0) + gradQ(:,1,i_0,j_0)*dist_to_ecen_down(1,1,i)&
                 + gradQ(:,2,i_0,j_0)*dist_to_ecen_down(2,1,i)
            else
              q1l = qold2(:,i_0,j_0)
              q1 = qold(:,i_0,j_0)
            end if
           call rotated_redist_fwave(q1l,q1,aux_up,aux_un, &
                x,y,amdq_wall,apdq_wall,ixy,wall_height)
                ! the bottom edge
            if (order .eq. 2) then
                 call riemann_solver(qold(:,i_0,j_0-2:j_0),qold(:,i_0,j_0-1:j_0+1),&
                   aux(1,i_0,j_0-2:j_0),aux(1,i_0,j_0-1:j_0+1),s,fwave,&
                   amdq_l,apdq_l,2,2,dt/dx,mthlim)
                   apdq_c = apdq_l(:,2)
            else
                  apdq_c = gp(:,i_0,j_0)
            end if
                  ! the top edge
              if (order .eq. 2) then
                  call riemann_solver(qold(:,i_0,j_0-1:j_0+1),qold(:,i_0,j_0:j_0+2),&
                    aux(1,i_0,j_0-1:j_0+1),aux(1,i_0,j_0:j_0+2),s,fwave, &
                     amdq_r,apdq_r,2,2,dt/dx,mthlim)
                   amdq_c = amdq_r(:,2)
              else
                    amdq_c = gm(:,i_0,j_0+1)
              end if
            if (area_sunder(i)>0.5d0) then
              cm_avg2 = qold(:,i_0,j_0)

               cm_avg2 = cm_avg2 - dt/dx * (lengths_sunder(ibot,i)) * &
                 (1/area_sunder(i)) * apdq_c

              cm_avg2 = cm_avg2 - dt/dx * (lengths_sunder(itop,i)) * &
                 (1/area_sunder(i)) * amdq_c
            ! the lateral edge
              if (order .eq. 2) then
               if (type_sunder(i) .eq. 7) then
                call riemann_solver(qold(:,ilat-1:i_0,j_0),qold(:,ilat:i_0+1,j_0),&
                  aux(1,ilat-1:i_0,j_0),aux(1,ilat:i_0+1,j_0),s,fwave,&
                  amdq,apdq,1,2,dt/dx,mthlim)
                 apdq_c = apdq(:,2)
               else
                 call riemann_solver(qold(:,i_0-1:ilat,j_0),qold(:,i_0:ilat+1,j_0),&
                   aux(1,i_0-1:ilat,j_0),aux(1,i_0:ilat+1,j_0),s,fwave,&
                   amdq,apdq,1,2,dt/dx,mthlim)
                  apdq_c = amdq(:,2)
                end if
              else
                if (type_sunder(i).eq.2) then
                  apdq_c = fm(:,i_0+1,j_0)
                else
                  apdq_c = fp(:,i_0,j_0)
                end if
              end if
              cm_avg2 = cm_avg2 - dt/dx * (1/area_sunder(i))*apdq_c
            ! the barrier edge
               cm_avg2 = cm_avg2 - dt/dx *lengths_sunder(1,i)*(1/(area_sunder(i)))&
                      *apdq_wall

            ! update the corresponding affected cells:
              qnew(:,i_0,j_0) = cm_avg2
            else
              cm_avg2 = (area_supper(i)*qold(:,i_0,j_0) + qold(:,ilat,j_0)) &
               / ( 1+ area_sunder(i))
               ! barrier edge:
               cm_avg2 = cm_avg2 - dt/dx *lengths_sunder(1,i)*(1/(1+ area_sunder(i)))&
                      *apdq_wall
              ! bottom edge:
              cm_avg2 = cm_avg2 - dt/dx * (lengths_sunder(ibot,i)) * &
                (1/( 1+ area_sunder(i)) ) * apdq_c
              ! top edge:
               cm_avg2 = cm_avg2 - dt/dx * (lengths_sunder(itop,i)) * &
                  (1/( 1+ area_sunder(i)) ) * amdq_c
              ! lateral top edge:
              if (order .eq. 2) then
                call riemann_solver(qold(:,ilat,j_0-1:j_0+1),qold(:,ilat,j_0:j_0+2),&
                  aux(1,ilat,j_0-1:j_0+1),aux(1,ilat,j_0:j_0+2),s,fwave,amdq,apdq,2,&
                  2,dt/dx,mthlim)
                amdq_c = amdq(:,2)
              else
                amdq_c = gm(:,ilat,j_0+1)
              end if
              cm_avg2 = cm_avg2 - dt/dx* (1/( 1+ area_sunder(i)) )*amdq_c
              ! lateral lateral edge:
              if (order .eq. 2) then
                if (type_sunder(i) .eq. 7) then
                  call riemann_solver(qold(:,ilat3-1:ilat,j_0),qold(:,ilat3:ilat+1,j_0),&
                    aux(1,ilat3-1:ilat,j_0),aux(1,ilat3:ilat+1,j_0),s,fwave,amdq,apdq,1,&
                    2,dt/dx,mthlim)
                   apdq_c = apdq(:,2)
                else
                  call riemann_solver(qold(:,ilat-1:ilat3,j_0),qold(:,ilat:ilat3+1,j_0),&
                    aux(1,ilat-1:ilat3,j_0),aux(1,ilat:ilat3+1,j_0),s,fwave,amdq,apdq,1,&
                    2,dt/dx,mthlim)
                   apdq_c = amdq(:,2)
                end if
              else
                if (type_sunder(i) .eq. 2) then
                  apdq_c = fm(:,ilat3,j_0)
                else
                  apdq_c = fp(:,ilat,j_0)
                end if
              end if
              cm_avg2 = cm_avg2 - dt/dx*(1/(1+area_sunder(i)))*apdq_c
              !lateral bottom edge :
              if (order .eq. 2) then
                call riemann_solver(qold(:,ilat,j_0-2:j_0),qold(:,ilat,j_0-1:j_0+1),&
                  aux(1,ilat,j_0-2:j_0),aux(1,ilat,j_0-1:j_0+1),s,fwave,amdq,apdq,2,&
                  2,dt/dx,mthlim)
                  apdq_c = apdq(:,2)
              else
                apdq_c = gp(:,ilat,j_0)
              end if
              cm_avg2 = cm_avg2 - dt/dx*(1/(1+area_sunder(i)))*apdq_c

              ! update the corresponding affected cells:
              qnew(:,i_0,j_0) =  cm_avg2!(area_sunder(i))/(area_sunder(i)+1)*cm_avg2
              qnew(:,ilat,j_0)= cm_avg2!(1)/(area_sunder(i)+1)*cm_avg2

              ! do the same for lower cell?
            end if
          end select
        end do

    end subroutine


    subroutine rotated_redist_fwave(ql,qr,aux_up,aux_un,x,&
          y,amdq_wall,apdq_wall,ixy,wall_height)

      ! ixy = 2 means the barrier is positively sloped
      ! ixy = 1 means the barrier is negatively sloped (important for rotating directions)
      implicit none
      real(8) :: ql(3),qr(3),aux_up,aux_un,wall_height
      real(8) :: amdq_wall(3), apdq_wall(3),hL,hR,huL,hvL,huR,hvR
      real(8) :: wave_wall(3,3), s_wall(3),hstarL,hstarR,ustarL,ustarR
      real(8) :: x,y
      integer :: ixy
      logical :: L2R,R2L

      real(8) :: n_vec(2),n_vec2(2),t_vec(2)

      amdq_wall = 0.d0
      apdq_wall = 0.d0
      n_vec = -(/-(y),x/)
      n_vec2 = (/(y),x/) ! for turning back to original
      n_vec = n_vec/(sqrt(x**2+y**2))
      n_vec2 = n_vec2/(sqrt(x**2+y**2))
      t_vec = -(/x,(y)/)
      t_vec = t_vec/(sqrt(x**2+y**2))

      call rotate_state(ql,ql,n_vec,t_vec)
      call rotate_state(qr,qr,n_vec,t_vec)

      hL = ql(1)
      hR = qr(1)
      huL = ql(2)
      huR = qr(2)
      hvL= ql(3)
      hvR = qr(3)

      call barrier_passing(hL,hR,huL,huR,aux_up,aux_un,wall_height,&
                L2R,R2L,hstarL,hstarR,ustarL,ustarR)
        call redistribute_fwave(1,ql,qr,aux_up,aux_un,wall_height,&
           1,wave_wall,s_wall,amdq_wall,apdq_wall,3,&
           3,L2R,R2L)

       if (ixy.eq.1) then
        call rotate_state(amdq_wall,amdq_wall,t_vec,n_vec2)
        call rotate_state(apdq_wall,apdq_wall,t_vec,n_vec2)
      else
        call rotate_state(amdq_wall,amdq_wall,-t_vec,-n_vec2)
        call rotate_state(apdq_wall,apdq_wall,-t_vec,-n_vec2)
      end if


    end subroutine



    subroutine grad_calc(ii,jj,mx,my,dx,dy,N_cells,gradQ,gradQ2,qold,qold2,type_supper,&
      dist_for_grad_up,dist_for_grad_down,cen_grid_down,cen_grid_up)
      ! calculates gradients over the grid, only for cut cells and their neighboring full cells
      ! the rest are zeros (sparse array)
      implicit none
      integer :: ii(:), jj(:), type_supper(:),N_cells,mx,my
      real(8) :: delQ1(3,3),delQ2(3,4),delR1(2,3),delR2(2,4),delRT1(3,2),delRT2(4,2)
      real(8) :: gradQ(3,2,-1:mx+2,-1:my+2),gradQ2(3,2,-1:mx+2,-1:my+2)
      real(8) :: qold(3,-1:mx+2,-1:my+2),qold2(3,-1:mx+2,-1:my+2)
      real(8) :: dist_for_grad_up(2,4,-1:mx+2,-1:my+2)
      real(8) :: dist_for_grad_down(2,4,-1:mx+2,-1:my+2)
      real(8) :: cen_grid_up(2,-1:mx+2,-1:my+2), cen_grid_down(2,-1:mx+2,-1:my+2)
      real(8) :: temp(2,2),temp1(3,2),temp2(4,2),dx,dy
      integer :: i, j, k, n , i_0,j_0


      gradQ2(:,:,ii(1)-1,jj(1)) = nabla_Q_normal(qold2,cen_grid_up,ii(1)-1,jj(1)&
            ,mx,my,dx,dy)
      gradQ(:,:,ii(1)-1,jj(1)) = nabla_Q_normal(qold,cen_grid_down,ii(1)-1,jj(1)&
            , mx,my,dx,dy)
      gradQ2(:,:,ii(N_cells)+1,jj(N_cells)) = nabla_Q_normal(qold2,cen_grid_up ,&
         ii(N_cells)+1,jj(N_cells),mx,my,dx,dy)
      gradQ(:,:,ii(N_cells)+1,jj(N_cells)) = nabla_Q_normal(qold,cen_grid_down ,&
         ii(N_cells)+1,jj(N_cells),mx,my,dx,dy)
      do i = 1, N_cells
        i_0 = ii(i)
        j_0 = jj(i)
        ! add the grad for neighboring normal full sized cells accordingly to types
        ! gradQ2(:,:,i_0,j_0+1) = nabla_Q_normal(qold2,i_0,j_0+1,mx,my,dx,dy)  ! if above
        ! gradQ(:,:,i_0,j_0-1) = nabla_Q_normal(qold,i_0,j_0-1,mx,my,dx,dy)  ! if below
        select case(type_supper(i))
        case(4)
          delQ1(:,1) = qold2(:,i_0-1,j_0) - qold2(:,i_0,j_0)
          delQ1(:,2) = qold2(:,i_0-1,j_0+1) - qold2(:,i_0,j_0)
          delQ1(:,3) = qold2(:,i_0,j_0+1) - qold2(:,i_0,j_0)
          delR1(:,:) = dist_for_grad_up(:,1:3,i_0,j_0)
          gradQ2(:,:,i_0,j_0) = nabla_Q(delQ1,delR1,1,dx)

          delQ2(:,1) = qold(:,i_0-1,j_0) - qold(:,i_0,j_0)
          delQ2(:,2) = qold(:,i_0,j_0+1) - qold(:,i_0,j_0)
          delQ2(:,3) = qold(:,i_0+1,j_0) - qold(:,i_0,j_0)
          delQ2(:,4) = qold(:,i_0,j_0-1) - qold(:,i_0,j_0)
          delR2(:,:) = dist_for_grad_down(:,1:4,i_0,j_0)
          gradQ(:,:,i_0,j_0) = nabla_Q(delQ2,delR2,2,dx)
          gradQ(:,:,i_0,j_0-1) = nabla_Q_normal(qold,cen_grid_down, &
                i_0,j_0-1,mx,my,dx,dy)

        case(3,6)
          gradQ2(:,:,i_0,j_0+1) = nabla_Q_normal(qold2,cen_grid_up, &
               i_0,j_0+1,mx,my,dx,dy)
          gradQ2(:,:,i_0,j_0+2) = nabla_Q_normal(qold2,cen_grid_up,&
              i_0,j_0+2,mx,my,dx,dy)
          gradQ(:,:,i_0,j_0-1) = nabla_Q_normal(qold,cen_grid_down,&
                i_0,j_0-1,mx,my,dx,dy)
          gradQ(:,:,i_0,j_0-2) = nabla_Q_normal(qold,cen_grid_down,&
              i_0,j_0-2,mx,my,dx,dy)
          delQ1(:,1) = qold2(:,i_0-1,j_0) - qold2(:,i_0,j_0)
          delQ1(:,2) = qold2(:,i_0,j_0+1) -qold2(:,i_0,j_0)
          delQ1(:,3) = qold2(:,i_0+1,j_0) -qold2(:,i_0,j_0)
          delR1(:,:) = dist_for_grad_up(:,1:3,i_0,j_0)
          gradQ2(:,:,i_0,j_0) = nabla_Q(delQ1,delR1,1,dx)
          delQ1(:,1) = qold(:,i_0-1,j_0) -qold(:,i_0,j_0)
          delQ1(:,2) = qold(:,i_0+1,j_0) - qold(:,i_0,j_0)
          delQ1(:,3) = qold(:,i_0,j_0-1) - qold(:,i_0,j_0)
          delR1(:,:) = dist_for_grad_down(:,1:3,i_0,j_0)
          gradQ(:,:,i_0,j_0) = nabla_Q(delQ1,delR1,1,dx)
        case(1)
          gradQ2(:,:,i_0,j_0+1) = nabla_Q_normal(qold2,cen_grid_up,&
               i_0,j_0+1,mx,my,dx,dy)  ! if above
          delQ2(:,1) = qold2(:,i_0-1,j_0) - qold2(:,i_0,j_0)
          delQ2(:,2) = qold2(:,i_0,j_0+1) - qold2(:,i_0,j_0)
          delQ2(:,3) = qold2(:,i_0+1,j_0) - qold2(:,i_0,j_0)
          delQ2(:,4) = qold2(:,i_0,j_0-1) - qold2(:,i_0,j_0)
          delR2(:,:) = dist_for_grad_up(:,1:4,i_0,j_0)
          gradQ2(:,:,i_0,j_0) = nabla_Q(delQ2,delR2,2,dx)
          delQ1(:,1) = qold(:,i_0+1,j_0)-qold(:,i_0,j_0)
          delQ1(:,2) = qold(:,i_0+1,j_0-1)-qold(:,i_0,j_0)
          delQ1(:,3) = qold(:,i_0,j_0-1) -qold(:,i_0,j_0)
          delR1(:,:) = dist_for_grad_down(:,1:3,i_0,j_0)
          gradQ(:,:,i_0,j_0) = nabla_Q(delQ1,delR1,1,dx)
        case(5)
          delQ1(:,1) = qold2(:,i_0,j_0+1) - qold2(:,i_0,j_0)
          delQ1(:,2) = qold2(:,i_0+1,j_0+1) - qold2(:,i_0,j_0)
          delQ1(:,3) = qold2(:,i_0+1,j_0) - qold2(:,i_0,j_0)
          delR1(:,:) = dist_for_grad_up(:,1:3,i_0,j_0)
          gradQ2(:,:,i_0,j_0) = nabla_Q(delQ1,delR1,1,dx)
          delQ2(:,1) = qold(:,i_0-1,j_0) - qold(:,i_0,j_0)
          delQ2(:,2) = qold(:,i_0,j_0+1) - qold(:,i_0,j_0)
          delQ2(:,3) = qold(:,i_0+1,j_0) - qold(:,i_0,j_0)
          delQ2(:,4) = qold(:,i_0,j_0-1) - qold(:,i_0,j_0)
          delR2(:,:) = dist_for_grad_down(:,1:4,i_0,j_0)
          gradQ(:,:,i_0,j_0) = nabla_Q(delQ2,delR2,2,dx)
          gradQ(:,:,i_0,j_0-1) = nabla_Q_normal(qold,cen_grid_down,&
              i_0,j_0-1,mx,my,dx,dy)

        case(7)
          gradQ2(:,:,i_0+1,j_0) = nabla_Q_normal(qold2,cen_grid_up,i_0+1,&
               j_0,mx,my,dx,dy)
          gradQ(:,:,i_0-1,j_0) = nabla_Q_normal(qold,cen_grid_down,i_0-1,&
                j_0,mx,my,dx,dy)
          delQ1(:,1) = qold2(:,i_0,j_0-1) - qold2(:,i_0,j_0)
          delQ1(:,2) = qold2(:,i_0+1,j_0) - qold2(:,i_0,j_0)
          delQ1(:,3) = qold2(:,i_0,j_0+1) - qold2(:,i_0,j_0)
          delR1(:,:) = dist_for_grad_up(:,1:3,i_0,j_0)
          gradQ2(:,:,i_0,j_0) = nabla_Q(delQ1,delR1,1,dx)
          delQ1(:,1) = qold(:,i_0,j_0-1) -qold(:,i_0,j_0)
          delQ1(:,2) = qold(:,i_0-1,j_0) - qold(:,i_0,j_0)
          delQ1(:,3) = qold(:,i_0,j_0+1) - qold(:,i_0,j_0)
          delR1(:,:) = dist_for_grad_down(:,1:3,i_0,j_0)
          gradQ(:,:,i_0,j_0) = nabla_Q(delQ1,delR1,1,dx)
        case(8)
          gradQ2(:,:,i_0,j_0+1) = nabla_Q_normal(qold2,cen_grid_up,&
              i_0,j_0+1,mx,my,dx,dy)  ! if above
          delQ2(:,1) = qold2(:,i_0-1,j_0) - qold2(:,i_0,j_0)
          delQ2(:,2) = qold2(:,i_0,j_0+1) - qold2(:,i_0,j_0)
          delQ2(:,3) = qold2(:,i_0+1,j_0) - qold2(:,i_0,j_0)
          delQ2(:,4) = qold2(:,i_0,j_0-1) - qold2(:,i_0,j_0)
          delR2(:,:) = dist_for_grad_up(:,1:4,i_0,j_0)
          gradQ2(:,:,i_0,j_0) = nabla_Q(delQ2,delR2,2,dx)
          delQ1(:,1) = qold(:,i_0,j_0-1)-qold(:,i_0,j_0)
          delQ1(:,2) = qold(:,i_0-1,j_0-1)-qold(:,i_0,j_0)
          delQ1(:,3) = qold(:,i_0-1,j_0) -qold(:,i_0,j_0)
          delR1(:,:) = dist_for_grad_down(:,1:3,i_0,j_0)
          gradQ(:,:,i_0,j_0) = nabla_Q(delQ1,delR1,1,dx)
        case(2)
          gradQ(:,:,i_0+1,j_0) = nabla_Q_normal(qold,cen_grid_up,i_0+1,&
               j_0,mx,my,dx,dy)
          gradQ2(:,:,i_0-1,j_0) = nabla_Q_normal(qold2,cen_grid_down,i_0-1,&
                j_0,mx,my,dx,dy)
          delQ1(:,1) = qold(:,i_0,j_0-1) - qold(:,i_0,j_0)
          delQ1(:,2) = qold(:,i_0+1,j_0) - qold(:,i_0,j_0)
          delQ1(:,3) = qold(:,i_0,j_0+1) - qold(:,i_0,j_0)
          delR1(:,:) = dist_for_grad_down(:,1:3,i_0,j_0)
          gradQ(:,:,i_0,j_0) = nabla_Q(delQ1,delR1,1,dx)
          delQ1(:,1) = qold2(:,i_0,j_0-1) -qold2(:,i_0,j_0)
          delQ1(:,2) = qold2(:,i_0-1,j_0) - qold2(:,i_0,j_0)
          delQ1(:,3) = qold2(:,i_0,j_0+1) - qold2(:,i_0,j_0)
          delR1(:,:) = dist_for_grad_up(:,1:3,i_0,j_0)
          gradQ2(:,:,i_0,j_0) = nabla_Q(delQ1,delR1,1,dx)

        end select
      enddo

    end subroutine

    function nabla_Q(delQ,delR,ixy,dx) result(nabla)
      implicit none
      real(8) :: delQ(:,:) , delR(:,:)
      real(8) :: nabla(3,2),dx
      real(8) :: temp(2,2), temp1(3,2),temp2(4,2)
      real(8) :: temp3(2,3),temp4(2,3),W3(3,3),W4(4,4)
      integer :: ixy,k,i,ix,iy

     if (ixy .eq. 1) then
       temp1 = transpose(delR)
      temp = matmul(delR,temp1)
      temp = matinv2(temp)
      temp3 = matmul(delR,transpose(delQ))
      nabla = transpose(matmul(temp,temp3))
    else
      temp2 = transpose(delR)
      temp = matmul(delR,temp2)
      temp = matinv2(temp)
      temp4 = matmul(delR,transpose(delQ))
      nabla = transpose(matmul(temp,temp4))
    end if
  end function nabla_Q

    function nabla_Q_normal(qold,cen_grid,i,j,mx,my,dx,dy) result(nabla)
      implicit none
      integer :: i ,j ,mx,my,k
      real(8) :: qold(3,-1:mx+2,-1:my+2),dx,dy
      real(8) :: cen_grid(2,-1:mx+2,-1:my+2)
      real(8) :: delQ(3,4), delR(2,4) , nabla(3,2),W(4,4)

      delQ(:,1) = qold(:,i-1,j) - qold(:,i,j)
      delQ(:,2) = qold(:,i,j+1) - qold(:,i,j)
      delQ(:,3) = qold(:,i+1,j) - qold(:,i,j)
      delQ(:,4) =qold(:,i,j-1) - qold(:,i,j)
      delR(:,1) = cen_grid(:,i-1,j) -cen_grid(:,i,j)
      delR(:,2) = cen_grid(:,i,j+1) - cen_grid(:,i,j)
      delR(:,3) = cen_grid(:,i+1,j) - cen_grid(:,i,j)
      delR(:,4) =cen_grid(:,i,j-1) - cen_grid(:,i,j)


      nabla = nabla_Q(delQ,delR,2,dx)
    end function nabla_Q_normal

    function matinv2(A) result(B)
      !! Performs a direct calculation of the inverse of a 2×2 matrix.
      real(8), intent(in) :: A(2,2)   !! Matrix
      real(8)             :: B(2,2)   !! Inverse matrix
      real(8)             :: detinv,det

      ! Calculate the inverse determinant of the matrix
      det = (A(1,1)*A(2,2) - A(1,2)*A(2,1))
      detinv = 1/det
      if (det .eq. 0.d0)then
        print*, "DETERMINANTN IS ZREOOOOOO"
      end if
      ! Calculate the inverse of the matrix
      B(1,1) = +detinv * A(2,2)
      B(2,1) = -detinv * A(2,1)
      B(1,2) = -detinv * A(1,2)
      B(2,2) = +detinv * A(1,1)
    end function matinv2

    function matinv3(A) result(B)
    !! Performs a direct calculation of the inverse of a 3×3 matrix.
    real(8), intent(in) :: A(3,3)   !! Matrix
    real(8)             :: B(3,3)   !! Inverse matrix
    real(8)             :: detinv

    ! Calculate the inverse determinant of the matrix
    detinv = 1/(A(1,1)*A(2,2)*A(3,3) - A(1,1)*A(2,3)*A(3,2)&
              - A(1,2)*A(2,1)*A(3,3) + A(1,2)*A(2,3)*A(3,1)&
              + A(1,3)*A(2,1)*A(3,2) - A(1,3)*A(2,2)*A(3,1))

    ! Calculate the inverse of the matrix
    B(1,1) = +detinv * (A(2,2)*A(3,3) - A(2,3)*A(3,2))
    B(2,1) = -detinv * (A(2,1)*A(3,3) - A(2,3)*A(3,1))
    B(3,1) = +detinv * (A(2,1)*A(3,2) - A(2,2)*A(3,1))
    B(1,2) = -detinv * (A(1,2)*A(3,3) - A(1,3)*A(3,2))
    B(2,2) = +detinv * (A(1,1)*A(3,3) - A(1,3)*A(3,1))
    B(3,2) = -detinv * (A(1,1)*A(3,2) - A(1,2)*A(3,1))
    B(1,3) = +detinv * (A(1,2)*A(2,3) - A(1,3)*A(2,2))
    B(2,3) = -detinv * (A(1,1)*A(2,3) - A(1,3)*A(2,1))
    B(3,3) = +detinv * (A(1,1)*A(2,2) - A(1,2)*A(2,1))
  end function


    subroutine rotate_state(q,q_rot,n_vec,t_vec)
      ! n_vec is the normal direction unit vector
      ! t_vec is the transverse direction unit vector, OG to n_vec
      ! q is the Cartesian coordinate aligned state vec
      ! q_rot is the rotated state vec
      implicit none
      real(8) :: q(3),q_rot(3),n_vec(2),t_vec(2)
      real(8) :: vel(2)

      q_rot(1) = q(1)
      vel = q(2:3)
      q_rot(2) = vel(1)*n_vec(1) + vel(2)*n_vec(2)
      q_rot(3) = vel(1)*t_vec(1) + vel(2)*t_vec(2)
    end subroutine











end module
