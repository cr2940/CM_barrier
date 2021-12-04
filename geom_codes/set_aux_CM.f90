program check_aux_SRD

use aux_module_CM

implicit none
integer, parameter :: mx=200,my=200,mbc=2,order=1
real(8) :: x_0,y_0,x_e,y_e,dx,dy,xlower,ylower,intersections(2,mx*4)
real(8) :: xupper,yupper
integer :: ii(mx*4),jj(mx*4),N,type_supper(mx*4),type_sunder(mx*4),m
real(8) :: lengths_supper(5,4*mx), lengths_sunder(5,4*mx)
real(8) :: area_supper(4*mx), area_sunder(4*mx)
integer :: all_undercells_i(mx*4), all_undercells_j(mx*4)
integer :: unS_cells_i(mx*4), unS_cells_j(mx*4)
integer :: upS_cells_i(mx*4), upS_cells_j(mx*4) ! THESE will tell you who the neighbor is for small cells
integer :: all_uppercells_i(mx*4), all_uppercells_j(mx*4),k_count_un,k_count_up
integer :: N_ij_un(1-mbc:mx+mbc,1-mbc:my+mbc),N_ij_up(1-mbc:mx+mbc,1-mbc:my+mbc)
integer :: num_frags_d(mx*4),index_frags_d(2,4,mx*4)
integer :: num_frags_u(mx*4),index_frags_u(2,4,mx*4)
real(8) :: hbox_areas_u(mx*4), area_frags_u(4,mx*4)
real(8) :: hbox_areas_d(mx*4), area_frags_d(4,mx*4)
real(8) :: up_area_ij(1-mbc:mx+mbc,1-mbc:my+mbc),un_area_ij(1-mbc:mx+mbc,1-mbc:my+mbc)
real(8) :: ecen_un_y(2,3,4*mx), ecen_up_y(2,3,4*mx)
real(8) :: ecen_un_x(2,2,4*mx), ecen_up_x(2,2,4*mx)
real(8) :: cen_grid_up(2,4*mx), cen_grid_down(2,4*mx)
real(8) :: dist_to_ecen_up(2,6,4*mx), dist_to_ecen_down(2,6,4*mx) ! the delta r vectors you will multiply to grad to get 2nd ord reconstructed vals, some have 6 MC edges, some 4
real(8) :: dist_for_grad_up(2,4,4*mx) ! the R^* in notes, stencil for grad approximation, some have 3 some 4
real(8) :: dist_for_grad_down(2,4,4*mx)
x_0 = 0.0d0
x_e = 1.d0!0.995d0
y_0 = 0.3d0!393d0!0.005d0
y_e = 0.653333d0!9999d0!6533d0!39d0!0.776d0!1.d0
dx = 1.d0/200.d0!0666666d0/6.d0!/4.d0
dy = 1.d0/200.d0!0.01d0!0666666d0/6.d0!/4.d0
xlower=0.d0
xupper = 1.d0
ylower=0.d0
yupper = 1.d0
call cut_cells_find(x_0,y_0,x_e,y_e,dx,dy,xlower,ylower,mx,my,ii,jj,intersections,N)
call small_cells_geom(ii(1:N),jj(1:N),intersections(:,1:N+1),mx,my,dx,dy,xlower,ylower,&
N,type_supper,type_sunder,lengths_supper,lengths_sunder,area_supper,area_sunder,&
  dist_to_ecen_up(:,:,1:N),dist_to_ecen_down(:,:,1:N),cen_grid_up,cen_grid_down,&
    dist_for_grad_up,dist_for_grad_down,ecen_un_y,ecen_up_y,&
    ecen_up_x,ecen_un_x)
call SRD_undercells(N,ii(1:N),jj(1:N),mx,my,area_sunder(1:N),x_0,y_0,x_e,y_e,dx,dy,&
  unS_cells_i,unS_cells_j,N_ij_un,all_undercells_i,all_undercells_j,k_count_un)
call SRD_uppercells(N,ii(1:N),jj(1:N),mx,my,area_supper(1:N),x_0,y_0,x_e,y_e,dx,dy,&
    upS_cells_i,upS_cells_j,N_ij_up,all_uppercells_i,all_uppercells_j,k_count_up)
call area_cells(ii(1:N),jj(1:N),mx,my,mbc,area_sunder(1:N),area_supper(1:N),up_area_ij,un_area_ij,.false.)

! try writing the data into a text file:
open (unit=1,file = "./small_cells_data.txt")
write (1,*) N
write (1,*) ii(1:N)
write (1,*) jj(1:N)
write (1,*) intersections(:,1:N+1)
write (1,*) type_supper(1:N)
write (1,*) type_sunder(1:N)
write (1,*) lengths_supper(:,1:N)
write (1,*) lengths_sunder(:,1:N)
write (1,*) unS_cells_i(1:N)
write (1,*) unS_cells_j(1:N)
write (1,*) upS_cells_i(1:N)
write (1,*) upS_cells_j(1:N)
write (1,*) area_supper(1:N)
write (1,*) area_sunder(1:N)
write (1,*) xlower
write (1,*) ylower
write (1,*) xupper
write (1,*) yupper
write (1,*) x_0
write (1,*) y_0
write (1,*) x_e
write (1,*) y_e
write (1,*) up_area_ij
write (1,*) un_area_ij
if (order .eq. 2) then
write (1,*) dist_to_ecen_up(:,:,1:N)
write (1,*) dist_to_ecen_down(:,:,1:N)
write (1,*) ecen_un_x(:,:,1:N)
write (1,*) ecen_un_y(:,:,1:N)
write (1,*) dist_for_grad_up
write (1,*) dist_for_grad_down
end if
ENDFILE 1
REWIND 1
close (1,status="keep")

! read(1,*) M
! read(1,*) ii(1:N)
print*, "The number of cut cells is: ", N
print*, "Cut cells i-index is given by:", ii(1:N)
print*, "Cut cells j-index is given by:", jj(1:N)
print*, "the intersections x are :", intersections(1,1:N+1)
print*, "the intersections y are :", intersections(2,1:N+1)
print*, "TYPES:UP",type_supper(1:N)
print*, "TYPES:UN",type_sunder(1:N)
print*, "Length ratio to dx are: UP", lengths_supper(:,1:N)
print*, "Length ratio to dx are: UN", lengths_sunder(:,1:N)
print*, "AREAs : UP", up_area_ij(ii(3),jj(3))
print*, "AREAs ; UN", un_area_ij(ii(3),jj(3))
print*, "TOTAL AREA :", area_supper(1:N) + area_sunder(1:N)
end program
