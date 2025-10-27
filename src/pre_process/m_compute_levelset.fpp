!>
!! @file m_compute_levelset.fpp
!! @brief Contains module m_compute_levelset

#:include 'macros.fpp'

!> @brief This module is used to handle all operations related to immersed
!!              boundary methods (IBMs)
module m_compute_levelset

    use m_derived_types        !< Definitions of the derived types

    use m_global_parameters    !< Definitions of the global parameters

    use m_mpi_proxy            !< Message passing interface (MPI) module proxy

    use m_helper_basic         !< Functions to compare floating point numbers

    implicit none

    private; public :: s_cylinder_levelset, s_circle_levelset, &
 s_airfoil_levelset, &
 s_3D_airfoil_levelset, &
 s_rectangle_levelset, &
 s_cuboid_levelset, &
 s_sphere_levelset, &
 s_triangle_levelset

contains

    pure subroutine s_circle_levelset(ib_patch_id, levelset, levelset_norm)

        type(levelset_field), intent(INOUT), optional :: levelset
        type(levelset_norm_field), intent(INOUT), optional :: levelset_norm
        integer, intent(IN) :: ib_patch_id

        real(wp) :: radius, dist
        real(wp) :: x_centroid, y_centroid
        real(wp), dimension(3) :: dist_vec

        integer :: i, j !< Loop index variables

        radius = patch_ib(ib_patch_id)%radius
        x_centroid = patch_ib(ib_patch_id)%x_centroid
        y_centroid = patch_ib(ib_patch_id)%y_centroid

        do i = 0, m
            do j = 0, n

                dist_vec(1) = x_cc(i) - x_centroid
                dist_vec(2) = y_cc(j) - y_centroid
                dist_vec(3) = 0
                dist = sqrt(sum(dist_vec**2))
                levelset%sf(i, j, 0, ib_patch_id) = dist - radius
                if (f_approx_equal(dist, 0._wp)) then
                    levelset_norm%sf(i, j, 0, ib_patch_id, :) = 0
                else
                    levelset_norm%sf(i, j, 0, ib_patch_id, :) = &
                        dist_vec(:)/dist
                end if

            end do
        end do

    end subroutine s_circle_levelset

    pure subroutine s_airfoil_levelset(ib_patch_id, levelset, levelset_norm)

        type(levelset_field), intent(INOUT), optional :: levelset
        type(levelset_norm_field), intent(INOUT), optional :: levelset_norm
        integer, intent(IN) :: ib_patch_id

        real(wp) :: dist, global_dist
        integer :: global_id
        real(wp) :: x_centroid, y_centroid, x_act, y_act, theta
        real(wp), dimension(3) :: dist_vec

        integer :: i, j, k !< Loop index variables

        x_centroid = patch_ib(ib_patch_id)%x_centroid
        y_centroid = patch_ib(ib_patch_id)%y_centroid
        theta = pi*patch_ib(ib_patch_id)%theta/180._wp

        do i = 0, m
            do j = 0, n

                if (.not. f_is_default(patch_ib(ib_patch_id)%theta)) then
                    x_act = (x_cc(i) - x_centroid)*cos(theta) - (y_cc(j) - y_centroid)*sin(theta) + x_centroid
                    y_act = (x_cc(i) - x_centroid)*sin(theta) + (y_cc(j) - y_centroid)*cos(theta) + y_centroid
                else
                    x_act = x_cc(i)
                    y_act = y_cc(j)
                end if

                if (y_act >= y_centroid) then
                    do k = 1, Np
                        dist_vec(1) = x_cc(i) - airfoil_grid_u(k)%x
                        dist_vec(2) = y_cc(j) - airfoil_grid_u(k)%y
                        dist_vec(3) = 0
                        dist = sqrt(sum(dist_vec**2))
                        if (k == 1) then
                            global_dist = dist
                            global_id = k
                        else
                            if (dist < global_dist) then
                                global_dist = dist
                                global_id = k
                            end if
                        end if
                    end do
                    dist_vec(1) = x_cc(i) - airfoil_grid_u(global_id)%x
                    dist_vec(2) = y_cc(j) - airfoil_grid_u(global_id)%y
                    dist_vec(3) = 0
                    dist = global_dist
                else
                    do k = 1, Np
                        dist_vec(1) = x_cc(i) - airfoil_grid_l(k)%x
                        dist_vec(2) = y_cc(j) - airfoil_grid_l(k)%y
                        dist_vec(3) = 0
                        dist = sqrt(sum(dist_vec**2))
                        if (k == 1) then
                            global_dist = dist
                            global_id = k
                        else
                            if (dist < global_dist) then
                                global_dist = dist
                                global_id = k
                            end if
                        end if
                    end do
                    dist_vec(1) = x_cc(i) - airfoil_grid_l(global_id)%x
                    dist_vec(2) = y_cc(j) - airfoil_grid_l(global_id)%y
                    dist_vec(3) = 0
                    dist = global_dist
                end if

                levelset%sf(i, j, 0, ib_patch_id) = dist
                if (f_approx_equal(dist, 0._wp)) then
                    levelset_norm%sf(i, j, 0, ib_patch_id, :) = 0
                else
                    levelset_norm%sf(i, j, 0, ib_patch_id, :) = &
                        dist_vec(:)/dist
                end if

            end do
        end do

    end subroutine s_airfoil_levelset

    pure subroutine s_3D_airfoil_levelset(ib_patch_id, levelset, levelset_norm)

        type(levelset_field), intent(INOUT), optional :: levelset
        type(levelset_norm_field), intent(INOUT), optional :: levelset_norm
        integer, intent(IN) :: ib_patch_id

        real(wp) :: dist, dist_surf, dist_side, global_dist
        integer :: global_id
        real(wp) :: x_centroid, y_centroid, z_centroid, lz, z_max, z_min, x_act, y_act, theta
        real(wp), dimension(3) :: dist_vec

        real(wp) :: length_z

        integer :: i, j, k, l !< Loop index variables

        x_centroid = patch_ib(ib_patch_id)%x_centroid
        y_centroid = patch_ib(ib_patch_id)%y_centroid
        z_centroid = patch_ib(ib_patch_id)%z_centroid
        lz = patch_ib(ib_patch_id)%length_z
        theta = pi*patch_ib(ib_patch_id)%theta/180._wp

        z_max = z_centroid + lz/2
        z_min = z_centroid - lz/2

        do l = 0, p
            do j = 0, n
                do i = 0, m

                    if (.not. f_is_default(patch_ib(ib_patch_id)%theta)) then
                        x_act = (x_cc(i) - x_centroid)*cos(theta) - (y_cc(j) - y_centroid)*sin(theta) + x_centroid
                        y_act = (x_cc(i) - x_centroid)*sin(theta) + (y_cc(j) - y_centroid)*cos(theta) + y_centroid
                    else
                        x_act = x_cc(i)
                        y_act = y_cc(j)
                    end if

                    if (y_act >= y_centroid) then
                        do k = 1, Np
                            dist_vec(1) = x_cc(i) - airfoil_grid_u(k)%x
                            dist_vec(2) = y_cc(j) - airfoil_grid_u(k)%y
                            dist_vec(3) = 0
                            dist_surf = sqrt(sum(dist_vec**2))
                            if (k == 1) then
                                global_dist = dist_surf
                                global_id = k
                            else
                                if (dist_surf < global_dist) then
                                    global_dist = dist_surf
                                    global_id = k
                                end if
                            end if
                        end do
                        dist_vec(1) = x_cc(i) - airfoil_grid_u(global_id)%x
                        dist_vec(2) = y_cc(j) - airfoil_grid_u(global_id)%y
                        dist_vec(3) = 0
                        dist_surf = global_dist
                    else
                        do k = 1, Np
                            dist_vec(1) = x_cc(i) - airfoil_grid_l(k)%x
                            dist_vec(2) = y_cc(j) - airfoil_grid_l(k)%y
                            dist_vec(3) = 0
                            dist_surf = sqrt(sum(dist_vec**2))
                            if (k == 1) then
                                global_dist = dist_surf
                                global_id = k
                            else
                                if (dist_surf < global_dist) then
                                    global_dist = dist_surf
                                    global_id = k
                                end if
                            end if
                        end do
                        dist_vec(1) = x_cc(i) - airfoil_grid_l(global_id)%x
                        dist_vec(2) = y_cc(j) - airfoil_grid_l(global_id)%y
                        dist_vec(3) = 0
                        dist_surf = global_dist
                    end if

                    dist_side = min(abs(z_cc(l) - z_min), abs(z_max - z_cc(l)))

                    if (dist_side < dist_surf) then
                        levelset%sf(i, j, l, ib_patch_id) = dist_side
                        if (f_approx_equal(dist_side, abs(z_cc(l) - z_min))) then
                            levelset_norm%sf(i, j, l, ib_patch_id, :) = (/0, 0, -1/)
                        else
                            levelset_norm%sf(i, j, l, ib_patch_id, :) = (/0, 0, 1/)
                        end if
                    else
                        levelset%sf(i, j, l, ib_patch_id) = dist_surf
                        if (f_approx_equal(dist, 0._wp)) then
                            levelset_norm%sf(i, j, l, ib_patch_id, :) = 0
                        else
                            levelset_norm%sf(i, j, l, ib_patch_id, :) = &
                                dist_vec(:)/dist_surf
                        end if
                    end if

                end do
            end do
        end do

    end subroutine s_3D_airfoil_levelset

    !>  Initialize IBM module
    pure subroutine s_rectangle_levelset(ib_patch_id, levelset, levelset_norm)

        type(levelset_field), intent(INOUT), optional :: levelset
        type(levelset_norm_field), intent(INOUT), optional :: levelset_norm

        integer, intent(in) :: ib_patch_id
        real(wp) :: top_right(2), bottom_left(2)
        real(wp) :: x, y, min_dist
        real(wp) :: side_dists(4)

        real(wp) :: x_centroid, y_centroid
        real(wp) :: length_x, length_y

        integer :: i, j, k !< Loop index variables
        integer :: idx !< Shortest path direction indicator

        length_x = patch_ib(ib_patch_id)%length_x
        length_y = patch_ib(ib_patch_id)%length_y
        x_centroid = patch_ib(ib_patch_id)%x_centroid
        y_centroid = patch_ib(ib_patch_id)%y_centroid

        top_right(1) = x_centroid + length_x/2
        top_right(2) = y_centroid + length_y/2

        bottom_left(1) = x_centroid - length_x/2
        bottom_left(2) = y_centroid - length_y/2

        do i = 0, m
            do j = 0, n

                x = x_cc(i)
                y = y_cc(j)

                if ((x > bottom_left(1) .and. x < top_right(1)) .or. &
                    (y > bottom_left(2) .and. y < top_right(2))) then

                    side_dists(1) = bottom_left(1) - x
                    side_dists(2) = x - top_right(1)
                    side_dists(3) = bottom_left(2) - y
                    side_dists(4) = y - top_right(2)
                    min_dist = initial_distance_buffer
                    idx = 1

                    do k = 1, 4
                        if (abs(side_dists(k)) < abs(min_dist)) then
                            idx = k
                            min_dist = side_dists(idx)
                        end if
                    end do

                    if (idx == 1) then
                        levelset%sf(i, j, 0, ib_patch_id) = side_dists(1)
                        if (f_approx_equal(side_dists(1), 0._wp)) then
                            levelset_norm%sf(i, j, 0, ib_patch_id, 1) = 0._wp
                        else
                            levelset_norm%sf(i, j, 0, ib_patch_id, 1) = side_dists(1)/ &
                                                                        abs(side_dists(1))
                        end if

                    else if (idx == 2) then
                        levelset%sf(i, j, 0, ib_patch_id) = side_dists(2)
                        if (f_approx_equal(side_dists(2), 0._wp)) then
                            levelset_norm%sf(i, j, 0, ib_patch_id, 1) = 0._wp
                        else
                            levelset_norm%sf(i, j, 0, ib_patch_id, 1) = -side_dists(2)/ &
                                                                        abs(side_dists(2))
                        end if

                    else if (idx == 3) then
                        levelset%sf(i, j, 0, ib_patch_id) = side_dists(3)
                        if (f_approx_equal(side_dists(3), 0._wp)) then

                            levelset_norm%sf(i, j, 0, ib_patch_id, 2) = 0._wp
                        else
                            levelset_norm%sf(i, j, 0, ib_patch_id, 2) = side_dists(3)/ &
                                                                        abs(side_dists(3))
                        end if

                    else if (idx == 4) then
                        levelset%sf(i, j, 0, ib_patch_id) = side_dists(4)
                        if (f_approx_equal(side_dists(4), 0._wp)) then

                            levelset_norm%sf(i, j, 0, ib_patch_id, 2) = 0._wp
                        else
                            levelset_norm%sf(i, j, 0, ib_patch_id, 2) = -side_dists(4)/ &
                                                                        abs(side_dists(4))
                        end if

                    end if

                end if

            end do
        end do

    end subroutine s_rectangle_levelset

    pure subroutine s_triangle_levelset(ib_patch_id, levelset, levelset_norm)
        ! use mod_precision
        implicit none

        type(levelset_field), intent(INOUT), optional :: levelset
        type(levelset_norm_field), intent(INOUT), optional :: levelset_norm
        integer, intent(in) :: ib_patch_id

        ! Triangle vertices
        real(wp) :: x1, y1, x2, y2, x3, y3
        real(wp) :: a(3), b(3), c(3), norm(3)
        real(wp) :: x, y, d(3), min_dist
        real(wp) :: length_x, length_y, x_corner, y_corner
        integer :: i, j, k, idx

        ! Fetch triangle data from patch
        length_x = patch_ib(ib_patch_id)%length_x
        length_y = patch_ib(ib_patch_id)%length_y
        x_corner = patch_ib(ib_patch_id)%x_centroid
        y_corner = patch_ib(ib_patch_id)%y_centroid

        !> Bottom point
        x1 = x_corner
        y1 = y_corner - length_y
        !> Leading point
        x2 = x_corner - length_x
        y2 = y_corner
        !> Top point
        x3 = x_corner
        y3 = y_corner + length_y

        ! Compute line coefficients for edges
        a(1) = y2 - y3;  b(1) = x3 - x2;  c(1) = x2*y3 - x3*y2
        a(2) = y3 - y1;  b(2) = x1 - x3;  c(2) = x3*y1 - x1*y3
        a(3) = y1 - y2;  b(3) = x2 - x1;  c(3) = x1*y2 - x2*y1

        ! Precompute edge normal magnitudes
        norm = sqrt(a**2 + b**2)

        ! Loop over grid
        do j = 0, n
            do i = 0, m
                x = x_cc(i)
                y = y_cc(j)

                ! Signed distances to each edge
                do k = 1, 3
                    d(k) = (a(k)*x + b(k)*y + c(k)) / norm(k)
                end do

                ! Pick smallest absolute distance
                min_dist = abs(d(1))
                idx = 1
                do k = 2, 3
                    if (abs(d(k)) < abs(min_dist)) then
                        min_dist = d(k)
                        idx = k
                    end if
                end do

                ! Assign levelset signed distance
                levelset%sf(i, j, 0, ib_patch_id) = d(idx)

                ! Assign corresponding normal components
                if (present(levelset_norm)) then
                    if (norm(idx) > 0._wp) then
                        levelset_norm%sf(i, j, 0, ib_patch_id, 1) = a(idx) / norm(idx)
                        levelset_norm%sf(i, j, 0, ib_patch_id, 2) = b(idx) / norm(idx)
                    else
                        levelset_norm%sf(i, j, 0, ib_patch_id, 1) = 0._wp
                        levelset_norm%sf(i, j, 0, ib_patch_id, 2) = 0._wp
                    end if
                end if

            end do
        end do

    end subroutine s_triangle_levelset

        
        
    pure subroutine s_cuboid_levelset(ib_patch_id, levelset, levelset_norm)

        type(levelset_field), intent(INOUT), optional :: levelset
        type(levelset_norm_field), intent(INOUT), optional :: levelset_norm

        integer, intent(IN) :: ib_patch_id
        real(wp) :: Right, Left, Bottom, Top, Front, Back
        real(wp) :: x, y, z, min_dist
        real(wp) :: side_dists(6)

        real(wp) :: x_centroid, y_centroid, z_centroid
        real(wp) :: length_x, length_y, length_z

        integer :: i, j, k !< Loop index variables

        length_x = patch_ib(ib_patch_id)%length_x
        length_y = patch_ib(ib_patch_id)%length_y
        length_z = patch_ib(ib_patch_id)%length_z

        x_centroid = patch_ib(ib_patch_id)%x_centroid
        y_centroid = patch_ib(ib_patch_id)%y_centroid
        z_centroid = patch_ib(ib_patch_id)%z_centroid

        Right = x_centroid + length_x/2
        Left = x_centroid - length_x/2

        Top = y_centroid + length_y/2
        Bottom = y_centroid - length_y/2

        Front = z_centroid + length_z/2
        Back = z_centroid - length_z/2

        do i = 0, m
            do j = 0, n
                do k = 0, p

                    x = x_cc(i)
                    y = y_cc(j)
                    z = z_cc(k)

                    if ((x > Left .and. x < Right) .or. &
                        (y > Bottom .and. y < Top) .or. &
                        (z > Back .and. z < Front)) then

                        side_dists(1) = Left - x
                        side_dists(2) = x - Right
                        side_dists(3) = Bottom - y
                        side_dists(4) = y - Top
                        side_dists(5) = Back - z
                        side_dists(6) = z - Front

                        min_dist = minval(abs(side_dists))

                        if (f_approx_equal(min_dist, abs(side_dists(1)))) then
                            levelset%sf(i, j, k, ib_patch_id) = side_dists(1)
                            if (f_approx_equal(side_dists(1), 0._wp)) then
                                levelset_norm%sf(i, j, k, ib_patch_id, 1) = 0._wp
                            else
                                levelset_norm%sf(i, j, k, ib_patch_id, 1) = side_dists(1)/ &
                                                                            abs(side_dists(1))
                            end if

                        else if (f_approx_equal(min_dist, abs(side_dists(2)))) then
                            levelset%sf(i, j, k, ib_patch_id) = side_dists(2)
                            if (f_approx_equal(side_dists(2), 0._wp)) then
                                levelset_norm%sf(i, j, k, ib_patch_id, 1) = 0._wp
                            else
                                levelset_norm%sf(i, j, k, ib_patch_id, 1) = -side_dists(2)/ &
                                                                            abs(side_dists(2))
                            end if

                        else if (f_approx_equal(min_dist, abs(side_dists(3)))) then
                            levelset%sf(i, j, k, ib_patch_id) = side_dists(3)
                            if (f_approx_equal(side_dists(3), 0._wp)) then
                                levelset_norm%sf(i, j, k, ib_patch_id, 2) = 0._wp
                            else
                                levelset_norm%sf(i, j, k, ib_patch_id, 2) = side_dists(3)/ &
                                                                            abs(side_dists(3))
                            end if

                        else if (f_approx_equal(min_dist, abs(side_dists(4)))) then
                            levelset%sf(i, j, k, ib_patch_id) = side_dists(4)
                            if (f_approx_equal(side_dists(4), 0._wp)) then
                                levelset_norm%sf(i, j, k, ib_patch_id, 2) = 0._wp
                            else
                                levelset_norm%sf(i, j, k, ib_patch_id, 2) = -side_dists(4)/ &
                                                                            abs(side_dists(4))
                            end if

                        else if (f_approx_equal(min_dist, abs(side_dists(5)))) then
                            levelset%sf(i, j, k, ib_patch_id) = side_dists(5)
                            if (f_approx_equal(side_dists(5), 0._wp)) then
                                levelset_norm%sf(i, j, k, ib_patch_id, 3) = 0._wp
                            else
                                levelset_norm%sf(i, j, k, ib_patch_id, 3) = side_dists(5)/ &
                                                                            abs(side_dists(5))
                            end if

                        else if (f_approx_equal(min_dist, abs(side_dists(6)))) then
                            levelset%sf(i, j, k, ib_patch_id) = side_dists(6)
                            if (f_approx_equal(side_dists(6), 0._wp)) then
                                levelset_norm%sf(i, j, k, ib_patch_id, 3) = 0._wp
                            else
                                levelset_norm%sf(i, j, k, ib_patch_id, 3) = -side_dists(6)/ &
                                                                            abs(side_dists(6))
                            end if

                        end if

                    end if

                end do
            end do
        end do

    end subroutine s_cuboid_levelset

    pure subroutine s_sphere_levelset(ib_patch_id, levelset, levelset_norm)

        type(levelset_field), intent(INOUT), optional :: levelset
        type(levelset_norm_field), intent(INOUT), optional :: levelset_norm
        integer, intent(IN) :: ib_patch_id

        real(wp) :: radius, dist
        real(wp) :: x_centroid, y_centroid, z_centroid
        real(wp), dimension(3) :: dist_vec

        integer :: i, j, k !< Loop index variables

        radius = patch_ib(ib_patch_id)%radius
        x_centroid = patch_ib(ib_patch_id)%x_centroid
        y_centroid = patch_ib(ib_patch_id)%y_centroid
        z_centroid = patch_ib(ib_patch_id)%z_centroid

        do i = 0, m
            do j = 0, n
                do k = 0, p
                    dist_vec(1) = x_cc(i) - x_centroid
                    dist_vec(2) = y_cc(j) - y_centroid
                    dist_vec(3) = z_cc(k) - z_centroid
                    dist = sqrt(sum(dist_vec**2))
                    levelset%sf(i, j, k, ib_patch_id) = dist - radius
                    if (f_approx_equal(dist, 0._wp)) then
                        levelset_norm%sf(i, j, k, ib_patch_id, :) = (/1, 0, 0/)
                    else
                        levelset_norm%sf(i, j, k, ib_patch_id, :) = &
                            dist_vec(:)/dist
                    end if
                end do
            end do
        end do

    end subroutine s_sphere_levelset

    pure subroutine s_cylinder_levelset(ib_patch_id, levelset, levelset_norm)

        type(levelset_field), intent(INOUT), optional :: levelset
        type(levelset_norm_field), intent(INOUT), optional :: levelset_norm
        integer, intent(IN) :: ib_patch_id

        real(wp) :: radius
        real(wp) :: x_centroid, y_centroid, z_centroid
        real(wp) :: length_x, length_y, length_z
        real(wp), dimension(3) :: pos_vec, centroid_vec, dist_sides_vec, dist_surface_vec
        real(wp) :: dist_side, dist_surface, side_pos
        type(bounds_info) :: boundary
        integer :: i, j, k !< Loop index variables

        radius = patch_ib(ib_patch_id)%radius
        x_centroid = patch_ib(ib_patch_id)%x_centroid
        y_centroid = patch_ib(ib_patch_id)%y_centroid
        z_centroid = patch_ib(ib_patch_id)%z_centroid
        length_x = patch_ib(ib_patch_id)%length_x
        length_y = patch_ib(ib_patch_id)%length_y
        length_z = patch_ib(ib_patch_id)%length_z

        if (.not. f_approx_equal(length_x, 0._wp)) then
            boundary%beg = x_centroid - 0.5_wp*length_x
            boundary%end = x_centroid + 0.5_wp*length_x
            dist_sides_vec = (/1, 0, 0/)
            dist_surface_vec = (/0, 1, 1/)
        else if (.not. f_approx_equal(length_y, 0._wp)) then
            boundary%beg = y_centroid - 0.5_wp*length_y
            boundary%end = y_centroid + 0.5_wp*length_y
            dist_sides_vec = (/0, 1, 0/)
            dist_surface_vec = (/1, 0, 1/)
        else if (.not. f_approx_equal(length_z, 0._wp)) then
            boundary%beg = z_centroid - 0.5_wp*length_z
            boundary%end = z_centroid + 0.5_wp*length_z
            dist_sides_vec = (/0, 0, 1/)
            dist_surface_vec = (/1, 1, 0/)
        end if

        do i = 0, m
            do j = 0, n
                do k = 0, p
                    pos_vec = [x_cc(i), y_cc(j), z_cc(k)]
                    centroid_vec = [x_centroid, y_centroid, z_centroid]
                    side_pos = dot_product(pos_vec, dist_sides_vec)
                    dist_side = min(abs(side_pos - boundary%beg), &
                                    abs(boundary%end - side_pos))

                    dist_surface = norm2((pos_vec - centroid_vec)*dist_surface_vec) &
                                   - radius

                    if (dist_side < abs(dist_surface)) then
                        levelset%sf(i, j, k, ib_patch_id) = -dist_side
                        if (f_approx_equal(dist_side, abs(side_pos - boundary%beg))) then
                            levelset_norm%sf(i, j, k, ib_patch_id, :) = -dist_sides_vec
                        else
                            levelset_norm%sf(i, j, k, ib_patch_id, :) = dist_sides_vec
                        end if
                    else
                        levelset%sf(i, j, k, ib_patch_id) = dist_surface

                        levelset_norm%sf(i, j, k, ib_patch_id, :) = &
                            (pos_vec - centroid_vec)*dist_surface_vec
                        levelset_norm%sf(i, j, k, ib_patch_id, :) = &
                            levelset_norm%sf(i, j, k, ib_patch_id, :)/ &
                            norm2(levelset_norm%sf(i, j, k, ib_patch_id, :))
                    end if
                end do
            end do
        end do

    end subroutine s_cylinder_levelset

end module m_compute_levelset
