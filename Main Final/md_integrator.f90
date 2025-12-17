module md_integrator
    use globals
    use ziggurat
    implicit none
contains

  !==========================================================
  !inicio posiciones
  !==========================================================
subroutine init_positions()
    implicit none
    integer :: i, j
    do j = 1, N
        do i = 1, 3
            r(i,j) = uni() * L
        end do
    end do

    ! Inicio de posiciones de particulas en paredes
    do j = N, N+N_plus
        do i = 1, 3
            if ( (N_plus/2 + N) > j) then! Partículas de abajo
                if (i == 3) then
                    r(i,j) = 1 * uni()
                else
                    r(i,j) = uni() * L
                end if
            else                  ! Partículas de arriba
                if (i == 3) then
                    r(i,j) = L - 1 * uni()
                else
                    r(i,j) = uni() * L
                end if
            end if
        end do
    end do

    print *, "Posiciones iniciales generadas. L =", L
end subroutine init_positions

  !==========================================================
  !inicio velocidades
  !==========================================================
subroutine init_velocities(Ttarget)
    implicit none
    real(kind=8), intent(in) :: Ttarget
    integer :: i, j, Nf 
    real(kind=8) :: sigma_v, vcm(3), K, Tinst, alpha, v2

    sigma_v = sqrt(Ttarget)

    do j = 1, N
        do i = 1, 3
            v(i,j) = rnor() * sigma_v
        end do
    end do

    !quitar velocidad del centro de masa
    vcm = 0d0
    do j = 1, N
        vcm = vcm + v(:,j)
    end do
    vcm = vcm / real(N,8)
    do j = 1, N
        v(:,j) = v(:,j) - vcm
    end do

    !reescalar para la T seteada
        K = 0.0d0
    do i = 1, N
        v2 = v(1,i)**2 + v(2,i)**2 + v(3,i)**2
        K = K + 0.5d0 * v2
    end do
    !--------------------------
    ! Temperatura
    Nf = 3*N - 3
    Tinst = 2d0*K / real(Nf,8)
    
    if (Tinst > 0d0) then
        alpha = sqrt(Ttarget / Tinst)
        v = alpha * v
    end if

    print *, "Velocidades iniciales listas. T =", Ttarget
end subroutine init_velocities

  !==========================================================
  !condiciones periodicas de contorno
  !==========================================================

subroutine apply_pbc_positions()
  use globals
  implicit none
  integer :: i, j
  
  do j = 1, N+N_plus
        do i = 1, 3
            ! Mover partículas que están fuera por la derecha
            if (r(i,j) > L) r(i,j) = MOD(r(i,j), L)
            ! Mover partículas que están fuera por la izquierda  
            if (r(i,j) < 0.0d0) r(i,j) = MOD(r(i,j), L) + L

        end do
        if ( j > N .and. r(3,j) > 1 .and. r(3,j) < L-1) r(3,j) = 1 
  end do

end subroutine apply_pbc_positions

  !==========================================================
  !Verlet. actualizo posiciones
  !==========================================================
subroutine verlet_positions(deltaT)
    implicit none
    real(kind=8), intent(in) :: deltaT
    integer :: j
    do j = 1, N
        r(:,j) = r(:,j) + v(:,j)*deltaT + 0.5d0*f(:,j)*deltaT*deltaT
    end do
end subroutine verlet_positions

  !==========================================================
  !Verlet. actualizo velocidades
  !==========================================================
subroutine verlet_velocities(deltaT, f_old)
    implicit none
    real(kind=8), intent(in) :: deltaT
    real(kind=8), intent(in) :: f_old(3,N+N_plus)
    integer :: j
    do j = 1, N
        v(:,j) = v(:,j) + 0.5d0*(f_old(:,j) + f(:,j))*deltaT
    end do
end subroutine verlet_velocities

  !=========================================================================================
  !medición de parámetros de interes. Energía cinetica, presion, temperatura y energía total
  !=========================================================================================

subroutine measurements(K, Tinst, P, Etot)
    use globals
    implicit none

    real(kind=8), intent(out) :: K, Tinst, P, Etot
    real(kind=8) :: v2, volumen
    integer :: i, Nf

   
    ! Energia cinética
    K = 0.0d0
    do i = 1, N
        v2 = v(1,i)**2 + v(2,i)**2 + v(3,i)**2
        K = K + 0.5d0 * v2
    end do
    
    ! Temperatura
    Nf = 3*N - 3
    Tinst = 2d0*K / real(Nf,8)
    
    ! Presión
    volumen = L*L*L
    P = (N / volumen) * Tinst + W_inst / (3d0 * volumen)

    ! Energía total
    Etot = K + E_pot

end subroutine measurements


  !==========================================================
  !!preparar muestreo g(r)
  !==========================================================
subroutine init_gr(nbins, rmax_in)        !bins= cantidad de cajones del histograma
    implicit none
    integer, intent(in) :: nbins
    real(kind=8), intent(in) :: rmax_in

    ngr_bins = nbins
    rmax = rmax_in
    dr = rmax / real(ngr_bins,8)

    if (allocated(gr_hist)) deallocate(gr_hist)
    allocate(gr_hist(0:ngr_bins-1))
    gr_hist = 0.0d0
    gr_count_frames = 0
end subroutine init_gr

  !==========================================================
  !!conteo para el histograma de g(r)
  !==========================================================
subroutine accumulate_gr()
    implicit none
    integer :: i, j, k
    real(kind=8) :: deltar(3), rij2, rij
    integer :: bin
    real(kind=8) :: rij_val

    ! contamos pares i<j
    do j = 1, N-1
        do i = j+1, N
            deltar = r(:,i) - r(:,j)
            ! mínima imagen
            do k = 1,2
                if (deltar(k) >  L/2d0) deltar(k) = deltar(k) - L
                if (deltar(k) < -L/2d0) deltar(k) = deltar(k) + L
            end do
            rij2 = deltar(1)**2 + deltar(2)**2 + deltar(3)**2
            rij = sqrt(rij2)
            if (rij < rmax) then
                bin = int(rij / dr)
                if (bin >= 0 .and. bin < ngr_bins) gr_hist(bin) = gr_hist(bin) + 2.0d0
                ! multiplicamos por 2 para contar i and j contributions (standard)
            end if
        end do
    end do

    gr_count_frames = gr_count_frames + 1
end subroutine accumulate_gr


!==========================================================
!salida .dat del g(r)
!==========================================================
subroutine normalize_and_write_gr(filename)
    implicit none
    character(len=*), intent(in) :: filename
    integer :: bin
    real(kind=8) :: rho, V, norm, r_mid, shell_vol
    integer :: count_int

    V = L**3
    rho = real(N,8) / V

    open(unit=90, file='output/gr.dat', status='unknown')

    write(90,'(A)') '# r   g(r)   counts'

    do bin = 0, ngr_bins-1
        r_mid = ( (real(bin,8) + 0.5d0) * dr )
        shell_vol = 4d0 * 3.141592653589793d0 * r_mid*r_mid * dr

        count_int = int(gr_hist(bin))     ! <--- convertir a integer

        if (gr_count_frames > 0) then
           norm = real(gr_count_frames,8) * rho * shell_vol * real(N,8)
           if (norm > 0d0) then
              write(90,*) r_mid, gr_hist(bin)/norm, count_int
           else
              write(90,*) r_mid, 0d0, count_int
           end if
        else
           write(90, *) r_mid, 0d0, count_int
        end if
    end do

    close(90)
end subroutine normalize_and_write_gr

!==========================================================
! Calcular perfil de velocidad (Versión Estática 100% limpia)
!==========================================================
subroutine velocity_prof_simple(mode)
    use globals
    implicit none
    integer, intent(in) :: mode
    
    integer :: j, capa
    real(kind=8) :: dz, z_actual, promedio_vx
    real(kind=8) :: densidad            !!!! AGREGADO
  

    ! Nota: Usamos 'n_capas' que ya viene definido como parameter en globals.f90

    select case (mode)
    
    ! --- MODO 0: INICIALIZACIÓN ---
    case(0)
        v_sum_x   = 0.0d0
        n_count_z = 0.0d0
        
        print *, " * Perfil de velocidad inicializado (", n_capas, " capas)."

    ! --- MODO 1: ACUMULACIÓN ---
    case(1)
        dz = L / real(n_capas, 8)

        do j = 1, N
            z_actual = r(3,j) 
            capa = int(z_actual / dz) + 1

            ! Blindaje de bordes
            if (capa < 1) capa = 1
            if (capa > n_capas) capa = n_capas
            
            v_sum_x(capa)   = v_sum_x(capa) + v(1,j)
            n_count_z(capa) = n_count_z(capa) + 1.0d0
        end do

    ! --- MODO 2: ESCRITURA ---
    case(2)
        print *, "--> Guardando perfil..." 
        
        dz   = L / real(n_capas, 8)
                         
        
        open(unit=74, file='output/perfil_velocidad.dat', status='unknown')
        open(unit=75, file='output/perfil_densidad.dat',  status='unknown')  !!!! AGREGADO

        write(74,*) '# Z_centro    Vx_promedio'
        write(75,*) '# Z_centro    densidad'   !!!! AGREGADO

        do j = 1, n_capas
            if (n_count_z(j) > 0.0d0) then
                promedio_vx = v_sum_x(j) / n_count_z(j)
                densidad    = n_count_z(j) / ( L * L   * dz)   !!!! AGREGADO
            else
                promedio_vx = 0.0d0
                densidad    = 0.0d0                         !!!! AGREGADO
            end if

            write(74, *) (dble(j)-0.5d0)*dz, promedio_vx, n_count_z(j)
            write(75, *) (dble(j)-0.5d0)*dz, densidad, n_count_z(j)       !!!! AGREGADO
        end do
        
        close(74)
        close(75)                         !!!! AGREGADO
        print *, "--> Archivos guardados." 
        
    end select

end subroutine velocity_prof_simple

end module md_integrator



