module globals
  implicit none

  integer :: N, N_plus
  real(kind=8), allocatable :: r(:,:), v(:,:), f(:,:)
  real(kind=8) :: L, sigma, epsilon, gama

  ! Potential / cut
  real(kind=8) :: E_pot
  real(kind=8) :: rc, rc2, V_shift !parametros de corte (se actualizan en force)

  ! Virial / presion
  real(kind=8) :: W_inst   ! virial instantáneo para presión

  ! g(r) histogram (use ngr_bins consistently)
  integer :: ngr_bins
  real(kind=8), allocatable :: gr_hist(:)
  real(kind=8) :: dr, rmax
  integer :: gr_count_frames

  ! Simulation parameters
  real(kind=8) :: deltaT
  real(kind=8) :: deltaTMinimizacion
  integer      :: pasosT, tMinimizacion, muestreo, muestreoMinimizacion
  real(kind=8) :: Ttarget
  real(kind=8) :: pump

  ! constants
  real(kind=8), parameter :: kB = 1.0d0

  ! Fluid-wall interaction parameters
  real(kind=8) :: a_wall(2)  ! Coefficients for bottom (2) and top (1) walls
  real(kind=8) :: sigma_wall ! Interaction parameter for particle-wall
  integer :: a_type          ! Type of particle (only one type in this case)

  ! --- Variables para el Perfil de Velocidad ---
  integer, parameter :: n_capas = 50 ! Cantidad de capas (bins) en Z
  real(kind=8) :: v_sum_x(n_capas)   ! Suma de velocidades Vx en cada capa
  real(kind=8) :: n_count_z(n_capas) ! Contador de partículas en cada capa
  ! --------------------------------------------------


end module globals


