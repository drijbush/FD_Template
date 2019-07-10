
Program testit

  Use numbers_module        , Only : wp
  Use FD_Laplacian_3d_module, Only : FD_Laplacian_3D

  Implicit None
  
  Type( FD_Laplacian_3d ) :: FD
  
  Real( wp ), Dimension( :, :, : ), Allocatable :: grid
  Real( wp ), Dimension( :, :, : ), Allocatable :: laplacian
  Real( wp ), Dimension( :, :, : ), Allocatable :: fd_laplacian

  Real( wp ), Dimension( 1:3, 1:3 ) :: grid_vecs

  Real( wp ), Dimension( 1:3 ) :: r, ri
  
  Real( wp ) :: alpha, alpha_sq
  Real( wp ) :: x, y, z
  Real( wp ) :: norm
  Real( wp ) :: arg
  Real( wp ) :: gauss, gauss_x2, gauss_y2, gauss_z2

  Integer, Dimension( 1:3 ) :: ng, nt
  
  Integer :: order
  Integer :: i3, i2, i1

  Write( *, * ) 'Grid vecs?'
  Read ( *, * ) grid_vecs

  Write( *, * ) 'ri ?'
  Read ( *, * ) ri

  Write( *, * ) 'alpha ?'
  Read ( *, * ) alpha
  alpha_sq = alpha * alpha
  
  Write( *, * ) 'order ?'
  Read ( *, * ) order

  Write( *, * ) 'ng ?'
  Read ( *, * ) ng
  
  nt = ng + order / 2

  Allocate(         grid( -nt( 1 ):nt( 1 ), -nt( 2 ):nt( 2 ), -nt( 3 ):nt( 3 ) ) )
  Allocate(    laplacian( -nt( 1 ):nt( 1 ), -nt( 2 ):nt( 2 ), -nt( 3 ):nt( 3 ) ) )
  Allocate( fd_laplacian( -ng( 1 ):ng( 1 ), -ng( 2 ):ng( 2 ), -ng( 3 ):ng( 3 ) ) )

  norm = ( alpha / Sqrt( 3.1415926535897932384626433832795_wp ) ) ** 3
  Do i3 = -nt( 3 ), nt( 3 )
     Do i2 = -nt( 2 ), nt( 2 )
        Do i1 = -nt( 1 ), nt( 1 )
           r = i1 * grid_vecs( :, 1 ) + i2 * grid_vecs( :, 2 ) + i3 * grid_vecs( :, 3 ) - ri
           arg = alpha_sq * Dot_Product( r, r )
           gauss = norm * Exp( - arg )
           grid( i1, i2, i3 ) = gauss
           x = r( 1 )
           y = r( 2 )
           z = r( 3 )
           gauss_x2 = gauss * ( - 2.0_wp * alpha_sq + 4.0 * alpha_sq * alpha_sq * x * x )
           gauss_y2 = gauss * ( - 2.0_wp * alpha_sq + 4.0 * alpha_sq * alpha_sq * y * y )
           gauss_z2 = gauss * ( - 2.0_wp * alpha_sq + 4.0 * alpha_sq * alpha_sq * z * z )
           laplacian( i1, i2, i3 ) = gauss_x2 + gauss_y2 + gauss_z2
        End Do
     End Do
  End Do

  Call FD%init( order, grid_vecs )
  Call FD%apply( -nt( 1 ), -nt( 2 ), -nt( 3 ), -ng( 1 ), -ng( 2 ), -ng( 3 ), &
       -ng( 1 ), -ng( 2 ), -ng( 3 ), ng( 1 ), ng( 2 ), ng( 3 ), &
       grid, fd_laplacian )

  Write( *, * ) Maxval( Abs( fd_laplacian - &
       laplacian( -ng( 1 ):ng( 1 ), -ng( 2 ):ng( 2 ), -ng( 3 ):ng( 3 ) ) ) )
  Write( *, * ) Maxval( Abs( laplacian( -ng( 1 ):ng( 1 ), -ng( 2 ):ng( 2 ), -ng( 3 ):ng( 3 ) ) ) )
  
End Program testit
