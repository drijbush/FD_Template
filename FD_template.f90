
Module numbers_module

  Implicit None

  Integer, Parameter, Public :: wp = Selected_real_kind( 12, 70 )

  Private

End Module numbers_module

Module grid_vectors_module

  Use numbers_module, Only : wp

  Implicit None

  Type, Public :: grid_vectors
     Real( wp ), Dimension( :, : ), Allocatable, Private :: dir_vecs
     Real( wp ),                                 Private :: volume
     Real( wp ), Dimension( :, : ), Allocatable, Private :: inv_vecs
   Contains
     Procedure, Public :: set_dir_vecs
     Procedure, Public :: get_dir_vec
     Procedure, Public :: get_volume
     Procedure, Public :: get_inv_vec
  End type grid_vectors
  
  Private

Contains

  Subroutine set_dir_vecs( g, vecs )

    Class( grid_vectors )        , Intent( InOut ) :: g
    Real( wp ), Dimension( :, : ), Intent( In    ) :: vecs

    Real( wp ) :: V
    
    g%dir_vecs = vecs
    Allocate( g%inv_vecs, Mold = g%dir_vecs )
    
    Call givens_invert( g%dir_vecs, g%inv_vecs, V )
    g%volume = Abs( V )

  End Subroutine set_dir_vecs

  Pure Function get_dir_vec( g, which ) Result( v )

    Real( wp ), Dimension( : ), Allocatable :: v
    
    Class( grid_vectors ), Intent( In ) :: g
    Integer              , Intent( In ) :: which

    v = g%dir_vecs( :, which )
    
  End Function get_dir_vec

  Pure Function get_inv_vec( g, which ) Result( v )

    Real( wp ), Dimension( : ), Allocatable :: v
    
    Class( grid_vectors ), Intent( In ) :: g
    Integer              , Intent( In ) :: which

    v = g%inv_vecs( which, : )
    
  End Function get_inv_vec

  Pure Function get_volume( g ) Result( V )

    Real( wp ) :: V
    
    Class( grid_vectors ), Intent( In ) :: g

    V = g%volume

  End Function get_volume

  Subroutine givens_invert( A, B, det )

    ! Invert the matrix B, returning the inverse in B
    ! and the determinant of A in det

    ! NOT FOR USE ON LARGE MATRICES

    ! Really designed as a robust, not totally disastorous, way
    ! to find inverse lattice vectors and the volume of the unit cell
    ! ( i.e. the determinant )

    ! The method used is to QR factorise ( in fact LQ factorise ),
    ! invert the triangular matrix and the from that form the
    ! inverse of the original matrix. The determinant is simply the
    ! product of the diagonal elements of the diagonal matrix.
    ! Method chosen as easy to implement and nicely numerically stable,
    ! especially for the small matrices need here.

    Real( wp ), Dimension( 1:, 1: ), Intent( In    )           :: A
    Real( wp ), Dimension( 1:, 1: ), Intent(   Out )           :: B
    Real( wp )                     , Intent(   Out ), Optional :: det

    Real( wp ), Dimension( :, : ), Allocatable :: Q, L

    Real( wp ) :: theta, c, s
    Real( wp ) :: ann, anm
    Real( wp ) :: Lim, Lin
    Real( wp ) :: Qmi, Qni
    Real( wp ) :: mv

    Integer :: dim
    Integer :: m, n
    Integer :: i

    ! Size of the problem
    dim = Size( A, dim = 1 )

    ! The factors of the matrix A
    Allocate( L( 1:dim, 1:dim ) )
    Allocate( Q( 1:dim, 1:dim ) )

    ! Factorise A by Givens rotations
    L = A
    Q = 0.0_wp
    Do i = 1, dim
       Q( i, i ) = 1.0_wp
    End Do
    LQ_factorise: Do m = 2, dim
       Do n = 1, m - 1
          ! For element mn generate the rotation which
          ! will zweo that element
          ann = L( n, n )
          anm = L( n, m )
          theta = atan2( - anm, ann )
          c = Cos( theta )
          s = Sin( theta )
          ! And apply that rotation to the original matrix
          Do i = 1, Dim
             Lim = L( i, m )
             Lin = L( i, n )
             L( i, n ) = c * Lin - s * Lim
             L( i, m ) = c * Lim + s * Lin
          End Do
          ! And then apply the rotation to the orthognal matrix
          Do i = 1, Dim
             Qmi = Q( m, i )
             Qni = Q( n, i )
             Q( n, i ) = c * Qni - s * Qmi 
             Q( m, i ) = c * Qmi + s * Qni
          End Do
       End Do
    End Do LQ_factorise

    ! Now invert the triangular matrix, returning the inverse in B
    B(  1, 1  ) = 1.0_wp / L( 1, 1 )
    B(  1, 2: ) = 0.0_wp
    Do m = 2, dim
       Do n = 1, m
          mv = 0.0_wp
          Do i = 1, m - 1
             mv = mv - L( m, i ) * B( i, n )
          End Do
          If( m == n ) Then
             mv = mv + 1.0_wp
          End If
          mv = mv / L( m, m )
          B( m, n ) = mv
       End Do
       B( m, n: ) = 0.0_wp
    End Do

    ! The inverse of A is now trivially formed
    B = Matmul( Transpose( Q ), B )

    ! And generate the determinant
    If( Present( det ) ) Then
       det = 1.0_wp
       Do i = 1, dim
          det = det * L( i, i )
       End Do
    End If

    Deallocate( Q )
    Deallocate( L )

  End Subroutine givens_invert
  
End Module grid_vectors_module

Module FD_template_module

  Use grid_vectors_module, Only : grid_vectors
  
  Implicit None

  Integer, Parameter, Public :: wp = Selected_real_kind( 12, 70 )

  Private

  Type, Abstract, Extends( grid_vectors ), Public :: FD_template
     Integer                                   , Private :: max_deriv
     Integer                                   , Private :: order
     Real( wp ), Dimension( :, : ), Allocatable, Private :: weights
   Contains
     Procedure, Public :: FD_init
     Procedure, Public :: set_order
     Procedure, Public :: get_order
     Procedure, Public :: get_max_deriv
     Procedure, Public :: get_weight
  End type FD_template

Contains

  Subroutine FD_init( FD, max_deriv, order, vecs )

    Class( FD_template )         , Intent(   Out ) :: FD
    Integer                      , Intent( In    ) :: max_deriv
    Integer                      , Intent( In    ) :: order
    Real( wp ), Dimension( :, : ), Intent( In    ) :: vecs

    Call FD%set_order( max_deriv, order )
    Call FD%set_dir_vecs( vecs )

  End Subroutine FD_init
  
  Subroutine set_order( FD, max_deriv, order )

    Class( FD_template ), Intent( InOut ) :: FD
    Integer             , Intent( In    ) :: max_deriv
    Integer             , Intent( In    ) :: order

    Integer :: i
    
    FD%max_deriv = max_deriv
    FD%order     = order

    Call weights( 0.0_wp, [ ( Real( i, wp ), i = -FD%order, FD%order ) ], &
         2 * FD%order, FD%max_deriv, FD%weights )

  End Subroutine set_order

  Pure Function get_order( FD ) Result( order )

    Integer :: order

    Class( FD_template ), Intent( In ) :: FD
    
    order = FD%order

  End Function get_order

  Pure Function get_max_deriv( FD ) Result( max_deriv )

    Integer :: max_deriv

    Class( FD_template ), Intent( In ) :: FD
    
    max_deriv = FD%max_deriv

  End Function get_max_deriv

  Pure Function get_weight( FD, deriv  ) Result( weight )

    Real( wp ), Allocatable, Dimension( : ) :: weight

    Class( FD_template ), Intent( In ) :: FD
    Integer             , Intent( In ) :: deriv

    Allocate( weight( -FD%order:FD%order ) )
    weight = FD%weights( :, deriv )

  End Function get_weight

  Subroutine weights (z,x,n,m,c)
    !---
    !--- Input Parameters
    !--- z location where approximations are to be accurate,
    !--- x(0:nd) grid point locations, found in x(0:n)
    !--- n one less than total number of grid points; n must
    !--- not exceed the parameter nd below,
    !--- nd dimension of x- and c-arrays in calling program
    !--- x(0:nd) and c(0:nd,0:m), respectively,
    !--- m highest derivative for which weights are sought,
    !--- Output Parameter
    !--- c(0:nd,0:m) weights at grid locations x(0:n) for derivatives
    !--- of order 0:m, found in c(0:n,0:m)
    !---
    Implicit None
    Real( wp ), Intent( In ) :: z
    Real( wp ), Dimension( 0: ), Intent( In ) :: x
    Integer, Intent( In ) :: n
    Integer, Intent( In ) :: m
    Real( wp ), Dimension( :, : ), Allocatable, Intent( Out ) :: c

    Real( wp ) :: c1, c2, c3, c4, c5
    Integer :: mn
    Integer :: i, j, k

    Allocate( c( 0:n, 0:m ) )

    c1 = 1.0_wp
    c4 = x(0)-z
    c = 0.0_wp
    c(0,0) = 1.0_wp
    Do i=1,n
       mn = Min(i,m)
       c2 = 1.0_wp
       c5 = c4
       c4 = x(i)-z
       Do j=0,i-1
          c3 = x(i)-x(j)
          c2 = c2*c3
          If (j.Eq.i-1) Then
             Do k=mn,1,-1
                c(i,k) = c1*(k*c(i-1,k-1)-c5*c(i-1,k))/c2
             End Do
             c(i,0) = -c1*c5*c(i-1,0)/c2
          Endif
          Do k=mn,1,-1
             c(j,k) = (c4*c(j,k)-k*c(j,k-1))/c3
          End Do
          c(j,0) = c4*c(j,0)/c3
       End Do
       c1 = c2
    End Do

  End Subroutine weights

End Module FD_template_module

Module FD_Laplacian_3d_module

  Use numbers_module    , Only : wp
  Use FD_template_module, Only : FD_template
  
  Implicit None

  Integer, Parameter :: n_cache = ( 2 ** 18 ) / ( 8 ) ! Number of reals in cache

  Integer, Parameter, Public :: XX = 1
  Integer, Parameter, Public :: XY = 2
  Integer, Parameter, Public :: XZ = 3
  Integer, Parameter, Public :: YY = 4
  Integer, Parameter, Public :: YZ = 5
  Integer, Parameter, Public :: ZZ = 6

  Type, Extends( FD_template ), Public :: FD_Laplacian_3d
     Integer   , Dimension( 1:3 ), Private :: nc_block
     Real( wp ), Dimension( 1:6 ), Private :: deriv_weights
   Contains
     Procedure, Public :: init
     Procedure, Public :: reset_vecs
     Procedure, Public :: apply
  End type FD_Laplacian_3d

  Private

Contains

  Subroutine init( FD, order, vecs )

    Class( FD_Laplacian_3d )     , Intent(   Out ) :: FD
    Integer                      , Intent( In    ) :: order
    Real( wp ), Dimension( :, : ), Intent( In    ) :: vecs

    Call FD%FD_init( 2, order / 2, vecs )

    FD%nc_block = 1
    Do While( usage( FD%nc_block, 2 * FD%get_order() ) < n_cache )
       FD%nc_block = FD%nc_block + 1
    End Do
    FD%nc_block( 1 ) = FD%nc_block( 1 ) + 1
    If( usage( FD%nc_block, 2 * FD%get_order() ) >= n_cache ) Then
       FD%nc_block( 1 ) = FD%nc_block( 1 ) - 1
    Else
       FD%nc_block( 2 ) = FD%nc_block( 2 ) + 1
       If( usage( FD%nc_block, 2 * FD%get_order() ) >= n_cache ) Then
          FD%nc_block( 2 ) = FD%nc_block( 2 ) - 1
       End If
    End If

    Call FD%reset_vecs( vecs )
    
  End Subroutine init

  Subroutine reset_vecs( FD, vecs )

    Class( FD_Laplacian_3d )     , Intent( InOut ) :: FD
    Real( wp ), Dimension( :, : ), Intent( In    ) :: vecs

    Integer :: this
    Integer :: i, j
    
!!$    Call FD%FD_init( 2, FD%get_order(), vecs )
    
    this = 0
    Do i = 1, 3
       Do j = i, 3
          this = this + 1
          FD%deriv_weights( this ) = Dot_Product( FD%get_inv_vec( i ), FD%get_inv_vec( j ) )
          If( i /= j ) FD%deriv_weights( this ) = 2.0_wp * FD%deriv_weights( this )
       End Do
    End Do
    
  End Subroutine reset_vecs

  Subroutine apply( FD, l1g, l2g, l3g, l1l, l2l, l3l, s1, s2, s3, f1, f2, f3, grid, Laplacian )

    Class( FD_Laplacian_3d )                 , Intent( In    ) :: FD
    Integer                                  , Intent( In    ) :: l1g ! lower bounds grid
    Integer                                  , Intent( In    ) :: l2g
    Integer                                  , Intent( In    ) :: l3g
    Integer                                  , Intent( In    ) :: l1l ! lower bounds laplacian
    Integer                                  , Intent( In    ) :: l2l
    Integer                                  , Intent( In    ) :: l3l
    Integer                                  , Intent( In    ) :: s1 ! start point for calculation
    Integer                                  , Intent( In    ) :: s2
    Integer                                  , Intent( In    ) :: s3
    Integer                                  , Intent( In    ) :: f1 ! final point for calculation
    Integer                                  , Intent( In    ) :: f2
    Integer                                  , Intent( In    ) :: f3
    Real( wp ), Dimension( l1g:, l2g:, l3g: ), Intent( In    ) :: grid ! Thing to be differentiated
    Real( wp ), Dimension( l1l:, l2l:, l3l: ), Intent(   Out ) :: laplacian

    Real( wp ), Dimension( : ), Allocatable :: w1, w2
    
    Integer :: order
    
    Integer :: i_block_3, i_block_2, i_block_1

    order = FD%get_order()

    w1 = FD%get_weight( 1 )
    w2 = FD%get_weight( 2 )

    !$omp do collapse( 3 )
    Do i_block_3 = s3, f3, FD%nc_block( 3 )
       Do i_block_2 = s2, f2, FD%nc_block( 2 )
          Do i_block_1 = s1, f1, FD%nc_block( 1 )
             Call apply_block( [ i_block_1, i_block_2, i_block_3 ], &
                  [ l1g, l2g, l3g ], [l1l, l2l, l3l ], FD%nc_block, [ f1, f2, f3 ], &
                  order, w1, w2, FD%deriv_weights, grid, laplacian )
          End Do
       End Do
    End Do
    !$omp end do

  End Subroutine apply

  Pure Subroutine apply_block( s, lg, ll, nb, f, order, w1, w2, deriv_weights, grid, laplacian )

    ! Apply the FD laplacian operator to part of the grid. An effort has been made
    ! ro make sure the required data is in cache
    
    Integer, Dimension( 1:3 )       , Intent( In    ) :: s   ! Where to start
    Integer, Dimension( 1:3 )       , Intent( In    ) :: lg  ! Lower bound of grid array
    Integer, Dimension( 1:3 )       , Intent( In    ) :: ll  ! Lower bound of laplacian array
    Integer, Dimension( 1:3 )       , Intent( In    ) :: nb  ! Cache blocking factors
    Integer, Dimension( 1:3 )       , Intent( In    ) :: f   ! Where to finish
    Integer                         , Intent( In    ) :: order ! Order of the FD approximation
    Real( wp ), Dimension( -order: ), Intent( In    ) :: w1    ! FD Weights for first derivs
    Real( wp ), Dimension( -order: ), Intent( In    ) :: w2    ! FD Weights for second derivs
    Real( wp ), Dimension( 1:6     ), Intent( In    ) :: deriv_weights ! See below
    Real( wp ), Dimension( lg( 1 ):, lg( 2 ):, lg( 3 ): ), Intent( In    ) :: grid ! The source
    Real( wp ), Dimension( ll( 1 ):, ll( 2 ):, ll( 3 ): ), Intent( InOut ) :: laplacian ! The result

    ! Deriv_weights: As we do NOT assume the grid is orthogonal our FD laplacian is of the form
    ! d_xx * del_xx + d_xy * del_xy + d_xz * del_xz + d_yy * del_yy + d_yz * del_yz + d_zz * del_zz
    ! as we must finite difference along the directions of the grid. deriv_weights are the d_xx, d_xy
    ! etc. coefficients in this expression

    Real( wp ), Parameter :: orthog_tol = 1.0e-14_wp
    
    Real( wp ) :: st1, st2, st3
    Real( wp ) :: fac1, fac2, fac3
    Real( wp ) :: st12, st13, st23
    Real( wp ) :: fac12, fac13, fac23

    Integer :: i3, i2, i1
    Integer :: it, it1, it2, it3

    ! Order the loops for the various terms so that the inner loop is stride 1 

    ! First do the xx, yy, zz terms. Assume the weights of these are always non-zero
    Do i3 = s( 3 ), Min( s( 3 ) + nb( 3 ) - 1, f( 3  ) )
       Do i2 = s( 2 ), Min( s( 2 ) + nb( 2 ) - 1, f( 2 ) )
          Do i1 = s( 1 ), Min( s( 1 ) + nb( 1 ) - 1, f( 1 ) )

             ! FD x^2, y^2, z^2 at grid point 
             laplacian( i1, i2, i3 ) =                           w2( 0 ) * deriv_weights( 1 ) * grid( i1, i2, i3 )
             laplacian( i1, i2, i3 ) = laplacian( i1, i2, i3 ) + w2( 0 ) * deriv_weights( 4 ) * grid( i1, i2, i3 )
             laplacian( i1, i2, i3 ) = laplacian( i1, i2, i3 ) + w2( 0 ) * deriv_weights( 6 ) * grid( i1, i2, i3 )

             ! FD x^2
             Do it = 1, order

                fac1 = w2( it ) * deriv_weights( 1 )
                st1 = ( grid( i1 + it, i2     , i3      ) + grid( i1 - it, i2     , i3       ) )
                laplacian( i1, i2, i3 ) = laplacian( i1, i2, i3 ) + fac1 * st1

             End Do

          End Do
       End Do
    End Do

    ! FD y^2, z^2
    Do i3 = s( 3 ), Min( s( 3 ) + nb( 3 ) - 1, f( 3 ) )
       Do i2 = s( 2 ), Min( s( 2 ) + nb( 2 ) - 1, f( 2 ) )
          Do it = 1, order
             Do i1 = s( 1 ), Min( s( 1 ) + nb( 1 ) - 1, f( 1 ) )

                fac2 = w2( it ) * deriv_weights( 4 )
                st2 = ( grid( i1     , i2 + it, i3      ) + grid( i1     , i2 - it, i3       ) )
                laplacian( i1, i2, i3 ) = laplacian( i1, i2, i3 ) + fac2 * st2

                fac3 = w2( it ) * deriv_weights( 6 )
                st3 = ( grid( i1     , i2     , i3 + it ) + grid( i1     , i2     , i3 - it  ) )
                laplacian( i1, i2, i3 ) = laplacian( i1, i2, i3 ) + fac3 * st3

             End Do
          End Do
       End Do
    End Do

    ! xy
    If( Abs( deriv_weights( 2 ) ) > orthog_tol ) Then
       Do i3 = s( 3 ), Min( s( 3 ) + nb( 3 ) - 1, f( 3 ) )
          Do i2 = s( 2 ), Min( s( 2 ) + nb( 2 ) - 1, f( 2 ) )
             Do it2 = 1, order
                Do i1 = s( 1 ), Min( s( 1 ) + nb( 1 ) - 1, f( 1 ) )
                   ! first derivs have zero weight at the grid point for central differences
                   Do it1 = 1, order
                      fac12 = w1( it1 ) * w1( it2 ) * deriv_weights( 2 )
                      st12 = grid( i1 + it1, i2 + it2, i3 ) - grid( i1 - it1, i2 + it2, i3 ) - &
                           ( grid( i1 + it1, i2 - it2, i3 ) - grid( i1 - it1, i2 - it2, i3 ) )
                      laplacian( i1, i2, i3 ) = laplacian( i1, i2, i3 ) + fac12 * st12
                   End Do
                End Do
             End Do
          End Do
       End Do
    End If

    ! xz
    If( Abs( deriv_weights( 3 ) ) > orthog_tol ) Then
       Do i3 = s( 3 ), Min( s( 3 ) + nb( 3 ) - 1, f( 3 ) )
          Do it3 = 1, order
             Do i2 = s( 2 ), Min( s( 2 ) + nb( 2 ) - 1, f( 2 ) )
                Do i1 = s( 1 ), Min( s( 1 ) + nb( 1 ) - 1, f( 1 ) )
                   ! first derivs have zero weight at the grid point for central differences
                   Do it1 = 1, order
                      fac13 = w1( it1 ) * w1( it3 ) * deriv_weights( 3 )
                      st13 = grid( i1 + it1, i2, i3 + it3 ) - grid( i1 - it1, i2, i3 + it3 ) - &
                           ( grid( i1 + it1, i2, i3 - it3 ) - grid( i1 - it1, i2, i3 - it3 ) )
                      laplacian( i1, i2, i3 ) = laplacian( i1, i2, i3 ) + fac13 * st13
                   End Do
                End Do
             End Do
          End Do
       End Do
    End If

    ! yz
    If( Abs( deriv_weights( 5 ) ) > orthog_tol ) Then
       Do i3 = s( 3 ), Min( s( 3 ) + nb( 3 ) - 1, f( 3 ) )
          Do it3 = 1, order
             Do i2 = s( 2 ), Min( s( 2 ) + nb( 2 ) - 1, f( 2 ) )
                ! first derivs have zero weight at the grid point for central differences
                Do it2 = 1, order
                   Do i1 = s( 1 ), Min( s( 1 ) + nb( 1 ) - 1, f( 1 ) )
                      fac23 = w1( it2 ) * w1( it3 ) * deriv_weights( 5 )
                      st23 = grid( i1, i2 + it2, i3 + it3 ) - grid( i1, i2 - it2, i3 + it3 )  - &
                           ( grid( i1, i2 + it2, i3 - it3 ) - grid( i1, i2 - it2, i3 - it3 ) )
                      laplacian( i1, i2, i3 ) = laplacian( i1, i2, i3 ) + fac23 * st23
                   End Do
                End Do
             End Do
          End Do
       End Do
    End If
       
  End Subroutine apply_block
  
  Pure Function usage( nc_block, accuracy ) Result( reals )
    
    Integer :: reals
    
    Integer, Dimension( 1:3 ), Intent( In ) :: nc_block
    Integer,                   Intent( In ) :: accuracy
    
    Integer :: grid
    Integer :: fd
    
    grid = Product( nc_block )
    fd   = Product( nc_block + 1 + accuracy )
    
    reals = grid + fd
    
  End Function usage
  
End Module FD_Laplacian_3d_module

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
