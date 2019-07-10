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
