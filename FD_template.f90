
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
     Procedure, Public :: set_vecs
  End type grid_vectors
  
  Private

Contains

  Subroutine set_vecs( g, vecs )

    Class( grid_vectors )        , Intent( InOut ) :: g
    Real( wp ), Dimension( :, : ), Intent( In    ) :: vecs

    Real( wp ) :: V
    
    g%dir_vecs = vecs

    Call givens_invert( g%dir_vecs, g%inv_vecs, V )
    g%volume = Abs( V )
    
  End Subroutine set_vecs

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

  Type, Abstract, Extends( grid_vectors ) :: FD_template
     Integer                                , Private :: max_deriv
     Integer                                , Private :: order
     Real( wp ), Dimension( : ), Allocatable, Private :: weights
   Contains
     Procedure, Public :: init
     Procedure, Public :: set_order
  End type FD_template

Contains

  Subroutine init( FD, max_deriv, order, vecs )

    Class( FD_template )         , Intent(   Out ) :: FD
    Integer                      , Intent( In    ) :: max_deriv
    Integer                      , Intent( In    ) :: order
    Real( wp ), Dimension( :, : ), Intent( In    ) :: vecs

    Call FD%set_order( max_deriv, order )
    Call FD%set_vecs( vecs )

  End Subroutine init
  
  Subroutine set_order( FD, max_deriv, order )

    Class( FD_template ), Intent( InOut ) :: FD
    Integer             , Intent( In    ) :: max_deriv
    Integer             , Intent( In    ) :: order

    FD%max_deriv = max_deriv
    FD%order     = order

  End Subroutine set_order

End Module FD_template_module
