


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

