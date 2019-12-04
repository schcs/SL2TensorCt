## The following file contains some functions that compute elements with
## the Lie (super)algebra sl(1|2) ⊗ C[t]
## Written by Lucas Calixto and Csaba Schneider

# This function is to objectify a list into a SL2TensorCt element

# the following function can be used to define an element of SL(1|2)[t]
# gap> SLTensorCtElement( [-1,[1,2,1],2,[2,3,2]] );
# -x(1,2)⊗ t^1+2*x(2,3)⊗ t^2

SLTensorCtElement := function( list )
    
    local length;
    
    if list = [] then
        list := [[],[]];
    fi;
    
    length := Length( list );
    
    if Length( list ) >= 1 and IsInt( list[1] ) then
    
        return Objectify( NewType( SLTensorCtFamily,  IsSLTensorCtElement and 
                       IsSLTensorCtElementRep ),
                       [ List( [1,3..length-1], x->list[x] ), 
                         List( [2,4..length], x->list[x] )]);
    else
        return Objectify( NewType( SLTensorCtFamily,  IsSLTensorCtElement and 
                       IsSLTensorCtElementRep ), list );
    fi;
end;

# the following function is to define a basis element of sl(1|2)[t] 
# gap> SLTensorCtBasisElement ([1,2,3] );    
# x(1,2)⊗ t^3


SLTensorCtBasisElement := function( list )
    local length;
    
    list := ShallowCopy( list );
    
    return Objectify( NewType( SLTensorCtBasisFamily,  
                   IsSLTensorCtBasisElement and 
                       IsSLTensorCtBasisElementRep ), list );
end;

CoefficientList := function( el )
    return el![1];
end;

MonomialList := function( el )
    return List( el![2], x->SLTensorCtBasisElement( x ));
end;

MonomialComponents := function( mon )
    return [mon![1],mon![2],mon![3]];
end;

# the following collects an element of sl(1|2)[t]


CollectSLTensorCtElement := function( el )
    local i, coeffs, monoms, list, m, pos;
    
    coeffs := CoefficientList( el );
    monoms := List( MonomialList( el ), MonomialComponents );
    
    list := [[],[]];
    
    for i in [1..Length( monoms )] do
        m := monoms[i];
        pos := Position( list[2], m );
        if pos = fail then
            Add( list[1], coeffs[i] );
            Add( list[2], monoms[i] );
        elif list[1][pos]+coeffs[i] <> 0 then
            list[1][pos] := list[1][pos]+coeffs[i];
        else
            RemoveElmList( list[1], pos );
            RemoveElmList( list[2], pos );
         fi;
    od;
        
    return SLTensorCtElement( list );
end;

InstallMethod( \+,
        "For elements of sl2 tensor C[t]",
        [ IsSLTensorCtElement and IsSLTensorCtElementRep, 
          IsSLTensorCtElement and IsSLTensorCtElementRep ],
        function( x, y )
    
    local coeffs, monoms;
    
    coeffs := []; monoms := [];
    Append( coeffs, CoefficientList( x ));
    Append( coeffs, CoefficientList( y ));
    
    Append( monoms, MonomialList( x ));
    Append( monoms, MonomialList( y ));
    
    return CollectSLTensorCtElement( SLTensorCtElement( [coeffs,monoms] ));

end );

InstallMethod( \=,
        "For elements of sl2 tensor C[t]",
        [ IsSLTensorCtElement and IsSLTensorCtElementRep, 
          IsSLTensorCtElement and IsSLTensorCtElementRep ],
        function( x, y )
    
    local el;
    
    el := x-y;
    return Length( CoefficientList( el )) = 0; end );
    
InstallMethod( \=,
        "For basis elements of sl2 tensor C[t]",
        [ IsSLTensorCtBasisElement and IsSLTensorCtBasisElementRep, 
          IsSLTensorCtBasisElement and IsSLTensorCtBasisElementRep ],
        function( x, y )
    
    return MonomialComponents( x ) = MonomialComponents( y );
    
end );

InstallMethod( \<, 
        "For basis elements of sl2 tensor C[t]",
        [ IsSLTensorCtBasisElement and IsSLTensorCtBasisElementRep, 
          IsSLTensorCtBasisElement and IsSLTensorCtBasisElementRep ],
        function( x, y )
    
    local xc, yc, l, xp, yp;
    
    # this list defines the order for the basis elements of sl(1|2)    
    l := [[1,2],[3,1],[3,2],[1,1],[3,3],[2,1],[1,3],[2,3]];
    
    # in this order, x tensor t^i < y tensor t^j iff x < y or x = y and j > i 
    
    xc := MonomialComponents( x ); yc := MonomialComponents( y );
    xp := Position( l, xc{[1,2]} ); yp := Position( l, yc{[1,2]});
    
    if xp < yp then 
    	return true; 
    elif xp > yp then 
        return false;
    elif xp = yp then
        return xc[3] > yc[3];
    fi;	                                         
end );

RandomSLTensorCtElement := function()
    
    local coeffs, mons;
    
    coeffs := List( [1..8], x->Random( Rationals ));
    mons := [[1,2,Random( [1..10] )], 
             [1,3,Random( [1..10] )], 
             [2,1,Random( [1..10] )], 
             [2,3,Random( [1..10] )], 
             [3,1,Random( [1..10] )], 
             [3,2,Random( [1..10] )],
             [1,1,Random( [1..10] )], 
             [3,3,Random( [1..10] )]];
    return SLTensorCtElement( [coeffs,mons] );
end;

Parity := function( mon )
    local vec, p0;
    
    vec := MonomialComponents( mon );
    
    if vec{[1,2]} in [[1,2],[1,3],[2,1],[3,1]] then
        p0 := 1;
    elif vec{[1,2]} in [[2,3],[3,2],[1,1],[3,3]] then
        p0 := 0;
    else
        Error( "Illegal basis element" );
    fi;
    
    return p0;
end;

InstallMethod( \*,
        "Scalar multiplication for basis elements of sl2 tensor C[t]",
        [ IsRat, IsSLTensorCtBasisElement and IsSLTensorCtBasisElementRep ],
        function( x, y )
    
    return SLTensorCtElement( [ [ x ], [ MonomialComponents( y )]]);
end );



InstallMethod( \*,
        "For basis elements of sl2 tensor C[t]",
        [ IsSLTensorCtBasisElement and IsSLTensorCtBasisElementRep, 
          IsSLTensorCtBasisElement and IsSLTensorCtBasisElementRep ],
        function( x, y )
    local table, elem, pos, prod, tpow, i, j, xc, yc;
    
    table := [[1,2,2,3,1,[1,3]],
              [1,2,2,1,1,[1,1]],
              [1,2,3,1,1,[3,2]],
              [1,3,2,1,1,[2,3]],
              [1,3,3,1,1,[3,3]],
              [1,3,3,2,1,[1,2]],
              [2,3,3,1,1,[2,1]],
              [2,3,3,2,1,[1,1],-1,[3,3]],
              [1,1,1,3,1,[1,3]],
              [1,1,2,3,1,[2,3]],
              [3,1,1,1,1,[3,1]],
              [3,2,1,1,1,[3,2]],
              [3,3,1,2,1,[1,2]],
              [2,3,3,3,1,[2,3]],
              [2,1,3,3,1,[2,1]],
              [3,3,3,2,1,[3,2]]];
    
    xc := MonomialComponents( x ); yc := MonomialComponents( y );
    elem := [xc[1],xc[2],yc[1],yc[2]];
    
    prod := [];
    
    pos := Position( List( table, x->[x[1],x[2],x[3],x[4]] ), elem );
    if IsInt( pos ) then 
        prod := table[pos]{[5..Length( table[pos])]};
    else
        elem := [yc[1],yc[2],xc[1],xc[2]];
        pos := Position( List( table, x->[x[1],x[2],x[3],x[4]] ), elem );
        if IsInt( pos ) then
            prod := table[pos]{[5..Length( table[pos])]};
        fi;
        if Parity( x )*Parity( y ) = 0 then
            for i in [1,3..Length( prod )-1] do
                prod[i] := -prod[i];
            od;
        fi;
    fi;
    
    if prod = [] then
        return SLTensorCtElement( [[],[]] );
    fi;
    
    tpow := xc[3]+yc[3];
    
    for i in [2,4..Length( prod )] do
        prod[i][3] := tpow;
    od;
    
    return SLTensorCtElement( [ List( [1,3..Length( prod )-1], x->prod[x]),
                   List( [2,4..Length( prod )], x->prod[x])]);
end );

InstallMethod( ZeroOp,
        "For elements of sl2 tensor C[t]",
        [ IsSLTensorCtElement and IsSLTensorCtElementRep ],
        x -> SLTensorCtElement( [[],[]] ));

InstallMethod( AdditiveInverseOp,
        "For elements of sl2 tensor C[t]",
        [ IsSLTensorCtElement and IsSLTensorCtElementRep ],
        x -> (-1)*x );


InstallMethod( \*,
        "scalar multiple of elements of sl2 tensor C[t]",
        [ IsRat, 
          IsSLTensorCtElement and IsSLTensorCtElementRep ],
        function( a, x )
    
    local coeffsx, monomsx, res, i;
    
    coeffsx := CoefficientList( x );
    monomsx := List( MonomialList( x ), MonomialComponents );
    
    if coeffsx = [] then
        return SLTensorCtElement( [[],[]] );
    fi;
    
    coeffsx := a*coeffsx;
    
    return SLTensorCtElement( [coeffsx,monomsx] );
end );

InstallMethod( \*,
        "product of two of elements of sl2 tensor C[t]",
        [ IsSLTensorCtElement and IsSLTensorCtElementRep,
          IsSLTensorCtElement and IsSLTensorCtElementRep ],
        function( x, y )
    
    local coeffsx, monomsx, coeffsy, monomsy, coeff, prod, i, j, p, 
          tpow;
    
    coeffsx := CoefficientList( x );
    monomsx := MonomialList( x );
    coeffsy := CoefficientList( y );
    monomsy := MonomialList( y );
    
    prod := 0*x;
    
    for i in [1..Length( coeffsx )] do
        for j in [1..Length( coeffsy )] do
            coeff := coeffsx[i]*coeffsy[j];
            p := coeff*(monomsx[i]*monomsy[j]);
            prod := prod + p;
        od;
    od;
    
    return prod;
end );
