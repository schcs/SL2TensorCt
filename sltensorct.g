## The following file contains some functions that compute elements with
## the Lie (super)algebra sl(1|2) ⊗ C[t]

# This function is to objectify a list into a SL2TensorCt element

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

SLTensorCtBasisElement := function( list )
    
    local length;
    
    list := ShallowCopy( list );
    
    
    return Objectify( NewType( SLTensorCtBasisFamily,  
                   IsSLTensorCtBasisElement and 
                       IsSLTensorCtBasisElementRep ), list );
end;

CollectSLTensorCtElement := function( el )
    local i, coeffs, monoms, list, m, pos;
    
    coeffs := el![1];
    monoms := el![2];
    
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
    
    local coeffs, monoms, el;
    
    coeffs := []; monoms := [];
    Append( coeffs, x![1] );
    Append( coeffs, y![1] );
    
    Append( monoms, x![2] );
    Append( monoms, y![2] );
    
    el := SLTensorCtElement( [coeffs,monoms] );
    el := CollectSLTensorCtElement( el );
    return el;

end );

InstallMethod( \=,
        "For elements of sl2 tensor C[t]",
        [ IsSLTensorCtElement and IsSLTensorCtElementRep, 
          IsSLTensorCtElement and IsSLTensorCtElementRep ],
        function( x, y )
    
    local el;
    
    el := x-y;
    return Length( el![1] ) = 0; end );
    
InstallMethod( \=,
        "For basis elements of sl2 tensor C[t]",
        [ IsSLTensorCtBasisElement and IsSLTensorCtBasisElementRep, 
          IsSLTensorCtBasisElement and IsSLTensorCtBasisElementRep ],
        function( x, y )
    
    return x![1] = y![1] and x![2] = y![2] and x![3] = y![3];
    
end );

InstallMethod( \<, 
        "For basis elements of sl2 tensor C[t]",
        [ IsSLTensorCtBasisElement and IsSLTensorCtBasisElementRep, 
          IsSLTensorCtBasisElement and IsSLTensorCtBasisElementRep ],
        function( x, y )
    if x![3] < y![3] then 
        return true; 
    elif x![2] < y![2] then
        return true; 
    elif x![1] < y![1] then
        return true;
    else
        return false;
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
    
    vec := [mon![1], mon![2], mon![3]];
    
    if vec{[1,2]} in [[1,2],[1,3],[2,1],[3,1]] then
        p0 := 1;
    elif vec{[1,2]} in [[2,3],[3,2],[1,1],[3,3]] then
        p0 := 0;
    else
        Error( "Illegal basis element" );
    fi;
    
    return (p0 + vec[3]) mod 2;
end;

InstallMethod( \*,
        "For basis elements of sl2 tensor C[t]",
        [ IsSLTensorCtBasisElement and IsSLTensorCtBasisElementRep, 
          IsSLTensorCtBasisElement and IsSLTensorCtBasisElementRep ],
        function( x, y )
    local table, elem, pos, prod, tpow, i, j;
    
    table := [[1,2,2,3,1,[1,3]],
              [1,2,2,1,1,[1,1]],
              [1,2,3,1,1,[3,2]],
              [1,3,2,1,1,[2,3]],
              [1,3,3,1,1,[3,3]],
              [1,3,3,2,1,[1,2]],
              [2,3,3,1,1,[1,2]],
              [2,3,3,2,1,[1,1],-1,[3,3]],
              [1,1,1,3,1,[1,3]],
              [1,1,2,3,1,[2,3]],
              [3,1,1,1,1,[3,1]],
              [3,2,1,1,1,[3,2]],
              [3,3,1,2,1,[1,2]],
              [2,3,3,3,1,[2,3]],
              [2,1,3,3,1,[2,1]],
              [3,3,3,2,1,[3,2]]];
              
              
    elem := [x![1],x![2],y![1],y![2]];
    
    prod := [];
    
    pos := Position( List( table, x->[x[1],x[2],x[3],x[4]] ), elem );
    if IsInt( pos ) then 
        prod := table[pos]{[5..Length( table[pos])]};
    else
        elem := [y![1],y![2],x![1],x![2]];
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
    
    tpow := x![3]+y![3];
    
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
    
    coeffsx := x![1];
    monomsx := x![2];
    
    if coeffsx = [] then
        return SLTensorCtElement( [[],[]] );
    fi;
    
    coeffsx := a*coeffsx;
    
    return SLTensorCtElement( [coeffsx,monomsx] );
end );

InstallMethod( \*,
        "scalar multiple of elements of sl2 tensor C[t]",
        [ IsSLTensorCtElement and IsSLTensorCtElementRep,
          IsSLTensorCtElement and IsSLTensorCtElementRep ],
        function( x, y )
    
    local coeffsx, monomsx, coeffsy, monomsy, coeff, prod, i, j, p, 
          monomx, monomy, tpow;
    
    coeffsx := x![1];
    monomsx := x![2];
    coeffsy := y![1];
    monomsy := y![2];
    
    prod := 0*x;
    
    for i in [1..Length( coeffsx )] do
        for j in [1..Length( coeffsy )] do
            coeff := coeffsx[i]*coeffsy[j];
            monomx := SLTensorCtBasisElement( monomsx[i] ); 
            monomy := SLTensorCtBasisElement( monomsy[j] );
            p := coeff*(monomx*monomy);
            prod := prod + p;
        od;
    od;
    
    return prod;
end );