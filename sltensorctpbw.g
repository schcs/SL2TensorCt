## The following file contains some functions that compute elements with
## the Lie (super)algebra sl(1|2) ⊗ C[t]


# This function is to objectify a list into a SL2TensorCt element

SLTensorCtPBWElement := function( list )
    
    local length;
    
    if list = [] then
        list := [[],[]];
    fi;
    
    length := Length( list );
    
    if Length( list ) >= 1 and IsInt( list[1] ) then
    
        return Objectify( NewType( SLTensorCtPBWFamily,  
                       IsSLTensorCtPBWElement and 
                       IsSLTensorCtPBWElementRep ),
                       [ List( [1,3..length-1], x->list[x] ), 
                         List( [2,4..length], x->list[x] )]);
    else
        return Objectify( NewType( SLTensorCtPBWFamily,  
                       IsSLTensorCtPBWElement and 
                       IsSLTensorCtPBWElementRep ), list );
    fi;
end;

SLTensorCtPBWMonomial := function( list )
    
    local length;
        
    return Objectify( NewType( SLTensorCtPBWMonomialFamily,  
                       IsSLTensorCtPBWMonomial and 
                       IsSLTensorCtPBWMonomialRep ), list );
end;

SplitSLTensorCtPBWMonomial := function( mon )
    local list, i;
    
    list := [];
    i := 1;
    while IsBound( mon![i] ) do
        Add( list, [mon![i][1], mon![i][2], mon![i][3]] );
        i := i+1;
    od;
    
    return list;
end;

CollectSLTensorCtPBWElement := function( el )
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
    return SLTensorCtPBWElement( list );
end;

InstallMethod( \+,
        "For PBW elements of sl2 tensor C[t]",
        [ IsSLTensorCtPBWElement and IsSLTensorCtPBWElementRep, 
          IsSLTensorCtPBWElement and IsSLTensorCtPBWElementRep ],
        function( x, y )
    
    local coeffs, monoms, el;
    
    coeffs := []; monoms := [];
    Append( coeffs, x![1] );
    Append( coeffs, y![1] );
    
    Append( monoms, x![2] );
    Append( monoms, y![2] );
    
    el := SLTensorCtPBWElement( [coeffs,monoms] );
    
    el := CollectSLTensorCtPBWElement( el );
    return el;

end );

InstallMethod( \=,
        "For PBW elements of sl2 tensor C[t]",
        [ IsSLTensorCtPBWElement and IsSLTensorCtPBWElementRep, 
          IsSLTensorCtPBWElement and IsSLTensorCtPBWElementRep ],
        function( x, y )
    
    local el;
    
    el := x-y;
    return Length( el![1] ) = 0; end );
      

InstallMethod( ZeroOp,
        "For PBW elements of sl2 tensor C[t]",
        [ IsSLTensorCtPBWElement and IsSLTensorCtPBWElementRep ],
        x -> SLTensorCtPBWElement( [[],[]] ));

InstallMethod( AdditiveInverseOp,
        "For PBW elements of sl2 tensor C[t]",
        [ IsSLTensorCtPBWElement and IsSLTensorCtPBWElementRep ],
        x -> (-1)*x );

InstallMethod( \*,
        "scalar multiple of PBW monomial of sl2 tensor C[t]",
        [ IsRat, 
          IsSLTensorCtPBWMonomial and IsSLTensorCtPBWMonomialRep ],
        function( a, x )
    
    return SLTensorCtPBWElement( [[a],[SplitSLTensorCtPBWMonomial( x )]] );
end );

InstallMethod( \*,
        "scalar multiple of PBW elements of sl2 tensor C[t]",
        [ IsRat, 
          IsSLTensorCtPBWElement and IsSLTensorCtPBWElementRep ],
        function( a, x )
    
    local coeffsx, monomsx, res, i;
   
    coeffsx := x![1];
    monomsx := x![2];
    
    if coeffsx = [] then
        return SLTensorCtPBWElement( [[],[]] );
    fi;
    
    res := [];
    
    for i in [1..Length( coeffsx )] do
        Append( res, [ a*coeffsx[i], monomsx[i] ] );
    od;
    
    return SLTensorCtPBWElement( res );
end );

InstallMethod( \*,
        "product PBW monomials of sl2 tensor C[t]",
        [ IsSLTensorCtPBWMonomial and IsSLTensorCtPBWMonomialRep, 
          IsSLTensorCtPBWMonomial and IsSLTensorCtPBWMonomialRep ],
        function( x, y )
    
    local xl, yl;
    
    xl := SplitSLTensorCtPBWMonomial( x );
    yl := SplitSLTensorCtPBWMonomial( y );
    Append( xl, yl );
    
    return SLTensorCtPBWMonomial( xl );
end );

InstallMethod( \*,
        "product PBW monomials of sl2 tensor C[t]",
        [ IsSLTensorCtPBWElement and IsSLTensorCtPBWElement, 
          IsSLTensorCtPBWElement and IsSLTensorCtPBWElementRep ],
        function( x, y )
    
    local coeffsx, coeffsy, monomsx, monomsy, res, i, j, coeff, monom;
    
    coeffsx := x![1];
    coeffsy := y![1];
    monomsx := x![2];
    monomsy := y![2];
    res := [[],[]];
    
    for i in [1..Length( coeffsx )] do
        for j in [1..Length( coeffsy )] do
            coeff := coeffsx[i]*coeffsy[j];
            monom := monomsx[i];
            Append( monom, monomsy[j] );
            Add( res[1], coeff );
            Add( res[2], monom );
        od;
    od;
    
    return SLTensorCtPBWElement( res );
    
end );

IsCollectedSLTensorCtPBWMonomial := function( mon )
    local els, i;
    
    els := List( SplitSLTensorCtPBWMonomial( mon ), SLTensorCtBasisElement );
    
    for i in [1..Length( els )-1] do
        if els[i] > els[i+1] then 
            return i;
        fi;
    od;
    
    return true;
end;

CollectSLTensorCtPBWMonomialAtPosition := function( mon, pos )
   local monlist, list1, list2, x1, y1, px, py, z, mon1, mon2, y1x1, zz, e1, e2; 
   
    monlist := SplitSLTensorCtPBWMonomial( mon );
    list1 := monlist{[1..pos-1]};
    list2 := monlist{[pos+2..Length( monlist )]};
    
    x1 := SLTensorCtBasisElement( monlist[pos] );
    y1 := SLTensorCtBasisElement( monlist[pos+1] );
    px := Parity( x1 ); py := Parity( y1 );
    
    z := x1*y1;
    
    mon1 := SLTensorCtPBWMonomial( list1 );
    mon2 := SLTensorCtPBWMonomial( list2 );
    y1x1 := SLTensorCtPBWMonomial( [monlist[pos+1],monlist[pos]] );
    zz := SLTensorCtPBWElement( [ z![1], List( [z![2], x->[x]])]);
    
    return (1*mon1)*(1*y1x1)*(1*mon2)-(px*py)*((1*mon1)*zz*(1*mon2));
end;

    
    
    
    
    