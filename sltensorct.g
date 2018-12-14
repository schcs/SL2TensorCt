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
    Error();
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
    
    res := [];
    
    for i in [1..Length( coeffsx )] do
        Append( res, [ a*coeffsx[i], monomsx[i] ] );
    od;
    
    return SLTensorCtElement( res );
end );

InstallMethod( \*,
        "scalar multiple of elements of sl2 tensor C[t]",
        [ IsSLTensorCtElement and IsSLTensorCtElementRep,
          IsSLTensorCtElement and IsSLTensorCtElementRep ],
        function( x, y )
    
    local coeffsx, monomsx, coeffsy, monomsy, res, coeff, elem, prod, i, j,
          tpow;
    
    coeffsx := x![1];
    monomsx := x![2];
    coeffsy := y![1];
    monomsy := y![2];
    
    res := [[],[]];
    
    prod := [];
    coeff := [];
    
    for i in [1..Length( coeffsx )] do
        for j in [1..Length( coeffsy )] do
            coeff := [coeffsx[i]*coeffsy[j]];
            elem := [monomsx[i][1],monomsx[i][2],monomsy[j][1],monomsy[j][2]];
            tpow := monomsx[i][3]+monomsy[j][3];
            if elem = [1,2,2,3] then
                prod := [[1,3,tpow]];
            elif elem = [2,3,1,2] then
                prod := [[1,3,tpow]];
                coeff[1] := -coeff[1];
            elif elem = [1,2,2,1] or elem = [2,1,1,2] then
                prod := [[1,1,tpow]];
            elif elem = [1,2,3,1] or elem = [3,1,1,2] then
                prod := [[3,2,tpow]];
            elif elem = [1,3,2,1] or elem = [2,1,1,3] then
                prod := [[2,3,tpow]];
            elif elem = [1,3,3,1] or elem = [3,1,1,3] then
                prod := [[3,3,tpow]];
            elif elem = [1,3,3,2] then
                prod := [[1,2,tpow]];
            elif elem = [3,2,1,3] then
                prod := [[1,2,tpow]];
                coeff[1] := -coeff[1];
            elif elem = [2,3,3,1] then
                prod := [[1,2,tpow]];
            elif elem = [3,1,2,3] then
                prod := [[1,2,tpow]];
                coeff[1] := -coeff[1];
            elif elem = [2,3,3,2] then 
                prod := [[1,1,tpow],[3,3,tpow]];
                Add( coeff, -coeff[1] );
            elif elem = [3,2,2,3] then
                prod := [[1,1,tpow],[3,3,tpow]];
                coeff[1] := -coeff[1];
                Add( coeff, -coeff[1] );
            elif elem = [1,1,1,3] then
                prod := [[1,3,tpow]];
            elif elem = [1,3,1,1] then
                prod := [[1,3,tpow]];
                coeff[1] := -coeff[1];
            elif elem = [1,1,2,3] then
                prod := [[2,3,tpow]];
            elif elem = [2,3,1,1] then
                prod := [[2,3,tpow]];
                coeff[1] := -coeff[1];
            elif elem = [1,1,3,1] then
                prod := [[3,1,tpow]];
                coeff[1] := -coeff[1];
            elif elem = [3,1,1,1] then
                prod := [[3,1,tpow]];
            elif elem = [1,1,3,2] then
                prod := [[3,2,tpow]];
                coeff[1] := -coeff[1];
            elif elem = [3,2,1,1] then
                prod := [[3,2,tpow]];
            elif elem = [3,3,1,2] then
                prod := [[1,2,tpow]];
            elif elem = [1,2,3,3] then
                prod := [[1,2,tpow]];
                coeff[1] := -coeff[1];
            elif elem = [3,3,2,3] then
                prod := [[2,3,tpow]];
                coeff[1] := -coeff[1];
            elif elem = [2,3,3,3] then
                prod := [[2,3,tpow]];
            elif elem = [3,3,2,1] then
                prod := [[2,1,tpow]];
                coeff[1] := -coeff[1];
            elif elem = [2,1,3,3] then
                prod := [[2,1,tpow]];
            elif elem = [3,3,3,2] then
                prod := [[3,2,tpow]];
            elif elem = [3,2,3,3] then
                prod := [[3,2,tpow]];
                coeff[1] := -coeff[1];
            fi;
            
            if prod <> [] then
                Append( res[1], coeff );
                Append( res[2], prod );
            fi;
            
        od;
    od;
    return CollectSLTensorCtElement( SLTensorCtElement( res ));
end );