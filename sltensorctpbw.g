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
        "scalar multiple of PBW elements of sl2 tensor C[t]",
        [ IsRat, 
          IsSLTensorCtPBWElement and IsSLTensorCtPBWElementRep ],
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
    
    return SLTensorCtPBWElement( res );
end );
