## The following file contains some functions that compute elements with
## the Lie (super)algebra sl(1|2) ⊗ C[t]


# This function is to objectify a list into a SL2TensorCt element

SLTensorCtElement := function( list )
    
    local length;
    
    length := Length( list );
    
    if IsInt( list[1] ) then
    
        return Objectify( NewType( SLTensorCtFamily,  IsSLTensorCtElement and 
                       IsSLTensorCtElementRep ),
                       [ List( [1,3..length-1], x->list[x] ), 
                         List( [2,4..length], x->list[x] )]);
    else
        
        return Objectify( NewType( SLTensorCtFamily,  IsSLTensorCtElement and 
                       IsSLTensorCtElementRep ), list );
    fi;
end;

InstallMethod( PrintObj,
        "For elements of sl2 tensor C[t]",
        [ IsSLTensorCtElement and IsSLTensorCtElementRep  ],
        function( x )
    
    local coeffs, monoms, l, i;
        
    coeffs := x![1];
    monoms := x![2];
    l := Length( coeffs );
    
    if l = 0 then
        Print( "0" );
        return;
    fi;
    
    Print( coeffs[1], "*x(", monoms[1][1],",",monoms[1][2],")⊗ ",
           "t^",monoms[1][3] );

    
    for i in [2..l] do
        if coeffs[i] > 0 then
            Print( "+" );
        else
            Print( "-" );
        fi;
        
        Print( AbsInt( coeffs[i] ), "*x(", monoms[i][1],",",monoms[i][2],")⊗ ", "t", 
               "^",monoms[i][3] );
    od;

end );

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


    
