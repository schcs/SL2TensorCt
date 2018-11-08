## The following file contains some functions that compute elements with
## the Lie (super)algebra sl(1|2) ⊗ C[t]


# This function is to objectify a list into a SL2TensorCt element

SLTensorCtElement := function( list )
    
    local length;
    
    length := Length( list );
    
    return Objectify( NewType( SLTensorCtFamily,  IsSLTensorCtElement and 
                   IsSLTensorCtElementRep ),
                   [ List( [1,3..length-1], x->list[x] ), 
                     List( [2,4..length], x->list[x] )]);
end;

InstallMethod( PrintObj,
        "For elements of sl2 tensor C[t]",
        [ IsSLTensorCtElement and IsSLTensorCtElementRep  ],
        function( x )
    
    local coeffs, monoms, l, i;
    coeffs := x![1];
    monoms := x![2];
    l := Length( coeffs );
    
    for i in [1..l-1] do
        Print( coeffs[i], "*x(", monoms[i][1],",",monoms[i][2],")⊗ ", "t", 
               "^",monoms[i][3],"+" );
    od;
    
    Print( coeffs[l], "*x(", monoms[l][1],",",monoms[l][2],")⊗ ",
           "t^",monoms[l][3] );
end );

InstallMethod( \+,
        "For elements of sl2 tensor C[t]",
        [ IsSLTensorCtElement and IsSLTensorCtElementRep, 
          IsSLTensorCtElement and IsSLTensorCtElementRep ],
        function( x, y )
    
    local coeffsx, coeffsy, monomsx, monomsy, res, i, pos, coeff;
    
    coeffsx := x![1];
    coeffsy := y![1];
    
    monomsx := x![2];
    monomsy := y![2];
    
    res := [];
    
    for i in [1..Length( coeffsx )] do
        if monomsx[i] in monomsy then
            pos := Position( monomsy, monomsx[i] );
            coeff := coeffsx[i] + coeffsy[pos];
        else
            coeff := coeffsx[i];
        fi;
        
        Append( res, [ coeff, monomsx[i] ] );
    od;
    
    return SLTensorCtElement( res );
end );

InstallMethod( \*,
        "scalar multiple of elements of sl2 tensor C[t]",
        [ IsRat, 
          IsSLTensorCtElement and IsSLTensorCtElementRep ],
        function( a, x )
    
    local coeffsx, monomsx, res, i;
    
    coeffsx := x![1];
    monomsx := x![2];
    
    res := [];
    
    for i in [1..Length( coeffsx )] do
        Append( res, [ a*coeffsx[i], monomsx[i] ] );
    od;
    
    return SLTensorCtElement( res );
end );


    
