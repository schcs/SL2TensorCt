## The following file contains some functions that compute elements with
## the Lie (super)algebra sl(1|2) ⊗ C[t]


# This function is to objectify a list into a SL2TensorCt element
# gap> SLTensorCtPBWElement( [-1,[[1,2,3],[2,3,1]],3,[[2,1,2],[3,2,-1]]] );
# -[x(1,2)⊗ t^3][x(2,3)⊗ t^1]+3*[x(2,1)⊗ t^2][x(3,2)⊗ t^-1]

SLTensorCtPBWElement := function( list )
    
    local length;
    
    if list = [] then
        list := [[],[]];
    fi;
    
    length := Length( list );
    
    if Length( list ) >= 1 and not IsList( list[1] ) then
    
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

ListCoefficientsSLTensorPBWElement := function( el )

	return el![1];
end;
	

ListMonomialsSLTensorCtPBWElement := function( el )
	
	local mons;
	
	mons := StructuralCopy( el![2] );	
	return List( mons, SLTensorCtPBWMonomial );
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
    
InstallMethod( \=,
        "For PBW elements of sl2 tensor C[t]",
        [ IsSLTensorCtPBWMonomial and IsSLTensorCtPBWMonomialRep, 
          IsSLTensorCtPBWMonomial and IsSLTensorCtPBWMonomialRep ],
        function( x, y )
    
    local el;
    
    el := 1*x-1*y;
    return Length( el![1] ) = 0; end );


InstallMethod( ZeroOp,
        "For PBW elements of sl2 tensor C[t]",
        [ IsSLTensorCtPBWElement and IsSLTensorCtPBWElementRep ],
        x -> SLTensorCtPBWElement( [[],[]] ));

InstallMethod( ZeroOp,
        "For PBW elements of sl2 tensor C[t]",
        [ IsSLTensorCtPBWMonomial and IsSLTensorCtPBWMonomialRep ],
        x -> SLTensorCtPBWElement( [[],[]] ));

InstallMethod( AdditiveInverseOp,
        "For PBW elements of sl2 tensor C[t]",
        [ IsSLTensorCtPBWElement and IsSLTensorCtPBWElementRep ],
        
        function( x )
        return (-1)*x; end );

InstallMethod( \*,
        "scalar multiple of PBW monomial of sl2 tensor C[t]",
        [ IsRat, 
          IsSLTensorCtPBWMonomial and IsSLTensorCtPBWMonomialRep ],
        function( a, x )
    if a = 0 then
        return SLTensorCtPBWElement( [[],[]] );
    else
        return SLTensorCtPBWElement( [[a],[SplitSLTensorCtPBWMonomial( x )]] );
    fi; 
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

IsCollectedSLTensorCtPBWMonomial := function( mon )
    local els, i;
    
    els := List( SplitSLTensorCtPBWMonomial( mon ), SLTensorCtBasisElement );
    
    for i in [1..Length( els )-1] do
        if els[i] > els[i+1] then 
            return i;
        elif els[i] = els[i+1] and Parity( els[i] ) = 1 then
            return i;
        fi;
    od;
    
    return true;
end;

ProductSLTensorCtPBWElements := function( x, y )
            
    local coeffsx, coeffsy, monomsx, monomsy, res, i, j, coeff, monom;
    
    coeffsx := StructuralCopy( x![1] );
    coeffsy := StructuralCopy( y![1] );
    monomsx := StructuralCopy( x![2] );
    monomsy := StructuralCopy( y![2] );
    res := [[],[]];
    for i in [1..Length( coeffsx )] do
        for j in [1..Length( coeffsy )] do
            coeff := coeffsx[i]*coeffsy[j];
            monom := StructuralCopy( monomsx[i] );
            Append( monom, monomsy[j] );
            Add( res[1], coeff );
            Add( res[2], monom );
        od;
    od;
    
    res := CollectSLTensorCtPBWElement( SLTensorCtPBWElement( res ));
    return res;
end;

# IsCollectedSLTensorCtPBWMonomial(m);
# 1
# CollectSLTensorCtPBWMonomialAtPosition(m,1);
# [x(3,1)⊗ t^3][x(1,2)⊗ t^2]-[x(3,2)⊗ t^5] (the result should be - this)

CollectSLTensorCtPBWMonomialAtPosition := function( mon, pos )
    local monlist, list1, list2, x1, y1, px, py, z, mon1, mon2, 
          y1x1, zz, e1, e2, newel; 
   
    monlist := SplitSLTensorCtPBWMonomial( mon );
    list1 := monlist{[1..pos-1]};
    list2 := monlist{[pos+2..Length( monlist )]};
    
    x1 := SLTensorCtBasisElement( monlist[pos] );
    y1 := SLTensorCtBasisElement( monlist[pos+1] );
    px := Parity( x1 ); py := Parity( y1 );
    
    if x1 = y1 then
        return 0*mon;
    fi;
    
    z := y1*x1;
    
    mon1 := 1*SLTensorCtPBWMonomial( list1 );
    mon2 := 1*SLTensorCtPBWMonomial( list2 );
    y1x1 := 1*SLTensorCtPBWMonomial( [monlist[pos+1],monlist[pos]] );
    zz := SLTensorCtPBWElement( [ z![1], List( z![2], x->[x])]);    
    
    newel := (-1)^(px*py)*(
           ProductSLTensorCtPBWElements(        
           ProductSLTensorCtPBWElements( mon1, y1x1 ), mon2 )-
           ProductSLTensorCtPBWElements(        
           ProductSLTensorCtPBWElements( mon1, zz ), mon2 ));
    
    return CollectSLTensorCtPBWElement( newel );
    
end;



CollectSLTensorCtPBWElement2 := function( el )

	local newel, coeffs, mons, i, iscollected, c, count;
		
	newel := StructuralCopy( el );	
 
        count := 0;
        repeat
	    iscollected := true;
		coeffs := ListCoefficientsSLTensorPBWElement( newel );
		mons := ListMonomialsSLTensorCtPBWElement( newel );
	
		for i in [1..Length(mons)] do
			c := IsCollectedSLTensorCtPBWMonomial( mons[i] );
			if c <> true then
		           iscollected := false;
                           newel := newel + coeffs[i]*((-1)*mons[i]+
                          CollectSLTensorCtPBWMonomialAtPosition( mons[i], c ));
                          #Print( mons[i], "\n", c, "\n", newel, "\n\n" ); 
                           break; 
			fi;
		od;
                newel := CollectSLTensorCtPBWElement( newel );
                count := count+1;
        until iscollected;
        Print( count, "\n" );    
        if false in List( ListMonomialsSLTensorCtPBWElement( newel ), 
                   IsCollectedSLTensorCtPBWMonomial ) then
            Error( "element is not properly collected" );
        fi;
			
	return newel;
end;
	
	
InstallMethod( \*,
        "product PBW elements of sl2 tensor C[t]",
        [ IsSLTensorCtPBWElement and IsSLTensorCtPBWElement, 
          IsSLTensorCtPBWElement and IsSLTensorCtPBWElementRep ],
        function( x, y )
    
    return CollectSLTensorCtPBWElement2( ProductSLTensorCtPBWElements( x, y ));
end );

    
TestCollection := function( l0 )
        
    local els, el, l1, l2, mons, coefs, i, e1, mon, w1, w2, u, w;
    
    els := [];
    
    for i in [1..3] do
        l1 := Random( [1..l0] );
        l2 := Random( [1..l0] );
        
        mons := List( [1..l1], x->List( [1..l2], x->Random( 
                        [[1,1],[1,2],[1,3],[2,1],[2,3],[3,1],[3,2],[3,3]] )));
        coefs := List( [1..l1], x->Random( Rationals ));
        
        for e1 in mons do
            for mon in e1 do
                Add( mon, Random( [1..5] ));
            od;
        od;
        els[i] := Sum( List( [1..Length(coefs)], 
                          x->coefs[x]*SLTensorCtPBWMonomial( mons[x] )));
    od;
    
    if (els[1]*els[2])*els[3] <> els[1]*(els[2]*els[3]) then
        Error( 101 );
    fi;
    
    return true;
end;
