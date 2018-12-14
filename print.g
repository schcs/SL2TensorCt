MonomString := function( monom ) 
    
    if monom{[1,2]} = [1,1] then
        return Concatenation( "h(1,2)⊗ t^", String( monom[3] ));
    elif monom{[1,2]} = [3,3] then
        return Concatenation( "h(1,3)⊗ t^", String( monom[3] )); 
    else
        return Concatenation( "x(", String( monom[1] ), ",", 
                       String( monom[2] ), ")⊗ t^", String( monom[3] ));;
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
    
    if AbsInt( coeffs[1] ) <> 1 then 
        Print( coeffs[1], "*" );
    fi;
    
    if coeffs[1] = -1 then
        Print( "-" );
        fi;
        
    Print( MonomString( monoms[1] ));
    
    for i in [2..l] do
        if coeffs[i] > 0 then
            Print( "+" );
        else
            Print( "-" );
        fi;
        
        if AbsInt( coeffs[i] ) <> 1 then 
            Print( AbsInt( coeffs[i] ), "*" );
        fi;
        
        Print( MonomString( monoms[i] ));
    od;

end );

InstallMethod( PrintObj,
        "For PBW elements of sl2 tensor C[t]",
        [ IsSLTensorCtPBWElement and IsSLTensorCtPBWElementRep  ],
        function( x )
    
    local coeffs, monoms, l, i, j;
        
    coeffs := x![1];
    monoms := x![2];
    l := Length( coeffs );
    
    if l = 0 then
        Print( "0" );
        return;
    fi;
    
    if AbsInt( coeffs[1] ) <> 1 then 
        Print( coeffs[1], "*" );
    fi;
    
    if coeffs[1] = -1 then
        Print( "-" );
        fi;
        
    for i in [1..Length( monoms[1] )] do 
        Print( "[", MonomString( monoms[1][i] ), "]" );
    od;
    
    for i in [2..l] do
        if coeffs[i] > 0 then
            Print( "+" );
        else
            Print( "-" );
        fi;
        
        if AbsInt( coeffs[i] ) <> 1 then 
            Print( AbsInt( coeffs[i] ), "*" );
        fi;
        for j in [1..Length( monoms[1] )] do 
            Print( "[", MonomString( monoms[i][j] ), "]" );
        od;
    od;

end );
