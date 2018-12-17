DeclareCategory( "IsSLTensorCtBasisElement", IsVector );
DeclareRepresentation( "IsSLTensorCtBasisElementRep",
        IsPositionalObjectRep, [1]);
BindGlobal( "SLTensorCtBasisFamily", NewFamily( "SLTensorCtBasisFamily" ));

DeclareCategory( "IsSLTensorCtElement", IsVector );
DeclareRepresentation( "IsSLTensorCtElementRep", IsPositionalObjectRep, [1,2]);
BindGlobal( "SLTensorCtFamily", NewFamily( "SLTensorCtFamily" ));

DeclareCategory( "IsSLTensorCtPBWElement", IsVector );
DeclareRepresentation( "IsSLTensorCtPBWElementRep",
        IsPositionalObjectRep, [1,2]);
BindGlobal( "SLTensorCtPBWFamily", NewFamily( "SLTensorCtPBWFamily" ));

DeclareCategory( "IsSLTensorCtPBWMonomial", IsVector );
DeclareRepresentation( "IsSLTensorCtPBWMonomialRep",
        IsPositionalObjectRep, [1]);
BindGlobal( "SLTensorCtPBWMonomialFamily",
        NewFamily( "SLTensorCtPBWMonomialFamily" ));

