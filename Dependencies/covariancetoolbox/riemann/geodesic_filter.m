function filteredCovset = geodesic_filter(CovSet,C,W)

% passage dans le plan tangent
S = Tangent_space(CovSet,C);

%filtering
Out = (W*((W'*W)\W'))*S;

% passage dans la vari�t�
filteredCovset = UnTangent_space(Out,C);
