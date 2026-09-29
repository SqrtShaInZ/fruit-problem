mkc := function(n)
    return EllipticCurve([0, 4*n^2 + 12*n - 3, 0, 32*(n + 3), 0]);
end function;
mkc2 := function(n)
    return EllipticCurve([0, 4*n^2 + 12*n - 3, 0, -128*n - 384, -512*n^3 - 3072*n^2 - 4224*n + 1152]);
end function;
mkc3 := function(n)
    return EllipticCurve([0, 4*n^2 + 12*n - 3, 0, -320*n^2 - 1248*n - 1104, -1024*n^4 - 7168*n^3 - 18944*n^2 - 24576*n - 15040]);
end function;
mkc6 := function(n)
    return EllipticCurve([0, 4*n^2 + 12*n - 3, 0, -640*n^3 - 6080*n^2 - 18528*n - 18384, -2048*n^5 - 39936*n^4 - 281088*n^3 - 935936*n^2 - 1503744*n - 941248]);
end function;

SquarefreePart := function(n)
    if n eq 0 then return 0; end if;
    sign := n lt 0 select -1 else 1;
    n := AbsoluteValue(n);
    f := Factorization(n);
    sqfree := 1;
    for p in f do
        if p[2] mod 2 eq 1 then
            sqfree *:= p[1];
        end if;
    end for;
    return sign * sqfree;
end function;

GetCurve := function(n, iso)
    case iso:
        when 1: return mkc(n);
        when 2: return mkc2(n);
        when 3: return mkc3(n);
        when 6: return mkc6(n);
        else error "Invalid isogenous configuration requested.";
    end case;
end function;

SolveIsogeny2 := function(E_domain, E_codomain, P)
    f2 := DivisionPolynomial(E_domain, 2);
    roots := Roots(f2, RationalField());
    Poly<x> := PolynomialRing(RationalField());
    for r in roots do
        kernel_poly := x - r[1];
        E_cod, phi := IsogenyFromKernel(E_domain, kernel_poly);
        is_iso, iso := IsIsomorphic(E_cod, E_codomain);
        if is_iso then
            P_cod := P @@ iso;
            return DualIsogeny(phi)(P_cod);
        end if;
    end for;
    error "SolveIsogeny2: No matching 2-isogeny found between the curves.";
end function;

SolveIsogeny3 := function(E_domain, E_codomain, P)
    f3 := DivisionPolynomial(E_domain, 3);
    roots := Roots(f3, RationalField());
    Poly<x> := PolynomialRing(RationalField());
    for r in roots do
        kernel_poly := x - r[1];
        E_cod, phi := IsogenyFromKernel(E_domain, kernel_poly);
        is_iso, iso := IsIsomorphic(E_cod, E_codomain);
        if is_iso then
            P_cod := P @@ iso;
            return DualIsogeny(phi)(P_cod);
        end if;
    end for;
    error "SolveIsogeny3: No matching 3-isogeny found between the curves.";
end function;

map2 := function(n, P)
    return SolveIsogeny2(mkc(n), mkc2(n), P);
end function;

map3 := function(n, P)
    return SolveIsogeny3(mkc(n), mkc3(n), P);
end function;

map6 := function(n, P)
    P_mkc2 := SolveIsogeny3(mkc2(n), mkc6(n), P);
    return SolveIsogeny2(mkc(n), mkc2(n), P_mkc2);
end function;

MyTwoDescent := function(E : RemoveGens:=[])
    S2, S2map := TwoSelmerGroup(E : RemoveTorsion, RemoveGens:=RemoveGens);
    S2elts := [s : s in S2 | not IsIdentity(s)];
    TD := [];
    TDmaps := [* *];
    for s in S2elts do 
        C2, C2toE := TwoCover(s @@ S2map : E:=E);
        Append(~TD, C2);
        Append(~TDmaps, C2toE);
    end for;
    n := Ngens(S2);
    inds := [Index(S2elts, S2.i) : i in [1..n]];
    CTmat := ZeroMatrix(GF(2), n, n);
    for i in [1..n], j in [1..i-1] do 
        ii := inds[i];
        jj := inds[j];
        CTij := CasselsTatePairing(TD[ii], TD[jj]);
        CTmat[i,j] := CTij;  
        CTmat[j,i] := CTij;
    end for; 
    printf "The Cassels-Tate pairing on Sel(2,E)/E[2] is\n%o\n", CTmat;
    error if not IsEven(Rank(CTmat)), "CasselsTatePairing failed the parity check";
    kernel := [v : v in Kernel(CTmat) | not IsZero(v)];
    kernel_inds := [Index(S2elts, S2!Eltseq(v)) : v in kernel];
    S2_4 := [S2elts[ii] : ii in kernel_inds];
    TDCT := [TD[ii] : ii in kernel_inds];
    return TDCT;
end function;

MyThreeDescent_dedup := function(Crv3, map3)
    P<x1,x2,x3> := PolynomialRing(Rationals(), 3);
    // the four "diagonal sign" substitutions (identity + three)
    subs := [[x1,x2,x3], [x1,-x2,-x3], [-x1,x2,-x3], [-x1,-x2,x3]];
    homs := [hom<P -> P | s> : s in subs];

    kept_pols := [P | ];
    kept_mods := [];          // GenusOneModel cache for the fallback
    Crv3_real := [];
    map3_real := [];

    for i in [1..#Crv3] do
        f := P ! DefiningPolynomial(Crv3[i]);
        dup := false;
        // fast test: equal to a kept cover up to a sign substitution and overall sign
        for g in kept_pols do
            for h in homs do
                if f eq h(g) or f eq -h(g) then dup := true; break; end if;
            end for;
            if dup then break; end if;
        end for;
        // slow fallback: general equivalence over Q (only for survivors)
        if not dup then
            M := GenusOneModel(f);
            for M2 in kept_mods do
                if IsEquivalent(M, M2) then dup := true; break; end if;
            end for;
        end if;
        if not dup then
            Append(~kept_pols, f);
            Append(~kept_mods, GenusOneModel(f));
            Append(~Crv3_real, Crv3[i]);
            map3_real := map3_real cat [map3[i]];
        end if;
    end for;
    return Crv3_real, map3_real;
end function;

MyThreeDescent := function(number, isogenous : NoMap:=false, HighRank:=false)
    if isogenous mod 3 eq 0 then
        return ThreeDescent(GetCurve(number, isogenous));
    end if;
    verb := GetVerbose("Minimisation");
    SetVerbose("Minimisation", 0);
    E := GetCurve(number, isogenous);
    E_phi := GetCurve(number, isogenous * 3);
    if NoMap and not HighRank then
        _, _, Crv3 := pIsogenyDescent(E_phi, 3);
        if #Crv3 gt 0 then
            SetVerbose("Minimisation", verb);
            return Crv3;
        end if;
    end if;
    Crv3, map3, Crv3_phi, map3_phi, isog := ThreeIsogenyDescent(E);
    map3 := [m * DualIsogeny(isog) : m in map3];
    for i in [1..#Crv3_phi] do
        try
            C3, m3 := pIsogenyDescent(Crv3_phi[i], E, E_phi);
            Crv3 := Crv3 cat C3;
            map3 := map3 cat [m * map3_phi[i] : m in m3];
        catch e
            ;
        end try;
    end for;
    SetVerbose("Minimisation", verb);
    if #Crv3 gt 0 then return MyThreeDescent_dedup(Crv3, map3); end if;
    return ThreeDescent(E);
end function;

ComputeGeneratorTS := function(number, isogenous, reg : TwoPowerDescOnly:=false, NoEightDesc:=false, HyperE:=[], descent_no:=0, Crv3:=[], mapA:=[])
    E := GetCurve(number, isogenous);
    if descent_no eq 0 then
        if (TwoPowerDescOnly and reg lt 135) or (not TwoPowerDescOnly and reg lt 120) then 
            descent_no := 4;
        elif reg lt 270 and not TwoPowerDescOnly then 
            descent_no := 6;
        elif reg lt 480 and not NoEightDesc then
            descent_no := 8;
        else // still have to do 12-descent since the point is so high
            descent_no := 12;
        end if;
        printf "Automatically selected descent depth: %o-descent\n", descent_no;
    end if;
    if #HyperE eq 0 then
        HyperE := MyTwoDescent(E);
    end if;
    case descent_no:
        when 4:
            Crvs4 := [];
            if #HyperE eq 0 then
                Crvs4 := FourDescent(E);
            else
                for HE in HyperE do
                    Crv4 := FourDescent(HE : RemoveTorsion := true);
                    Crvs4 := Crvs4 cat Crv4;
                end for;
            end if;
            index := 1;
            Ps := [];
            bound := Max(10^5, Round(10^(reg / 20 + 2)));
            while #Ps eq 0 do
                if index gt #Crvs4 then
                    Crvs4 := [];
                    if #HyperE eq 0 then
                        Crvs4 := FourDescent(E);
                    else
                        for HE in HyperE do
                            Crv4 := FourDescent(HE : RemoveTorsion := true);
                            Crvs4 := Crvs4 cat Crv4;
                        end for;
                    end if;
                    bound := Min(bound * 10, 10^10);
                    index := 1;
                end if;
                Ps := PointsQI(Crvs4[index], bound : OnlyOne := true);
                index := index + 1;
            end while;
            index := index - 1;
            A, m := AssociatedEllipticCurve(Crvs4[index] : E := E);
            P := Saturation([m(Ps[1])], 1000 : TorsionFree := true)[1];
        when 6:
            if #Crv3 eq 0 or #mapA eq 0 then
                Crv3, mapA := MyThreeDescent(number, isogenous);
            end if;
            P6 := [];
            index_2 := 1;
            index_3 := 1;
            bound := Max(10^2, Round(10^(reg / 25 + 1)));
            flag := false;
            while #P6 eq 0 do
                if index_2 gt #HyperE then
                    bound := bound * 10;
                    index_2 := 1;
                    index_3 := 1;
                end if;
                while #P6 eq 0 do
                    if index_3 gt #Crv3 then break; end if;
                    Crv, map6to3 := SixDescent(HyperE[index_2], Crv3[index_3]);
                    P6 := PointSearch(Crv, bound : OnlyOne := true);
                    index_3 := index_3 + 1;
                end while;
                if #P6 eq 0 and index_3 gt #Crv3 then
                    index_2 := index_2 + 1;
                    index_3 := 1;
                end if;
            end while;
            winning_index := index_3 - 1;
            P3_internal := map6to3(Domain(map6to3) ! P6[1]);
            winning_map := mapA[winning_index];
            comps := Components(winning_map);
            if #comps le 1 then
                P := Saturation([winning_map(P3_internal)], 1000 : TorsionFree := true)[1];
            else
                P3_fixed := Domain(comps[1]) ! Eltseq(P3_internal);
                P_intermediate := comps[1](P3_fixed);
                P_intermediate_fixed := Domain(comps[2]) ! Eltseq(P_intermediate);
                PE := comps[2](P_intermediate_fixed);
                P := Saturation([PE], 1000 : TorsionFree := true)[1];
            end if;
        when 8:
            index := 1;
            Crvs4 := [];
            if #HyperE eq 0 then
                Crvs4 := FourDescent(E);
            else
                for HE in HyperE do
                    Crv4 := FourDescent(HE : RemoveTorsion := true);
                    Crvs4 := Crvs4 cat Crv4;
                end for;
            end if;
            Ps := [];
            print "Doing preliminary PointsQI's...\n";
            while #Ps eq 0 do
                if index gt #Crvs4 then break; end if;
                Ps := PointsQI(Crvs4[index], Min(10^8, Max(10^5, Round(10^(reg / 400 + 2)))) : OnlyOne := true);
                index := index + 1;
            end while;
            if #Ps gt 0 then
                A, m := AssociatedEllipticCurve(Crvs4[index - 1] : E := E);
                P := Saturation([m(Ps[1])], 1000 : TorsionFree := true)[1];
            else
                Crvs8, maps8 := EightDescent(Crvs4[1]);
                for i in [2..#Crvs4] do
                    Crv8, map8 := EightDescent(Crvs4[i]);
                    Crvs8 := Crvs8 cat Crv8;
                    maps8 := maps8 cat map8;
                end for;
                bound := Max(10^2, Round(10^(reg / 30 + 1)));
                index := 1;
                P8 := [];
                while #P8 eq 0 do
                    if index gt #Crvs8 then break; end if;
                    P8 := PointSearch(Crvs8[index], Round(bound^(3/11)) : OnlyOne := true);
                    if #P8 eq 0 then
                        P8 := PointSearch(Crvs8[index], bound : OnlyOne := true);
                    end if;
                    index := index + 1;
                end while;
                if #P8 gt 0 then
                    index := index - 1;
                    P4 := maps8[index](P8[1]);
                    A, m := AssociatedEllipticCurve(Codomain(maps8[index]) : E := E);
                    P := Saturation([m(P4)], 1000 : TorsionFree := true)[1];
                else
                    return [0];
                end if;
            end if;
        when 12:
            if #Crv3 eq 0 then
                Crv3 := MyThreeDescent(number, isogenous : NoMap);
            end if;
            index := 1;
            Crvs4 := [];
            if #HyperE eq 0 then
                Crvs4 := FourDescent(E);
            else
                for HE in HyperE do
                    Crv4 := FourDescent(HE : RemoveTorsion := true);
                    Crvs4 := Crvs4 cat Crv4;
                end for;
            end if;
            Ps := [];
            print "Doing preliminary PointsQI's...\n";
            while #Ps eq 0 do
                if index gt #Crvs4 then break; end if;
                Ps := PointsQI(Crvs4[index], Min(10^8, Max(10^5, Round(10^(reg / 400 + 2)))) : OnlyOne := true);
                index := index + 1;
            end while;
            if #Ps gt 0 then
                A, mapA := AssociatedEllipticCurve(Crvs4[index - 1] : E := E);
                P := Saturation([mapA(Ps[1])], 1000 : TorsionFree := true)[1];
            else
                Crvs12 := [];
                maps12 := [];
                for C3 in Crv3 do
                    // avoid computing fiber product of two covers that produce different generators
                    if #PointSearch(C3, 10^4) gt 0 then
                        continue;
                    end if;
                    for C4 in Crvs4 do
                        Crv12, map12 := TwelveDescent(C3, C4);
                        Crvs12 := Crvs12 cat Crv12;
                        maps12 := maps12 cat map12;
                    end for;
                end for;
                P12 := [];
                index := 1;
                bound := Max(10^2, Round(10^(reg / 60 + 1)));
                while #P12 eq 0 do
                    while #P12 eq 0 and index le #Crvs12 do
                        P12 := PointSearch(Crvs12[index], bound : OnlyOne);
                        index := index + 1;
                    end while;
                    if #P12 gt 0 then
                        break;
                    end if;
                    if index gt #Crvs12 then
                        index := 1;
                    end if;
                    bound := bound * 10;
                end while;
                index := index - 1;
                P4 := maps12[index](P12[1]);
                A, m := AssociatedEllipticCurve(Codomain(maps12[index]) : E := E);
                P := Saturation([m(P4)], 1000 : TorsionFree := true)[1];
            end if;
    end case;
    if isogenous eq 1 then
        P_orig := P;
    elif isogenous eq 2 then
        P_orig := map2(number, P);
    elif isogenous eq 3 then
        P_orig := map3(number, P);
    elif isogenous eq 6 then
        P_orig := map6(number, P);
    end if;
    E_orig := mkc(number);
    P_final := Saturation([E_orig ! P_orig], 1000 : TorsionFree := true)[1];
    print "-----------------------------------------";
    print "True Canonical Height on E_", number, ":", CanonicalHeight(P_final);
    print "Ratio against BSD-predicted regulator:", RealField(15)!(CanonicalHeight(P) / reg);
    print "-----------------------------------------";
    return Eltseq(P_final)[1..2];
end function;
Sha4Order := function(N)
    p := Valuation(N, 2);
    m := N div 2^p;
    q := Ilog(2, m + 1);
    error if 2^q ne m + 1,
        Sprintf("#FourDescent = %o is not of the form 2^p*(2^q-1)", N);
    error if q - 1 gt p or IsOdd(q - 1 - p),
        Sprintf("(p,q) = (%o,%o) violates Cassels-Tate constraints", p, q);
    return 2^(p + q - 1);
end function;
TSSize := function(number, isogenous)
    Sel4_size := #FourDescent(GetCurve(number, isogenous));
    TS_order := Sha4Order(Sel4_size);
    Crv3, mapA := MyThreeDescent(number, isogenous);
    if isogenous mod 3 eq 0 then
        TS_order := TS_order * (2 * #Crv3 + 1) / 3;
    else
        if #Crv3 gt 4 then TS_order := TS_order * 9; end if;
    end if;
    return TS_order, Crv3, mapA;
end function;

function BSDEasyTermsQ(E : Precision := 6)
    p_prec := Max(Precision, 20);
    dsc := Integers()!Discriminant(E);
    c_inf := (dsc gt 0) select 2 else 1;
    omegas := c_inf * RealPeriod(E : Precision := p_prec);
    loc := LocalInformation(E);
    tam := &*[ Integers() | l[4] : l in loc ];
    om := &*[ Rationals() | l[1]^((Valuation(dsc, l[1]) - l[2]) div 12) : l in loc ];
    tors := #TorsionSubgroup(E);
    return (tam * om * omegas) / (tors^2);
end function;

ComputeGeneratorTrial := function(number, isogenous)
    E := GetCurve(number, isogenous);
    HyperE := MyTwoDescent(E)[1];
    Crv4 := FourDescent(HyperE : RemoveTorsion);
    for C4 in Crv4 do
        P4 := PointsQI(C4, 10^7 : OnlyOne);
        if #P4 gt 0 then
            A, m := AssociatedEllipticCurve(C4 : E := E);
            P := Saturation([m(P4[1])], 1000 : TorsionFree := true)[1];
            return [P];
        end if;
    end for;
    Crv3, mapA := MyThreeDescent(number, isogenous);
    P6 := [];
    index_2 := 1;
    index_3 := 1;
    flag := false;
    while #P6 eq 0 do
        if index_3 gt #Crv3 then break; end if;
        Crv, map6to3 := SixDescent(HyperE, Crv3[index_3]);
        P6 := PointSearch(Crv, 10^8 : OnlyOne := true);
        index_3 := index_3 + 1;
    end while;
    if #P6 gt 0 then
        winning_index := index_3 - 1;
        P3_internal := map6to3(Domain(map6to3) ! P6[1]);
        winning_map := mapA[winning_index];
        comps := Components(winning_map);
        if #comps le 1 then
            P := Saturation([winning_map(P3_internal)], 1000 : TorsionFree := true)[1];
        else
            P3_fixed := Domain(comps[1]) ! Eltseq(P3_internal);
            P_intermediate := comps[1](P3_fixed);
            P_intermediate_fixed := Domain(comps[2]) ! Eltseq(P_intermediate);
            PE := comps[2](P_intermediate_fixed);
            P := Saturation([PE], 1000 : TorsionFree := true)[1];
        end if;
        return [P];
    end if;
    return [];
end function;

ComputeGenerator := function(number, isogenous, descent_no : TwoPowerDescOnly:=true, NoEightDesc:=false, MaxReg:=0)
    if RootNumber(mkc(number)) eq 1 then
        error "Rank of curve must be 1";
    end if;
    if isogenous eq 0 then
        verb_3desc := GetVerbose("ThreeDescent");
        verb_selmer := GetVerbose("Selmer");
        SetVerbose("ThreeDescent", 0);
        SetVerbose("Selmer", 0);
        reg1 := 1 / BSDEasyTermsQ(mkc(number)) / TSSize(number, 1);
        reg6 := 1 / BSDEasyTermsQ(mkc6(number)) / TSSize(number, 6);
        SetVerbose("ThreeDescent", verb_3desc);
        SetVerbose("Selmer", verb_selmer);
        ratio := reg1 / reg6;
        frac := BestApproximation(RealField(10)!ratio, 100);
        printf "Ratio = %o\n", ratio;
        isogenous := SquarefreePart(Numerator(frac));
        printf "Automatically determined isogenous target curve index: %o\n", isogenous;
    end if;
    E := GetCurve(number, isogenous);
    P := ComputeGeneratorTrial(number, isogenous);
    if #P gt 0 then
        P := P[1];
        if isogenous eq 1 then
            P_orig := P;
        elif isogenous eq 2 then
            P_orig := map2(number, P);
        elif isogenous eq 3 then
            P_orig := map3(number, P);
        elif isogenous eq 6 then
            P_orig := map6(number, P);
        end if;
        E_orig := mkc(number);
        P_final := Saturation([E_orig ! P_orig], 1000 : TorsionFree := true)[1];
        print "-----------------------------------------";
        print "True Canonical Height on E_", number, ":", CanonicalHeight(P_final);
        print "-----------------------------------------";
        return Eltseq(P_final)[1..2];
    end if;
    rnk, lead := AnalyticRank(E : Precision := 6);
    if rnk gt 1 then
        error "Rank of curve must be 1";
    end if;
    TS_order, Crv3, mapA := TSSize(number, isogenous);
    reg := lead / BSDEasyTermsQ(E) / TS_order;
    printf "Regulator is %o based on BSD formula\n", reg;
    printf "Tate-Shafarevich group has order %o\n", TS_order;
    if MaxReg ne 0 and reg gt MaxReg then
        printf "Regulator too large; aborting\n";
        return [reg, TS_order];
    end if;
    HyperE := MyTwoDescent(E);
    return ComputeGeneratorTS(number, isogenous, reg : TwoPowerDescOnly:=(Valuation(TS_order, 3) gt 0), NoEightDesc:=NoEightDesc, \
        HyperE:=HyperE, descent_no:=descent_no, Crv3:=Crv3, mapA:=mapA);
end function;

SaturateGeneratorFull := function(number, isogenous, gens)
    if isogenous eq 1 then
        gens_orig := gens;
    elif isogenous eq 2 then
        gens_orig := [map2(number, P) : P in gens];
    elif isogenous eq 3 then
        gens_orig := [map3(number, P) : P in gens];
    elif isogenous eq 6 then
        gens_orig := [map6(number, P) : P in gens];
    end if;
    E_orig := mkc(number);
    gens_orig := Saturation(gens_orig, 1000 : TorsionFree := true);
    print "-----------------------------------------";
    print "Verification status:", [(P in E_orig) : P in gens_orig];
    print "Regulator of E_", number, ":", Determinant(HeightPairingMatrix(gens_orig));
    print "-----------------------------------------";
    return [[Eltseq(P)[1], Eltseq(P)[2]] : P in gens_orig];
end function;

forward ComputeGeneratorFull;

ComputeGeneratorFull := function(number, isogenous, rank : known_gens:=[], TwoPowerDescOnly:=true, NoEightDesc:=false)
    if isogenous eq 1 then
        E := mkc(number);
    elif isogenous eq 2 then
        E := mkc2(number);
    elif isogenous eq 3 then
        E := mkc3(number);
    elif isogenous eq 6 then
        E := mkc6(number);
    else
        error "Invalid isogenous configuration requested.";
    end if;
    if isogenous mod 3 eq 0 then
        NoEightDesc:=true;
    end if;
    if #known_gens gt 0 then
        known_gens := [E!P : P in known_gens];
        gens := Saturation(known_gens, 1000 : TorsionFree := true);
    else
        gens := [];
    end if;
    if #gens eq rank then
        return SaturateGeneratorFull(number, isogenous, gens);
    end if;
    HyperE := MyTwoDescent(E : RemoveGens := gens);
    for HE in HyperE do
        Ps := RationalPoints(HE : Bound := 162755);
        if #Ps gt 0 then
            A, mapA := AssociatedEllipticCurve(HE : E := E);
            P := mapA(Ps[1]);
            printf "Found rational point %o\n", Eltseq(P)[1..2];
            Append(~gens, P);
            return ComputeGeneratorFull(number, isogenous, rank : known_gens:=gens, TwoPowerDescOnly:=TwoPowerDescOnly, NoEightDesc:=NoEightDesc);
        end if;
    end for;
    Crv4 := [];
    lim_4desc := 3;
    if isogenous mod 3 eq 0 then lim_4desc := 1; end if;
    for a in [1..lim_4desc] do // try three times
        HyperE := MyTwoDescent(E : RemoveGens := gens);
        Crv4 := [];
        for HE in HyperE do
            Crv4 := Crv4 cat FourDescent(HE : RemoveTorsion := true, RemoveGensEC := gens);
        end for;
        for C4 in Crv4 do
            Ps := PointsQI(C4, 10^7 : OnlyOne);
            if #Ps gt 0 then
                A, mapA := AssociatedEllipticCurve(C4 : E := E);
                P := mapA(Ps[1]);
                printf "Found rational point %o\n", Eltseq(P)[1..2];
                Append(~gens, P);
                return ComputeGeneratorFull(number, isogenous, rank : known_gens:=gens, TwoPowerDescOnly:=TwoPowerDescOnly, NoEightDesc:=NoEightDesc);
            end if;
        end for;
    end for;
    Crvs4 := [];
    if not NoEightDesc then
        Crvs8, maps8 := EightDescent(Crv4[1]);
        if #Crvs8 gt 0 then
            Crvs4 := [Crv4[1]];
        end if;
        for i in [2..#Crv4] do
            Crv8, map8 := EightDescent(Crv4[i]);
            Crvs8 := Crvs8 cat Crv8;
            maps8 := maps8 cat map8;
            if #Crv8 gt 0 then
                Append(~Crvs4, Crv4[i]);
            end if;
        end for;
        P8 := [];
        for bound in [10^11, 10^15] do
            index := 1;
            while #P8 eq 0 do
                if index gt #Crvs8 then break; end if;
                P8 := PointSearch(Crvs8[index], bound : OnlyOne := true);
                index := index + 1;
            end while;
            if #P8 gt 0 then
                index := index - 1;
                P4 := maps8[index](P8[1]);
                A, mapA := AssociatedEllipticCurve(Codomain(maps8[index]) : E := E);
                gens := Saturation(gens cat [mapA(P4)], 1000 : TorsionFree := true);
                return ComputeGeneratorFull(number, isogenous, rank : known_gens:=gens, TwoPowerDescOnly:=TwoPowerDescOnly, NoEightDesc:=NoEightDesc);
            end if;
        end for;
    else
        Crvs4 := Crv4;
    end if;
    // 4-descent fails to find a point
    Crvs3, maps3 := MyThreeDescent(number, isogenous : HighRank);
    Crvs6 := [];
    maps6 := [];
    parent_C3 := [];
    index := 1;
    for Crv3 in Crvs3 do
        // avoid computing fiber product of two covers that produce different generators
        if #PointSearch(Crv3, 10^4) gt 0 then
            index := index + 1;
            continue;
        end if;
        for HE in HyperE do
            Crv6, map6 := SixDescent(HE, Crv3);
            Crvs6 := Crvs6 cat [Crv6];
            maps6 := maps6 cat [map6];
            Append(~parent_C3, index);
        end for;
        index := index + 1;
    end for;
    P6 := [];
    index := 1;
    for bound in [10^8, 10^11] do
        while #P6 eq 0 and index le #Crvs6 do
            P6 := PointSearch(Crvs6[index], bound : OnlyOne);
            index := index + 1;
        end while;
        if #P6 gt 0 then
            break;
        end if;
        if index gt #Crvs6 then
            index := 1;
        end if;
    end for;
    if #P6 gt 0 then
        index := index - 1;
        P3_internal := maps6[index](P6[1]);
        comps := Components(maps3[parent_C3[index]]);
        if #comps le 1 then
            PE := maps3[parent_C3[index]](P3_internal);
            printf "Found rational point %o\n", Eltseq(PE)[1..2];
            gens := Saturation(gens cat [PE], 1000 : TorsionFree := true);
        else
            P3_fixed := Domain(comps[1]) ! Eltseq(P3_internal);
            P_intermediate := comps[1](P3_fixed);
            P_intermediate_fixed := Domain(comps[2]) ! Eltseq(P_intermediate);
            PE := comps[2](P_intermediate_fixed);
            printf "Found rational point %o\n", Eltseq(PE)[1..2];
            gens := Saturation(gens cat [PE], 1000 : TorsionFree := true);
        end if;
        return ComputeGeneratorFull(number, isogenous, rank : known_gens:=gens, TwoPowerDescOnly:=TwoPowerDescOnly, NoEightDesc:=NoEightDesc);
    end if;
    Crvs12 := [];
    maps12 := [];
    for Crv3 in Crvs3 do
        // avoid computing fiber product of two covers that produce different generators
        if #PointSearch(Crv3, 10^4) gt 0 then
            continue;
        end if;
        for C4 in Crvs4 do
            Crv12, map12 := TwelveDescent(Crv3, C4);
            Crvs12 := Crvs12 cat Crv12;
            maps12 := maps12 cat map12;
        end for;
    end for;
    P12 := [];
    index := 1;
    for bound in [10^7, 10^13, 10^18, 10^23] do
        while #P12 eq 0 and index le #Crvs12 do
            P12 := PointSearch(Crvs12[index], bound : OnlyOne);
            index := index + 1;
        end while;
        if #P12 gt 0 then
            break;
        end if;
        if index gt #Crvs12 then
            index := 1;
        end if;
    end for;
    index := index - 1;
    P4 := maps12[index](P12[1]);
    A, mapA := AssociatedEllipticCurve(Codomain(maps12[index]) : E := E);
    P := mapA(P4);
    printf "Found rational point %o\n", Eltseq(P)[1..2];
    gens := Saturation(gens cat [P], 1000 : TorsionFree := true);
    return ComputeGeneratorFull(number, isogenous, rank : known_gens:=gens, TwoPowerDescOnly:=TwoPowerDescOnly, NoEightDesc:=NoEightDesc);
end function;

SetMemoryLimit(2^31);
SetClassGroupBounds("GRH");
SetVerbose("ThreeDescent", 2);
SetVerbose("PointSearch", 2);
SetVerbose("TwelveDescent", 1);
SetVerbose("FourDescent", 1);
SetVerbose("EightDescent", 2);
SetVerbose("Heegner", 1);
SetVerbose("NineDescent", 2);
SetVerbose("QISearch", 1);
SetVerbose("Selmer", 2);
SetVerbose("Conic", 2);
SetVerbose("Minimisation", 1);
SetColumns(0);
SetDefaultRealField(RealField(1000));
