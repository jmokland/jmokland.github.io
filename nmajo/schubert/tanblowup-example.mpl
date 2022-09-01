> HomeLib := "/home/nmajo/schubertlib"; libname := libname, HomeLib;
                          "/home/nmajo/schubertlib"
             "/usr/local/maple13/lib", "/home/nmajo/schubertlib"
> with(SF);
[Par, add_basis, char2sf, conjugate, dominate, dual_basis, evalsf, hooks,

  itensor, jt_matrix, nextPar, omega, plethysm, scalar, sf2char, skew, stdeg,

  subPar, theta, toe, toh, top, tos, varset, zee]
> with(schubert);
[&-!, &-*, &/, &@, &^*, End, Grass, Hom, POINT, Proj, Symm, adams,

  additivebasis, betti, blowup, blowuppoints, bundle, bundlesection, chern,

  chi, codimension, compose, curve, determinant, dimension, division, down,

  dual, grass, grobnerbasis, insertedge, integral, integral2, koszul,

  lowershriek, lowerstar, monomials, monomialvalues, morphism, multiplepoint,

  normalbundle, normalform, o, porteous, porteous2, productvariety, proj,

  rank, schur, schurfunctor, schurfunctor2, segre, setvariety, sheaf, strip,

  symm, tangentbundle, tensor, todd, toddclass, toricvariety, totalspace,

  twist, up, upperstar, variety, verifyduality, wedge, where, whichcone, wproj

  ]
> tanblowup := proc (f) local n, lambda, codim, c, i, j, alpha, z, directimage, totalchern, class, classes; if `not`(assigned(f[source_][tangentbundle_])) then ERROR(` - source variety needs a tangent bundle`) end if; if `not`(assigned(f[target_][tangentbundle_])) then ERROR(` - target variety needs a tangent bundle`) end if; blowup(f); setvariety(f[source_]); Proj(n[f], dual(normalbundle(f)), lambda[f], all); totalspace(n[f], f[source_], tan); codim := -f[dimension_]; for i from 0 to codim do if dimension(f[source_]) < i then c[i] := 0 else c[i] := chern(i, normalbundle(f)) end if end do; i := 'i'; j := 'j'; sum(t^(codim-i)*c[codim-i], i = 0 .. codim)-(1-t*lambda[f])*(sum((1+t*lambda[f])^j*c[codim-j]*t^(codim-j), j = 0 .. codim)); expand(%); alpha[f] := expand(%/(t*lambda[f])); setvariety(f[source_]); chern(f[source_][tangentbundle_]); z[f] := expand(alpha[f]*%); morphism(jj[f], Tn[f], B || f, [f[upperstardata_], E || f = -lambda[f]]); directimage := 0; for i to nops(z[f]) do directimage := directimage+t*lowerstar(jj[f], op(z[f])[i]) end do; z[f] := expand(directimage); setvariety(f[target_]); z[f] := z[f]+chern(f[target_][tangentbundle_]); series(z[f], t, dimension(f[target_])+1); totalchern[f] := convert(%, polynom); for i to dimension(f[target_]) do class[f][i] := coeff(totalchern[f], t, i) end do; classes[f] := [seq(class[f][i], i = 1 .. dimension(f[target_]))]; setvariety(B || f); sheaf(dimension(f[target_]), classes[f]); B || f[tangentbundle_] := %; `currentvariety_ is ` || currentvariety_ || `, tangent bundle computed,
> DIM is ` || DIM end proc;
> NULL;
> k := 3; n := 6;
                                      3
                                      6
> grass(k, n, c, all);
                       currentvariety_ is Gc, DIM is 9
> F := dual(determinant(Qc));
               1  2   2   1  3   3   1   4   4    1   5   5    1   6   6
    1 - c1 t + - t  c1  - - t  c1  + -- t  c1  - --- t  c1  + --- t  c1
               2          6          24          120          720

          1    7   7     1    8   8     1     9   9
       - ---- t  c1  + ----- t  c1  - ------ t  c1
         5040          40320          362880
> m := 6; G := expand(m*F);
                                      6
                     2   2    3   3   1  4   4   1   5   5    1   6   6
     6 - 6 c1 t + 3 t  c1  - t  c1  + - t  c1  - -- t  c1  + --- t  c1
                                      4          20          120

           1   7   7    1    8   8     1    9   9
        - --- t  c1  + ---- t  c1  - ----- t  c1
          840          6720          60480
> G := tensor(F, tensor(F, F));
                 9  2   2   9  3   3   27  4   4   81  5   5   81  6   6
    1 - 3 c1 t + - t  c1  - - t  c1  + -- t  c1  - -- t  c1  + -- t  c1
                 2          2          8           40          80

         243  7   7   729   8   8   243   9   9
       - --- t  c1  + ---- t  c1  - ---- t  c1
         560          4480          4480
> bundlesection(A, m*F);
                       currentvariety_ is A, DIM is 3
> tanblowup(iA);
              currentvariety_ is BiA, tangent bundle computed,

                DIM is 9
> H := dual(tangentbundle(BiA));
                        1 /     2                   2\  2   1 /     3
9 + (-6 c1 + 5 EiA) t + - \7 EiA  + 12 c1 EiA + 2 c1 / t  + - \-6 c1
                        2                                   6

                                              2        3         2    \  3
   + 18 c1 c2 - 18 c3 + 762 c2 EiA + 18 c1 EiA  + 5 EiA  - 363 c1  EiA/ t  +

  1  /    4                   2        4                             3
  -- |2 c1  + 24 c1 c3 - 24 c2  + 7 EiA  - 5976 c2 c1 EiA + 24 c1 EiA
  24 \

            2    2          3                  2   48552       \  4    1  /
   - 1160 c1  EiA  + 1856 c1  EiA + 2392 c2 EiA  + ----- c3 EiA| t  + --- \
                                                     5         /      120
                       5                      2           2           3
-73758 c3 c1 EiA - 6 c1  + 30 c3 c2 - 30 c3 c1  - 30 c1 c2  + 30 c2 c1

          5                 2              2                4          2    3
   + 5 EiA  - 9840 c2 c1 EiA  + 24310 c2 c1  EiA + 30 c1 EiA  - 3400 c1  EiA

            3    2          4                  3         2                   2
   + 2540 c1  EiA  - 5795 c1  EiA + 6920 c2 EiA  - 150 c2  EiA + 20496 c3 EiA

  \  5    1  /                      2           3        2   2        3
  / t  + --- |-360 c3 c1 c2 + 360 c3  + 60 c3 c1  - 60 c2  c1  - 30 c2
         720 \

         6         6                  3               2    2
   + 2 c1  - 23 EiA  - 17712 c2 c1 EiA  - 169800 c2 c1  EiA

                3              2          316368          2
   - 71232 c2 c1  EiA + 8748 c2  c1 EiA - ------ c3 c1 EiA
                                            5

     1422792      2                                  5          2    4
   + ------- c3 c1  EiA - 8964 c3 c2 EiA - 144 c1 EiA  - 8664 c1  EiA
        5

            3    3           4    2           5                   4
   + 2472 c1  EiA  + 41859 c1  EiA  + 13980 c1  EiA + 16608 c2 EiA

              2    2   247968       3\  6    1   /               2   2
   + 192972 c2  EiA  + ------ c3 EiA | t  + ---- |-2633400 c1 EiA  c2
                         5           /      5040 \

               3    2                2    3              2       2
   + 2098166 c1  EiA  c2 - 1015210 c1  EiA  c2 - 58898 c1  EiA c2

              4                      4              6              3
   + 170870 c1  EiA c2 - 18732 c1 EiA  c2 - 1176 EiA  c1 - 10843 c2  EiA

              5    2            4    3           6                5
   - 397348 c1  EiA  + 249732 c1  EiA  - 28707 c1  EiA + 39508 EiA  c2

              5   2           4   3           5        2   3             3
   - 18242 EiA  c1  - 3640 EiA  c1  + 42 c2 c1  - 84 c2  c1  + 1092 c1 c2

               2           2           4          7   873936
   - 1302 c1 c3  - 42 c3 c2  - 42 c3 c1  + 1050 c1  + ------ c3 c2 c1 EiA
                                                        5

            7   408408          3   10526334      2    2   3907064      3
   - 409 EiA  - ------ c3 c1 EiA  - -------- c3 c1  EiA  - ------- c3 c1  EiA
                  5                    5                      5

     21459144          2   1004892       4   502446   2                 2    3
   + -------- c3 c2 EiA  + ------- c3 EiA  - ------ c3  EiA + 1061564 c2  EiA
        5                     5                5

                  2\  7     1   /               4   2              4   4
   + 1386 c3 c2 c1 | t  + ----- |-4154944 c2 EiA  c1  + 1012848 EiA  c1
                   /      40320 \

                4   2             6    2             5    3   411616   7
   + 4238576 EiA  c2  + 2005368 c1  EiA  - 1586624 c1  EiA  + ------ c1  EiA
                                                                7

              6               5   3            6   2   22800    7
   + 87456 EiA  c2 - 34096 EiA  c1  - 32232 EiA  c1  - ----- EiA  c1
                                                         7

             3    2          3   2         2   4         4          2
   - 90304 c2  EiA  - 3352 c2  c1  - 112 c2  c1  - 504 c2  - 1344 c3  c2

            2   2            5   50674   8   13991    8
   + 5064 c3  c1  + 112 c3 c1  - ----- c1  - ----- EiA
                                   7           7

     293356704             2   5129184         2       3445152          4
   - --------- c3 c2 c1 EiA  - ------- c3 c2 c1  EiA - ------- c3 c1 EiA
         5                        5                       5

     49191552      2    3   120087664      3    2                4
   - -------- c3 c1  EiA  + --------- c3 c1  EiA  + 1755120 c3 c1  EiA
        5                       5

     99393504          3   522864      2       4103088   2
   + -------- c3 c2 EiA  - ------ c3 c2  EiA + ------- c3  c1 EiA
        5                    5                    5

     4074144       5   584539536   2    2                2                 3
   + ------- c3 EiA  + --------- c3  EiA  + 1008 c3 c1 c2  - 46672 c3 c2 c1
        5                 25

                    3   2              2    2   2              4    2
   - 11088096 c1 EiA  c2  + 18750048 c1  EiA  c2  - 12710928 c1  EiA  c2

               3    3                     3            3       2
   + 8668864 c1  EiA  c2 + 65520 c1 EiA c2  + 227280 c1  EiA c2

              5          525216    5      \  8     1    /            4   5
   - 357872 c1  EiA c2 - ------ EiA  c2 c1| t  + ------ |-4649256 EiA  c1
                           7              /      362880 \

             4                  3    3            6    3   51901938   7    2
   - 13455 c2  EiA + 48703740 c2  EiA  + 246345 c1  EiA  - -------- c1  EiA
                                                              7

                 5   2             7                 5   4            6   3
   + 13767372 EiA  c2  + 173178 EiA  c2 + 3400920 EiA  c1  - 93582 EiA  c1

     357507    7   2   55350    8      46981    9           2   3
   - ------ EiA  c1  - ----- EiA  c1 - ----- EiA  - 11988 c3  c1
       7                 7               7

               3           6   1004684364             3
   + 1404 c3 c2  - 54 c3 c1  - ---------- c3 c2 c1 EiA
                                   5

     2068910262         2    2                   3                   2
   + ---------- c3 c2 c1  EiA  + 3499938 c3 c2 c1  EiA + 370710 c3 c2  c1 EiA
         5

                      5                 2    4   396071064      3    3
   - 3194424 c3 c1 EiA  - 32929416 c3 c1  EiA  + --------- c3 c1  EiA
                                                     5

     705695922      4    2   17266662      5       355252068          4
   - --------- c3 c1  EiA  - -------- c3 c1  EiA + --------- c3 c2 EiA
         5                      5                      5

     21303324      2    2   8037900756   2       2   17582427   2   2
   - -------- c3 c2  EiA  - ---------- c3  c1 EiA  - -------- c3  c1  EiA
        5                       25                      5

     1509192   2          13087926       6   2228243472   2    3   180534   9
   + ------- c3  c2 EiA + -------- c3 EiA  + ---------- c3  EiA  + ------ c1
        5                    5                   25                  7

            3          2                   2   2                  4
   - 1638 c3  + 3402 c3  c1 c2 - 1944 c3 c2  c1  + 291870 c3 c2 c1

     854253   8            7            5   2           3   3             4
   - ------ c1  EiA + 54 c1  c2 - 162 c1  c2  + 10080 c1  c2  - 1404 c1 c2
       7

                 6                3   2                2   4
   + 679266 c2 c1  EiA - 168570 c2  c1  EiA - 656640 c2  c1  EiA

                2    3   2             4    3                    2   3
   - 12028644 c1  EiA  c2  - 4096638 c1  EiA  c2 + 3831192 c1 EiA  c2

                3    2   2              5    2      93120228    5      2
   - 94086702 c1  EiA  c2  + 54490302 c1  EiA  c2 - -------- EiA  c2 c1
                                                       7

     2338254    6                        4   2                  4   3\  9
   - ------- EiA  c2 c1 - 35824572 c1 EiA  c2  + 27110088 c2 EiA  c1 | t
        7                                                            /
> J := schurfunctor(conjugate([6]), H);
Warning,  computation interrupted
> chi(J);
                           -888949306755221155684
> conjugate([2]);
                                   [1, 1]
> seq(chi(wedge(j, H)), j = DIM-1 .. DIM);
                                  1555, -1
> seq(chi(schurfunctor2(conjugate([j]), H)), j = 0 .. 1);
                                  1, -1555
> op(chi);
proc()  ...  end;
> toddclass(BiA);
                         toddclass computed for BiA
> op(BiA);

> chern(4, J);
                                                                 4
65521428238233487690161291482329623100079355163920444883097760 c1

   + 124629811876917661547210104844475 c1 c3

                                       2
   + 1958468472351292663135312220865 c2

                                                          2
   - 39241047611096159189215899871955638688279990460 c2 c1

                                                                      4
   + 1677348562898765230018302933258624816194742656034453805660500 EiA

   + 31392838088875145999667304273979757988945662990 c2 c1 EiA

                                                                          3
   - 16773485628987733921562060412160601896574154634693123128054300 c1 EiA

                                                                      2    2
   + 62900571108704179575392928699383221154852180467765823970424260 c1  EiA

                                                                       3
   - 104834285181173713723819944097222934476414286774814674702162130 c1  EiA

                                                          2
   - 6278567617774660082039199104432238562176057600 c2 EiA

   - 49851924750773441360593806883650 c3 EiA
> m := 3; k := 3; n := 5; grass(k, n, c, all); F := dual(determinant(Qc)); bundlesection(A, m*F); tanblowup(iA); H := dual(tangentbundle(BiA)); seq(chi(wedge(j, H)), j = 0 .. DIM);
                                      3
                                      3
                                      5
                       currentvariety_ is Gc, DIM is 6
               1  2   2   1  3   3   1   4   4    1   5   5    1   6   6
    1 - c1 t + - t  c1  - - t  c1  + -- t  c1  - --- t  c1  + --- t  c1
               2          6          24          120          720
                       currentvariety_ is A, DIM is 3
              currentvariety_ is BiA, tangent bundle computed,

                DIM is 6
                        1 /             2               2\  2   1 /
6 + (-5 c1 + 2 EiA) t + - \6 c1 EiA + c1  + 2 c2 + 4 EiA / t  + - \210 c2 EiA
                        2                                       6

          3           2       3                            2    \  3   1  /  4
   + 2 EiA  + 9 c1 EiA  - 5 c1  + 15 c1 c2 - 15 c3 - 117 c1  EiA/ t  + -- \c1
                                                                       24

                      2          2                             3
   + 84 c1 c3 + 102 c2  + 4 c2 c1  - 1308 c2 c1 EiA + 64 c1 EiA

           2    2         3                 2              \  4    1  /
   - 234 c1  EiA  + 476 c1  EiA + 720 c2 EiA  + 2544 c3 EiA/ t  + --- \
                                                                  120
                       5                         2            2           3
-16015 c3 c1 EiA - 5 c1  + 1645 c3 c2 - 425 c3 c1  - 825 c1 c2  + 25 c2 c1

          5                 2             2                 4         2    3
   - 8 EiA  - 3235 c2 c1 EiA  + 4310 c2 c1  EiA + 130 c1 EiA  - 715 c1  EiA

           3    2          4                  3         2                  2\
   + 410 c1  EiA  - 1215 c1  EiA + 1790 c2 EiA  + 600 c2  EiA + 8790 c3 EiA /

   5    1  /                        2             3          2   2          4
  t  + --- \-9792 c3 c1 c2 + 5529 c3  + 1302 c3 c1  + 2427 c2  c1  + 6 c2 c1
       720

           3     6         6                 3             2    2
   + 314 c2  + c1  - 20 EiA  - 5784 c2 c1 EiA  - 8802 c2 c1  EiA

                3              2                         2              2
   - 10578 c2 c1  EiA - 1452 c2  c1 EiA - 41478 c3 c1 EiA  + 50916 c3 c1  EiA

                                5          2    4         3    3
   + 5724 c3 c2 EiA + 276 c1 EiA  - 1041 c1  EiA  + 906 c1  EiA

            4    2          5                  4           2    2
   + 4239 c1  EiA  + 2466 c1  EiA + 3456 c2 EiA  + 16380 c2  EiA

                 3\  6
   + 21156 c3 EiA / t
                        1, 53, -221, 336, -221, 53, 1
> m := 4; k := 2; n := 5; grass(k, n, c, all); F := dual(determinant(Qc)); bundlesection(A, m*F); tanblowup(iA); H := dual(tangentbundle(BiA)); seq(chi(wedge(j, H)), j = 0 .. DIM);
                                      4
                                      2
                                      5
                       currentvariety_ is Gc, DIM is 6
               1  2   2   1  3   3   1   4   4    1   5   5    1   6   6
    1 - c1 t + - t  c1  - - t  c1  + -- t  c1  - --- t  c1  + --- t  c1
               2          6          24          120          720
                       currentvariety_ is A, DIM is 2
              currentvariety_ is BiA, tangent bundle computed,

                DIM is 6
                        1 /    2               2           \  2
6 + (-5 c1 + 3 EiA) t + - \3 c1  - 2 c2 + 5 EiA  + 8 c1 EiA/ t
                        2

     1 /     3            2        3   855                           2    \  3
   + - |-5 c1  + 12 c1 EiA  + 3 EiA  + --- c2 EiA + 15 c1 c2 - 159 c1  EiA| t  +
     6 \                                2                                 /

  1  /    4        2           2        4                             3
  -- \3 c1  - 10 c2  - 40 c2 c1  - 7 EiA  - 2866 c2 c1 EiA - 32 c1 EiA
  24

           2    2         3                 2\  4    1  /     5            2
   - 356 c1  EiA  + 704 c1  EiA + 770 c2 EiA / t  + --- |-5 c1  - 775 c1 c2
                                                    120 \

              3         5   3965          2              2                 4
   + 250 c2 c1  - 37 EiA  - ---- c2 c1 EiA  + 10400 c2 c1  EiA - 200 c1 EiA
                             2

           2    3         3    2          4       25       3   1515   2    \
   - 435 c1  EiA  + 620 c1  EiA  - 2025 c1  EiA - -- c2 EiA  - ---- c2  EiA|
                                                  2             2          /

   5    1  /       2   2            4          3       6          6
  t  + --- |4305 c2  c1  - 762 c2 c1  - 3626 c2  + 3 c1  - 103 EiA
       720 \

                   3              2    2              3              2
   + 8163 c2 c1 EiA  - 45471 c2 c1  EiA  - 28176 c2 c1  EiA + 5802 c2  c1 EiA

               5         2    4          3    3          4    2          5
   - 672 c1 EiA  - 198 c1  EiA  - 1024 c1  EiA  + 7941 c1  EiA  + 4632 c1  EiA

                4   238551   2    2\  6
   - 4065 c2 EiA  + ------ c2  EiA | t
                      4            /
                       1, -52, 158, -209, 158, -52, 1
> m := 4; k := 3; n := 6; grass(k, n, c, all); F := dual(determinant(Qc)); bundlesection(A, m*F); tanblowup(iA); H := dual(tangentbundle(BiA)); seq(chi(wedge(j, H)), j = 0 .. DIM);
                                      4
                                      3
                                      6
                       currentvariety_ is Gc, DIM is 9
               1  2   2   1  3   3   1   4   4    1   5   5    1   6   6
    1 - c1 t + - t  c1  - - t  c1  + -- t  c1  - --- t  c1  + --- t  c1
               2          6          24          120          720

          1    7   7     1    8   8     1     9   9
       - ---- t  c1  + ----- t  c1  - ------ t  c1
         5040          40320          362880
                       currentvariety_ is A, DIM is 5
              currentvariety_ is BiA, tangent bundle computed,

                DIM is 9
                        1 /                2       2\  2
9 + (-6 c1 + 3 EiA) t + - \8 c1 EiA + 5 EiA  + 2 c1 / t
                        2

     1 /     3                           2            3            2\  3   1
   + - \-6 c1  + 18 c1 c2 - 18 c3 + 12 c1  EiA + 3 EiA  + 12 c1 EiA / t  + --
     6                                                                     24

  /      4                   2        4                             3
  \-10 c1  + 24 c1 c3 - 24 c2  - 7 EiA  + 1890 c2 c1 EiA - 32 c1 EiA

          2    2         3                  \  4    1  /
   - 48 c1  EiA  - 752 c1  EiA - 1890 c3 EiA/ t  + --- |22059 c3 c1 EiA
                                                   120 \

          5                        2            2           3         5
   + 84 c1  + 930 c3 c2 - 180 c3 c1  - 780 c1 c2  + 30 c2 c1  - 37 EiA

     1575          2              2                 4         2    3
   + ---- c2 c1 EiA  - 13515 c2 c1  EiA - 200 c1 EiA  - 440 c1  EiA
      2

           3    2          4            2       1575       2\  5    1  /
   - 780 c1  EiA  + 3895 c1  EiA - 90 c2  EiA - ---- c3 EiA | t  + --- |
                                                 2          /      720 \
                         2             3          2   2          3         6
-12168 c3 c1 c2 + 5688 c3  + 1140 c3 c1  + 5340 c2  c1  - 4470 c2  - 340 c1

            6                  3             2    2              3
   - 103 EiA  - 11718 c2 c1 EiA  + 7326 c2 c1  EiA  + 52065 c2 c1  EiA

            2          278838          2   566649      2
   - 2952 c2  c1 EiA - ------ c3 c1 EiA  - ------ c3 c1  EiA - 10020 c3 c2 EiA
                         5                   5

               5          2    4          3    3         4    2
   - 672 c1 EiA  - 1824 c1  EiA  + 1808 c1  EiA  + 858 c1  EiA

             5             2    2               3\  6    1   /
   - 12012 c1  EiA - 252 c2  EiA  + 11718 c3 EiA | t  + ---- |
                                                 /      5040 \
            2   2   95277   3    2               2    3              2       2
-8316 c1 EiA  c2  - ----- c1  EiA  c2 + 102606 c1  EiA  c2 + 33621 c1  EiA c2
                      2

              4                      4              6              3
   - 147714 c1  EiA c2 - 49833 c1 EiA  c2 - 1680 EiA  c1 - 15645 c2  EiA

           5    2           4    3           6               5   2
   - 322 c1  EiA  - 15113 c1  EiA  + 28595 c1  EiA - 5488 EiA  c1

             4   3           5           2   3              3              2
   + 9016 EiA  c1  - 84 c2 c1  - 20034 c2  c1  + 32172 c1 c2  - 38598 c1 c3

             2             4         7                               7
   - 42 c3 c2  - 3906 c3 c1  + 918 c1  + 17262 c3 c2 c1 EiA - 221 EiA

     1480794          3   4350381      2    2   1879542      3
   - ------- c3 c1 EiA  + ------- c3 c1  EiA  + ------- c3 c1  EiA
        5                   10                     5

                     2                4           2       1070790   2    3
   + 326774 c3 c2 EiA  + 128205 c3 EiA  + 39375 c3  EiA + ------- c2  EiA
                                                            11

                    2\  7     1   /             4   2            4   4
   + 140994 c3 c2 c1 | t  + ----- |257568 c2 EiA  c1  - 28744 EiA  c1
                     /      40320 \

     3681840    4   2            6    2           5    3           7
   + ------- EiA  c2  + 142952 c1  EiA  + 64528 c1  EiA  - 53440 c1  EiA
       11

              5   3            6   2           7              3    2
   + 26496 EiA  c1  - 12672 EiA  c1  - 3328 EiA  c1 - 41720 c2  EiA

              3   2            6           2   4         4          2
   - 123032 c2  c1  + 864 c2 c1  + 52832 c2  c1  - 504 c2  - 9984 c3  c2

              2   2             5          8          8
   + 141480 c3  c1  + 9664 c3 c1  - 2002 c1  - 383 EiA

                         2   2251432         2       6783744          4
   - 4138450 c3 c2 c1 EiA  + ------- c3 c2 c1  EiA - ------- c3 c1 EiA
                                5                       5

     7774248      2    3   4643696      3    2   4740968      4
   + ------- c3 c1  EiA  - ------- c3 c1  EiA  - ------- c3 c1  EiA
        5                     5                     5

                      3              2       1680792   2
   + 1144112 c3 c2 EiA  + 22536 c3 c2  EiA - ------- c3  c1 EiA
                                                5

                  5            2    2                 2                  3
   + 410832 c3 EiA  + 952665 c3  EiA  + 16848 c3 c1 c2  - 742416 c3 c2 c1

     7579104       3   2            2    2   2            4    2
   - ------- c1 EiA  c2  + 929065 c1  EiA  c2  - 568176 c1  EiA  c2
       11

     3070000   3    3                     3            3       2
   - ------- c1  EiA  c2 + 63784 c1 EiA c2  - 166496 c1  EiA c2
        7

              5                    5      \  8     1    /          4   5
   + 335448 c1  EiA c2 - 147876 EiA  c2 c1| t  + ------ |201870 EiA  c1
                                          /      362880 \

            4               3    3           6    3             7    2
   - 1215 c2  EiA - 87165 c2  EiA  + 64803 c1  EiA  - 1323072 c1  EiA

     9681066    5   2           5   4            6   3            7   2
   + ------- EiA  c2  - 6219 EiA  c1  + 66288 EiA  c1  - 21312 EiA  c1
       11

             8             9            2   3              3              6
   - 4896 EiA  c1 - 501 EiA  - 374976 c3  c1  + 31644 c3 c2  - 17982 c3 c1

     14693643             3   194252958         2    2
   - -------- c3 c2 c1 EiA  + --------- c3 c2 c1  EiA
        2                         5

     15461451         3       5762259      2          15572196          5
   - -------- c3 c2 c1  EiA + ------- c3 c2  c1 EiA - -------- c3 c1 EiA
        5                        5                       5

     26609229      2    4                3    3   64746801      4    2
   + -------- c3 c1  EiA  - 4118652 c3 c1  EiA  - -------- c3 c1  EiA
        5                                            10

     9641097      5                        4              2    2
   + ------- c3 c1  EiA + 3051162 c3 c2 EiA  + 84672 c3 c2  EiA
        5

     184789701   2       2   8593722   2   2                2
   - --------- c3  c1 EiA  + ------- c3  c1  EiA - 704349 c3  c2 EiA
        10                      5

                   6   2692467   2    3          9           3
   + 1021734 c3 EiA  + ------- c3  EiA  + 3450 c1  - 49590 c3
                          4

              2                     2   2                   4           8
   + 167346 c3  c1 c2 - 163296 c3 c2  c1  + 2420712 c3 c2 c1  + 53658 c1  EiA

            7               5   2            3   3              4
   - 3186 c1  c2 - 111672 c1  c2  + 351720 c1  c2  - 33264 c1 c2

                 6                3   2                2   4
   - 565731 c2 c1  EiA - 794025 c2  c1  EiA + 821682 c2  c1  EiA

     132693129   2    3   2   3311442   4    3                  2   3
   + --------- c1  EiA  c2  + ------- c1  EiA  c2 + 97218 c1 EiA  c2
        44                       7

     22913055   3    2   2   15140385   5    2                5      2
   - -------- c1  EiA  c2  + -------- c1  EiA  c2 + 489591 EiA  c2 c1
        2                       2

     671571    6         19685268       4   2   6950988       4   3\  9
   - ------ EiA  c2 c1 - -------- c1 EiA  c2  - ------- c2 EiA  c1 | t
       2                    11                     7               /
        1, -1168, 7928, -22533, 35970, -35970, 22533, -7928, 1168, -1
> m := 4; k := 4; n := 6; grass(k, n, c, all); F := dual(determinant(Qc)); bundlesection(A, m*F); tanblowup(iA); H := dual(tangentbundle(BiA)); seq(chi(wedge(j, H)), j = 0 .. DIM);
                                      4
                                      4
                                      6
                       currentvariety_ is Gc, DIM is 8
               1  2   2   1  3   3   1   4   4    1   5   5    1   6   6
    1 - c1 t + - t  c1  - - t  c1  + -- t  c1  - --- t  c1  + --- t  c1
               2          6          24          120          720

          1    7   7     1    8   8
       - ---- t  c1  + ----- t  c1
         5040          40320
                       currentvariety_ is A, DIM is 4
              currentvariety_ is BiA, tangent bundle computed,

                DIM is 8
                        1 /                       2\  2
8 + (-6 c1 + 3 EiA) t + - \8 c1 EiA + 4 c2 + 5 EiA / t
                        2

     1 /     3                           2            3            2\  3   1
   + - \-6 c1  + 18 c1 c2 - 18 c3 + 12 c1  EiA + 3 EiA  + 12 c1 EiA / t  + --
     6                                                                     24

  /                       2          2        4                            3
  \8 c4 + 16 c1 c3 - 28 c2  + 8 c2 c1  + 5 EiA  + 248 c2 c1 EiA + 16 c1 EiA

          2    2         3                  \  4    1  /
   + 24 c1  EiA  - 840 c1  EiA + 2438 c3 EiA/ t  + --- |-18885 c3 c1 EiA
                                                   120 \

         5                       2           2           3        5
   - 6 c1  - 720 c3 c2 - 30 c3 c1  - 30 c1 c2  + 30 c2 c1  + 8 EiA

                  2            2                 4         2    3
   + 670 c2 c1 EiA  - 350 c2 c1  EiA - 110 c1 EiA  - 530 c1  EiA

            3    2          4             2       1955       2
   - 1550 c1  EiA  + 4190 c1  EiA - 135 c2  EiA + ---- c3 EiA  + 180 c1 c4
                                                   2

                 \  5    1  /                      2           3        2   2
   + 12790 c4 EiA| t  + --- \5244 c3 c1 c2 + 174 c3  + 48 c3 c1  - 78 c2  c1
                 /      720

             4         3         6                 3            2    2
   + 12 c2 c1  + 124 c2  + 23 EiA  + 1452 c2 c1 EiA  + 180 c2 c1  EiA

              3             2                        2              2
   - 852 c2 c1  EiA + 204 c2  c1 EiA + 6102 c3 c1 EiA  + 67614 c3 c1  EiA

                                5          2    4          3    3
   + 4938 c3 c2 EiA - 420 c1 EiA  - 1752 c1  EiA  + 1172 c1  EiA

            4    2           5                  4         2    2
   + 3276 c1  EiA  - 11820 c1  EiA + 1560 c2 EiA  - 378 c2  EiA

                3                       2              2                  \  6
   - 1617 c3 EiA  - 252 c4 c2 - 72 c4 c1  - 1092 c4 EiA  - 91176 c4 c1 EiA/ t  +

   1   /          2   2         3    2            2    3            2       2
  ---- |-14 c1 EiA  c2  - 798 c1  EiA  c2 + 854 c1  EiA  c2 + 133 c1  EiA c2
  5040 \

           4                     4              6            3
   + 406 c1  EiA c2 - 8456 c1 EiA  c2 - 1134 EiA  c1 + 378 c2  EiA

            5    2           4    3           6               5
   - 5530 c1  EiA  - 13286 c1  EiA  + 25242 c1  EiA + 5460 EiA  c2

             5   2           4   3           5        2   3           3
   - 5320 EiA  c1  + 7280 EiA  c1  + 42 c2 c1  - 84 c2  c1  + 42 c1 c2

                          3           2             2           4
   + 42 c4 c3 - 2730 c4 c1  - 42 c1 c3  - 2142 c3 c2  - 42 c3 c1

                    2       7                              7
   + 81676 c4 c1 EiA  - 6 c1  - 14637 c3 c2 c1 EiA + 52 EiA

                    3              2    2               3
   + 30779 c3 c1 EiA  - 58863 c3 c1  EiA  - 157983 c3 c1  EiA

                   2   34111       4           2             2    3
   - 3731 c3 c2 EiA  - ----- c3 EiA  - 25368 c3  EiA - 819 c2  EiA
                         2

                                 3               2
   + 336 c4 c1 c2 - 199598 c4 EiA  + 322364 c4 c1  EiA + 34846 c4 c2 EiA

                   2\  7     1   /            4   2            4   4
   - 35924 c3 c2 c1 | t  + ----- |39360 c2 EiA  c1  - 43760 EiA  c1
                    /      40320 \

             4   2            6    2           5    3           7
   + 2576 EiA  c2  + 188960 c1  EiA  + 28816 c1  EiA  - 45968 c1  EiA

              6               5   3            6   2           7
   + 14560 EiA  c2 + 18704 EiA  c1  - 12432 EiA  c1  - 2320 EiA  c1

            3    2         3   2           6         2   4         4
   + 1008 c2  EiA  + 480 c2  c1  + 16 c2 c1  - 152 c2  c1  - 508 c2

            2            2   2           5          2             2
   + 5872 c3  c2 + 480 c3  c1  + 96 c3 c1  - 1128 c4  + 1808 c4 c2

                4          8                      2                 2
   + 13584 c4 c1  + 109 EiA  + 355152 c3 c2 c1 EiA  - 84424 c3 c2 c1  EiA

                     4               2    3               3    2
   + 129512 c3 c1 EiA  - 114032 c3 c1  EiA  - 840576 c3 c1  EiA

                 4                      3              2
   + 284296 c3 c1  EiA - 18040 c3 c2 EiA  - 31148 c3 c2  EiA

              2                      5             2    2
   + 201584 c3  c1 EiA - 42664 c3 EiA  + 1477193 c3  EiA  - 960 c4 c1 c3

                 2                                      4                    3
   - 336 c4 c2 c1  - 121696 c4 c2 c1 EiA - 722304 c4 EiA  + 1449760 c4 c1 EiA

                 2    2               3                     2
   - 432192 c4 c1  EiA  - 787760 c4 c1  EiA - 5488 c4 c2 EiA

                                     2                  3              3   2
   - 112600 c4 c3 EiA + 8160 c3 c1 c2  + 162592 c3 c2 c1  - 1168 c1 EiA  c2

             2    2   2            4    2             3    3
   + 16152 c1  EiA  c2  - 111664 c1  EiA  c2 + 7408 c1  EiA  c2

                   3          3       2           5          88648    5      \
   - 2336 c1 EiA c2  + 7920 c1  EiA c2  + 10640 c1  EiA c2 - ----- EiA  c2 c1|
                                                               3             /

   8
  t
              1, -376, 2006, -4522, 5777, -4522, 2006, -376, 1
> m := 5; k := 2; n := 6; grass(k, n, c, all); F := dual(determinant(Qc)); bundlesection(A, m*F); tanblowup(iA); H := dual(tangentbundle(BiA)); seq(chi(wedge(j, H)), j = 0 .. DIM);
                                      5
                                      2
                                      6
                       currentvariety_ is Gc, DIM is 8
               1  2   2   1  3   3   1   4   4    1   5   5    1   6   6
    1 - c1 t + - t  c1  - - t  c1  + -- t  c1  - --- t  c1  + --- t  c1
               2          6          24          120          720

          1    7   7     1    8   8
       - ---- t  c1  + ----- t  c1
         5040          40320
                       currentvariety_ is A, DIM is 3
              currentvariety_ is BiA, tangent bundle computed,

                DIM is 8
                       1 /                 2       2       \  2
8 + (4 EiA - 6 c1) t + - \10 c1 EiA + 6 EiA  + 4 c1  - 4 c2/ t
                       2

     1 /4158                   2       3                   3         2    \  3
   + - |---- c2 EiA + 15 c1 EiA  - 6 c1  + 18 c1 c2 + 4 EiA  - 282 c1  EiA| t  +
     6 \ 5                                                                /

  1  /    4        2          2        4                             3
  -- |4 c1  - 20 c2  - 8 c2 c1  + 6 EiA  - 4200 c2 c1 EiA + 20 c1 EiA
  24 \

           2    2          3       11256       2\  4    1  /     5           2
   - 774 c1  EiA  + 1520 c1  EiA + ----- c2 EiA | t  + --- \-6 c1  - 90 c1 c2
                                     5          /      120

              3         5                 2              2                 4
   + 110 c2 c1  + 24 EiA  - 4900 c2 c1 EiA  + 10154 c2 c1  EiA + 125 c1 EiA

            2    3          3    2          4                  3
   - 1945 c1  EiA  + 2000 c1  EiA  - 5150 c1  EiA + 6146 c2 EiA

            2    \  5    1  /       2   2            4         3       6
   - 2912 c2  EiA/ t  + --- |2766 c2  c1  - 588 c2 c1  - 940 c2  + 4 c1
                        720 \

           6                3   771348      2    2              3
   + 96 EiA  + 756 c2 c1 EiA  - ------ c2 c1  EiA  - 15372 c2 c1  EiA
                                  5

     111126   2                    5          2    4          3    3
   + ------ c2  c1 EiA + 600 c1 EiA  - 5091 c1  EiA  + 1930 c1  EiA
       5

             4    2           5       93324       4   5582598   2    2\  6
   + 23538 c1  EiA  + 13566 c1  EiA + ----- c2 EiA  + ------- c2  EiA | t  +
                                        5               25            /

   1   /  10008747       2   2             3    2      4008676   2    3
  ---- |- -------- c1 EiA  c2  + 1427594 c1  EiA  c2 - ------- c1  EiA  c2
  5040 \     5                                            5

             2       2           4                      4              6
   - 59724 c1  EiA c2  + 12824 c1  EiA c2 + 27440 c1 EiA  c2 + 2100 EiA  c1

     50666   3                5    2            4    3           6
   - ----- c2  EiA - 248997 c1  EiA  + 128338 c1  EiA  - 30380 c1  EiA
       5

     268912    5               5   2           4   3             5
   + ------ EiA  c2 - 12558 EiA  c1  + 1925 EiA  c1  + 2282 c2 c1
       5

             2   3              3       7          7   26841808   2    3\  7
   - 18620 c2  c1  + 31794 c1 c2  - 6 c1  + 284 EiA  + -------- c2  EiA | t  +
                                                          25            /

    1   /  14610608       4   2             4   4   98512512    4   2
  ----- |- -------- c2 EiA  c1  + 442962 EiA  c1  + -------- EiA  c2
  40320 \     5                                        25

               6    2            5    3           7       719712    6
   + 1360840 c1  EiA  - 900512 c1  EiA  + 60672 c1  EiA + ------ EiA  c2
                                                            5

              5   3            6   2           7      32964736   3    2
   - 14000 EiA  c1  - 29908 EiA  c1  + 5800 EiA  c1 - -------- c2  EiA
                                                         25

              3   2             6           2   4          4       8
   - 201696 c2  c1  - 6416 c2 c1  + 70392 c2  c1  + 2316 c2  + 4 c1

            8   34247072       3   2   233436672   2    2   2
   + 686 EiA  - -------- c1 EiA  c2  + --------- c1  EiA  c2
                   5                      25

     34868128   4    2                3    3                     3
   - -------- c1  EiA  c2 + 5158096 c1  EiA  c2 + 97088 c1 EiA c2
        5

     607608   3       2   26368   5          728352    5      \  8
   + ------ c1  EiA c2  - ----- c1  EiA c2 + ------ EiA  c2 c1| t
       5                    5                  5              /
               1, 395, -1655, 2916, -3311, 2916, -1655, 395, 1
> m := 6; k := 3; n := 6; grass(k, n, c, all); F := dual(determinant(Qc)); bundlesection(A, m*F); tanblowup(iA); H := dual(tangentbundle(BiA)); seq(chi(wedge(j, H)), j = 1 .. DIM);
                                      6
                                      3
                                      6
                       currentvariety_ is Gc, DIM is 9
               1  2   2   1  3   3   1   4   4    1   5   5    1   6   6
    1 - c1 t + - t  c1  - - t  c1  + -- t  c1  - --- t  c1  + --- t  c1
               2          6          24          120          720

          1    7   7     1    8   8     1     9   9
       - ---- t  c1  + ----- t  c1  - ------ t  c1
         5040          40320          362880
                       currentvariety_ is A, DIM is 3
              currentvariety_ is BiA, tangent bundle computed,

                DIM is 9
                        1 /     2                   2\  2   1 /     3
9 + (-6 c1 + 5 EiA) t + - \7 EiA  + 12 c1 EiA + 2 c1 / t  + - \-6 c1
                        2                                   6

                                              2        3         2    \  3
   + 18 c1 c2 - 18 c3 + 762 c2 EiA + 18 c1 EiA  + 5 EiA  - 363 c1  EiA/ t  +

  1  /    4                   2        4                             3
  -- |2 c1  + 24 c1 c3 - 24 c2  + 7 EiA  - 5976 c2 c1 EiA + 24 c1 EiA
  24 \

            2    2          3                  2   48552       \  4    1  /
   - 1160 c1  EiA  + 1856 c1  EiA + 2392 c2 EiA  + ----- c3 EiA| t  + --- \
                                                     5         /      120
                       5                      2           2           3
-73758 c3 c1 EiA - 6 c1  + 30 c3 c2 - 30 c3 c1  - 30 c1 c2  + 30 c2 c1

          5                 2              2                4          2    3
   + 5 EiA  - 9840 c2 c1 EiA  + 24310 c2 c1  EiA + 30 c1 EiA  - 3400 c1  EiA

            3    2          4                  3         2                   2
   + 2540 c1  EiA  - 5795 c1  EiA + 6920 c2 EiA  - 150 c2  EiA + 20496 c3 EiA

  \  5    1  /                      2           3        2   2        3
  / t  + --- |-360 c3 c1 c2 + 360 c3  + 60 c3 c1  - 60 c2  c1  - 30 c2
         720 \

         6         6                  3               2    2
   + 2 c1  - 23 EiA  - 17712 c2 c1 EiA  - 169800 c2 c1  EiA

                3              2          316368          2
   - 71232 c2 c1  EiA + 8748 c2  c1 EiA - ------ c3 c1 EiA
                                            5

     1422792      2                                  5          2    4
   + ------- c3 c1  EiA - 8964 c3 c2 EiA - 144 c1 EiA  - 8664 c1  EiA
        5

            3    3           4    2           5                   4
   + 2472 c1  EiA  + 41859 c1  EiA  + 13980 c1  EiA + 16608 c2 EiA

              2    2   247968       3\  6    1   /               2   2
   + 192972 c2  EiA  + ------ c3 EiA | t  + ---- |-2633400 c1 EiA  c2
                         5           /      5040 \

               3    2                2    3              2       2
   + 2098166 c1  EiA  c2 - 1015210 c1  EiA  c2 - 58898 c1  EiA c2

              4                      4              6              3
   + 170870 c1  EiA c2 - 18732 c1 EiA  c2 - 1176 EiA  c1 - 10843 c2  EiA

              5    2            4    3           6                5
   - 397348 c1  EiA  + 249732 c1  EiA  - 28707 c1  EiA + 39508 EiA  c2

              5   2           4   3           5        2   3             3
   - 18242 EiA  c1  - 3640 EiA  c1  + 42 c2 c1  - 84 c2  c1  + 1092 c1 c2

               2           2           4          7   873936
   - 1302 c1 c3  - 42 c3 c2  - 42 c3 c1  + 1050 c1  + ------ c3 c2 c1 EiA
                                                        5

            7   408408          3   10526334      2    2   3907064      3
   - 409 EiA  - ------ c3 c1 EiA  - -------- c3 c1  EiA  - ------- c3 c1  EiA
                  5                    5                      5

     21459144          2   1004892       4   502446   2                 2    3
   + -------- c3 c2 EiA  + ------- c3 EiA  - ------ c3  EiA + 1061564 c2  EiA
        5                     5                5

                  2\  7     1   /               4   2              4   4
   + 1386 c3 c2 c1 | t  + ----- |-4154944 c2 EiA  c1  + 1012848 EiA  c1
                   /      40320 \

                4   2             6    2             5    3   411616   7
   + 4238576 EiA  c2  + 2005368 c1  EiA  - 1586624 c1  EiA  + ------ c1  EiA
                                                                7

              6               5   3            6   2   22800    7
   + 87456 EiA  c2 - 34096 EiA  c1  - 32232 EiA  c1  - ----- EiA  c1
                                                         7

             3    2          3   2         2   4         4          2
   - 90304 c2  EiA  - 3352 c2  c1  - 112 c2  c1  - 504 c2  - 1344 c3  c2

            2   2            5   50674   8   13991    8
   + 5064 c3  c1  + 112 c3 c1  - ----- c1  - ----- EiA
                                   7           7

     293356704             2   5129184         2       3445152          4
   - --------- c3 c2 c1 EiA  - ------- c3 c2 c1  EiA - ------- c3 c1 EiA
         5                        5                       5

     49191552      2    3   120087664      3    2                4
   - -------- c3 c1  EiA  + --------- c3 c1  EiA  + 1755120 c3 c1  EiA
        5                       5

     99393504          3   522864      2       4103088   2
   + -------- c3 c2 EiA  - ------ c3 c2  EiA + ------- c3  c1 EiA
        5                    5                    5

     4074144       5   584539536   2    2                2                 3
   + ------- c3 EiA  + --------- c3  EiA  + 1008 c3 c1 c2  - 46672 c3 c2 c1
        5                 25

                    3   2              2    2   2              4    2
   - 11088096 c1 EiA  c2  + 18750048 c1  EiA  c2  - 12710928 c1  EiA  c2

               3    3                     3            3       2
   + 8668864 c1  EiA  c2 + 65520 c1 EiA c2  + 227280 c1  EiA c2

              5          525216    5      \  8     1    /            4   5
   - 357872 c1  EiA c2 - ------ EiA  c2 c1| t  + ------ |-4649256 EiA  c1
                           7              /      362880 \

             4                  3    3            6    3   51901938   7    2
   - 13455 c2  EiA + 48703740 c2  EiA  + 246345 c1  EiA  - -------- c1  EiA
                                                              7

                 5   2             7                 5   4            6   3
   + 13767372 EiA  c2  + 173178 EiA  c2 + 3400920 EiA  c1  - 93582 EiA  c1

     357507    7   2   55350    8      46981    9           2   3
   - ------ EiA  c1  - ----- EiA  c1 - ----- EiA  - 11988 c3  c1
       7                 7               7

               3           6   1004684364             3
   + 1404 c3 c2  - 54 c3 c1  - ---------- c3 c2 c1 EiA
                                   5

     2068910262         2    2                   3                   2
   + ---------- c3 c2 c1  EiA  + 3499938 c3 c2 c1  EiA + 370710 c3 c2  c1 EiA
         5

                      5                 2    4   396071064      3    3
   - 3194424 c3 c1 EiA  - 32929416 c3 c1  EiA  + --------- c3 c1  EiA
                                                     5

     705695922      4    2   17266662      5       355252068          4
   - --------- c3 c1  EiA  - -------- c3 c1  EiA + --------- c3 c2 EiA
         5                      5                      5

     21303324      2    2   8037900756   2       2   17582427   2   2
   - -------- c3 c2  EiA  - ---------- c3  c1 EiA  - -------- c3  c1  EiA
        5                       25                      5

     1509192   2          13087926       6   2228243472   2    3   180534   9
   + ------- c3  c2 EiA + -------- c3 EiA  + ---------- c3  EiA  + ------ c1
        5                    5                   25                  7

            3          2                   2   2                  4
   - 1638 c3  + 3402 c3  c1 c2 - 1944 c3 c2  c1  + 291870 c3 c2 c1

     854253   8            7            5   2           3   3             4
   - ------ c1  EiA + 54 c1  c2 - 162 c1  c2  + 10080 c1  c2  - 1404 c1 c2
       7

                 6                3   2                2   4
   + 679266 c2 c1  EiA - 168570 c2  c1  EiA - 656640 c2  c1  EiA

                2    3   2             4    3                    2   3
   - 12028644 c1  EiA  c2  - 4096638 c1  EiA  c2 + 3831192 c1 EiA  c2

                3    2   2              5    2      93120228    5      2
   - 94086702 c1  EiA  c2  + 54490302 c1  EiA  c2 - -------- EiA  c2 c1
                                                       7

     2338254    6                        4   2                  4   3\  9
   - ------- EiA  c2 c1 - 35824572 c1 EiA  c2  + 27110088 c2 EiA  c1 | t
        7                                                            /
         -1555, 6590, -11625, 13179, -13179, 11625, -6590, 1555, -1
> H; op(wedge);
                       1 /                 2       2       \  2
8 + (4 EiA - 6 c1) t + - \10 c1 EiA + 6 EiA  + 4 c1  - 4 c2/ t
                       2

     1 /4158                   2       3                   3         2    \  3
   + - |---- c2 EiA + 15 c1 EiA  - 6 c1  + 18 c1 c2 + 4 EiA  - 282 c1  EiA| t  +
     6 \ 5                                                                /

  1  /    4        2          2        4                             3
  -- |4 c1  - 20 c2  - 8 c2 c1  + 6 EiA  - 4200 c2 c1 EiA + 20 c1 EiA
  24 \

           2    2          3       11256       2\  4    1  /     5           2
   - 774 c1  EiA  + 1520 c1  EiA + ----- c2 EiA | t  + --- \-6 c1  - 90 c1 c2
                                     5          /      120

              3         5                 2              2                 4
   + 110 c2 c1  + 24 EiA  - 4900 c2 c1 EiA  + 10154 c2 c1  EiA + 125 c1 EiA

            2    3          3    2          4                  3
   - 1945 c1  EiA  + 2000 c1  EiA  - 5150 c1  EiA + 6146 c2 EiA

            2    \  5    1  /       2   2            4         3       6
   - 2912 c2  EiA/ t  + --- |2766 c2  c1  - 588 c2 c1  - 940 c2  + 4 c1
                        720 \

           6                3   771348      2    2              3
   + 96 EiA  + 756 c2 c1 EiA  - ------ c2 c1  EiA  - 15372 c2 c1  EiA
                                  5

     111126   2                    5          2    4          3    3
   + ------ c2  c1 EiA + 600 c1 EiA  - 5091 c1  EiA  + 1930 c1  EiA
       5

             4    2           5       93324       4   5582598   2    2\  6
   + 23538 c1  EiA  + 13566 c1  EiA + ----- c2 EiA  + ------- c2  EiA | t  +
                                        5               25            /

   1   /  10008747       2   2             3    2      4008676   2    3
  ---- |- -------- c1 EiA  c2  + 1427594 c1  EiA  c2 - ------- c1  EiA  c2
  5040 \     5                                            5

             2       2           4                      4              6
   - 59724 c1  EiA c2  + 12824 c1  EiA c2 + 27440 c1 EiA  c2 + 2100 EiA  c1

     50666   3                5    2            4    3           6
   - ----- c2  EiA - 248997 c1  EiA  + 128338 c1  EiA  - 30380 c1  EiA
       5

     268912    5               5   2           4   3             5
   + ------ EiA  c2 - 12558 EiA  c1  + 1925 EiA  c1  + 2282 c2 c1
       5

             2   3              3       7          7   26841808   2    3\  7
   - 18620 c2  c1  + 31794 c1 c2  - 6 c1  + 284 EiA  + -------- c2  EiA | t  +
                                                          25            /

    1   /  14610608       4   2             4   4   98512512    4   2
  ----- |- -------- c2 EiA  c1  + 442962 EiA  c1  + -------- EiA  c2
  40320 \     5                                        25

               6    2            5    3           7       719712    6
   + 1360840 c1  EiA  - 900512 c1  EiA  + 60672 c1  EiA + ------ EiA  c2
                                                            5

              5   3            6   2           7      32964736   3    2
   - 14000 EiA  c1  - 29908 EiA  c1  + 5800 EiA  c1 - -------- c2  EiA
                                                         25

              3   2             6           2   4          4       8
   - 201696 c2  c1  - 6416 c2 c1  + 70392 c2  c1  + 2316 c2  + 4 c1

            8   34247072       3   2   233436672   2    2   2
   + 686 EiA  - -------- c1 EiA  c2  + --------- c1  EiA  c2
                   5                      25

     34868128   4    2                3    3                     3
   - -------- c1  EiA  c2 + 5158096 c1  EiA  c2 + 97088 c1 EiA c2
        5

     607608   3       2   26368   5          728352    5      \  8
   + ------ c1  EiA c2  - ----- c1  EiA c2 + ------ EiA  c2 c1| t
       5                    5                  5              /
proc(p, A)  ...  end;
> wedge := proc (p, A) local r, m; option remember; if nargs = 0 then RETURN(`usage: wedge(p,A)`) end if; r := rank(A); if p < 0 or r < p then 0 elif p = 0 then 1 elif p = 1 then A elif r < 2*p then tensor(determinant(A), dual(wedge(r-p, A))) else `schubert/strip`(add((-1)^(p-m+1)*tensor(wedge(m, A), subs(t = t*(p-m), A)), m = 0 .. p-1)/p) end if end proc;
proc(p, A)  ...  end;
> wedge(2, H);
     43  4    4   2789  7    6      27007  7    4   3   252243  8    4   4
28 + -- t  EiA  - ---- t  EiA  c1 - ----- t  EiA  c1  + ------ t  EiA  c1
     6            180                360                 1120

     668021  8    4   2   481  7    5   2   4091  8    7
   + ------ t  EiA  c2  + --- t  EiA  c1  - ---- t  EiA  c1
      450                 15                360

     11575  8    6   2   514  8    5   3   1211  7   6       5342  8   7
   + ----- t  EiA  c1  - --- t  EiA  c1  + ---- t  c1  EiA - ---- t  c1  EiA
      336                 9                 18               105

     8527  7   4    3   9591  7   5    2   12662  8   5    3
   + ---- t  c1  EiA  - ---- t  c1  EiA  - ----- t  c1  EiA
      30                 20                 35

     481667  8   6    2      2              3       2   697  6       5
   + ------ t  c1  EiA  + 6 t  c1 EiA + 12 t  c1 EiA  - --- t  c1 EiA
      1260                                              60

         3             4      2   8  5      3      5      2   101  6      4
   + 24 t  c1 c2 - 22 t  c2 c1  + - t  c2 c1  + 5 t  c1 c2  + --- t  c2 c1
                                  3                           10

     421  6   2   2   266972   2  7    3   18181   3  7
   - --- t  c2  c1  + ------ c2  t  EiA  + ----- c2  t  EiA
      5                125                  450

     168076   3  8    2                    2   12617  8    8   167  7    7
   - ------ c2  t  EiA  - 42 c1 t - 12 c2 t  - ----- t  EiA  - --- t  EiA
      1125                                     10080           90

         2   2       3   3      4   4      4   2   13  5   5    6   6
   + 30 t  c1  - 16 t  c1  + 8 t  c1  + 2 t  c2  - -- t  c1  + t  c1
                                                   5

         6   3   4148  4             3587  5          2   11369  5      2
   + 33 t  c2  - ---- t  c2 c1 EiA - ---- t  c2 c1 EiA  + ----- t  c2 c1  EiA
                  5                   15                   15

          6          3   105361  6      2    2   22487  6      3
   - 205 t  c2 c1 EiA  - ------ t  c2 c1  EiA  - ----- t  c2 c1  EiA
                           75                     30

     24509  6   2          231071  8   3   2   873  8   2   4   5717  8      6
   + ----- t  c2  c1 EiA + ------ t  c2  c1  - --- t  c2  c1  + ---- t  c2 c1
      150                   420                 5               420

     661  8   4   59   8   8        3   2       2732  3              4       3
   - --- t  c2  + --- t  c1  - 210 t  c1  EiA + ---- t  c2 EiA + 21 t  c1 EiA
     168          840                            5

     369  4   2    2        4   3       2742  4       2   13  5       4
   - --- t  c1  EiA  + 288 t  c1  EiA + ---- t  c2 EiA  + -- t  c1 EiA
      2                                  5                3

     263  5   2    3   161  5   3    2        5   4       5699  5       3
   - --- t  c1  EiA  + --- t  c1  EiA  - 135 t  c1  EiA + ---- t  c2 EiA
      2                 2                                  15

     432  5   2       867  6   2    4   71  6   3    3   1492  6   4    2
   - --- t  c2  EiA - --- t  c1  EiA  + -- t  c1  EiA  + ---- t  c1  EiA
      5               40                6                 5

     102  6   5       6599  6       4   702881  6   2    2   7   7   7
   - --- t  c1  EiA - ---- t  c2 EiA  + ------ t  c2  EiA  - -- t  c1
      5               150                375                 30

                    2    2   44  3    3   7  5    5   257  6    6
   + 28 EiA t + 26 t  EiA  + -- t  EiA  + - t  EiA  - --- t  EiA
                             3            5           180

     81581  7    4         575257  8    4      2   101557  8    5
   - ----- t  EiA  c1 c2 - ------ t  EiA  c2 c1  - ------ t  EiA  c2 c1
      900                   1400                    1575

     116237  7   4          38969  7   2       2   415679  8   5
   + ------ t  c1  EiA c2 - ----- t  c1  EiA c2  - ------ t  c1  EiA c2
      180                    150                    1050

     3329777  8   3       2   552568  8          3   325297  7   2    3
   + ------- t  c1  EiA c2  - ------ t  c1 EiA c2  - ------ t  c1  EiA  c2
      12600                    1575                   300

     2214427  7   3    2      1599119  8   3    3      30576031  8   4    2
   + ------- t  c1  EiA  c2 + ------- t  c1  EiA  c2 - -------- t  c1  EiA  c2
       900                     1050                     12600

     11437121  8   2    2   2   5269999  7       2   2
   + -------- t  c1  EiA  c2  - ------- t  c1 EiA  c2
       3000                      1500

     88930607  8       3   2   10393     7   3   32461  7   2   3
   - -------- t  c1 EiA  c2  - ----- c1 t  c2  + ----- t  c2  c1
      31500                     30                180

     796  7      5   21911  7    5      50152  8    6
   - --- t  c2 c1  - ----- t  EiA  c2 - ----- t  EiA  c2
     45               90                 225
>

toddclass

1+(3*c1-(5/2)*EiA)*t+((53/12)*c1^2-8*c1*EiA+(17/6)*EiA^2)*t^2+((17/4)*c1^3-(301/24)*c1^2*EiA+(39/4)*c1*EiA^2-(15/8)*EiA^3)*t^3+((721/240)*c1^4-(1097/90)*c1^3*EiA+(1151/72)*c1^2*EiA^2-(211/30)*c1*EiA^3+(137/180)*EiA^4+(1/120)*c1*c3-(1/120)*c2^2-(83/40)*c2*c1*EiA+(299/360)*c2*EiA^2+(2023/600)*c3*EiA)*t^4+((12113/1200)*c3*c1*EiA+(133/80)*c1^5+(1/40)*c3*c1^2-(1/40)*c1*c2^2-(1/6)*EiA^5+(1843/240)*c2*c1*EiA^2-(249/40)*c2*c1^2*EiA+(16/5)*c1*EiA^4-(17009/1440)*c1^2*EiA^3+(2177/144)*c1^3*EiA^2-(3697/480)*c1^4*EiA-(299/144)*c2*EiA^3+(1/48)*c2^2*EiA-(2023/240)*c3*EiA^2)*t^5+((1/504)*c3*c1*c2-(1/504)*c3^2+(1103/30240)*c3*c1^3-(1103/30240)*c2^2*c1^2+(1/6048)*c2^3+(4535/6048)*c1^6+(1/6048)*EiA^6-(62627/5040)*c2*c1*EiA^3+(641213/30240)*c2*c1^2*EiA^2-(7579/864)*c2*c1^3*EiA+(31/1680)*c2^2*c1*EiA-(134069/5040)*c3*c1*EiA^2+(95447/7200)*c3*c1^2*EiA+(83/1680)*c3*c2*EiA-(863/1008)*c1*EiA^5+(40583/7560)*c1^2*EiA^4-(6371/630)*c1^3*EiA^3+(436133/60480)*c1^4*EiA^2-(14561/5040)*c1^5*EiA+(11399/5040)*c2*EiA^4-(8219/7560)*c2^2*EiA^2+(33407/3600)*c3*EiA^3)*t^6+(-(2269/720)*c1*EiA^2*c2^2+(379571/12096)*c1^3*EiA^2*c2-(282481/8640)*c1^2*EiA^3*c2-(2477/60480)*c1^2*EiA*c2^2-(11003/1440)*c1^4*EiA*c2+(231233/20160)*c1*EiA^4*c2+(89/1344)*EiA^6*c1-(5/12096)*c2^3*EiA-(817/630)*c1^5*EiA^2-(4499/3780)*c1^4*EiA^3-(11981/60480)*c1^6*EiA-(16069/12096)*EiA^5*c2-(8423/6048)*EiA^5*c1^2+(7403/2160)*EiA^4*c1^3-(347/10080)*c2^2*c1^3+(1/2016)*c1*c2^3-(1/168)*c1*c3^2+(347/10080)*c3*c1^4+(2857/10080)*c1^7+(361/2520)*c3*c2*c1*EiA+(265/24192)*EiA^7+(3141563/100800)*c3*c1*EiA^3-(750719/20160)*c3*c1^2*EiA^2+(2879119/302400)*c3*c1^3*EiA-(83/672)*c3*c2*EiA^2-(5413/960)*c3*EiA^4+(5/1008)*c3^2*EiA+(32351/12096)*c2^2*EiA^3+(1/168)*c3*c2*c1^2)*t^7+((26225777/907200)*c2*EiA^4*c1^2-(3079609/907200)*EiA^4*c1^4-(338209/151200)*EiA^4*c2^2-(1478693/302400)*c1^6*EiA^2+(76603/9450)*c1^5*EiA^3+(822839/1587600)*c1^7*EiA+(38489/100800)*EiA^6*c2-(67717/604800)*EiA^5*c1^3+(385669/3628800)*EiA^6*c1^2+(86861/2116800)*EiA^7*c1-(179/11340)*c2^3*EiA^2+(199/518400)*c2^3*c1^2-(21397/907200)*c2^2*c1^4-(1/57600)*c2^4-(1/7200)*c3^2*c2-(397/48384)*c3^2*c1^2+(21397/907200)*c3*c1^5+(769873/8467200)*c1^8-(7771/3175200)*EiA^8-(20324299/1512000)*c3*c2*c1*EiA^2+(13277/168000)*c3*c2*c1^2*EiA-(31153523/1512000)*c3*c1*EiA^4+(67108399/1512000)*c3*c1^2*EiA^3-(31616863/1296000)*c3*c1^3*EiA^2+(1655483/504000)*c3*c1^4*EiA+(1887959/378000)*c3*c2*EiA^3-(7843/201600)*c3*c2^2*EiA+(129803/1008000)*c3^2*c1*EiA+(2836471/1512000)*c3*EiA^5+(12239183/1512000)*c3^2*EiA^2+(1/28800)*c3*c1*c2^2+(2383/604800)*c3*c2*c1^3+(840857/151200)*c1*EiA^3*c2^2-(10181/28800)*c1^2*EiA^2*c2^2+(6229313/226800)*c1^4*EiA^2*c2-(14618671/302400)*c1^3*EiA^3*c2+(13753/604800)*c1*EiA*c2^3-(18251/201600)*c1^3*EiA*c2^2-(2057837/453600)*c1^5*EiA*c2-(13589371/2116800)*EiA^5*c2*c1)*t^8+(-(2384411/241920)*EiA^4*c1^5+(1/23040)*c2^4*EiA+(1951/48384)*c2^3*EiA^3+(5413381/518400)*c1^6*EiA^3-(47191393/12700800)*c1^7*EiA^2+(13879/362880)*EiA^5*c2^2+(1283/120960)*EiA^7*c2+(23712043/7257600)*EiA^5*c1^4-(133057/403200)*EiA^6*c1^3+(1186351/25401600)*EiA^7*c1^2-(20963/1411200)*EiA^8*c1+(17/188160)*EiA^9-(109/16128)*c3^2*c1^3+(151/12096)*c3*c1^6+(144162511/3024000)*c3*c2*c1*EiA^3-(5680679/144000)*c3*c2*c1^2*EiA^2-(1044859/6048000)*c3*c2*c1^3*EiA-(47093/403200)*c3*c2^2*c1*EiA+(2654299/336000)*c3*c1*EiA^5-(167541041/6048000)*c3*c1^2*EiA^4+(62880029/3024000)*c3*c1^3*EiA^3+(6949139/3024000)*c3*c1^4*EiA^2-(5225149/9072000)*c3*c1^5*EiA-(7396211/604800)*c3*c2*EiA^4+(7843/80640)*c3*c2^2*EiA^2+(48382717/2016000)*c3^2*c1*EiA^2+(4381033/12096000)*c3^2*c1^2*EiA+(1/2880)*c3^2*c2*EiA-(5823/22400)*c3*EiA^6-(4081811/201600)*c3^2*EiA^3+(66457/2822400)*c1^9-(1/2400)*c3^2*c1*c2+(1/9600)*c3*c2^2*c1^2-(1217/201600)*c3*c2*c1^4+(1331479/3386880)*c1^8*EiA-(151/12096)*c1^5*c2^2-(407/1209600)*c1^3*c2^3-(1/19200)*c1*c2^4-(35761/18900)*c2*c1^6*EiA+(515143/7257600)*c2^3*c1^2*EiA-(67481/907200)*c2^2*c1^4*EiA-(5322839/907200)*c1^2*EiA^3*c2^2-(10183969/259200)*c1^4*EiA^3*c2-(8653/80640)*c1*EiA^2*c2^3+(1098641/134400)*c1^3*EiA^2*c2^2+(2850427/226800)*c1^5*EiA^2*c2-(127474757/8467200)*EiA^5*c2*c1^2+(2786573/1411200)*EiA^6*c2*c1-(474203/604800)*c1*EiA^4*c2^2+(146288659/3628800)*c2*EiA^4*c1^3)*t^9

