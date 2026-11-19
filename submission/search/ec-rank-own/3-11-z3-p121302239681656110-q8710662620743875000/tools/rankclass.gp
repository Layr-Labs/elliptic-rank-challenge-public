\\ Owner stopgap ranker (generic over the isogeny class). Calibration + late stage.
\\ rankclass(E, eff, target): E minimal model. Returns [lo, hi, pts_on_E, ms].
rankclass(E, eff, target) = {
  my(t0 = getabstime(), iso = ellisomat(E), L = iso[1], n = #L, pts = [], hi = 10^6, r, lo);
  for(i = 1, n,
    my(Ei = ellinit(L[i][1]));
    r = ellrank(Ei, eff);
    hi = min(hi, r[2]);
    if(hi < target, return([0, hi, [], getabstime() - t0]));
    if(#r[4], pts = concat(pts, apply(P -> ellisogenyapply(L[i][3], P), r[4])));
  );
  pts = select(P -> #P == 2, pts);      \\ drop points at infinity ([0])
  pts = Set(pts);
  r = ellrank(E, eff, pts);
  [r[1], min(hi, r[2]), r[4], getabstime() - t0]
};
