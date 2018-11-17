
module Numeric.Optimization.Brent (
    FVal(..)
  , fmin
  ) where

{-
      an approximation  x  to the point where  f  attains a minimum  on
  the interval  (ax,bx)  is determined.


  input..

  ax    left endpoint of initial interval
  bx    right endpoint of initial interval
  f     function subprogram which evaluates  f(x)  for any  x
        in the interval  (ax,bx)
  tol   desired length of the interval of uncertainty of the final
        result ( .ge. 0.0d0)


  output..

  fmin  abcissa approximating the point where  f  attains a minimum


      the method used is a combination of  golden  section  search  and
  successive parabolic interpolation.  convergence is never much slower
  than  that  for  a  fibonacci search.  if  f  has a continuous second
  derivative which is positive at the minimum (which is not  at  ax  or
  bx),  then  convergence  is  superlinear, and usually of the order of
  about  1.324....
      the function  f  is never evaluated at two points closer together
  than  eps*abs(fmin) + (tol/3), where eps is  approximately the square
  root  of  the  relative  machine  precision.   if   f   is a unimodal
  function and the computed values of   f   are  always  unimodal  when
  separated by at least  eps*abs(x) + (tol/3), then  fmin  approximates
  the abcissa of the global minimum of  f  on the interval  ax,bx  with
  an error less than  3*eps*abs(fmin) + tol.  if   f   is not unimodal,
  then fmin may approximate a local, but perhaps non-global, minimum to
  the same accuracy.
      this function subprogram is a slightly modified  version  of  the
  algol  60 procedure  localmin  given in richard brent, algorithms for
  minimization without derivatives, prentice - hall, inc. (1973).
-}

--  eps is approximately the square root of the relative machine
--  precision.
eps :: Double
eps = sqrt . head . dropWhile (\x -> 1 + x > 1) $ iterate (/ 2) 1

--  c is the squared inverse of the golden ratio
c :: Double
c = 0.5 * (3 - sqrt 5.0)

data FVal = FVal { val :: {-# UNPACK #-} !Double, fval :: {-# UNPACK #-} !Double }
  deriving (Eq, Show, Read)

data Bounds = Bounds { bLeft :: {-# UNPACK #-} !Double, bRight :: {-# UNPACK #-} !Double }
  deriving (Eq, Show, Read)

bnddist, midpoint :: Bounds -> Double
bnddist (Bounds a b) = b - a
midpoint (Bounds a b) = a + 0.5 * (b - a)

data Section = Section {
    fvalV :: {-# UNPACK #-} !FVal
  , fvalW :: {-# UNPACK #-} !FVal
  , fvalX :: {-# UNPACK #-} !FVal
  , sectionBounds :: {-# UNPACK #-} !Bounds
  , dist :: {-# UNPACK #-} !Double
  }

-- translation of http://www.netlib.org/fmm/fmin.f
fmin :: Double -> Double -> (Double -> Double) -> Double -> FVal
fmin ax bx f tol = loop f tol $ Section fv fv fv (Bounds ax bx) 0.0
  where
    bnd = Bounds ax bx
    v = bLeft bnd + c * bnddist bnd
    fv = FVal v (f v)

loop :: (Double -> Double) -> Double -> Section -> FVal
loop f tol s
  | abs (val (fvalX s) - xm) <= (tol2 - 0.5 * bnddist (sectionBounds s)) = fvalX s
  | abs (dist s) < tol1 = loop f tol goldenSection
  | otherwise = loop f tol fitParabola
  where
    xm = midpoint (sectionBounds s)
    tol1 = eps * abs (val (fvalX s)) + tol / 3.0
    tol2 = 2.0 * tol1

    fitParabola
      | abs p >= abs (0.5 * q * e) || p <= q * (bLeft bnd - val fx) || p >= q * (bRight bnd - val fx) = goldenSection
      | otherwise = update s { dist = d } fu
      where
        fv = fvalV s
        fw = fvalW s
        fx = fvalX s
        bnd = sectionBounds s
        e = dist s

        r = (val fx - val fw) * (fval fx - fval fv)
        qq = (val fx - val fv) * (fval fx - fval fw)
        pp = (val fx - val fv) * qq - (fval fx - fval fw) * r
        q = 2.0 * (qq - r)
        p = if q > 0 then negate pp else pp
        aq = abs q
        --
        --  a parabolic interpolation step
        --
        d = p / aq
        u = val fx + d
        fu = bumpDiff fx (if u - bLeft bnd < tol2 || bRight bnd - u < tol2 then signum (xm - val fx) * tol1 else d)

    goldenSection = update s { dist = d } $ bumpDiff fx d
      where
        fx = fvalX s
        bnd = sectionBounds s

        ee = if val fx >= xm then bLeft bnd - val fx else bRight bnd - val fx
        d = c * ee

    bumpDiff fx d = let u = if abs d >= tol1 then val fx + d else val fx + signum d * tol1 in FVal u (f u)

    update (Section fv fw fx bnd e) fu
      | fval fu <= fval fx = Section fw fx fu bndLE e
      | fval fu <= fval fw || val fw == val fx = Section fw fu fx bndGT e
      | val fu <= fval fv || val fv == val fx || val fv == val fw = Section fu fw fx bndGT e
      | otherwise = Section fw fv fx bndGT e
      where
        bndLE = if val fu >= val fx then bnd { bLeft = val fx } else bnd { bRight = val fx }
        bndGT = if val fu <= val fx then bnd { bLeft = val fu } else bnd { bRight = val fu }
