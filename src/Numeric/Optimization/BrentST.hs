  {-# LANGUAGE NamedFieldPuns #-}

module Numeric.Optimization.BrentST (
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

data FVal = FVal {
    val :: {-# UNPACK #-} !Double
  , fval :: {-# UNPACK #-} !Double
  } deriving (Eq, Show, Read)

data BrentState = BrentState {
    fu :: {-# UNPACK #-} !FVal
  , fv :: {-# UNPACK #-} !FVal
  , fw :: {-# UNPACK #-} !FVal
  , fx :: {-# UNPACK #-} !FVal
  , a, b, d, e, xm, p, q, r, tol1, tol2 :: {-# UNPACK #-} !Double
  }

fmin :: Double -> (Double -> Double) -> Double -> Double -> FVal
fmin tol f ax bx = fx fstart
  where
    fstart = loop20 $ BrentState {
        fu = FVal 0 0
      , fv = FVal v' fx'
      , fw = FVal w' fx'
      , fx = FVal x' fx'
      , a = a'
      , b = b'
      , e = e'
      , d = 0, xm = 0, p = 0, q = 0, r = 0, tol1 = 0, tol2 = 0
      }
      where
        a' = ax
        b' = bx
        v' = a' + c * (b' - a')
        w' = v'
        x' = v'
        e' = 0
        fx' = f x'

    loop20 bs@BrentState{fv, fw, fx, a, b, d, e}
      | abs (val fx - xm') <= tol2' - 0.5 * (b - a) = bs { xm = xm', tol1 = tol1', tol2 = tol2' }
      | abs e <= tol1' = loop40 $ bs { xm = xm', tol1 = tol1', tol2 = tol2' }
      | otherwise = loop30 $ bs { e = e', p = p'', q = q''', r = r'', xm = xm', tol1 = tol1', tol2 = tol2' }
      where
        xm' = 0.5 * (a + b)
        tol1' = eps * abs (val fx) + tol / 3
        tol2' = 2 * tol1'
        r' = (val fx - val fw) * (fval fx - fval fv)
        q' = (val fx - val fv) * (fval fx - fval fw)
        p' = (val fx - val fv) * q' - (val fx - val fw) * r'
        q'' = 2 * (q' - r')
        p'' = if q'' > 0 then negate p' else p'
        q''' = abs q''
        r'' = e
        e' = d

    loop30 bs@BrentState{fu, fx, a, b, d, xm, p, q, r, tol1, tol2}
      | abs p >= abs (0.5 * q * r) = loop40 bs
      | p <= q * (a - val fx) = loop40 bs
      | p >= q * (b - val fx) = loop40 bs
      | otherwise = loop40 $ bs { fu = FVal u' (fval fu), d = d''' }
      where
        d' = p / q
        u' = val fx + d
        d'' = if (val fu - a) < tol2 then signum (xm - val fx) * tol1 else d'
        d''' = if (b - val fu) < tol2 then signum (xm - val fx) * tol1 else d''

    loop40 bs@BrentState{fx, a, b, e, xm} = loop50 $ bs { d = d', e = e'' }
      where
        e' = if val fx >= xm then a - val fx else e
        e'' = if val fx < xm then b - val fx else e'
        d' = c * e''

    loop50 bs@BrentState{fu, fw, fx, a, b, d, tol1}
      | fval fu' >= fval fx = loop60 $ bs { fu = fu' }
      | otherwise = loop20 $ bs { fu = fu', fv = fv', fw = fw', fx = fx', a = a', b = b' }
      where
        u' = if abs d >= tol1 then val fx + d else val fu
        u'' = if abs d < tol1 then val fx + signum d * tol1 else u'
        fu' = FVal u'' (f u'')
        a' = if val fu' > val fx then val fx else a
        b' = if val fu' < val fx then val fx else b
        fv' = fw
        fw' = fx
        fx' = fu'

    loop60 bs@BrentState{fu, fw, fv, fx, a, b}
      | fval fu <= fval fw || val fw == val fx = loop70 $ bs { a = a', b = b' }
      | (fval fu <= fval fv) || (val fv == val fx) || (val fv == val fw) = loop80 $ bs { a = a', b = b' }
      | otherwise = loop20 $ bs { a = a', b = b' }
      where
        a' = if val fu <  val fx then val fu else a
        b' = if val fu >= val fx then val fu else b

    loop70 bs@BrentState{fu, fw} = loop20 $ bs { fv = fw, fw = fu }

    loop80 bs@BrentState{fu} = loop20 $ bs { fv = fu }
