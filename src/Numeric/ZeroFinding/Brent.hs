  {-# LANGUAGE NamedFieldPuns #-}

module Numeric.ZeroFinding.Brent (
    FVal(..)
  , BrentError(..)
  , BrentResult
  , zeroin
  ) where

import Numeric.Optimization.Brent (FVal(..))

{-
      a zero of the function  f(x)  is computed in the interval ax,bx .

  input..

  ax     left endpoint of initial interval
  bx     right endpoint of initial interval
  f      function subprogram which evaluates f(x) for any x in
         the interval  ax,bx
  tol    desired length of the interval of uncertainty of the
         final result (.ge.0.)

  output..

  zeroin abscissa approximating a zero of  f  in the interval ax,bx

      it is assumed  that   f(ax)   and   f(bx)   have  opposite  signs
  this is checked, and an error message is printed if this is not
  satisfied.   zeroin  returns a zero  x  in the given interval
  ax,bx  to within a tolerance  4*macheps*abs(x)+tol, where macheps  is
  the  relative machine precision defined as the smallest representable
  number such that  1.+macheps .gt. 1.
      this function subprogram is a slightly  modified  translation  of
  the algol 60 procedure  zero  given in  richard brent, algorithms for
  minimization without derivatives, prentice-hall, inc. (1973).
  -}

--  eps is approximately the square root of the relative machine
--  precision.
eps :: Double
eps = sqrt . head . dropWhile (\x -> 1 + x > 1) $ iterate (/ 2) 1

data BrentError = NotOppositeSigns
  deriving (Eq, Show, Read)

type BrentResult = Either BrentError

data BrentState = BrentState {
    fa :: {-# UNPACK #-}!FVal
  , fb :: {-# UNPACK #-}!FVal
  , fc :: {-# UNPACK #-}!FVal
  , d, e, tol1, xm, p, q, r, s :: {-# UNPACK #-}!Double
  }

zeroin :: Double -> (Double -> Double) -> Double -> Double -> BrentResult FVal
zeroin tol f ax bx
  | signum fa0 == signum fb0 = Left NotOppositeSigns
  | otherwise = Right . fb . loop20 $ BrentState {
      fa = FVal ax fa0
    , fb = FVal bx fb0
    , fc = FVal 0 0
    , d = 0, e = 0, tol1 = eps + 1, xm = 0, p = 0, q = 0, r = 0, s = 0
    }
  where
    fa0 = f ax
    fb0 = f bx

    loop20 bs@BrentState{fa, fb} = loop30 $ bs { fc = fa, d = d', e = e'}
      where
        d' = val fb - val fa
        e' = d'

    loop30 bs@BrentState{fa, fb, fc}
      | abs (fval fc) >= abs (fval fb) = loop40 bs
      | otherwise = loop40 $ bs { fa = fb, fb = fa, fc = fb }

    loop40 bs@BrentState{fa, fb, fc, e, tol1}
      | abs xm' <= tol1 || fval fb == 0 = bs { xm = xm', tol1 = tol1' }
      | (abs e >= tol1') && (abs (fval fa) > abs (fval fb)) = loop50 $ bs { xm = xm', tol1 = tol1' }
      | otherwise = loop110 $ bs { d = xm', e = xm', xm = xm', tol1 = tol1' }
      where
        tol1' = 2.0 * eps * abs (val fb) + 0.5 * tol
        xm' = 0.5 * (val fc - val fb)

    loop50 bs@BrentState{fa, fb, fc, xm}
      | val fa /= val fc = loop60 $ bs { s = s' }
      | otherwise = loop70 $ bs { p = 2 * xm * s', q = 1 - s', s = s' }
      where s' = fval fb / fval fa

    loop60 bs@BrentState{fa, fb, fc, xm, s} = loop70 $ bs { q = q'', p = p', r = r' }
      where
        q' = fval fa / fval fc
        r' = fval fb / fval fc
        p' = s * (2 * xm * q' * (q' - r') - (val fb - val fa) * (r' - 1))
        q'' = (q' - 1) * (r' - 1) * (s - 1)

    loop70 bs@BrentState{p, q}
      | p <= 0 = loop80 bs
      | otherwise = loop90 $ bs { q = negate q }

    loop80 bs@BrentState{p} = loop90 $ bs { p = negate p }

    loop90 bs@BrentState{d, e, xm, tol1, p, q}
      | (2 * p > 3 * xm * q - abs (tol1 * q)) || (p >= abs (5 * s' * q)) = loop100 $ bs { e = e', s = s' }
      | otherwise = loop110 $ bs { d = d', e = e', s = s' }
      where
        s' = e
        e' = d
        d' = p / q

    loop100 bs@BrentState{xm} = loop110 $ bs { d = xm, e = xm }

    loop110 bs@BrentState{fb, d, tol1}
      | abs d <= tol1 = loop120 $ bs { fa = fa' }
      | otherwise = loop140 $ bs { fa = fa', fb = FVal (val fb + d) (fval fb) }
      where fa' = fb

    loop120 bs@BrentState{fb, tol1, xm}
      | xm <= 0 = loop130 bs
      | otherwise = loop140 $ bs { fb = FVal (val fb + tol1) (fval fb) }

    loop130 bs@BrentState{fb, tol1} = loop140 $ bs { fb = FVal (val fb - tol1) (fval fb) }

    loop140 bs@BrentState{fb, fc}
      | signum (fval fb') == signum (fval fc) = loop20 $ bs { fb = fb' }
      | otherwise = loop30 $ bs { fb = fb' }
      where fb' = FVal (val fb) (f (val fb))
