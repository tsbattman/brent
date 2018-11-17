
module Main where

import Criterion.Main

import Numeric.Optimization.Brent
import Numeric.ZeroFinding.Brent
import qualified Numeric.Optimization.BrentST as BT

parab, hump :: Double -> Double -> Double
parab a x = (x - a)^(2 :: Int)
hump a x = exp $ negate (x - a)^(2 :: Int)

brentOpt :: Benchmark
brentOpt = bgroup "brent optimization" [
    bench "parabola" $ whnf (fmin 0 1 (parab 0.5)) 0.0001
  , bench "hump" $ whnf (fmin (-2) 2 (hump 0)) 0.0001
  ]

brentOptST :: Benchmark
brentOptST = bgroup "brent st optimization" [
    bench "parabola" $ whnf (BT.fmin  0.0001 (parab 0.5) 0) 1
  , bench "hump" $ whnf (BT.fmin 0.0001 (hump 0) (-2)) 2
  ]

brentZero :: Benchmark
brentZero = bgroup "brent zerofinding" [
    bench "parabola" $ whnf (zeroin 0.0001 (parab 0.5) 0) 1
  , bench "hump" $ whnf (zeroin 0.0001 (hump 0) (-2)) 2
  , bench "sqrt" $ whnf (zeroin 0.0001 (\x -> x*x - 5) 0) 5
  ]

main :: IO ()
main = defaultMain [
    brentOpt
  , brentOptST
  , brentZero
  ]

