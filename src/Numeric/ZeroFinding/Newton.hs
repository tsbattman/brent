
module Numeric.ZeroFinding.Newton (
    NewtonCriteria(..)
  , NewtonError(..)
  , NewtonResult
  , newtonZero
  ) where

data NewtonCriteria = NewtonCriteria {
    newtonMaxIter :: Int
  , newtonEpsilon :: Double
  , newtonTolerance :: Double
  } deriving (Eq, Show, Read)

data NewtonError = NewtonMaxIter | NewtonStuck
  deriving (Eq, Show, Read)

type NewtonResult = Either NewtonError

newtonZero :: NewtonCriteria -> (Double -> Double) -> (Double -> Double) -> Double -> NewtonResult Double
newtonZero crit f f' = go 0
  where
    go n x0
      | n >= newtonMaxIter crit = Left NewtonMaxIter
      | y' <= newtonEpsilon crit = Left NewtonStuck
      | abs dx <= newtonTolerance crit * abs x1 = Right x1
      | otherwise = go (n + 1) x1
      where
        y = f x0
        y' = f' x0
        dx = y / y'
        x1 = x0 - dx
