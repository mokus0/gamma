module Math.Gamma where

class Floating a => Gamma a where
    gamma :: a -> a
    lnGamma :: a -> a
    lnFactorial :: Integral b => b -> a

instance Gamma Float where
instance Gamma Double where
