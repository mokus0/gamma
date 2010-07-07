{-# LANGUAGE ParallelListComp #-}
module Math.Gamma.Incomplete
    ( lowerGammaCF, pCF
    , upperGammaCF, qCF, qNeg
    , lowerGammaHypGeom, pHypGeom
    ) where

import {-# SOURCE #-}  Math.Gamma
import Math.ContinuedFraction
import Math.Sequence.Converge

-- |Special case of Kummer's confluent hypergeometric function, used
-- in lower gamma functions.
-- 
-- m_1_sp1 s z = M(1;s+1;z)
-- 
m_1_sp1 s z = converge . scanl (+) 0 . scanl (*) 1 $
    [z / x | x <- iterate (1+) (s+1)]

-- |Lower incomplete gamma function, computed using Kummer's confluent
-- hypergeometric function M(a;b;x).  Specifically, this uses the identity:
-- 
-- gamma(s,x) = x**s * exp (-x) / s * M(1; 1+s; x)
-- 
-- From Abramowitz & Stegun (6.5.12).
lowerGammaHypGeom 0 0 = 0/0
lowerGammaHypGeom s x = sign (exp (log (abs x) * s - x) / s * m_1_sp1 s x)
        where
            sign 
                | x < 0 = case properFraction s of
                    (sI, 0) | s < 0     -> const (0/0)
                            | even sI   -> id 
                            | otherwise -> negate
                    _                   -> const (0/0)
                | otherwise = id

-- |Regularized lower incomplete gamma function, computed using Kummer's
-- confluent hypergeometric function.  Uses same identity as 'lowerGammaHypGeom'.
pHypGeom 0 0 = 0/0
pHypGeom s x
    | s < 0
    = sin (pi*s) / (-pi)
    * exp (s * log x - x + lnGamma  (-s)) * m_1_sp1 s x

    | s == 0 || x == 0
    = 0

    | otherwise
    = exp (s * log x - x - lnGamma (s+1)) * m_1_sp1 s x

-- Only valid for infinite lists.  Only used in the following definitions of continued fractions.
interleave (x:xs) (y:ys) = x:y:interleave xs ys

-- |Continued fraction representtion of the lower incomplete gamma function.
lowerGammaCF s z = gcf 0
    [ (p,q)
    | p <- exp (s * log z - z)
        : interleave
            [negate spn * z | spn <- [s..]]
            [n * z   | n   <- [1..]]
    | q <- [s..]
    ]

-- |Continued fraction representtion of the upper incomplete gamma function.
upperGammaCF s z = gcf 0
    [ (p,q)
    | p <- exp (log z * s - z)
        : zipWith (*) [1..] (iterate (subtract 1) (s-1))
    | q <- [n + z - s | n <- [1,3..]]
    ]

-- A common subexpression appearing in both 'pCF' and 'qCF'.
pow_x_s_div_gamma_s_div_exp_x s x 
    | x > 0     = exp (log x * s - x - lnGamma s)
    | otherwise = x ** s / (exp x * gamma s)


-- |Continued fraction representtion of the regularized lower incomplete gamma function.
pCF s x = gcf 0
    [ (p,q)
    | p <- pow_x_s_div_gamma_s_div_exp_x s x
        : interleave
            [negate spn * x | spn <- [s..]]
            [n * x          | n   <- [1..]]
    | q <- [s..]
    ]

-- |Continued fraction representtion of the regularized upper incomplete gamma function.
qCF s x = gcf 0
    [ (p,q)
    | p <- pow_x_s_div_gamma_s_div_exp_x s x
        : zipWith (*) [1..] (iterate (subtract 1) (s-1))
    | q <- [n + x - s | n <- [1,3..]]
    ]

-- |Q(s,x) for real x < 0
qNeg s x = case properFraction s of
    (sI, 0) | s > 0 -> exp (-x) * sum (scanl (*) 1 [x / fromIntegral k | k <- [1 .. sI-1]])
    _               -> 0/0
