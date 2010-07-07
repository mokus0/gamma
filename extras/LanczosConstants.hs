-- This makes use of a not-yet-released matrix library.  It could be rewritten
-- to use any of the existing ones on hackage, but I don't know of any of them
-- that support matrices over arbitrary types - they are all focused on
-- efficiently packing the matrices and/or calling foreign libraries
-- (BLAS/GSL/etc.) and do not support any types other than Double, Float, and
-- Complex Double/Float.
-- 
-- I am keeping it around anyway and including this file in the source 
-- distribution, because with a very small amount of work an end-user 
-- could fill in the gaps and use this code to generate their own constants 
-- for lanczos gamma function approximations, which one may wish to do if 
-- they wanted to implement, say, a gamma function for a very high precision
-- floating point type.
--
-- Note that these really need to be run with significantly higher precision
-- than the target type or truncation error will make the results useless.
-- 
-- The algorithm implemented here is by Paul Godfrey, and is described in full
-- at http://www.numericana.com/answer/info/godfrey.htm (as of 21 June 2010).
module LanczosConstants where

import Math.Matrix
import Math.Matrix.Alias

cs g n = vectorToList (applyRat dbc f)
    where
        applyRat :: (Real t, Fractional t) => IMatrix Rational -> IVector t -> IVector t
        applyRat m v = fromRatVec (apply m (toRatVec v))
        fromRatVec :: (Vector v t, Fractional t) => IVector Rational -> v t
        fromRatVec = convertByV fromRational
        toRatVec :: (Vector v t, Real t) => v t -> IVector Rational
        toRatVec   = convertByV toRational
        
        dbc = dbcMat n
        f = fVec g n

dbcMat n = multRat d (multRat b c)
    where
        multRat :: (Real a, Matrix m1 a, Real b, Matrix m2 b) => m1 a -> m2 b -> IMatrix Rational
        multRat = multiplyWith sum (\d b -> toRational d * toRational b)
        
        d = dMat n
        b = bMat n
        c = cMat n

fVec :: (Floating b, Vector v b) => b -> Int -> v b
fVec g n = vector n f
    where
        f a = sqrt (2 / pi)
            * product [fromIntegral i - 0.5 | i <-[1..a]]
            * exp (a' + g + 0.5)
            / (a' + g + 0.5) ** (a' + 0.5)
            where a' = fromIntegral a

cMat :: Int -> IMatrix Rational
cMat n = matrix n n m
    where
        m 0 0 = 1/2
        m i j = fromInteger (c (2*i) (2*j))
        
        c 0 0 = 1
        c 1 1 = 1
        c i 0 = negate (c (i-2) 0)
        c i j
            | i == j    = 2 * c (i-1) (j-1)
            | i > j     = 2 * c (i-1) (j-1) - c (i-2) j
            | otherwise = 0

dMat :: Int -> IAlias Mat Integer
dMat n = AsDiag (IVec (ivector n dFunc)) 0
    where
        dFunc    0  = 1
        dFunc (i+1) = negate (factorial (2*i+2) `div` (2 * factorial i * factorial (i+1)))
        factorial n = product [1..toInteger n]

bMat :: Int -> IMatrix Integer
bMat n = matrixFromList bList
    where
        bList = take n . map (take n) $
            repeat 1 : 
            [ replicate i 0 ++ bicofs (negate (toInteger i*2))
            | i <- [1..]
            ]
            
        
        bFunc 0 _ = 1
        bFunc i j
            | i > j = 0
        bFunc i j = bicofs (toInteger (2 * j - 1)) !! i
        
        bicofs x = go x 1 1
            where
                go num denom x = x : go (num+signum num) (denom+signum denom) (x * num `div` denom)

-- 
-- p g k = sum [c (2*k+1) (2*a+1) * f a | a <- [0..k]]
--         where
--             k' = fromIntegral k
--             f a = 
{-# INLINE risingPowers #-}
risingPowers x = scanl1 (*) (iterate (1+) x)

{-# INLINE fallingPowers #-}
fallingPowers x = scanl1 (*) (iterate (subtract 1) x)
