{-# LANGUAGE ParallelListComp #-}
-- |Stirling's approximation to the gamma function
module Math.Gamma.Stirling (lnGammaStirling, cs, s, abs_s, terms) where

import qualified Data.Vector as V

-- |Convergent when Re(z) > 0.  The first argument is the c_n series to use 
-- ('cs' is an ineffecient but generic definition of the full infinite series.
-- Some precomputed finite prefix of 'cs' should be fed to this function, the 
-- length of which will determine the accuracy achieved.)
{-# INLINE lnGammaStirling #-}
lnGammaStirling cs z = (z - 0.5) * log z - z + 0.5 * log (2*pi) + sum [c / q | c <- cs | q <- risingPowers (z+1)]
    where

{-# INLINE risingPowers #-}
risingPowers x = scanl1 (*) (iterate (1+) x)

-- |The c_n series in the convergent version of stirling's approximation given
-- on wikipedia at
-- http://en.wikipedia.org/wiki/Stirling%27s_approximation#A_convergent_version_of_Stirling.27s_formula
-- as fetched on 11 June 2010.
cs :: (Fractional a, Ord a) => [a]
cs = map c [1..]

c :: (Fractional a, Ord a) => Int -> a
c n = 0.5 * recip n' * sum [k' * fromInteger (abs_s n k) / ((k' + 1) * (k' + 2)) | k <- [1..n], let k' = fromIntegral k]
    where n' = fromIntegral n

-- |The (signed) Stirling numbers of the first kind.
s :: Int -> Int -> Integer
s n k
    | n < 0     = error "s n k: n < 0"
    | k < 0     = error "s n k: k < 0"
    | k > n     = error "s n k: k > n"
    | otherwise = s n k
    
    where
        table = [V.generate (n+1) $ \k -> s n k | n <- [0..]]
        s 0 0 = 1
        s _ 0 = 0
        s n k 
            | n == k    = 1
            | otherwise = s (n-1) (k-1) - (toInteger n-1) * s (n-1) k
            where
                s n k = table !! n V.! k

-- |The (unsigned) Stirling numbers of the first kind.
abs_s :: Int -> Int -> Integer
abs_s n k
    | n < 0     = error "abs_s n k: n < 0"
    | k < 0     = error "abs_s n k: k < 0"
    | k > n     = error "abs_s n k: k > n"
    | otherwise = abs_s n k
    
    where
        table = [V.generate (n+1) $ \k -> abs_s n k | n <- [0..]]
        abs_s 0 0 = 1
        abs_s _ 0 = 0
        abs_s n k 
            | n == k    = 1
            | otherwise = abs_s (n-1) (k-1) + (toInteger n-1) * abs_s (n-1) k
            where
                abs_s n k = table !! n V.! k

-- Compute the number of terms required to achieve a given precision for a
-- given value of z.  The mamimum will typically (always?) be around 1, and 
-- seems to be more or less independent of the precision desired (though not 
-- of the machine epsilon - essentially, near zero I think this method is
-- extremely numerically unstable).
terms prec z = converge (eps z) (f z)
    where
        cs' = cs
        f z = scanl1 (+) [c / q | c <- cs' | q <- risingPowers (z+1)]
        -- (eps is 0 at z=0.86639115674955 and z=2.087930091329227)
        eps z = prec * abs ((z - 0.5) * log z - z + 0.5 * log (2*pi))
        converge eps xs = go 1 xs where go n (x:y:zs) | abs(x-y)<=eps = n | otherwise = go (n+1) (y:zs)

f z = scanl1 (+) [c / q | c <- cs | q <- risingPowers (z+1)]
