{-# LANGUAGE FlexibleInstances #-}
module Math.Gamma
    ( Gamma(..)
    , Factorial(..)
    , IncGamma(..)
    , beta
    ) where

import Math.Gamma.Lanczos
import Math.Gamma.Incomplete

import Data.Complex
import Data.List (sortBy)
import Data.Ord (comparing)
import qualified Data.Vector.Unboxed as V
import Math.ContinuedFraction
import Math.Sequence.Converge

-- |Gamma function.  Minimal definition is ether gamma or lnGamma.
class Floating a => Gamma a where
    -- |The gamma function:  gamma z == integral from 0 to infinity of
    -- @\t -> t**(z-1) * exp (negate t)@
    gamma :: a -> a
    gamma 0 = 0/0
    gamma z
        | z == abs z    = exp (lnGamma z)
        | otherwise     = pi / (sin (pi * z) * exp (lnGamma (1-z)))


    -- |Natural log of the gamma function
    lnGamma :: a -> a
    lnGamma z = log (gamma z)
    
    -- |Natural log of the factorial function
    lnFactorial :: Integral b => b -> a
    lnFactorial n = lnGamma (fromIntegral n+1)

instance Gamma Float where
    gamma = realToFrac . reflect (gammaLanczos g cs) . realToFrac
        where
            g :: Double
            g = pi
            cs = [1.0000000249904433,9.100643759042066,-4.3325519094475,
                  0.12502459858901147,1.1378929685052916e-4,-9.555011214455924e-5]
    
    lnGamma = realToFrac . reflectLn (lnGammaLanczos g cs) . realToFrac
        where
            g :: Double
            g = pi
            cs = [1.0000000249904433,9.100643759042066,-4.3325519094475,
                  0.12502459858901147,1.1378929685052916e-4,-9.555011214455924e-5]
    
    lnFactorial n
        | n' < 0                = error "lnFactorial n: n < 0"
        | n' < toInteger nFacs  = facs V.! fromIntegral n
        | otherwise             = lnGamma (fromIntegral n+1)
        where
            n' = toInteger n
            nFacs       = 2000 -- limited only by time and space
            facs        = V.map lnGamma (V.enumFromN 1 nFacs)

instance Gamma Double where
    gamma = reflect (gammaLanczos g cs)
        where
            g = 2*pi
            cs = [1.0000000000000002,311.60117505414695,-498.65119046033163,244.08472899875767,-38.67036462939322,1.3350899103585203,-1.8972831806242229e-3,-3.935368195357295e-7,2.592464641764731e-6,-3.2263565156368265e-6,2.5666169886566876e-6,-1.3737776806198937e-6,4.4551204024819644e-7,-6.576826592057796e-8]

    lnGamma = reflectLn (lnGammaLanczos g cs)
        where
            g = exp pi / pi
            cs = [1.0000000000000002,1002.5049417114732,-1999.6140446432912,1352.1626218340114,-360.6486475548049,33.344988357090685,-0.6637188712004668,5.16644552377916e-4,1.684651140163429e-7,-1.8148207145896904e-7,6.171532716135051e-8,-9.014004881476154e-9]

    lnFactorial n
        | n' < 0                = error "lnFactorial n: n < 0"
        | n' < toInteger nFacs  = facs V.! fromIntegral n
        | otherwise             = lnGamma (fromIntegral n+1)
        where
            n' = toInteger n
            nFacs       = 2000 -- limited only by time and space
            facs        = V.map lnGamma (V.enumFromN 1 nFacs)

complexDoubleToFloat :: Complex Double -> Complex Float
complexDoubleToFloat (a :+ b) = realToFrac a :+ realToFrac b
complexFloatToDouble :: Complex Float -> Complex Double
complexFloatToDouble (a :+ b) = realToFrac a :+ realToFrac b

instance Gamma (Complex Float) where
    gamma = complexDoubleToFloat . reflectC (gammaLanczos g cs) . complexFloatToDouble
        where
            g = pi
            cs = [1.0000000249904433,9.100643759042066,-4.3325519094475,
                  0.12502459858901147,1.1378929685052916e-4,-9.555011214455924e-5]
    
    lnGamma = complexDoubleToFloat . reflectLnC (lnGammaLanczos g cs) . complexFloatToDouble
        where
            g = pi
            cs = [1.0000000249904433,9.100643759042066,-4.3325519094475,
                  0.12502459858901147,1.1378929685052916e-4,-9.555011214455924e-5]

    
    lnFactorial n
        | n' < 0                = error "lnFactorial n: n < 0"
        | n' < toInteger nFacs  = facs V.! fromIntegral n
        | otherwise             = lnGamma (fromIntegral n+1)
        where
            n' = toInteger n
            nFacs       = 2000 -- limited only by time and space
            facs        = V.map lnGamma (V.enumFromN 1 nFacs)

instance Gamma (Complex Double) where
    gamma = reflectC (gammaLanczos g cs)
        where
            g = 2*pi
            cs = [1.0000000000000002,311.60117505414695,-498.65119046033163,244.08472899875767,-38.67036462939322,1.3350899103585203,-1.8972831806242229e-3,-3.935368195357295e-7,2.592464641764731e-6,-3.2263565156368265e-6,2.5666169886566876e-6,-1.3737776806198937e-6,4.4551204024819644e-7,-6.576826592057796e-8]

    lnGamma = reflectLnC (lnGammaLanczos g cs)
        where
            g = exp pi / pi
            cs = [1.0000000000000002,1002.5049417114732,-1999.6140446432912,1352.1626218340114,-360.6486475548049,33.344988357090685,-0.6637188712004668,5.16644552377916e-4,1.684651140163429e-7,-1.8148207145896904e-7,6.171532716135051e-8,-9.014004881476154e-9]

    lnFactorial n
        | n' < 0                = error "lnFactorial n: n < 0"
        | n' < toInteger nFacs  = facs V.! fromIntegral n
        | otherwise             = lnGamma (fromIntegral n+1)
        where
            n' = toInteger n
            nFacs       = 2000 -- limited only by time and space
            facs        = V.map lnGamma (V.enumFromN 1 nFacs)


-- |Incomplete gamma functions.  Minimal definition is either 'p' or 'q', preferably both.
class Gamma a => IncGamma a where
    -- |Lower gamma function: lowerGamma s x == integral from 0 to x of 
    -- @\t -> t**(s-1) * exp (negate t)@
    lowerGamma :: a -> a -> a
    lowerGamma s x = exp (lnLowerGamma s x)
    -- |Natural log of lower gamma function
    lnLowerGamma :: a -> a -> a 
    lnLowerGamma s x = lnGamma s + log (p s x)
    -- |Regularized lower incomplete gamma function: lowerGamma z / gamma z
    p :: a -> a -> a
    p s x = 1 - q s x
    
    -- |Upper gamma function: lowerGamma s x == integral from x to infinity of 
    -- @\t -> t**(s-1) * exp (negate t)@
    upperGamma :: a -> a -> a
    upperGamma s x = exp (lnUpperGamma s x)
    -- |Natural log of upper gamma function
    lnUpperGamma :: a -> a -> a
    lnUpperGamma s x = lnGamma s + log (q s x)
    -- |Regularized upper incomplete gamma function: upperGamma z / gamma z
    q :: a -> a -> a
    q s x = 1 - p s x

instance IncGamma Float where
    lowerGamma   s x = realToFrac $ (lowerGamma   :: Double -> Double -> Double) (realToFrac s) (realToFrac x)
    lnLowerGamma s x = realToFrac $ (lnLowerGamma :: Double -> Double -> Double) (realToFrac s) (realToFrac x)
    p s x = realToFrac $ (p :: Double -> Double -> Double) (realToFrac s) (realToFrac x)
    
    upperGamma   s x = realToFrac $ (upperGamma   :: Double -> Double -> Double) (realToFrac s) (realToFrac x)
    lnUpperGamma s x = realToFrac $ (lnUpperGamma :: Double -> Double -> Double) (realToFrac s) (realToFrac x)
    q s x = realToFrac $ (q :: Double -> Double -> Double) (realToFrac s) (realToFrac x)
instance IncGamma Double where
    lowerGamma = lowerGammaHypGeom
    p s x
        | x < 0     = 1 - qNeg s x
        | x == 0    = 0
        | x >= s+1  = 1 - q s x
        | otherwise = pHypGeom s x
    
    q s x
        | x < 0     = qNeg s x
        | x == 0    = 1
        | x < s+1   = 1 - p s x
        | otherwise =
            converge . concat
            $ modifiedLentz 1e-30 (qCF s x)

-- |Factorial function
class Num a => Factorial a where
    factorial :: Integral b => b -> a
    factorial = fromInteger . factorial

instance Factorial Integer where
    factorial n = product [1..toInteger n]

instance Factorial Float where
    factorial = realToFrac . (factorial :: Integral a => a -> Double)
instance Factorial Double where
    factorial n
        | n < 0         = error "factorial: n < 0"
        | n < nFacs     = facs V.! fromIntegral n
        | otherwise     = infinity
        where
            nFacs :: Num a => a
            nFacs       = 171 -- any more is pointless, everything beyond here is "Infinity"
            facs        = V.scanl (*) 1 (V.enumFromN 1 nFacs)
            infinity    = facs V.! nFacs

beta :: Gamma a => a -> a -> a
beta z w = exp (lnGamma z + lnGamma w - lnGamma (z+w))
