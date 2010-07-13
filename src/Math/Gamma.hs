{-# LANGUAGE FlexibleInstances, TemplateHaskell #-}
module Math.Gamma
    ( Gamma(..)
    , Factorial(..)
    , IncGamma(..)
    , beta
    ) where

import Math.Gamma.Lanczos
import Math.Gamma.Incomplete

import Data.Complex
import Data.List (sortBy, findIndex)
import Data.Ord (comparing)
import qualified Data.Vector.Unboxed as V
import Language.Haskell.TH (litE, Lit(IntegerL))
import Math.ContinuedFraction
import Math.Sequence.Converge

-- |Gamma function.  Minimal definition is ether 'gamma' or 'lnGamma'.
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

floatGammaInfCutoff :: Double
floatGammaInfCutoff = $( do
        let Just cutoff = findIndex isInfinite (scanl (*) (1::Float) [1..])
        litE (IntegerL (1 + toInteger cutoff))
    )

instance Gamma Float where
    gamma = realToFrac . gam . realToFrac
        where
            gam :: Double -> Double
            gam x 
                | x >= floatGammaInfCutoff  = 1/0
                | otherwise = case properFraction x of
                (n,0) | n < 1     -> 0/0
                      | otherwise -> factorial (n-1)
                _     | x < (-20) -> let s = pi / sin (pi * x)
                                      in signum s * exp (log (abs s) - lnGamma (1-x))
                      | otherwise -> reflect (gammaLanczos g cs) x
            
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

doubleGammaInfCutoff :: Double
doubleGammaInfCutoff = $( do
        let Just cutoff = findIndex isInfinite (scanl (*) (1::Double) [1..])
        litE (IntegerL (1 + toInteger cutoff))
    )

instance Gamma Double where
    gamma x 
        | x >= doubleGammaInfCutoff  = 1/0
        | otherwise = case properFraction x of
        (n,0) | n < 1     -> 0/0
              | otherwise -> factorial (n-1)
        _     | x < (-50) -> let s = pi / sin (pi * x)
                              in signum s * exp (log (abs s) - lnGamma (1-x))
              | otherwise -> reflect (gammaLanczos g cs) x
        where
            g = 2*pi
            cs = [9.99999999999985803099845245112243194630609100328127e-1,3.11601175054147181416018805783854078966968054820574e2,-4.98651190460363940327267435809982353665915898566851e2,2.44084728999768763967887353973264446808195804052725e2,-3.86703646430741942964006117268082794993800219289745e1,1.33509001013705483157804869520994460219627829058483e0,-1.89772218995656823839469050681023636391151679063349e-3,8.47526461434914856885138928027675324645466808374844e-7,2.59715567376858005948233185894865034384724173599575e-7,-2.71664378506075157744515910362825103696629060955620e-7,6.15111480613629944231662042128164100055404439541778e-8]
            -- cs = [1.0000000000000002,311.60117505414695,-498.65119046033163,244.08472899875767,-38.67036462939322,1.3350899103585203,-1.8972831806242229e-3,-3.935368195357295e-7,2.592464641764731e-6,-3.2263565156368265e-6,2.5666169886566876e-6,-1.3737776806198937e-6,4.4551204024819644e-7,-6.576826592057796e-8]

            -- g = 9
            -- cs = [1.00000000000000017466330156623527316911425637266631e0,5.71640018827434137913574603541678950268182092108649e3,-1.48153042676841390904407305200012978786350039702570e4,1.42914927765747855402510960116666846723365741870550e4,-6.34816021764145881328945494441494058092852514365618e3,1.30160828605832187410470515341489270916658800187122e3,-1.08176705351436963467921849912127524570831371271716e2,2.60569650561175582772877838762204314602400379296189e0,-7.42345251020141615152744454168760537276562314755790e-3,5.38413643250956406296099893482669586535250446433247e-8,-4.02353314126823637206733624633774281557006605066199e-9]

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
    gamma = complexDoubleToFloat . gamma . complexFloatToDouble
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
    factorial n
        | n < 0     = error "factorial: n < 0"
        | otherwise = product [1..toInteger n]

instance Factorial Float where
    factorial = realToFrac . (factorial :: Integral a => a -> Double)
instance Factorial Double where
    factorial n
        | n < 0         = 0/0
        | n < nFacs     = facs V.! fromIntegral n
        | otherwise     = infinity
        where
            nFacs :: Num a => a
            nFacs       = 171 -- any more is pointless, everything beyond here is "Infinity"
            facs        = V.scanl (*) 1 (V.enumFromN 1 nFacs)
            infinity    = facs V.! nFacs

-- |The beta function: @beta z w@ == @gamma z * gamma w / gamma (z+w)@
beta :: Gamma a => a -> a -> a
beta z w = exp (lnGamma z + lnGamma w - lnGamma (z+w))
