{-# LANGUAGE ForeignFunctionInterface, ImplicitParams #-}
module GammaTests where

import Control.Applicative
import Data.Complex
import Math.Gamma
import Test.Framework (testGroup, Test)
import Test.Framework.Providers.QuickCheck2 (testProperty)
import Test.QuickCheck

foreign import ccall "math.h lgamma" lgamma :: Double -> Double
foreign import ccall "math.h tgamma" tgamma :: Double -> Double

-- instance (RealFloat a, Arbitrary a) => Arbitrary (Complex a) where
--     arbitrary = liftA2 (:+) arbitrary arbitrary

eps :: RealFloat a => a
eps = eps'
    where
        eps' = encodeFloat 1 (1 - floatDigits eps')

infix 4 ~=
x ~= y 
    =  absErr <= ?eps
    || absErr <= ?eps * min (?mag x) (?mag y)
    where absErr = ?mag (x-y)


isSane x = all (\f -> not (f x)) [isNaN, isInfinite, isDenormalized]

tests = 
    [ testGroup "gamma"
        [ testGroup "Float"             (realGammaTests    (gamma  :: Float  -> Float ))
        , testGroup "Double"            (realGammaTests    (gamma  :: Double -> Double))
        , testGroup "Double (tgamma)"   (realGammaTests    (tgamma :: Double -> Double))
        , testGroup "Complex Float"     (complexGammaTests (gamma  :: Complex Float  -> Complex Float ))
        , testGroup "Complex Double"    (complexGammaTests (gamma  :: Complex Double -> Complex Double))
        ]
    , testGroup "lnGamma"
        [ testGroup "Float"             (realLogGammaTests     gamma (lnGamma :: Float          -> Float         ))
        , testGroup "Double"            (realLogGammaTests     gamma (lnGamma :: Double         -> Double        ))
        , testGroup "Double (lgamma)"   (realLogGammaTests    tgamma (lgamma  :: Double         -> Double        ))
        , testGroup "Complex Float"     (complexLogGammaTests  gamma (lnGamma :: Complex Float  -> Complex Float ))
        , testGroup "Complex Double"    (complexLogGammaTests  gamma (lnGamma :: Complex Double -> Complex Double))
        ]
    , testGroup "lnFactorial"
        [ testGroup "Float"          (logFactorialTests abs       (eps :: Float))
        , testGroup "Double"         (logFactorialTests abs       (eps :: Double))
        , testGroup "Complex Float"  (logFactorialTests magnitude (eps :: Float))
        , testGroup "Complex Double" (logFactorialTests magnitude (eps :: Double))
        ]
    ]

realGammaTests gamma = 
    let ?mag = abs
        ?complex = False
     in gammaTests gamma id (const 0) ++
        [ testProperty "between factorials" $ \(Positive x) -> 
            let gam x = fromInteger (product [1..x-1])
                gamma_x = gamma x `asTypeOf` eps
             in x > 2 && isSane gamma_x
                ==> gam (floor x) <= gamma_x && gamma_x <= gam (ceiling x)
        , testProperty "agrees with C tgamma" $ \x ->
            let a = gamma x
                b = realToFrac (tgamma (realToFrac x))
             in isSane a ==> 
                let ?eps = 512*eps in a ~= b
        ]

complexGammaTests gamma = 
    let ?mag = magnitude
        ?eps = eps
        ?complex = True
     in gammaTests gamma realPart imagPart ++
        [ testProperty "conjugate" $ \x -> 
            let gam = gamma x
             in isSane (magnitude gam) ==> conjugate gam ~= gamma (conjugate x)
        , testProperty "real argument" $ \(Positive x) ->
            let z = x :+ 0
                gam = Math.Gamma.gamma x
             in isSane gam ==> (gam :+ 0) ~= gamma z
        ]

gammaTests gamma real imag =
    [ testProperty "increment arg" $ \x ->
        let a = gamma x
            b = gamma (x + 1)
            margin  | ?complex  = 8
                    | otherwise = 4
         in all (isSane . ?mag) [a,b,recip a, recip b]
            ==> ?mag (a - b/x) <= margin * (max 2 (1 + recip (?mag x))) * eps * ?mag a
             || ?mag (a*x - b) <= margin * (max 2 (1 +        ?mag x))  * eps * ?mag b
    , testProperty "reflect" $ \x ->
        ?mag x > 0 ==>
        let a = gamma x
            b = gamma (1 - x)
            c = sin (pi * x) / pi
            margin  | ?complex  = 16
                    | otherwise = 256
         in all (isSane . ?mag) [a,b,c]
            -- There may be tighter bounds to be found but I haven't
            -- been able to derive them yet.
            ==> ?mag (a*b*c-1) <= margin * eps * (1 + ?mag c * (1 + ?mag (a+b)))
             || ?mag (a*b*c-1) <= margin * eps * (?mag (a*b) + ?mag (a*c) + ?mag (b*c))
    ]

logGammaTests gamma lnGamma real imag =
    [ testProperty "increment arg" $ \x ->
        let gam = lnGamma (x+1)
         in real x > 0 && isSane (?mag gam)
            ==> 
            let ?eps = 16 * eps
             in gam - log x ~= lnGamma x
                || gam ~= log x + lnGamma x
    , testProperty "reflect" $ \x ->
        ?mag x > 0 ==>
        let a = lnGamma x
            b = lnGamma (1 - x)
            c = log pi - c';    c' = log (sin (pi * x))
         in all (isSane . ?mag) [a,b,c] 
            ==> 
            let ?eps = 512 * eps
             in     a + b ~= c
                 || a ~= c - b
                 || b ~= c - a
                 || a - b - c' ~= log pi
    , testProperty "agrees with log . gamma" $ \x ->
        let a = log b'    ; a' = exp b
            b = lnGamma x ; b' = gamma x
            margin  | ?complex   = 1024
                    | otherwise  = 16
         in   (real x > 1 || abs (imag x) > 1)
            && all (isSane . ?mag) [a,b,a',b'] ==>
            let ?eps = margin * eps
             in a ~= b || a' ~= b'
    ]

realLogGammaTests gamma lnGamma = 
    let ?mag = abs
        ?complex = False
     in logGammaTests gamma lnGamma id (const 0) ++
        [ testProperty "between factorials" $ \(Positive x) -> 
            let gam x = sum $ map (log.fromInteger) [1..x-1]
                gamma_x = lnGamma x `asTypeOf` eps
             in x > 2 && isSane gamma_x
                ==> gam (floor x) <= gamma_x && gamma_x <= gam (ceiling x)
        , let ?eps = 2 * eps
           in testProperty "agrees with C lgamma" $ \(NonNegative x) ->
            let a = lnGamma x
                b = realToFrac (lgamma (realToFrac x))
             in isSane a ==> a ~= b
        ]

complexLogGammaTests gamma lnGamma = 
    let ?mag = magnitude
        ?complex = True
     in logGammaTests gamma lnGamma realPart imagPart ++
        [ testProperty "real argument" $ \(Positive x) ->
            let z = x :+ 0
                gam = Math.Gamma.lnGamma x
             in isSane gam ==> 
                let ?eps = 8 * eps
                 in (gam :+ 0) ~= lnGamma z
        ]

logFactorialTests abs eps =
    [ testProperty "agrees with lnGamma" $ \x ->
        let gam = lnGamma (fromInteger x + 1)
         in isSane (abs gam)
            ==> gam ~= lnFactorial x
    ]
    where
        infix 4 ~=
        x ~= y 
            =  absErr <= eps
            || absErr <= eps * min (abs x) (abs y)
            where absErr = abs (x-y)
