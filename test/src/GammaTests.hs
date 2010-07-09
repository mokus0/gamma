{-# LANGUAGE ImplicitParams #-}
module GammaTests where

import Reference

import Control.Applicative
import Data.Complex
import Math.Gamma
import Test.Framework (testGroup, Test)
import Test.Framework.Providers.QuickCheck2 (testProperty)
import Test.QuickCheck

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
    [ testGroup "Float"
        [ testGroup "Native"  (realTests gamma (lnGamma :: Float -> Float))
        , testGroup "Complex" (complexTests gamma (lnGamma :: Complex Float -> Complex Float))
        ]
    , testGroup "Double"
        [ testGroup "Native"  (realTests    gamma (lnGamma :: Double -> Double))
        , testGroup "FFI"     (realTests    tgamma lgamma)
        , testGroup "Complex" (complexTests gamma (lnGamma :: Complex Double -> Complex Double))
        ]
    ]

realTests gamma lnGamma = 
    let ?mag = abs
        ?eps = eps
        ?complex = False
     in [ testGroup "gamma"       (realGammaTests    gamma)
        , testGroup "lnGamma"     (realLogGammaTests gamma lnGamma)
        , testGroup "lnFactorial" (logFactorialTests lnGamma lnFactorial)
        ]

complexTests gamma lnGamma = 
    let ?mag = magnitude
        ?eps = eps
        ?complex = True
     in [ testGroup "gamma"   (complexGammaTests    gamma)
        , testGroup "lnGamma" (complexLogGammaTests gamma lnGamma)
        , testGroup "lnFactorial" (logFactorialTests lnGamma lnFactorial)
        ]

realGammaTests gamma = 
    gammaTests gamma id (const 0) ++
    [ testProperty "between factorials" $ \(Positive x) -> 
        let gam x = fromInteger (product [1..x-1])
            gamma_x = gamma x `asTypeOf` eps
         in x > 2 && isSane gamma_x
            ==> gam (floor x) <= gamma_x && gamma_x <= gam (ceiling x)
    , testProperty "agrees with factorial" $ \(Positive x) ->
        let gam x = fromInteger (product [1..x-1])
            gamma_x = gamma (fromInteger x)
         in isSane gamma_x ==> 
            let ?eps = 16*eps in gam x ~= gamma_x
    , testProperty "agrees with C tgamma" $ \x ->
        let a = gamma x
            b = realToFrac (tgamma (realToFrac x))
         in isSane a ==> 
            let ?eps = 512*eps in a ~= b
    , testProperty "agrees with reference implementation" $ \x ->
        let a = gamma x
         in isSane a ==> snd (err gamma x) <= 256*eps
    , testProperty "monotone when x>2" $ \(Positive x) ->
        let x' = x * (1+eps)
            a = gamma x
            b = gamma x'
         in (x > 2) && (x <= x) && all isSane [a,b] ==> a <= b
    ]

complexGammaTests gamma = 
    gammaTests gamma realPart imagPart ++
    [ testProperty "conjugate" $ \x -> 
        let gam = gamma x
         in isSane (magnitude gam) ==> conjugate gam ~= gamma (conjugate x)
    , testProperty "real argument" $ \(Positive x) ->
        let z = x :+ 0
            gam = Math.Gamma.gamma x
         in isSane gam ==> 
            let ?mag = abs; ?eps = 512 * eps
             in gam ~= realPart (gamma z)
    ]
    
gammaTests gamma real imag =
    [ testProperty "increment arg" $ \x ->
        let a = gamma x
            b = gamma (x + 1)
            margin  | ?complex  = 32
                    | otherwise = 32
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
            let ?eps = 32 * eps
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
             in let ?eps = 512 * eps
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

logFactorialTests lnGamma lnFactorial =
    [ testProperty "agrees with lnGamma" $ \x ->
        let gam = lnGamma (fromInteger x + 1)
         in isSane (?mag gam)
            ==> gam ~= lnFactorial x
    ]
