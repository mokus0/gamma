module GammaTests where

import Control.Applicative
import Data.Complex
import Math.Gamma
import Test.Framework (testGroup)
import Test.Framework.Providers.QuickCheck2 (testProperty)
import Test.QuickCheck

-- instance (RealFloat a, Arbitrary a) => Arbitrary (Complex a) where
--     arbitrary = liftA2 (:+) arbitrary arbitrary

eps :: RealFloat a => a
eps = eps'
    where
        eps' = encodeFloat 1 (1 - floatDigits eps')

isSane x = all (\f -> not (f x)) [isNaN, isInfinite, isDenormalized]

tests = 
    [ testGroup "gamma"
        [ testGroup "Float"          (realGammaTests    (1024  * eps :: Float))
        , testGroup "Double"         (realGammaTests    (512   * eps :: Double))
        , testGroup "Complex Float"  (complexGammaTests (16636 * eps :: Float))
        , testGroup "Complex Double" (complexGammaTests (16636 * eps :: Double))
        ]
    , testGroup "lnGamma"
        [ testGroup "Float"          (realLogGammaTests    (1024  * eps :: Float))
        , testGroup "Double"         (realLogGammaTests    (512   * eps :: Double))
        , testGroup "Complex Float"  (complexLogGammaTests (16636 * eps :: Float))
        , testGroup "Complex Double" (complexLogGammaTests (16636 * eps :: Double))
        ]
    , testGroup "lnFactorial"
        [ testGroup "Float"          (logFactorialTests abs       (eps :: Float))
        , testGroup "Double"         (logFactorialTests abs       (eps :: Double))
        , testGroup "Complex Float"  (logFactorialTests magnitude (eps :: Float))
        , testGroup "Complex Double" (logFactorialTests magnitude (eps :: Double))
        ]
    ]

realGammaTests eps = gammaTests abs id (const 0) eps ++
    [ testProperty "between factorials" $ \(Positive x) -> 
        let gam x = fromInteger (product [1..x-1])
            gamma_x = gamma x `asTypeOf` eps
         in x > 2 && isSane gamma_x
            ==> gam (floor x) <= gamma_x && gamma_x <= gam (ceiling x)
    ]

complexGammaTests eps = gammaTests magnitude realPart imagPart eps ++
    [ testProperty "conjugate" $ \x -> 
        let gam = gamma x
         in isSane (magnitude gam) ==> conjugate gam ~= gamma (conjugate x)
    , testProperty "real argument" $ \(Positive x) ->
        let z = x :+ 0
            gam = gamma x
         in isSane gam ==> (gam :+ 0) ~= gamma z
    ]
    where
        infix 4 ~=
        x ~= y 
            =  absErr <= eps
            || absErr <= eps * min (magnitude x) (magnitude y)
            where absErr = magnitude (x-y)

gammaTests abs real imag eps =
    [ testProperty "increment arg" $ \x ->
        let gam = gamma (x+1)
         in abs x > 0 && isSane (abs gam)
            ==> gam / x ~= gamma x
             || gam ~= x * gamma x
    , testProperty "reflect" $ \x ->
        abs x > 0 ==>
        let a = gamma x
            b = gamma (1 - x)
            c = pi / sin (pi * x)
         in all (isSane.abs) [a,b,c] 
            ==> a*b ~= c
             || a ~= c/b
             || b ~= c/a
             || a*b*sin(pi*x) ~= pi
    ]
    where
        infix 4 ~=
        x ~= y 
            =  absErr <= eps
            || absErr <= eps * min (abs x) (abs y)
            where absErr = abs (x-y)

logGammaTests abs real imag eps =
    [ testProperty "increment arg" $ \x ->
        let gam = lnGamma (x+1)
         in real x > 0 && isSane (abs gam)
            ==> gam - log x ~= lnGamma x
             || gam ~= log x + lnGamma x
    , testProperty "reflect" $ \x ->
        abs x > 0 ==>
        let a = lnGamma x
            b = lnGamma (1 - x)
            c = log (pi / sin (pi * x))
         in all (isSane.abs) [a,b,c] 
            ==> a + b ~= c
             || a ~= c - b
             || b ~= c - a
             || a - b - log (sin(pi*x)) ~= pi
    ]
    where
        infix 4 ~=
        x ~= y 
            =  absErr <= eps
            || absErr <= eps * min (abs x) (abs y)
            where absErr = abs (x-y)


realLogGammaTests eps = logGammaTests abs id (const 0) eps ++
    [ testProperty "between factorials" $ \(Positive x) -> 
        let gam x = sum $ map (log.fromInteger) [1..x-1]
            gamma_x = lnGamma x `asTypeOf` eps
         in x > 2 && isSane gamma_x
            ==> gam (floor x) <= gamma_x && gamma_x <= gam (ceiling x)
    ]

complexLogGammaTests eps = logGammaTests magnitude realPart imagPart eps ++
    [ testProperty "real argument" $ \(Positive x) ->
        let z = x :+ 0
            gam = lnGamma x
         in isSane gam ==> (gam :+ 0) ~= lnGamma z
    ]
    where
        infix 4 ~=
        x ~= y 
            =  absErr <= eps
            || absErr <= eps * min (magnitude x) (magnitude y)
            where absErr = magnitude (x-y)


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
