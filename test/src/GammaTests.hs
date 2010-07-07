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
        x ~= y = (errBy magnitude x y <= eps)

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

err a b = errBy abs a b
errBy abs a b 
    | isNaN relErr  = absErr
    | otherwise     = min absErr relErr
    where
        absErr = abs (a-b)
        relErr = absErr / max (abs a) (abs b)
