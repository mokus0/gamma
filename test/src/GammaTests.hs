module GammaTests where

import Control.Applicative
import Data.Complex
import Math.Gamma
import Test.Framework (testGroup)
import Test.Framework.Providers.QuickCheck2 (testProperty)
import Test.QuickCheck

instance (RealFloat a, Arbitrary a) => Arbitrary (Complex a) where
    arbitrary = liftA2 (:+) arbitrary arbitrary

eps :: RealFloat a => a
eps = eps'
    where
        eps' = encodeFloat 1 (1 - floatDigits eps')

isSane x = all (\f -> not (f x)) [isNaN, isInfinite, isDenormalized]

tests = 
    [ testGroup "gamma"
        [ testGroup "Float"          (realGammaTests    (256  * eps :: Float))
        , testGroup "Double"         (realGammaTests    (512  * eps :: Double))
        , testGroup "Complex Float"  (complexGammaTests (2048 * eps :: Float))
        , testGroup "Complex Double" (complexGammaTests (4096 * eps :: Double))
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
    [
    ]

gammaTests abs real imag eps =
    [ testProperty "increment arg" $ \x ->
        abs x > 1 ==>
        let gam = gamma (x+1)
         in isSane (abs gam) ==> gam ~= x * gamma x
    , testProperty "reflect" $ \x ->
        abs x > 1 ==>
        let a = gamma x
            b = gamma (1 - x)
            c = pi / sin (pi * x)
         in all (isSane.abs) [a,b,c] 
            ==> a*b ~= c
             || a*b*sin(pi*x) ~= pi
    ]
    where
        infix 4 ~=
        x ~= y = (errBy abs x y <= eps)

err a b = errBy abs a b
errBy abs a b 
    | isNaN relErr  = absErr
    | otherwise     = min absErr relErr
    where
        absErr = abs (a-b)
        relErr = absErr / max (abs a) (abs b)
