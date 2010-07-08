module IncGammaTests where

import GammaTests (eps, isSane)

import Data.Complex
import Math.Gamma
import Test.Framework (testGroup)
import Test.Framework.Providers.QuickCheck2 (testProperty)
import Test.QuickCheck

tests = 
    [ testGroup "incomplete gamma"
        [ testGroup "Float"  (incompleteGammaTests (128 * eps :: Float))
        , testGroup "Double" (incompleteGammaTests (512 * eps :: Double))
        ]
    ]

incompleteGammaTests eps =
    [ testProperty "lowerGamma + upperGamma" $ \s x ->
        let a = lowerGamma s x
            b = upperGamma s x
            c = gamma s
         in all isSane [a,b,c] 
            ==> a+b ~= c
             || a ~= c-b
             || b ~= c-a
    , testProperty "p + q" $ \s x ->
        let a = p s x
            b = q s x
         in all isSane [a,b] 
            ==> a+b ~= 1
             || a   ~= 1-b
             || b   ~= 1-a
    ]
    where
        infix 4 ~=
        x ~= y 
            =  absErr <= eps
            || absErr <= eps * min (abs x) (abs y)
            where absErr = abs (x-y)
