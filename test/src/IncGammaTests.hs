{-# LANGUAGE ImplicitParams #-}
module IncGammaTests where

import GammaTests (eps, isSane, (~=))

import Data.Complex
import Math.Gamma
import Test.Framework (testGroup)
import Test.Framework.Providers.QuickCheck2 (testProperty)
import Test.QuickCheck

tests = 
    [ testGroup "incomplete gamma"
        [ testGroup "Float"  (incompleteGammaTests (eps :: Float))
        , testGroup "Double" (incompleteGammaTests (eps :: Double))
        ]
    ]

incompleteGammaTests eps =
    let ?mag = abs
     in [ testProperty "lowerGamma + upperGamma" $ \s x ->
            let a = lowerGamma s x
                b = upperGamma s x
                c = gamma s
             in all isSane [a,b,c] ==> 
                let ?eps = 512 * eps
                 in  a+b ~= c
                    || a ~= c-b
                    || b ~= c-a
        , testProperty "p + q" $ \s x ->
            let a = p s x
                b = q s x
             in all isSane [a,b] ==>
                let ?eps = 256 * eps
                 in  a+b ~= 1
                    || a ~= 1-b
                    || b ~= 1-a
        , testGroup "upperGamma"
            [ testProperty "increment s" $ \s x ->
                let a = upperGamma (s+1) x
                    b = s * upperGamma s x
                    c = x ** s * exp (-x)
                 in all isSane [a,b,c] ==>
                    let ?eps = 2048*eps
                     in a ~= b+c || a-b ~= c || a-c ~= b
            , testProperty "upperGamma _ 0" $ \s ->
                let a = upperGamma s 0
                    b = gamma s
                 in all isSane [a,b] ==>
                    let ?eps = 512*eps
                     in a ~= b
            , testProperty "upperGamma 1 _" $ \x ->
                let ?eps = 256*eps
                 in x > 0 ==> upperGamma 1 x ~= exp (-x)
            ]
        , testGroup "lowerGamma"
            [ testProperty "increment s" $ \s x ->
                let a = lowerGamma (s+1) x
                    b = s * lowerGamma s x
                    c = x ** s * exp (-x)
                 in all isSane [a,b,c] ==>
                    let ?eps = 2048*eps
                     in a ~= b-c || a+c ~= b || b-a ~= c
            ]
        ]
