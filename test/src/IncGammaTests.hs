{-# LANGUAGE ImplicitParams #-}
module IncGammaTests where

import GammaTests (eps, isSane, (~=))

import Data.Complex
import Data.Number.Erf
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
     in [ testProperty "lowerGamma + upperGamma" $ \s (NonNegative x) ->
            let a = lowerGamma s x
                b = upperGamma s x
                c = gamma s
             in all isSane [a,b,c] ==> 
                let ?eps = 512 * eps
                 in  a+b ~= c
                    || a ~= c-b
                    || b ~= c-a
        , testProperty "p + q" $ \s (NonNegative x) ->
            let a = p s x
                b = q s x
             in all isSane [a,b] ==>
                let ?eps = 256 * eps
                 in  a+b ~= 1
                    || a ~= 1-b
                    || b ~= 1-a
        , testGroup "upperGamma"
            [ testProperty "increment s" $ \s (NonNegative x) ->
                let a = upperGamma (s+1) x
                    b = s * upperGamma s x
                    c = x ** s * exp (-x)
                 in all isSane [a,b,c] ==>
                    let ?eps = 1024*eps
                     in a ~= b+c || a-b ~= c || a-c ~= b
            , testProperty "x = 0" $ \s ->
                let a = upperGamma s 0
                    b = gamma s
                 in all isSane [a,b] ==>
                    let ?eps = 512*eps
                     in a ~= b
            , testProperty "s = 1, x > 0" $ \(Positive x) ->
                let ?eps = 256*eps
                 in x /= 0 ==> upperGamma 1 x ~= exp (-x)
--            , testProperty "s = 1, x < 0" $ \(Positive x) ->
--                let ?eps = 256*eps
--                 in x /= 0 ==> upperGamma 1 (-x) ~= exp x
            , testProperty "s = 0.5" $ \(Positive x) ->
                let ?eps = 128 * eps
                 in upperGamma 0.5 (x*x) / sqrt pi + erf x ~= 1
            ]
        , testGroup "lnUpperGamma"
            [ testProperty "s = 1, x > 0" $ \(Positive x) ->
                let a = lnUpperGamma 1 x
                 in let ?eps = 1024 * eps
                     in isSane a ==> a ~= (-x)
--            , testProperty "s = 1, x < 0" $ \(Positive x) ->
--                let a = lnUpperGamma 1 (-x)
--                 in let ?eps = 1024 * eps
--                     in isSane a ==> a ~= x
            , testProperty "agrees with upperGamma" $ \s (NonNegative x) ->
                let a = lnUpperGamma s x;   a' = log b
                    b = upperGamma s x;     b' = exp a
                    
                 in all isSane [a,a',b,b'] ==> 
                    let ?eps = 128*eps
                     in a ~= a' || b ~= b'
            ]
        , testGroup "lowerGamma"
            [ testProperty "increment s" $ \s (NonNegative x) ->
                let a = lowerGamma (s+1) x
                    b = s * lowerGamma s x
                    c = x ** s * exp (-x)
                 in not (s==0 || x==0) && all isSane [a,b,c] ==>
                    let ?eps = 1024*eps
                     in a ~= b-c || a+c ~= b || b-a ~= c
            , testProperty "s = 0.5" $ \(NonNegative x) ->
                let ?eps = 16 * eps
                 in lowerGamma 0.5 (x*x) / sqrt pi + erfc x ~= 1
            , testProperty "s = 1, x > 0" $ \(Positive x) ->
                let ?eps = 256*eps
                 in x /= 0 ==> lowerGamma 1 x + exp (-x) ~= 1
--            , testProperty "s = 1, x < 0" $ \(Positive x) ->
--                let ?eps = 256*eps
--                 in x /= 0 ==> lowerGamma 1 (-x) + exp x ~= 1
            ]
        ]
