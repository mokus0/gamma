#!/usr/bin/env runhaskell
module Main where

import qualified GammaTests (tests)
import qualified IncGammaTests (tests)

import Test.Framework (defaultMain, testGroup)

main = defaultMain 
    [ testGroup "Math.Gamma"            GammaTests.tests
    , testGroup "Math.Gamma.Incomplete" IncGammaTests.tests
    ]