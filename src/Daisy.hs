module Daisy

where

import GHC.Generics

data GWB a = GWB { g :: a, w :: a, b :: a} deriving (Eq, Functor, Generic, Show)

instance Foldable GWB where
  foldMap f (GWB a b c) = f a <> f b <> f c

instance Applicative GWB where
  pure a = GWB a a a
  (GWB fa fb fc) <*> (GWB a b c) = GWB (fa a) (fb b) (fc c)

data Config = Config {
  x0 :: GWB Double,
  albedo :: GWB Double,
  q :: Int,
  gamma :: Double,
  tOpt :: Double,
  tmax :: Int,
  dt :: Double,
  s :: Double,
  sigma :: Double
}

defaultConfig :: Config
defaultConfig = Config
  (GWB 1 0 0)
  (GWB 0.5 0.75 0.25)
  20
  0.3
  22.5
  200
  0.01
  917
  5.67e-8

steps :: Integral b => Config -> b
steps cfg = floor (fromIntegral (tmax cfg)/dt cfg )

type Time = Double

l :: Config -> Time -> Double
l cfg t = 0.6 + 1.2 * (t / fromIntegral (tmax cfg))

albedoP :: Config -> GWB Double -> Double
albedoP cfg alpha = sum $ (*) <$> alpha <*> albedo cfg

-- ((S*L(t)*(1-albedo_p(t, alpha_w, alpha_b, alpha_g)))/(sigma))**(0.25) - 273

tP :: Config -> Time -> GWB Double -> Double
tP cfg t alpha = (s cfg * l cfg t * (1 - albedoP cfg alpha) / sigma cfg) ** 0.25 - 273

tGWB :: Config -> Time -> GWB Double -> GWB Double
tGWB cfg t alpha = (\x -> fromIntegral (q cfg) * (albedoP cfg alpha - x) + tP cfg t alpha) <$> albedo cfg

beta :: Config -> Time -> GWB Double -> GWB Double
beta cfg t alpha = fmap (\x -> 1 - ((1/17.5)**2) * (tOpt cfg - x) ** 2) (tGWB cfg t alpha)

alphaDot :: Config -> Time -> GWB Double -> GWB Double
alphaDot cfg t alpha = GWB g' w' b'
  where
    w' = w alpha * (g alpha * w (beta cfg t alpha) - gamma cfg) + 0.001
    b' = b alpha * (g alpha * b (beta cfg t alpha) - gamma cfg) + 0.001
    g' = - w' - b'

data Results = Results { resTime :: Double, resAlpha :: GWB Double} deriving (Eq, Show, Generic)

runge :: Config -> Results -> Results
runge cfg (Results t gwb) = Results (t+1) gwbIntegrated
  where
    k = fmap (dt cfg *) (alphaDot cfg t gwb)
    k' = fmap (dt cfg *) (alphaDot cfg t ((\gwb' gwbk -> gwb' + 0.5 * gwbk) <$> gwb <*> k))
    k'' = fmap (dt cfg *) (alphaDot cfg t ((\gwb' gwbk -> gwb' + 0.5 * gwbk) <$> gwb <*> k'))
    k''' = fmap (dt cfg *) (alphaDot cfg t ((\gwb' gwbk -> gwb' + 1 * gwbk) <$> gwb <*> k''))
    gwbIntegrated = (\x0 x1 x2 x3 x4 -> x0 + (x1 + 2 * x2 + 2 * x3 + x4) / 6) <$> gwb <*> k <*> k' <*> k'' <*> k'''

tempP :: Config -> Results -> Double
tempP cfg x = (s cfg * l cfg (resTime x) * (1 - albedoP cfg (resAlpha x))/sigma cfg) ** 0.25 - 273

tempB :: Config -> Time -> Double
tempB cfg t = (s cfg * l cfg t * (1 - 0.5)/sigma cfg) ** 0.25 - 273

alphas :: Config -> [Results]
alphas cfg = iterate (runge cfg) (Results 0 (x0 cfg))
