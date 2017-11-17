import numpy as np
from bggrad import BgGrad


class TestBasic:
    params = {}

    def test_basic(self):
        np.random.seed(42)
        BgGrad.noise_variance = None
        bg = BgGrad(**self.params)
        img = 128 + 64 * np.random.rand(5, 5)
        img = img.astype(np.uint8)
        for i in range(100):
            noise = 10 * np.random.randn(5, 5)
            bg.segment(img + noise.astype(np.uint8))
        assert bg.segment(img).max() < 128
        assert bg.segment(img).max() < 128
        img[1,1] += 42
        assert bg.segment(img).max() > 128
        img[1,1] -= 42
        assert bg.segment(img).max() < 128


class TestExposureCheck(TestBasic):
    params = {'maxfcnt': 90, 'exposure_check': True}

class TestNoiseVariance(TestBasic):
    params = {'noise_variance': 100}

class TestDeshake(TestBasic):
    params = {'deshake': 1, 'noise_variance': 100}

class TestDeshake(TestBasic):
    params = {'deshake': 1, 'noise_variance': 100, 'blocksize': 2}
